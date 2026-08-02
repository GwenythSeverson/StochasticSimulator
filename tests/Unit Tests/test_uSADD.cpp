/*
 * Unit tests for uSADD, the uGEMM unary scaled adder.
 *
 * These are the structural checks: exact arithmetic, register sizing, the collapse property the
 * exhaustive fault campaigns depend on, and the fault-injection hooks. Accuracy sweeps live in
 * the functional test next door.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "../../modules/uSADD.hpp"

using namespace StochasticSimulator;

namespace {

// A stream of `length` bits carrying exactly `ones` ones, scattered at random.
std::vector<bool> scattered(std::size_t length, std::size_t ones, std::mt19937& rng) {
    std::vector<bool> s(length, false);
    std::vector<std::size_t> idx(length);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    for (std::size_t i = 0; i < ones; ++i) s[idx[i]] = true;
    return s;
}

std::size_t count_ones(const std::vector<bool>& s) {
    std::size_t n = 0;
    for (bool b : s) if (b) ++n;
    return n;
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// CYCLE-BY-CYCLE TRACE -- prints the datapath so the operation can be checked by hand
// ---------------------------------------------------------------------------------------------
// Not an assertion-light test: every line it prints is also checked. It exists so the unit can be
// verified against the uGEMM paper / UnarySim RTL by reading, rather than by trusting this file.
TEST(uSADDTest, PrintCycleByCycleTrace) {
    const std::vector<bool> a = {1, 0, 1, 1, 0, 0, 1, 1};   // 5 ones
    const std::vector<bool> b = {0, 1, 1, 0, 1, 0, 0, 1};   // 4 ones
    const int n = 2;

    UnaryScaledAdder unit(n);
    std::cout << "\n  uSADD trace, n = 2, 8 cycles.  a = 5/8 = 0.625, b = 4/8 = 0.500\n"
              << "  exact scaled sum (0.625 + 0.500)/2 = 0.5625\n\n"
              << "  cyc | a b | PC s | A before | A+s | emit | A after | running out\n"
              << "  ----+-----+------+----------+-----+------+---------+------------\n";

    std::size_t emitted = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const uint32_t before = unit.get_accumulator();
        const uint32_t s = static_cast<uint32_t>(a[i]) + static_cast<uint32_t>(b[i]);
        const bool out = unit.add(a[i], b[i]);
        if (out) ++emitted;

        // The three invariants of the datapath, checked on every single cycle.
        EXPECT_EQ(out, (before + s) >= static_cast<uint32_t>(n)) << "emit rule broke at cycle " << i;
        EXPECT_EQ(unit.get_accumulator(), out ? (before + s - n) : (before + s))
            << "drain rule broke at cycle " << i;
        EXPECT_LT(unit.get_accumulator(), static_cast<uint32_t>(n))
            << "residue left A outside [0, n) at cycle " << i;

        std::cout << "  " << std::setw(3) << i << " |  " << a[i] << " " << b[i]
                  << "  |  " << s << "   |    " << before
                  << "     |  " << (before + s) << "  |  " << (out ? "1" : "0")
                  << "   |    " << unit.get_accumulator()
                  << "    | " << emitted << "\n";
    }

    const std::size_t total_in = count_ones(a) + count_ones(b);
    std::cout << "\n  total input ones = " << total_in
              << ",  output ones = " << emitted
              << " = floor(" << total_in << "/2)\n"
              << "  residue left in A = " << unit.get_accumulator()
              << "  (the whole of the healthy unit's error)\n"
              << "  reported value " << emitted << "/8 = "
              << (static_cast<double>(emitted) / a.size())
              << " vs exact 0.5625, short by " << unit.get_accumulator() << "/(2*8)\n"
              << std::endl;

    EXPECT_EQ(emitted, total_in / 2);
    EXPECT_EQ(unit.get_accumulator(), total_in % 2);
}

// ---------------------------------------------------------------------------------------------
// SIZING
// ---------------------------------------------------------------------------------------------
TEST(uSADDTest, AccumulatorHoldsTheResidueAndNoMore) {
    // uGEMM takes the output from the CARRY of the accumulator, so the register stores only the
    // residue [0, n-1] -- the wide pre-drain sum lives in the adder, not in a flip-flop. The
    // register is therefore ceil(log2(n)) bits, NOT ceil(log2(2n)).
    for (unsigned n : {1u, 2u, 3u, 4u, 5u, 8u, 16u, 32u}) {
        const unsigned w = UnaryScaledAdder::width_for_inputs(n);
        EXPECT_GE((1u << w) - 1u, n - 1u)
            << "width " << w << " cannot hold the residue n-1 = " << (n - 1) << " for n=" << n;
        // ...and is the SMALLEST such width, so no flip-flop is wasted on the fault surface.
        if (w > 1) {
            EXPECT_LT((1u << (w - 1)) - 1u, n - 1u)
                << "width " << w << " is one bit wider than the residue needs for n=" << n;
        }
    }
    // The case the poster measures: a two-input scaled adder has a ONE-BIT accumulator.
    EXPECT_EQ(UnaryScaledAdder::width_for_inputs(2), 1u);
    EXPECT_EQ(UnaryScaledAdder(2).get_acc_width(), 1u);
}

TEST(uSADDTest, RejectsAnAccumulatorTooNarrowForTheResidue) {
    // n = 2 leaves a residue of at most 1, so one bit is exactly right and is accepted. n = 4
    // leaves up to 3, so one bit must be refused rather than silently dropping credit.
    EXPECT_NO_THROW(UnaryScaledAdder(2, 1));
    EXPECT_THROW(UnaryScaledAdder(4, 1), std::invalid_argument);
    EXPECT_NO_THROW(UnaryScaledAdder(4, 2));
    EXPECT_THROW(UnaryScaledAdder(0), std::invalid_argument);
    EXPECT_THROW(UnaryScaledAdder(UnaryScaledAdder::MAX_INPUTS + 1), std::invalid_argument);
}

TEST(uSADDTest, PcBusWidthHoldsZeroThroughN) {
    for (unsigned n : {1u, 2u, 3u, 4u, 7u, 8u, 9u, 16u}) {
        UnaryScaledAdder unit(n);
        EXPECT_GE((1u << unit.get_pc_width()) - 1u, n)
            << "PC bus cannot represent a full house of " << n << " ones";
    }
}

// ---------------------------------------------------------------------------------------------
// EXACT ARITHMETIC
// ---------------------------------------------------------------------------------------------
TEST(uSADDTest, RailsAreExact) {
    // All zeros in -> all zeros out. All ones in -> all ones out, since s = n every cycle drains
    // exactly one emission's worth of credit.
    std::vector<bool> zeros(64, false), ones(64, true);
    EXPECT_EQ(count_ones(usadd_stream(zeros, zeros)), 0u);
    EXPECT_EQ(count_ones(usadd_stream(ones, ones)), 64u);
    // One rail against the other: (1 + 0)/2 = 0.5 exactly, since the credit alternates.
    EXPECT_EQ(count_ones(usadd_stream(ones, zeros)), 32u);
}

TEST(uSADDTest, OutputIsFloorOfTotalOnesOverN) {
    // The defining identity. If this holds, everything else about the unit follows.
    std::mt19937 rng(20260801);
    for (int trial = 0; trial < 200; ++trial) {
        std::uniform_int_distribution<int> pick(0, 64);
        const std::vector<bool> a = scattered(64, pick(rng), rng);
        const std::vector<bool> b = scattered(64, pick(rng), rng);
        const std::size_t total = count_ones(a) + count_ones(b);
        EXPECT_EQ(count_ones(usadd_stream(a, b)), total / 2)
            << "output was not floor(total ones / 2)";
    }
}

TEST(uSADDTest, ResidueIsTheOnlyHealthyError) {
    // The healthy unit's only error is the unspent credit left in A at the end: strictly less
    // than n ones, and always a DEFICIT, never a surplus.
    std::mt19937 rng(7);
    UnaryScaledAdder unit(2);
    for (int trial = 0; trial < 100; ++trial) {
        unit.reset();
        std::uniform_int_distribution<int> pick(0, 128);
        const std::vector<bool> a = scattered(128, pick(rng), rng);
        const std::vector<bool> b = scattered(128, pick(rng), rng);
        for (std::size_t i = 0; i < 128; ++i) unit.add(a[i], b[i]);

        const std::size_t total = count_ones(a) + count_ones(b);
        const double exact = static_cast<double>(total) / 2.0;
        const double got = static_cast<double>(unit.get_emitted_ones());
        EXPECT_LE(got, exact);                 // never overshoots
        EXPECT_LT(exact - got, 2.0);           // and is short by less than n
        EXPECT_LT(unit.get_accumulator(), 2u); // residue lives in A, in [0, n)
    }
}

TEST(uSADDTest, ThreeAndFourInputScaling) {
    // The unit is general: n inputs scale by 1/n.
    std::vector<bool> all(60, true), none(60, false);
    // (1 + 1 + 0)/3 = 2/3 -> 40 of 60
    EXPECT_EQ(count_ones(usadd_stream({all, all, none})), 40u);
    // (1 + 0 + 0 + 0)/4 = 1/4 -> 15 of 60
    EXPECT_EQ(count_ones(usadd_stream({all, none, none, none})), 15u);
}

// ---------------------------------------------------------------------------------------------
// THE COLLAPSE -- the property every exhaustive fault campaign on this unit relies on
// ---------------------------------------------------------------------------------------------
TEST(uSADDTest, OutputDependsOnlyOnTheTotalOnesCount) {
    // Same ones-counts, wildly different arrangements, including adversarially correlated and
    // anti-correlated ones. If ANY of these disagree, the exhaustive enumeration in the poster
    // code is invalid.
    std::mt19937 rng(1337);
    const std::size_t N = 128;

    for (std::size_t na : {0u, 1u, 31u, 64u, 97u, 128u}) {
        for (std::size_t nb : {0u, 17u, 64u, 128u}) {
            std::size_t reference = 0;
            for (int shuffle = 0; shuffle < 12; ++shuffle) {
                std::vector<bool> a, b;
                if (shuffle == 0) {
                    // fully correlated: ones packed at the front of both
                    a.assign(N, false); b.assign(N, false);
                    for (std::size_t i = 0; i < na; ++i) a[i] = true;
                    for (std::size_t i = 0; i < nb; ++i) b[i] = true;
                } else if (shuffle == 1) {
                    // anti-correlated: a packed at the front, b at the back
                    a.assign(N, false); b.assign(N, false);
                    for (std::size_t i = 0; i < na; ++i) a[i] = true;
                    for (std::size_t i = 0; i < nb; ++i) b[N - 1 - i] = true;
                } else {
                    a = scattered(N, na, rng);
                    b = scattered(N, nb, rng);
                }
                const std::size_t got = count_ones(usadd_stream(a, b));
                if (shuffle == 0) reference = got;
                else EXPECT_EQ(got, reference)
                    << "arrangement changed the answer at na=" << na << " nb=" << nb
                    << " (shuffle " << shuffle << ") -- the collapse does not hold";
            }
            EXPECT_EQ(reference, (na + nb) / 2);
        }
    }
}

TEST(uSADDTest, NoCorrelationPenaltyUnlikeTheMux) {
    // The headline claim: feeding the SAME stream to both inputs is the worst case for a MUX
    // adder and a non-event here. (0.5 + 0.5)/2 = 0.5 whether or not the operands are the same
    // bits.
    std::mt19937 rng(99);
    const std::vector<bool> s = scattered(256, 128, rng);
    EXPECT_EQ(count_ones(usadd_stream(s, s)), 128u);          // perfectly correlated
    const std::vector<bool> t = scattered(256, 128, rng);
    EXPECT_EQ(count_ones(usadd_stream(s, t)), 128u);          // independent
    std::vector<bool> inv(256);
    for (std::size_t i = 0; i < 256; ++i) inv[i] = !s[i];
    EXPECT_EQ(count_ones(usadd_stream(s, inv)), 128u);        // perfectly anti-correlated
}

// ---------------------------------------------------------------------------------------------
// FAULT INJECTION HOOKS
// ---------------------------------------------------------------------------------------------
TEST(uSADDTest, AccumulatorFlipIsBoundedAndOneTime) {
    // THE uSADD RESULT. A holds a residue, not an operand, so an upset perturbs the running
    // credit once and the unit keeps integrating correctly afterwards. The total output moves by
    // about (perturbation / n) bits and then stops moving -- unlike uMUL, where a value-register
    // upset rescales every remaining cycle.
    std::mt19937 rng(4242);
    const std::size_t N = 256;
    const std::vector<bool> a = scattered(N, 128, rng);
    const std::vector<bool> b = scattered(N, 128, rng);
    const std::size_t clean = count_ones(usadd_stream(a, b));

    for (unsigned bit = 0; bit < UnaryScaledAdder(2).get_acc_width(); ++bit) {
        for (std::size_t strike : {std::size_t{0}, N / 4, N / 2, N - 1}) {
            UnaryScaledAdder unit(2);
            std::size_t ones = 0;
            for (std::size_t i = 0; i < N; ++i) {
                if (i == strike) unit.flip_accumulator_bit(bit);
                if (unit.add(a[i], b[i])) ++ones;
            }
            const long long delta = static_cast<long long>(ones) - static_cast<long long>(clean);
            // One flip of bit b injects at most 2^b credit, worth at most 2^b / n output bits.
            const long long bound = (1LL << bit) / 2 + 1;
            EXPECT_LE(std::llabs(delta), bound)
                << "accumulator bit " << bit << " struck at cycle " << strike
                << " moved the output by " << delta << ", past the bounded-damage claim";
        }
    }
}

TEST(uSADDTest, PcFlipBehavesLikeInjectedCredit) {
    // A PC glitch is a one-cycle bus event but NOT a one-cycle error: A integrates the spurious
    // credit, so it behaves exactly as if extra input ones had arrived that cycle.
    std::mt19937 rng(5150);
    const std::size_t N = 256;
    const std::vector<bool> a = scattered(N, 128, rng);
    const std::vector<bool> b = scattered(N, 128, rng);
    const std::size_t clean = count_ones(usadd_stream(a, b));

    UnaryScaledAdder unit(2);
    std::size_t ones = 0;
    for (std::size_t i = 0; i < N; ++i) {
        UnaryScaledAdder::CycleUpset up;
        if (i == 100) up.pc_flip = 0x1;   // +/- 1 unit of credit on the low PC wire
        if (unit.add({a[i], b[i]}, up)) ++ones;
    }
    const long long delta = static_cast<long long>(ones) - static_cast<long long>(clean);
    EXPECT_LE(std::llabs(delta), 1) << "one unit of spurious credit moved the output by " << delta;
}

TEST(uSADDTest, OutFlipIsExactlyOneBit) {
    std::mt19937 rng(3);
    const std::size_t N = 64;
    const std::vector<bool> a = scattered(N, 32, rng);
    const std::vector<bool> b = scattered(N, 32, rng);
    const std::size_t clean = count_ones(usadd_stream(a, b));

    UnaryScaledAdder unit(2);
    std::size_t ones = 0;
    for (std::size_t i = 0; i < N; ++i) {
        UnaryScaledAdder::CycleUpset up;
        up.out_flip = (i == 7);
        if (unit.add({a[i], b[i]}, up)) ++ones;
    }
    // The output wire is downstream of all state, so a glitch there is confined to its own bit.
    EXPECT_EQ(std::llabs(static_cast<long long>(ones) - static_cast<long long>(clean)), 1);
}

TEST(uSADDTest, FlipBitRangeIsChecked) {
    UnaryScaledAdder unit(2);
    EXPECT_NO_THROW(unit.flip_accumulator_bit(unit.get_acc_width() - 1));
    EXPECT_THROW(unit.flip_accumulator_bit(unit.get_acc_width()), std::out_of_range);
}

TEST(uSADDTest, ResetClearsEverything) {
    std::mt19937 rng(11);
    UnaryScaledAdder unit(2);
    const std::vector<bool> a = scattered(32, 21, rng), b = scattered(32, 13, rng);
    for (std::size_t i = 0; i < 32; ++i) unit.add(a[i], b[i]);
    EXPECT_GT(unit.get_cycles(), 0u);
    unit.reset();
    EXPECT_EQ(unit.get_accumulator(), 0u);
    EXPECT_EQ(unit.get_cycles(), 0u);
    EXPECT_EQ(unit.get_emitted_ones(), 0u);
}

// ---------------------------------------------------------------------------------------------
// INTERFACE GUARDS
// ---------------------------------------------------------------------------------------------
TEST(uSADDTest, RejectsMalformedInput) {
    UnaryScaledAdder unit(2);
    EXPECT_THROW(unit.add(std::vector<bool>{true}), std::invalid_argument);
    EXPECT_THROW(unit.add(std::vector<bool>{true, true, true}), std::invalid_argument);

    UnaryScaledAdder three(3);
    EXPECT_THROW(three.add(true, false), std::invalid_argument);  // pair overload on a 3-input unit

    std::vector<bool> shortStream(8, true), longStream(16, true);
    EXPECT_THROW(usadd_stream(shortStream, longStream), std::invalid_argument);
    EXPECT_THROW(usadd_stream(std::vector<bool>{}, std::vector<bool>{}), std::invalid_argument);
}
