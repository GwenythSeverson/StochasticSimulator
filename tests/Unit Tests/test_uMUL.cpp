#include <gtest/gtest.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdio>

#include "../../modules/uMUL.hpp"
#include "../../modules/multiplier.hpp"
#include "../../general_functions.hpp"
#include "../../bsg/lfsr.hpp"
#include "../../bsg/sng.hpp"
#include "../../bsg/sobol.hpp"

// Unit tests for the uGEMM Figure 3(a) unipolar uMUL -- the LOADED form, matching the synthesized
// uMUL_uni.sv: in_0 is a unary bitstream, in_1 is a binary number in a register.
//
// Two kinds of test in here, kept deliberately separate:
//   STRUCTURAL -- exact assertions about the circuit (what the enable does, what the register
//                 holds, what the AND gate can emit). These catch a rewiring bug.
//   ACCURACY   -- statistical, so they carry a tolerance. in_0 streams are built from an LFSR
//                 with an EXPLICIT seed; generate_valid_stream() seeds itself off the clock and
//                 would make these flaky. Tolerances come from the measured sweep at the bottom.
//
// Note what is NOT in here, and cannot be: a correlated-operand test. in_1 is a number, so there
// is no second stream for in_0 to correlate with. The failure mode the AND gate has simply does
// not exist in this form of the unit.

using namespace StochasticSimulator;

// =============================================================================================
// IS THE TOP BIT A SIGN BIT?  NO -- and this test exists to show why it looks like one.
// =============================================================================================
// The value register is an UNSIGNED magnitude in [0, 2^width - 1], read as value/2^width, and
// the comparator is a plain "value > random". There is no sign bit, no two's complement, and no
// negative operand anywhere in the unipolar unit.
//
// What makes bit 7 look special in the poster figure is the OPERAND, not the hardware:
//
//     a flip can only CLEAR a bit that was SET   -> the value falls
//     a flip can only SET a bit that was CLEAR   -> the value rises
//
// The poster uses p1 = 0.5, i.e. value = 128 = 1000_0000, which has exactly ONE bit set -- the
// MSB. So seven of the eight flips push the answer UP and only one pushes it DOWN, and the one
// that goes down is the heaviest bit. Change the operand and the pattern moves with it.
//
// The magnitude is always the same though: flipping bit b moves the answer by
//     |delta| = (2^b / 2^width) * p0
// so the MSB is worth p0/2 = 0.25 here REGARDLESS of the operand. Only the SIGN depends on the
// data. That is the honest statement, and it is what this test pins down.
// SUITE NAME IS DELIBERATELY *NOT* UnaryMultiplierTest. That name belongs to the TEST_F fixture
// further down, and googletest refuses to run a suite that mixes TEST and TEST_F -- it fails
// EVERY test in the suite with "All tests in the same test suite must use the same test fixture".
// This test and RegisterBitDamageCountsSobolPointsCrossed below need no fixture (they build their
// own multiplier per run), so they live in their own suite instead.
TEST(UnaryMultiplierFaultTest, TopBitIsMagnitudeNotSign) {
    constexpr unsigned WIDTH = 8;
    constexpr int N = 256;
    constexpr int ONES = 128;                 // p0 = 0.5, so 128 enabled cycles

    auto run = [&](uint64_t value) {
        UnaryMultiplier m(WIDTH, 1);
        m.load_value(value);
        int ones = 0;
        // in_0 packed at the front: position never matters, only the enabled count.
        for (int i = 0; i < N; ++i) if (m.multiply(i < ONES)) ++ones;
        return ones;
    };

    std::cout << "\n  uMUL single-flip damage vs the operand it is carrying (p0 = 0.5)\n"
              << "  value  p1      bits       flip b=7   flip b=6   flip b=5   set bits\n";

    for (uint64_t base : {uint64_t{128}, uint64_t{64}, uint64_t{192}, uint64_t{255},
                          uint64_t{1}, uint64_t{200}}) {
        const int clean = run(base);
        std::cout << "   " << std::setw(4) << base << "  "
                  << std::fixed << std::setprecision(3)
                  << (static_cast<double>(base) / 256.0) << "  ";
        for (int b = 7; b >= 0; --b) std::cout << ((base >> b) & 1);
        std::cout << " ";

        int setbits = 0;
        for (int b = 0; b < 8; ++b) if ((base >> b) & 1) ++setbits;

        for (int b : {7, 6, 5}) {
            const uint64_t flipped = base ^ (1ull << b);
            const double delta = static_cast<double>(run(flipped) - clean) / N;
            std::cout << std::setw(11) << std::showpos << std::setprecision(4) << delta
                      << std::noshowpos;

            // THE MAGNITUDE IS FIXED BY THE BIT, NOT THE DATA. Flipping bit b moves the value by
            // 2^b, which moves the answer by (2^b / 256) * p0, up to the comparator's 1/256
            // quantisation.
            const double expect_mag = (static_cast<double>(1u << b) / 256.0) * 0.5;
            EXPECT_NEAR(std::fabs(delta), expect_mag, 1.5 / N)
                << "bit " << b << " magnitude on operand " << base;

            // THE SIGN IS FIXED BY THE DATA. A set bit can only be cleared (down); a clear bit
            // can only be set (up). If this were a sign bit, bit 7 would go one way regardless.
            const bool was_set = ((base >> b) & 1) != 0;
            if (was_set) EXPECT_LT(delta, 0.0) << "clearing a set bit must lower the operand";
            else         EXPECT_GT(delta, 0.0) << "setting a clear bit must raise the operand";
        }
        std::cout << std::setw(9) << setbits << "\n";
    }
    std::cout << std::endl;
}

// ================================================================================================
// WHY THE PER-BIT DAMAGE LADDER IS NOT A CLEAN 2^b, AND EXACTLY WHEN IT IS
// ================================================================================================
//
// The poster's register panel plots |delta| for each register bit and calls the result a
// geometric ladder, 1x 1x 2x 4x 8x 16x 32x 64x in units of one stream bit. Bits 0 and 1 land on
// the SAME value, which looks like a broken ladder. It is not, and "output quantisation" is too
// vague an explanation. The exact law is:
//
//     |delta| in OUTPUT BITS  =  the number of consumed Sobol points lying in the interval
//                                between the two operand values
//
// because the output is out(v) = #{ j < n1 : v > sobol[j] }, so out(v2) - out(v1) counts exactly
// the points the operand swept past. Flipping bit b moves the operand between m and m + 2^b where
// m is v with bit b cleared -- a DYADIC-ALIGNED interval of width 2^b.
//
// NOW THE SOBOL PROPERTY DOES THE REST. The first n1 = 2^k points of a width-w Sobol sequence form
// a (0,k,1)-net in base 2: every dyadic interval of length 2^(w-k) holds exactly one point. So an
// aligned interval of width 2^b holds exactly
//
//     2^b / 2^(w-k) = 2^(b-w+k)  points,  PROVIDED b >= w-k
//
// and for b < w-k the interval is narrower than the net's resolution, so it holds 0 or 1 -- the
// ladder cannot resolve those bits and they all collapse onto the same rung.
//
//     THE THRESHOLD IS  b* = w - k = log2(2^w / n1) = log2(1/p0).
//
// AT p0 = 0.5 (w=8, k=7) b* = 1, so exactly ONE bit -- bit 0 -- is sub-resolution, and it shares
// bit 1's rung. That is the whole of the "anomaly", and it IS situational: at p0 = 0.25 the
// threshold moves to b* = 2 and bits 0 AND 1 both collapse. This test measures both.
//
// The 2^b/2^w * p0 formula quoted on the figure is therefore exact for b >= log2(1/p0) and an
// over-estimate below it. Nothing here depends on the operand VALUE -- only on the stream density.
TEST(UnaryMultiplierFaultTest, RegisterBitDamageCountsSobolPointsCrossed) {
    constexpr unsigned WIDTH = 8;
    constexpr int N = 256;
    constexpr int SCALE = 1 << WIDTH;         // 256 operand values, 0..255

    // One real gate run. in_0 packed at the front -- position never matters, only the count.
    auto run = [&](uint64_t value, int ones) {
        UnaryMultiplier m(WIDTH, 1);
        m.load_value(value);
        int k = 0;
        for (int i = 0; i < N; ++i) if (m.multiply(i < ones)) ++k;
        return k;
    };

    for (int n1 : {128, 64}) {
        const double p0 = static_cast<double>(n1) / N;

        // out[v] for every operand, measured. out[SCALE] is not expressible in an 8-bit register,
        // but it is definitionally n1: an operand of 256 exceeds every point in [0,255].
        std::vector<int> out(SCALE + 1);
        for (int v = 0; v < SCALE; ++v) out[v] = run(static_cast<uint64_t>(v), n1);
        out[SCALE] = n1;

        // RECOVER THE CONSUMED POINT SET BY DIFFERENCING, rather than reimplementing Sobol here.
        // out(v+1) - out(v) = #{ j < n1 : sobol[j] == v }, so this reads the generator's actual
        // behaviour out of real gate runs. A reimplementation could drift from the module; this
        // cannot.
        std::vector<int> consumed(SCALE);
        int total = 0;
        for (int v = 0; v < SCALE; ++v) {
            consumed[v] = out[v + 1] - out[v];
            EXPECT_GE(consumed[v], 0) << "output must be monotone in the operand";
            EXPECT_LE(consumed[v], 1) << "a Sobol point cannot be consumed twice";
            total += consumed[v];
        }
        ASSERT_EQ(total, n1) << "exactly n1 points should be consumed at p0 = " << p0;

        // ---- THE INTERVAL LAW, over every operand and every bit -------------------------------
        for (int v = 0; v < SCALE; ++v) {
            for (int b = 0; b < static_cast<int>(WIDTH); ++b) {
                const int m = v & ~(1 << b);            // interval base: v with bit b cleared
                int points = 0;
                for (int u = m; u < m + (1 << b); ++u) points += consumed[u];

                const int measured = std::abs(out[v ^ (1 << b)] - out[v]);
                ASSERT_EQ(measured, points)
                    << "damage from bit " << b << " on operand " << v << " at p0 = " << p0
                    << " is not the count of Sobol points crossed";
            }
        }

        // ---- THE NET PROPERTY, i.e. WHICH RUNGS ARE EXACT -------------------------------------
        // b* = log2(2^w / n1). At or above it the damage is a constant 2^(b-b*) for EVERY operand;
        // below it the bit is finer than the net can resolve.
        int bstar = 0;
        while ((n1 << bstar) < SCALE) ++bstar;

        std::cout << "\n  p0 = " << std::fixed << std::setprecision(3) << p0
                  << "   n1 = " << n1 << "   sub-resolution below bit b* = " << bstar << "\n"
                  << "   bit   ideal 2^b/2^w*p0   measured |delta|   in stream bits   exact?\n";

        for (int b = 0; b < static_cast<int>(WIDTH); ++b) {
            // Damage from bit b, over every operand. If the net property holds it is constant.
            int lo = N, hi = 0;
            for (int v = 0; v < SCALE; ++v) {
                const int d = std::abs(out[v ^ (1 << b)] - out[v]);
                lo = std::min(lo, d);  hi = std::max(hi, d);
            }
            const double ideal = (static_cast<double>(1 << b) / SCALE) * p0;

            std::cout << "   " << std::setw(3) << b << std::setw(16) << std::setprecision(6)
                      << ideal << std::setw(17)
                      << (static_cast<double>(hi) / N) << std::setw(16) << hi
                      << (b >= bstar ? "     yes\n" : "     no, floored\n");

            if (b >= bstar) {
                // AT OR ABOVE THE RESOLUTION LIMIT the damage is exact and operand-independent.
                EXPECT_EQ(lo, hi) << "bit " << b << " damage varies with the operand at p0 = " << p0;
                EXPECT_EQ(hi, (n1 << b) / SCALE)
                    << "bit " << b << " is not worth 2^b * p0 output bits at p0 = " << p0;
                EXPECT_NEAR(static_cast<double>(hi) / N, ideal, 1e-12)
                    << "the figure's (2^b/2^w)*p0 formula should be EXACT here";
            } else {
                // BELOW IT the bit cannot be resolved: at most one point sits in the interval, so
                // the rung is 0 or 1 output bits and the ideal formula OVER-estimates nothing --
                // it under-estimates, because a whole bit is the smallest step that exists.
                EXPECT_LE(hi, 1) << "a sub-resolution bit cannot move the answer by 2 bits";
                EXPECT_LT(ideal, 1.0 / N)
                    << "bit " << b << " should be finer than one output bit at p0 = " << p0;
            }
        }

        // THE POSTER'S HEADLINE, asserted: the top bit is worth half the operand's full scale, so
        // p0/2 of the answer, whatever the operand holds. 0.25 at p0 = 0.5.
        for (int v = 0; v < SCALE; ++v) {
            const double d = std::abs(out[v ^ (1 << (WIDTH - 1))] - out[v]) / static_cast<double>(N);
            EXPECT_NEAR(d, p0 / 2.0, 1e-12) << "top-bit damage on operand " << v;
        }
    }
    std::cout << std::endl;
}

class UnaryMultiplierTest : public ::testing::Test {
protected:
    static const std::size_t N = 1024;
    static const uint16_t SCALE = 1024;

    // Measured max error over the 361-point sweep below is 0.0035 for an LFSR in_0 and 0.0020 for
    // a Sobol in_0. This leaves several times that.
    static constexpr double TOL = 0.012;

    std::vector<bool> make_stream(double p, uint16_t seed) const {
        BitstreamGenerator bsg;
        FlexibleLFSR lfsr(StreamLength::Length_1024, seed);
        std::vector<bool> s;
        s.reserve(N);
        for (std::size_t i = 0; i < N; ++i) {
            s.push_back(bsg.generate_bit(p, lfsr.next(), SCALE));
        }
        return s;
    }

    // Two streams off ONE LFSR, thresholded differently: maximally correlated (SCC = 1). Only
    // needed to show what the plain AND gate does with them.
    void make_correlated_pair(double p0, double p1,
                              std::vector<bool>& a, std::vector<bool>& b) const {
        BitstreamGenerator bsg;
        FlexibleLFSR lfsr(StreamLength::Length_1024, 7);
        a.clear();
        b.clear();
        for (std::size_t i = 0; i < N; ++i) {
            uint16_t r = lfsr.next();
            a.push_back(bsg.generate_bit(p0, r, SCALE));
            b.push_back(bsg.generate_bit(p1, r, SCALE));
        }
    }

    double and_gate(const std::vector<bool>& a, const std::vector<bool>& b) const {
        Multiplier m;
        std::vector<bool> o;
        o.reserve(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) o.push_back(m.multiply(a[i], b[i]));
        return calculate_probability(o);
    }

    std::vector<bool> constant_stream(bool v) const { return std::vector<bool>(N, v); }
};

// =====================================================================================
// STRUCTURAL -- the circuit itself
// =====================================================================================

// The enable [en in Figure 3(a)] is in_0. G's RNG must step exactly once per in_0 one, and never
// on an in_0 zero. This is the single property that makes uMUL correlation-immune.
TEST_F(UnaryMultiplierTest, EnableFiresExactlyOnInput0Ones) {
    std::vector<bool> a = make_stream(0.3, 11);

    UnaryMultiplier umul(10);
    umul.load_probability(0.7);
    for (std::size_t i = 0; i < N; ++i) {
        uint64_t before = umul.get_rng_index();
        umul.multiply(a[i]);
        EXPECT_EQ(umul.get_rng_index() != before, a[i]) << "at cycle " << i;
    }
    EXPECT_EQ(umul.get_enabled_cycles(), count_ones(a));
}

// C is a register, not an estimator: loaded once, unchanged by running, and it survives reset()
// so one operand can be run against several in_0 streams.
TEST_F(UnaryMultiplierTest, LoadedValueIsExactAndSurvivesTheRun) {
    UnaryMultiplier umul(10);
    umul.load_probability(0.5);
    EXPECT_EQ(umul.get_value(), 512u);
    EXPECT_DOUBLE_EQ(umul.get_probability(), 0.5);

    std::vector<bool> a = make_stream(0.4, 5);
    for (std::size_t i = 0; i < N; ++i) umul.multiply(a[i]);
    EXPECT_EQ(umul.get_value(), 512u) << "running must not disturb the operand register";

    umul.reset();
    EXPECT_EQ(umul.get_value(), 512u) << "reset() clears the run, not the operand";
    EXPECT_EQ(umul.get_rng_index(), 0u);
    EXPECT_EQ(umul.get_enabled_cycles(), 0u);

    // p1 = 1 SATURATES to 2^width - 1. The register is width bits, exactly like `iB` in
    // uMUL_uni.sv, so 2^width does not fit and p1 = 1.0 is not an expressible operand.
    umul.load_probability(1.0);
    EXPECT_EQ(umul.get_value(), umul.get_entry() - 1);
    umul.load_probability(0.0);
    EXPECT_EQ(umul.get_value(), 0u);
}

// load_from_stream() is Figure 3(a)'s counter run to completion up front: it counts every one in
// the stream, exactly, then scales once.
TEST_F(UnaryMultiplierTest, LoadFromStreamCountsEveryOne) {
    std::vector<bool> b = make_stream(0.625, 9);
    uint64_t ones = count_ones(b);

    UnaryMultiplier umul(10);
    umul.load_from_stream(b);
    EXPECT_EQ(umul.get_value(), (ones << 10) / N);

    umul.load_from_stream(constant_stream(true));
    EXPECT_EQ(umul.get_value(), umul.get_entry() - 1) << "all-ones saturates; 2^width does not fit";
    umul.load_from_stream(constant_stream(false));
    EXPECT_EQ(umul.get_value(), 0u);
}

// The output is an AND with in_0, so it can never emit a one where in_0 had a zero.
TEST_F(UnaryMultiplierTest, OutputIsAlwaysSubsetOfInput0) {
    std::vector<bool> a = make_stream(0.35, 3);
    std::vector<bool> out = umul_stream(a, 0.8);

    ASSERT_EQ(out.size(), a.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        if (out[i]) EXPECT_TRUE(a[i]) << "emitted a one at cycle " << i << " where in_0 was zero";
    }
}

TEST_F(UnaryMultiplierTest, DeterministicAcrossReset) {
    std::vector<bool> a = make_stream(0.45, 31);

    UnaryMultiplier umul(10);
    umul.load_probability(0.55);
    std::vector<bool> first, second;
    for (std::size_t i = 0; i < N; ++i) first.push_back(umul.multiply(a[i]));
    umul.reset();
    for (std::size_t i = 0; i < N; ++i) second.push_back(umul.multiply(a[i]));
    EXPECT_EQ(first, second);
}

// width only has to cover the run without the Sobol sequence wrapping: ceil(log2(N)).
TEST_F(UnaryMultiplierTest, WidthForStreamLengthIsCeilLog2) {
    EXPECT_EQ(UnaryMultiplier::width_for_stream_length(256), 8u);
    EXPECT_EQ(UnaryMultiplier::width_for_stream_length(257), 9u);
    EXPECT_EQ(UnaryMultiplier::width_for_stream_length(1024), 10u);
    EXPECT_EQ(UnaryMultiplier::width_for_stream_length(4096), 12u);
    EXPECT_LE(UnaryMultiplier::width_for_stream_length(1ull << 40), UnaryMultiplier::MAX_WIDTH);

    UnaryMultiplier sized = UnaryMultiplier::for_stream_length(1024);
    EXPECT_EQ(sized.get_width(), 10u);
    EXPECT_EQ(sized.get_entry(), 1024u);
}

TEST_F(UnaryMultiplierTest, RejectsMalformedInput) {
    std::vector<bool> empty;
    EXPECT_THROW(umul_stream(empty, 0.5), std::invalid_argument);
    EXPECT_THROW(umul_stream(constant_stream(true), 1.5), std::invalid_argument);
    EXPECT_THROW(umul_stream(constant_stream(true), -0.1), std::invalid_argument);
    EXPECT_THROW(UnaryMultiplier(0), std::invalid_argument);
    EXPECT_THROW(UnaryMultiplier(UnaryMultiplier::MAX_WIDTH + 1), std::invalid_argument);
    EXPECT_THROW(UnaryMultiplier::width_for_stream_length(0), std::invalid_argument);

    UnaryMultiplier umul(10);
    EXPECT_NO_THROW(umul.load_value(1023));            // 2^width - 1 is the top of the register
    EXPECT_THROW(umul.load_value(1024), std::invalid_argument);  // 2^width does not fit in 10 bits
    EXPECT_THROW(umul.load_from_stream(empty), std::invalid_argument);
}

// =====================================================================================
// ACCURACY
// =====================================================================================

// The worked case: 896/1024 as the enable stream, 512/1024 loaded as the binary operand.
TEST_F(UnaryMultiplierTest, WorkedCase896Times512) {
    std::vector<bool> a = make_stream(896.0 / 1024.0, 101);

    UnaryMultiplier umul(10);
    umul.load_value(512);
    ASSERT_EQ(umul.get_value(), 512u);

    std::vector<bool> out;
    for (std::size_t i = 0; i < N; ++i) out.push_back(umul.multiply(a[i]));

    // G lands 512/1024 of the enabled cycles, and the enable fires on in_0's ones, so the output
    // density is the product. 896 never enters the comparator -- it only gates it.
    EXPECT_EQ(umul.get_enabled_cycles(), count_ones(a));
    EXPECT_NEAR(calculate_probability(out), (896.0 / 1024.0) * (512.0 / 1024.0), TOL);
}

// With in_1 loaded there is no estimator, so the ZERO rail is bit-exact. The ONE rail is not
// reachable at all: the value register is `width` bits, so the largest operand is 2^width - 1 and
// the densest stream G can generate is (2^width - 1)/2^width. That is the RTL's behaviour --
// uMUL_uni.sv's 8-bit iB tops out at 255/256 for exactly the same reason -- not a rounding slip.
TEST_F(UnaryMultiplierTest, ZeroRailIsExactAndOneRailSaturates) {
    std::vector<bool> ones = constant_stream(true);
    std::vector<bool> zeros = constant_stream(false);

    // Multiplying by zero, or enabling with zero, is bit-perfect.
    EXPECT_DOUBLE_EQ(umul_probability(ones, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(umul_probability(zeros, 1.0), 0.0);
    EXPECT_DOUBLE_EQ(umul_probability(zeros, 0.0), 0.0);

    // p1 = 1.0 saturates to 1023/1024. Exactly one cycle in the 1024 is lost: the one where the
    // Sobol draw is 1023 and the strict ">" fails against a value of 1023.
    EXPECT_DOUBLE_EQ(umul_probability(ones, 1.0), 1023.0 / 1024.0);
    EXPECT_EQ(count_ones(umul_stream(ones, 1.0)), N - 1);

    // WHERE the lost bit falls is exactly determined, not statistical. Sobol dimension 1 visits
    // the value 1023 once per period, at index 682, so p1 = 1.0 behaves as a perfect pass-through
    // until the enable has fired 683 times and drops precisely one bit thereafter.
    std::vector<bool> sparse = make_stream(0.375, 13);   // 384 ones: never reaches index 682
    ASSERT_LT(count_ones(sparse), 683u);
    EXPECT_EQ(umul_stream(sparse, 1.0), sparse) << "below 683 enabled cycles it is a pass-through";

    std::vector<bool> dense = make_stream(0.875, 13);    // 896 ones: passes index 682
    ASSERT_GT(count_ones(dense), 683u);
    EXPECT_EQ(count_ones(umul_stream(dense, 1.0)), count_ones(dense) - 1)
        << "past 683 enabled cycles exactly one bit is lost -- never more";
}

TEST_F(UnaryMultiplierTest, MultipliesAgainstAKnownOperand) {
    for (double p0 = 0.1; p0 <= 0.9; p0 += 0.2) {
        for (double p1 = 0.1; p1 <= 0.9; p1 += 0.2) {
            std::vector<bool> a = make_stream(p0, 101);
            EXPECT_NEAR(umul_probability(a, p1), p0 * p1, TOL)
                << "at p0 = " << p0 << ", p1 = " << p1;
        }
    }
}

// in_0's own generator does not matter -- the enable trick means the unit never has to reason
// about in_0's structure, only about which cycles it lets through.
TEST_F(UnaryMultiplierTest, InsensitiveToInput0StreamSource) {
    for (double p0 = 0.1; p0 <= 0.9; p0 += 0.2) {
        for (double p1 = 0.1; p1 <= 0.9; p1 += 0.2) {
            std::vector<bool> sobol_a = generate_sobol_stream(p0, N, 2);
            EXPECT_NEAR(umul_probability(sobol_a, p1), p0 * p1, TOL)
                << "sobol in_0, at p0 = " << p0 << ", p1 = " << p1;
        }
    }
}

// The plain AND gate collapses to min(p0, p1) when its two operands are correlated. uMUL cannot
// be made to do that, because there is no second stream to correlate with in_0 in the first place.
TEST_F(UnaryMultiplierTest, BeatsThePlainAndGateOnCorrelatedStreams) {
    std::vector<bool> a, b;
    make_correlated_pair(0.5, 0.5, a, b);

    EXPECT_NEAR(and_gate(a, b), 0.5, 0.05);                 // the AND gate reports min(p0, p1)
    EXPECT_NEAR(umul_probability(a, b), 0.25, TOL);         // uMUL counts b, then multiplies
    EXPECT_NEAR(umul_probability(a, 0.5), 0.25, TOL);       // same answer from the value directly
}

// =====================================================================================
// CHARACTERIZATION -- prints the numbers the header quotes
// =====================================================================================

TEST_F(UnaryMultiplierTest, CharacterizeAccuracy) {
    struct Acc {
        double sq = 0.0; double max = 0.0; int n = 0;
        void add(double e) { sq += e * e; if (e > max) max = e; ++n; }
        double rmse() const { return n ? std::sqrt(sq / n) : 0.0; }
    };

    Acc u_lfsr, u_sobol, a_ind, a_cor;

    for (int i = 1; i <= 19; ++i) {
        for (int j = 1; j <= 19; ++j) {
            double p0 = i * 0.05, p1 = j * 0.05;
            double expect = p0 * p1;

            std::vector<bool> a = make_stream(p0, 101);
            u_lfsr.add(std::fabs(umul_probability(a, p1) - expect));

            std::vector<bool> sa = generate_sobol_stream(p0, N, 2);
            u_sobol.add(std::fabs(umul_probability(sa, p1) - expect));

            std::vector<bool> b = make_stream(p1, 211);
            a_ind.add(std::fabs(and_gate(a, b) - expect));

            std::vector<bool> ca, cb;
            make_correlated_pair(p0, p1, ca, cb);
            a_cor.add(std::fabs(and_gate(ca, cb) - expect));
        }
    }

    std::printf("\n  N = 1024, width 10, 361-point p0 x p1 sweep over 0.05..0.95\n");
    std::printf("  %-34s %8s %8s\n", "", "RMSE", "max");
    std::printf("  %-34s %8.4f %8.4f\n", "  uMUL, in_0 from an LFSR", u_lfsr.rmse(), u_lfsr.max);
    std::printf("  %-34s %8.4f %8.4f\n", "  uMUL, in_0 from Sobol", u_sobol.rmse(), u_sobol.max);
    std::printf("  %-34s %8.4f %8.4f\n", "  AND gate, independent operands", a_ind.rmse(), a_ind.max);
    std::printf("  %-34s %8.4f %8.4f\n", "  AND gate, SCC = 1", a_cor.rmse(), a_cor.max);
    std::printf("\n");

    EXPECT_LT(u_lfsr.rmse(), a_cor.rmse() / 10.0);
}
