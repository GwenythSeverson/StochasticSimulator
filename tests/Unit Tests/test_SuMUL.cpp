#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "../../modules/SuMUL.hpp"
#include "../../modules/multiplier.hpp"
#include "../../general_functions.hpp"
#include "../../bsg/lfsr.hpp"
#include "../../bsg/sng.hpp"

// Unit tests for the uGEMM Figure 3(a) unipolar SuMUL.
//
// Two kinds of test in here, kept deliberately separate:
//   STRUCTURAL -- exact assertions about the circuit (what the enable does, what the AND gate can
//                 emit, that the counter is binary). These are the ones that catch a rewiring bug.
//   ACCURACY   -- statistical, so they carry a tolerance. Every stream is built from an LFSR with
//                 an EXPLICIT seed; generate_valid_stream() seeds itself off the clock and would
//                 make these flaky. Tolerances are set from a measured sweep, not guessed.

using namespace StochasticSimulator;

class SUnaryMultiplierTest : public ::testing::Test {
protected:
    static const std::size_t N = 1024;
    static const uint16_t SCALE = 1024;

    // Tolerance for a SuMUL result at N = 1024 on auto width. Measured max error over a 361-point
    // p0,p1 sweep was 0.016; this leaves roughly 2x headroom.
    static constexpr double TOL = 0.035;

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

    // Two streams off ONE LFSR, thresholded differently: maximally correlated (SCC = 1) but still
    // well mixed. This is the realistic failure mode -- operands sharing an RNG -- and the exact
    // case SuMUL was built for.
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

    // The old single-AND-gate multiplier, for side-by-side comparison.
    double and_gate(const std::vector<bool>& a, const std::vector<bool>& b) const {
        Multiplier m;
        std::vector<bool> o;
        o.reserve(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) {
            o.push_back(m.multiply(a[i], b[i]));
        }
        return calculate_probability(o);
    }
};

// =====================================================================================
// STRUCTURAL -- the circuit itself
// =====================================================================================

// The enable [en in Figure 3(a)] is in_0. G's RNG must step exactly once per in_0 one, and
// never on an in_0 zero. This is the single property that makes SuMUL correlation-immune, so if
// only one test survives, keep this one.
TEST_F(SUnaryMultiplierTest, EnableFiresExactlyOnInput0Ones) {
    std::vector<bool> a = make_stream(0.30, 7);
    std::vector<bool> b = make_stream(0.70, 113);

    SUnaryMultiplier sumul(6);
    for (std::size_t i = 0; i < N; ++i) {
        sumul.multiply(a[i], b[i]);
    }

    EXPECT_EQ(sumul.get_enabled_cycles(), count_ones(a));
    EXPECT_LT(sumul.get_enabled_cycles(), N);  // in_0 is not all ones, so some cycles must have held
}

// (2) ANDMUL: the output passes through an AND gate with in_0, so it can never carry a 1 on a
// cycle where in_0 was 0, no matter what the regenerated in_1 bit was.
TEST_F(SUnaryMultiplierTest, OutputIsAlwaysSubsetOfInput0) {
    std::vector<bool> a = make_stream(0.40, 21);
    std::vector<bool> b = make_stream(0.90, 305);
    std::vector<bool> out = sumul_stream(a, b);

    ASSERT_EQ(out.size(), a.size());
    for (std::size_t i = 0; i < N; ++i) {
        if (!a[i]) {
            EXPECT_FALSE(out[i]) << "output was 1 at cycle " << i << " where in_0 was 0";
        }
    }
}

// x * 0 is exactly 0 when the zero is on in_0 -- the AND gate forces it, no statistics involved.
// Contrast with EdgeProbabilities below, where the zero sits on in_1 and warm-up leaks a few ones.
TEST_F(SUnaryMultiplierTest, ZeroOnInput0GivesExactlyZero) {
    std::vector<bool> zeros(N, false);
    std::vector<bool> b = make_stream(0.75, 9);

    EXPECT_DOUBLE_EQ(sumul_probability(zeros, b), 0.0);
    EXPECT_EQ(count_ones(sumul_stream(zeros, b)), 0u);
}

// The unit is stateful, so reset() must genuinely return it to power-on -- otherwise a reused
// instance silently reports different numbers for the same inputs.
TEST_F(SUnaryMultiplierTest, DeterministicAcrossReset) {
    std::vector<bool> a = make_stream(0.35, 33);
    std::vector<bool> b = make_stream(0.65, 404);

    SUnaryMultiplier sumul(6);
    std::vector<bool> first, second;
    for (std::size_t i = 0; i < N; ++i) first.push_back(sumul.multiply(a[i], b[i]));

    sumul.reset();
    for (std::size_t i = 0; i < N; ++i) second.push_back(sumul.multiply(a[i], b[i]));

    EXPECT_EQ(first, second);
    EXPECT_EQ(sumul.get_enabled_cycles(), count_ones(a));  // not double-counted across the reset
}

// The counter window boots at 0,1,0,1,... -- half ones, i.e. the unit assumes p1 = 0.5 until it
// has seen real data. This is the source of the warm-up error, so pin it down.
TEST_F(SUnaryMultiplierTest, CounterStartsHalfFull) {
    SUnaryMultiplier sumul(8);
    EXPECT_EQ(sumul.get_entry(), 256u);
    EXPECT_EQ(sumul.get_counter_value(), 128u);
    EXPECT_DOUBLE_EQ(sumul.get_tracked_probability(), 0.5);
    EXPECT_EQ(sumul.get_enabled_cycles(), 0u);
}

// C carries a BINARY count, not a unary bit -- the thick line in Figure 3(a). Constant inputs
// drive it to its exact rails, which is an exact assertion rather than a statistical one.
TEST_F(SUnaryMultiplierTest, CounterIsBinaryAndReachesBothRails) {
    const uint64_t ENTRY = 64;  // width 6
    std::vector<bool> ones(N, true);
    std::vector<bool> zeros(N, false);

    SUnaryMultiplier high(6);
    for (std::size_t i = 0; i < ENTRY; ++i) high.multiply(true, true);
    EXPECT_EQ(high.get_counter_value(), ENTRY);          // window saturated with ones
    EXPECT_DOUBLE_EQ(high.get_tracked_probability(), 1.0);

    SUnaryMultiplier low(6);
    for (std::size_t i = 0; i < ENTRY; ++i) low.multiply(true, false);
    EXPECT_EQ(low.get_counter_value(), 0u);              // window flushed to zeros
    EXPECT_DOUBLE_EQ(low.get_tracked_probability(), 0.0);
}

// C shifts on EVERY cycle, enable or not -- it is tracking in_1's value, not producing output.
// If someone ever gates the counter along with the RNG, this catches it.
TEST_F(SUnaryMultiplierTest, CounterShiftsEvenWhileEnableIsLow) {
    SUnaryMultiplier sumul(6);
    for (std::size_t i = 0; i < 64; ++i) {
        sumul.multiply(false, true);  // in_0 held LOW the whole time
    }
    EXPECT_EQ(sumul.get_enabled_cycles(), 0u);   // RNG never stepped
    EXPECT_EQ(sumul.get_counter_value(), 64u);   // ...but the counter still filled
}

// in_0 is the enable and passes through untouched; in_1 gets regenerated. Both orderings compute
// p0 * p1, but they are not the same circuit and must not produce the same bits.
TEST_F(SUnaryMultiplierTest, OperandsAreNotInterchangeable) {
    std::vector<bool> a = make_stream(0.25, 55);
    std::vector<bool> b = make_stream(0.80, 601);

    std::vector<bool> forward = sumul_stream(a, b);
    std::vector<bool> reverse = sumul_stream(b, a);

    EXPECT_NE(forward, reverse) << "swapping the operands should change the circuit";

    double expect = calculate_probability(a) * calculate_probability(b);
    EXPECT_NEAR(calculate_probability(forward), expect, TOL);
    EXPECT_NEAR(calculate_probability(reverse), expect, TOL);
}

// =====================================================================================
// SIZING
// =====================================================================================

// width = 1 + floor(log2(N)/2), i.e. a window near 2 * sqrt(N), clamped to MAX_WIDTH.
TEST_F(SUnaryMultiplierTest, WidthForStreamLengthFollowsTheSqrtRule) {
    EXPECT_EQ(SUnaryMultiplier::width_for_stream_length(256), 5u);
    EXPECT_EQ(SUnaryMultiplier::width_for_stream_length(512), 5u);
    EXPECT_EQ(SUnaryMultiplier::width_for_stream_length(1024), 6u);
    EXPECT_EQ(SUnaryMultiplier::width_for_stream_length(4096), 7u);
    EXPECT_EQ(SUnaryMultiplier::width_for_stream_length(16384), 8u);

    // Never exceeds the flip-flop budget, however long the stream.
    EXPECT_LE(SUnaryMultiplier::width_for_stream_length(1ull << 40), SUnaryMultiplier::MAX_WIDTH);

    SUnaryMultiplier sized = SUnaryMultiplier::for_stream_length(1024);
    EXPECT_EQ(sized.get_width(), 6u);
    EXPECT_EQ(sized.get_entry(), 64u);
}

// AUTO_WIDTH must be exactly equivalent to passing the width the rule would have picked.
TEST_F(SUnaryMultiplierTest, AutoWidthMatchesTheExplicitWidth) {
    std::vector<bool> a = make_stream(0.45, 12);
    std::vector<bool> b = make_stream(0.55, 808);

    unsigned expected_width = SUnaryMultiplier::width_for_stream_length(N);
    EXPECT_EQ(sumul_stream(a, b), sumul_stream(a, b, expected_width));
}

TEST_F(SUnaryMultiplierTest, StreamLengthConstructorMatchesSobol) {
    SUnaryMultiplier sumul(StreamLength::Length_1024, 1);
    EXPECT_EQ(sumul.get_width(), 10u);
    EXPECT_EQ(sumul.get_entry(), 1024u);
}

// =====================================================================================
// ACCURACY
// =====================================================================================

TEST_F(SUnaryMultiplierTest, MultipliesIndependentStreams) {
    const double probs[] = {0.10, 0.25, 0.50, 0.75, 0.90};
    for (double p0 : probs) {
        for (double p1 : probs) {
            std::vector<bool> a = make_stream(p0, 7);
            std::vector<bool> b = make_stream(p1, 113);
            double expect = calculate_probability(a) * calculate_probability(b);
            EXPECT_NEAR(sumul_probability(a, b), expect, TOL)
                << "p0=" << p0 << " p1=" << p1;
        }
    }
}

// The headline claim of the paper. Both operands come off the same LFSR, so they are perfectly
// correlated -- and SuMUL still multiplies them.
TEST_F(SUnaryMultiplierTest, StaysAccurateOnMaximallyCorrelatedStreams) {
    const double probs[] = {0.20, 0.40, 0.60, 0.80};
    for (double p0 : probs) {
        for (double p1 : probs) {
            std::vector<bool> a, b;
            make_correlated_pair(p0, p1, a, b);
            double expect = calculate_probability(a) * calculate_probability(b);
            EXPECT_NEAR(sumul_probability(a, b), expect, TOL)
                << "p0=" << p0 << " p1=" << p1;
        }
    }
}

// Characterises the OLD unit so the regression is visible: a plain AND gate fed correlated
// streams does not multiply at all, it collapses to min(p0, p1). This test failing would mean
// the correlation problem went away and SuMUL is no longer needed -- worth knowing either way.
TEST_F(SUnaryMultiplierTest, PlainAndGateCollapsesToMinWhenCorrelated) {
    // 0.5 and 0.6 sit where the damage is near its worst: min is 0.50 against a true product of
    // 0.30, so the correlation error is ~0.20. Note the gap is min(p0,p1) - p0*p1, which peaks at
    // p0 = p1 = 0.5 and shrinks toward the corners -- pick the pair deliberately, not at random.
    std::vector<bool> a, b;
    make_correlated_pair(0.50, 0.60, a, b);

    double pa = calculate_probability(a);
    double pb = calculate_probability(b);
    double product = pa * pb;

    // Thresholding one RNG twice nests the streams, so the AND gate just hands back the smaller.
    EXPECT_NEAR(and_gate(a, b), std::fmin(pa, pb), 0.02);

    double and_err = std::fabs(and_gate(a, b) - product);
    double sumul_err = std::fabs(sumul_probability(a, b) - product);

    EXPECT_GT(and_err, 0.15);                              // the AND gate is badly wrong
    EXPECT_NEAR(sumul_probability(a, b), product, TOL);     // SuMUL is not
    EXPECT_GT(and_err, sumul_err * 5.0);                    // and by a wide margin
}

TEST_F(SUnaryMultiplierTest, BeatsThePlainAndGateOnCorrelatedStreams) {
    const double probs[] = {0.20, 0.40, 0.60, 0.80};
    for (double p0 : probs) {
        for (double p1 : probs) {
            std::vector<bool> a, b;
            make_correlated_pair(p0, p1, a, b);
            double expect = calculate_probability(a) * calculate_probability(b);

            double sumul_err = std::fabs(sumul_probability(a, b) - expect);
            double and_err = std::fabs(and_gate(a, b) - expect);

            EXPECT_LE(sumul_err, and_err) << "p0=" << p0 << " p1=" << p1;
        }
    }
}

// Warm-up is real and lives at the rails: the window boots believing p1 = 0.5, so a few wrong
// bits escape before it fills. These bounds document the cost rather than pretending it is zero.
TEST_F(SUnaryMultiplierTest, EdgeProbabilitiesPayTheWarmUpCost) {
    std::vector<bool> ones(N, true);
    std::vector<bool> zeros(N, false);

    // Not exact -- but must stay inside the warm-up budget.
    EXPECT_NEAR(sumul_probability(ones, ones), 1.0, 0.03);
    EXPECT_NEAR(sumul_probability(ones, zeros), 0.0, 0.03);
    EXPECT_LT(sumul_probability(ones, ones), 1.0);   // the leak is real, in this direction
    EXPECT_GT(sumul_probability(ones, zeros), 0.0);

    // 1 * p passes p straight through.
    std::vector<bool> b = make_stream(0.50, 77);
    EXPECT_NEAR(sumul_probability(ones, b), calculate_probability(b), TOL);
    EXPECT_NEAR(sumul_probability(b, ones), calculate_probability(b), TOL);
}

// =====================================================================================
// ERROR HANDLING
// =====================================================================================

TEST_F(SUnaryMultiplierTest, RejectsMalformedStreams) {
    std::vector<bool> two = {true, false};
    std::vector<bool> three = {true, false, true};
    std::vector<bool> empty;

    EXPECT_THROW(sumul_stream(two, three), std::invalid_argument);
    EXPECT_THROW(sumul_stream(empty, empty), std::invalid_argument);
    EXPECT_THROW(sumul_probability(two, three), std::invalid_argument);
}

TEST_F(SUnaryMultiplierTest, RejectsOutOfRangeWidth) {
    EXPECT_THROW(SUnaryMultiplier(0), std::invalid_argument);
    EXPECT_THROW(SUnaryMultiplier(SUnaryMultiplier::MAX_WIDTH + 1), std::invalid_argument);
    EXPECT_THROW(SUnaryMultiplier::width_for_stream_length(0), std::invalid_argument);

    // The (void) casts keep these from parsing as variable declarations.
    EXPECT_NO_THROW((void)SUnaryMultiplier(1));
    EXPECT_NO_THROW((void)SUnaryMultiplier(SUnaryMultiplier::MAX_WIDTH));
}
