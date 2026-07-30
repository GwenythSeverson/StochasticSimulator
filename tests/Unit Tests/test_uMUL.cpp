#include <gtest/gtest.h>
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

    // p1 = 1 must reach 2^width -- the only value that beats every random number.
    umul.load_probability(1.0);
    EXPECT_EQ(umul.get_value(), umul.get_entry());
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
    EXPECT_EQ(umul.get_value(), umul.get_entry());
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
    EXPECT_NO_THROW(umul.load_value(1024));            // 2^width is legal
    EXPECT_THROW(umul.load_value(1025), std::invalid_argument);
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

// With in_1 loaded there is no estimator, so the rails are exact rather than merely close.
TEST_F(UnaryMultiplierTest, RailsAreExact) {
    std::vector<bool> ones = constant_stream(true);
    std::vector<bool> zeros = constant_stream(false);

    EXPECT_DOUBLE_EQ(umul_probability(ones, 1.0), 1.0);
    EXPECT_DOUBLE_EQ(umul_probability(ones, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(umul_probability(zeros, 1.0), 0.0);
    EXPECT_DOUBLE_EQ(umul_probability(zeros, 0.0), 0.0);

    // p1 = 1 passes in_0 through untouched, bit for bit.
    std::vector<bool> a = make_stream(0.375, 13);
    EXPECT_EQ(umul_stream(a, 1.0), a);
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
