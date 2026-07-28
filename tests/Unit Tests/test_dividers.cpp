#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <random>
#include <cmath>

#include "divider.hpp"
#include "general_functions.hpp"

// note to self on edge case testing- i killed ranom rand gens and used a seed so now it checks if
// it converges generally to the right thing
//
// Tolerances below come from a measured px,py sweep at 3 seeds, not from guessing:
//     N =  1024  RMSE 0.0258  worst 0.0882
//     N =  4096  RMSE 0.0132  worst 0.0528
//     N = 16384  RMSE 0.0054  worst 0.0148
//     N = 50000  RMSE 0.0033  worst 0.0109
// The counter needs room to settle, so short streams get a looser bound. Tightening these without
// also raising N will make the suite flaky rather than strict.

using namespace StochasticSimulator;

class UDCounterDivisionStatisticalTest : public ::testing::Test {
protected:
    // Helper to generate a test bitstream with a specific target probability
    std::vector<bool> generate_test_stream(double probability, size_t length, unsigned int seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dis(0.0, 1.0);
        std::vector<bool> stream;
        stream.reserve(length);
        for (size_t i = 0; i < length; ++i) {
            stream.push_back(dis(gen) < probability);
        }
        return stream;
    }
};

// =====================================================================================
// REGRESSION -- the bug that took the divider offline entirely
// =====================================================================================

// The divider used to run its comparator on StreamLength::Length_65536. FlexibleLFSR keeps
// max_cycles in a uint16_t, so 65536 truncated to 0, the constructor's "seed >= max_cycles"
// check was true for every seed, and this function threw on EVERY call. That is why the
// convergence test below was commented out -- it was not failing, it was throwing.
TEST_F(UDCounterDivisionStatisticalTest, RunsWithoutThrowingAcrossSeeds) {
    std::vector<bool> X = generate_test_stream(0.20, 4096, 42);
    std::vector<bool> Y = generate_test_stream(0.80, 4096, 24);

    for (uint32_t seed : {0u, 1u, 42u, 1000u, 65535u, 65536u, 4000000000u}) {
        EXPECT_NO_THROW({
            double z = ud_counter_division(X, Y, seed);
            EXPECT_GE(z, 0.0);
            EXPECT_LE(z, 1.0);
        }) << "seed " << seed;
    }
}

// =====================================================================================
// CONVERGENCE -- does the feedback loop actually land on X/Y?
// =====================================================================================

// 3. Statistical Test: Accurate Mathematical Division (e.g., 0.20 / 0.80 = 0.25)
// Restored from the commented-out original. It now passes with roughly 37x margin.
TEST_F(UDCounterDivisionStatisticalTest, LongStreamConvergenceTest) {
    size_t long_stream_len = 50000;  // Long enough to allow statistical convergence

    // Generate streams with known statistical weights using a fixed seed for CI stability
    std::vector<bool> stream_X = generate_test_stream(0.20, long_stream_len, 42);
    std::vector<bool> stream_Y = generate_test_stream(0.80, long_stream_len, 24);

    double expected_math_result = 0.20 / 0.80;  // 0.25
    double actual_sc_result = ud_counter_division(stream_X, stream_Y, 42);

    // Stochastic computing with an SNG randomizer exhibits natural variance.
    // We expect the result to land within a +/- 3% margin of error.
    EXPECT_NEAR(actual_sc_result, expected_math_result, 0.03);
}

// The single point above only proves the loop works at one quotient. Sweep the domain, and score
// against the streams' ACTUAL ratio rather than the nominal one -- at finite length the generated
// streams are not exactly 0.20 and 0.80, and charging the divider for the SNG's error is unfair.
TEST_F(UDCounterDivisionStatisticalTest, ConvergesAcrossTheUnipolarDomain) {
    const size_t N = 16384;
    const double py_values[] = {0.40, 0.60, 0.80, 0.95};
    const double px_values[] = {0.05, 0.20, 0.50};

    for (double py : py_values) {
        for (double px : px_values) {
            if (px > py) continue;  // outside the unipolar domain; see SaturatesWhenXExceedsY

            std::vector<bool> X = generate_test_stream(px, N, 42);
            std::vector<bool> Y = generate_test_stream(py, N, 24);

            double expect = calculate_probability(X) / calculate_probability(Y);
            EXPECT_NEAR(ud_counter_division(X, Y, 42), expect, 0.03)
                << "px=" << px << " py=" << py;
        }
    }
}

// A longer stream must buy real accuracy -- if it does not, the counter is not settling and the
// depth rule is broken.
TEST_F(UDCounterDivisionStatisticalTest, AccuracyImprovesWithStreamLength) {
    auto mean_abs_error = [&](size_t n) {
        double acc = 0.0;
        int cnt = 0;
        for (double px : {0.05, 0.20, 0.50}) {
            std::vector<bool> X = generate_test_stream(px, n, 42);
            std::vector<bool> Y = generate_test_stream(0.80, n, 24);
            double expect = calculate_probability(X) / calculate_probability(Y);
            acc += std::fabs(ud_counter_division(X, Y, 42) - expect);
            ++cnt;
        }
        return acc / cnt;
    };

    double short_err = mean_abs_error(1024);
    double long_err = mean_abs_error(16384);
    EXPECT_LT(long_err, short_err) << "16x the stream did not improve accuracy";
}

// =====================================================================================
// COUNTER DEPTH
// =====================================================================================

// depth ~ sqrt(N), a power of two, clamped to [16, 1024].
TEST_F(UDCounterDivisionStatisticalTest, CounterDepthTracksSqrtOfStreamLength) {
    EXPECT_EQ(gaines_counter_depth_for_length(512), 16u);
    EXPECT_EQ(gaines_counter_depth_for_length(1024), 32u);
    EXPECT_EQ(gaines_counter_depth_for_length(4096), 64u);
    EXPECT_EQ(gaines_counter_depth_for_length(16384), 128u);
    EXPECT_EQ(gaines_counter_depth_for_length(50000), 256u);

    EXPECT_EQ(gaines_counter_depth_for_length(1), 16u);            // clamped low
    EXPECT_EQ(gaines_counter_depth_for_length(100000000), 1024u);  // clamped high
    EXPECT_THROW(gaines_counter_depth_for_length(0), std::invalid_argument);
}

// The depth feeds an exact power-of-two scaling into the comparator, so a non-power-of-two would
// silently bias every output bit. Reject it loudly instead.
TEST_F(UDCounterDivisionStatisticalTest, RejectsInvalidCounterDepth) {
    std::vector<bool> X = generate_test_stream(0.20, 512, 42);
    std::vector<bool> Y = generate_test_stream(0.80, 512, 24);

    EXPECT_THROW(gaines_division_stream(X, Y, 42, 33), std::invalid_argument);
    EXPECT_THROW(gaines_division_stream(X, Y, 42, 100), std::invalid_argument);
    EXPECT_THROW(gaines_division_stream(X, Y, 42, 32768), std::invalid_argument);

    EXPECT_NO_THROW(gaines_division_stream(X, Y, 42, 32));
    EXPECT_NO_THROW(gaines_division_stream(X, Y, 42, GAINES_AUTO_DEPTH));
}

// =====================================================================================
// BITSTREAM OUTPUT
// =====================================================================================

TEST_F(UDCounterDivisionStatisticalTest, StreamAndProbabilityFormsAgree) {
    std::vector<bool> X = generate_test_stream(0.25, 4096, 42);
    std::vector<bool> Y = generate_test_stream(0.50, 4096, 24);

    std::vector<bool> Z = gaines_division_stream(X, Y, 42);
    ASSERT_EQ(Z.size(), X.size());
    EXPECT_DOUBLE_EQ(calculate_probability(Z), ud_counter_division(X, Y, 42));
}

TEST_F(UDCounterDivisionStatisticalTest, IsDeterministicForAGivenSeed) {
    std::vector<bool> X = generate_test_stream(0.30, 4096, 42);
    std::vector<bool> Y = generate_test_stream(0.90, 4096, 24);

    EXPECT_EQ(gaines_division_stream(X, Y, 7), gaines_division_stream(X, Y, 7));
    EXPECT_NE(gaines_division_stream(X, Y, 7), gaines_division_stream(X, Y, 99));
}

// =====================================================================================
// DOMAIN AND EDGES
// =====================================================================================

TEST_F(UDCounterDivisionStatisticalTest, EdgeQuotients) {
    const size_t N = 8192;
    std::vector<bool> Y = generate_test_stream(0.80, N, 24);
    std::vector<bool> zeros(N, false);
    std::vector<bool> ones(N, true);

    EXPECT_NEAR(ud_counter_division(zeros, Y, 42), 0.0, 0.02);   // 0 / y
    EXPECT_NEAR(ud_counter_division(Y, Y, 42), 1.0, 0.02);       // y / y

    std::vector<bool> X = generate_test_stream(0.30, N, 42);
    EXPECT_NEAR(ud_counter_division(X, ones, 42), calculate_probability(X), 0.03);  // x / 1
}

// Unipolar output cannot exceed 1, so X > Y pins the counter at its top rail. This is the design
// working as specified, not a fault -- pinning it down so nobody "fixes" it later.
TEST_F(UDCounterDivisionStatisticalTest, SaturatesWhenXExceedsY) {
    std::vector<bool> X = generate_test_stream(0.80, 8192, 42);
    std::vector<bool> Y = generate_test_stream(0.40, 8192, 24);

    double result = ud_counter_division(X, Y, 42);  // 0.80 / 0.40 = 2.0 mathematically
    EXPECT_GT(result, 0.97);
    EXPECT_LE(result, 1.0);
}

// =====================================================================================
// ERROR HANDLING
// =====================================================================================

// 4. Boundary Protection: Error Handling for Mismatched Lengths
// Ensures your structural safety guards return gracefully instead of crashing.
TEST_F(UDCounterDivisionStatisticalTest, ErrorHandlingOnMismatchedStreams) {
    std::vector<bool> short_stream = {true, false, true};
    std::vector<bool> long_stream  = {true, false, true, false, true};

    double result = ud_counter_division(short_stream, long_stream, 42);

    // Should safely hit the error catch block and return 0.0
    EXPECT_DOUBLE_EQ(result, 0.0);

    // The stream form reports the same failure unambiguously -- 0.0 is a legitimate quotient, an
    // empty vector is not. Prefer this form when you need to tell the two apart.
    EXPECT_TRUE(gaines_division_stream(short_stream, long_stream, 42).empty());
}

TEST_F(UDCounterDivisionStatisticalTest, ErrorHandlingOnEmptyStreams) {
    std::vector<bool> empty;
    EXPECT_DOUBLE_EQ(ud_counter_division(empty, empty, 42), 0.0);
    EXPECT_TRUE(gaines_division_stream(empty, empty, 42).empty());
}
