#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <random>

// Forward declaration matching your HPP signature
namespace StochasticSimulator {
    double ud_counter_division(const std::vector<bool>& stream_X, const std::vector<bool>& stream_Y);
}

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

// 1. Test Extreme Saturated Bounds: Constant Increment (X=1, Y=0)
// Even with an RNG, if X is always 1 and Y is always 0, the counter pegs at MAX.
// Because counter == COUNTER_MAX, it is always greater than any RN drawn, forcing 1.0 output.
TEST_F(UDCounterDivisionStatisticalTest, SaturatedMaximumAlwaysOnes) {
    size_t stream_len = 1000;
    std::vector<bool> stream_X(stream_len, true);  // All 1s
    std::vector<bool> stream_Y(stream_len, false); // All 0s

    double result = ud_counter_division(stream_X, stream_Y);
    // Allows a tiny 1% margin for the counter to climb from 0 to 32- aka takes 31 cycles for it to climb to 32 and RNG pciks a number between 1 and 32 even on the first clocks so it cant technically be 1.0 all the way?
  EXPECT_NEAR(result, 1.0, 0.01); // STILL FAILS WHEN EXPECT NEAR AND NOT EXACT 1.0. counter logic?
}

// 2. Test Extreme Saturated Bounds: Constant Decrement (X=0, Y=1)
// If X is always 0, the counter pegs at 0. Since counter == 0, it can never be 
// strictly greater than an RN picked from [0..31]. Output must be 0.0.
TEST_F(UDCounterDivisionStatisticalTest, SaturatedMinimumAllZeros) {
    size_t stream_len = 1000;
    std::vector<bool> stream_X(stream_len, false); // All 0s
    std::vector<bool> stream_Y(stream_len, true);  // All 1s

    double result = ud_counter_division(stream_X, stream_Y);
    
    // Must be completely deterministic at the floor boundary
    EXPECT_DOUBLE_EQ(result, 0.0);
}

// 3. Statistical Test: Accurate Mathematical Division (e.g., 0.25 / 0.50 = 0.50)
// Instead of checking an exact fraction, we check if the RNG-driven loop 
// converges within an acceptable mathematical margin of error over a long stream.
TEST_F(UDCounterDivisionStatisticalTest, LongStreamConvergenceTest) {
    size_t long_stream_len = 50000; // Long enough to allow statistical convergence
    
    // Generate streams with known statistical weights using a fixed seed for CI stability
    std::vector<bool> stream_X = generate_test_stream(0.20, long_stream_len, 42);
    std::vector<bool> stream_Y = generate_test_stream(0.80, long_stream_len, 24);

    double expected_math_result = 0.20 / 0.80; // 0.25
    double actual_sc_result = ud_counter_division(stream_X, stream_Y);

    // Stochastic computing with an SNG randomizer exhibits natural variance.
    // We expect the result to land within a +/- 3% margin of error.
    EXPECT_NEAR(actual_sc_result, expected_math_result, 0.03);
}

// 4. Boundary Protection: Error Handling for Mismatched Lengths
// Ensures your structural safety guards return gracefully instead of crashing.
TEST_F(UDCounterDivisionStatisticalTest, ErrorHandlingOnMismatchedStreams) {
    std::vector<bool> short_stream = {true, false, true};
    std::vector<bool> long_stream  = {true, false, true, false, true};

    double result = ud_counter_division(short_stream, long_stream);
    
    // Should safely hit the error catch block and return 0.0
    EXPECT_DOUBLE_EQ(result, 0.0);
}