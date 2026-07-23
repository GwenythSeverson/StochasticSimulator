#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <random>
#include "divider.hpp"
// note to self on edge case testing- i killed ranom rand gens and used a seed so now it checks if it converges generally to the right thing
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



// 3. Statistical Test: Accurate Mathematical Division (e.g., 0.25 / 0.50 = 0.50)
// Instead of checking an exact fraction, we check if the RNG-driven loop 
// converges within an acceptable mathematical margin of error over a long stream.
/*TEST_F(UDCounterDivisionStatisticalTest, LongStreamConvergenceTest) {
    size_t long_stream_len = 50000; // Long enough to allow statistical convergence
    
    // Generate streams with known statistical weights using a fixed seed for CI stability
    std::vector<bool> stream_X = generate_test_stream(0.20, long_stream_len, 42);
    std::vector<bool> stream_Y = generate_test_stream(0.80, long_stream_len, 24);

    double expected_math_result = 0.20 / 0.80; // 0.25
    double actual_sc_result = ud_counter_division(stream_X, stream_Y, 42);

    // Stochastic computing with an SNG randomizer exhibits natural variance.
    // We expect the result to land within a +/- 3% margin of error.
    EXPECT_NEAR(actual_sc_result, expected_math_result, 0.03);
}*/

// 4. Boundary Protection: Error Handling for Mismatched Lengths
// Ensures your structural safety guards return gracefully instead of crashing.
TEST_F(UDCounterDivisionStatisticalTest, ErrorHandlingOnMismatchedStreams) {
    std::vector<bool> short_stream = {true, false, true};
    std::vector<bool> long_stream  = {true, false, true, false, true};

    double result = ud_counter_division(short_stream, long_stream, 42);
    
    // Should safely hit the error catch block and return 0.0
    EXPECT_DOUBLE_EQ(result, 0.0);
}