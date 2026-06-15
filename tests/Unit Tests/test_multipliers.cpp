#include <gtest/gtest.h>
#include "../../modules/multiplier.hpp"
#include "../../bsg/lfsr.hpp"
#include "../../bsg/sng.hpp"

using namespace StochasticSimulator;


// TEST 1: Hardware Truth Table Logic
// Checks if the stateless AND gate functions perfectly bit-by-bit.
TEST(MultiplierTest, StrictAndGateLogic) {
    Multiplier mul;

    // 0 x 0 = 0
    EXPECT_FALSE(mul.multiply(false, false));
    // 0 x 1 = 0
    EXPECT_FALSE(mul.multiply(false, true));
    // 1 x 0 = 0
    EXPECT_FALSE(mul.multiply(true, false));
    // 1 x 1 = 1
    EXPECT_TRUE(mul.multiply(true, true));
}