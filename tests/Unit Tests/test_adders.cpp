#include <gtest/gtest.h>
#include "../../modules/adder.hpp"
#include "../../bsg/lfsr.hpp"
#include "../../bsg/sng.hpp"

using namespace StochasticSimulator;

// TEST 1: Hardware Mux Logic 
// bool add_scaled_or_weighted(bool lhs, bool rhs, bool select_bit) const;
TEST(MuxScaledWeightedAdderTest, StrictMuxLogic){
    Adder adder;


    // Test case 1; select_bit is false (0)
    // Mux should strictly rute the Right hand side (rhs) value. 
    EXPECT_TRUE(adder.add_scaled_or_weighted(true, true, false)); // Should return true (rhs value 1)
    EXPECT_FALSE(adder.add_scaled_or_weighted(true, false, false)); // Should return false (rhs value 0)
    EXPECT_TRUE(adder.add_scaled_or_weighted(false, true, false)); // Should return true (rhs value 1)
    EXPECT_FALSE(adder.add_scaled_or_weighted(false, false, false)); // Should return false (rhs value 0)

    // Test case 2; select_bit is true (1)
    // Mux should strictly route the Left hand side (lhs) value.
    EXPECT_TRUE(adder.add_scaled_or_weighted(true, true, true)); // Should return true (lhs value 1)
    EXPECT_TRUE(adder.add_scaled_or_weighted(true, false, true)); // Should return true (lhs value 1)
    EXPECT_FALSE(adder.add_scaled_or_weighted(false, true, true)); // Should return false (lhs value 0)
    EXPECT_FALSE(adder.add_scaled_or_weighted(false, false, true)); // Should return false (lhs value 0)

}