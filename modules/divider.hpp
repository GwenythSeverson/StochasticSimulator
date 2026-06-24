#pragma once

#include "general_functions.hpp" 
#include <cstdint>

namespace StochasticSimulator {

/**
 * @brief Simulates a Stochastic Divider using an Up/Down Counter feedback loop.
 * * Computes Z = X / Y in the stochastic domain. For accurate unipolar division,
 * target_x should be <= target_y so that the result doesn't saturate at 1.0.
 * * @param stream_X The numerator bitstream.
 * @param stream_Y The denominator bitstream.
 * @return double The decoded output probability of the division.
 */
double ud_counter_division(const std::vector<bool>& stream_X, const std::vector<bool>& stream_Y);
} // namespace StochasticSimulator