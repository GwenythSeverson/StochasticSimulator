#pragma once
#include <vector>
#include "bsg/lfsr.hpp"
#include <string>

namespace StochasticSimulator {

inline bool terminate_early_flag = false;



// Determines whether an extreme probability can short-circuit the stream.
bool can_terminate_early(double confidence);

// Calculates Julie Hsiao's Zero Correlation Error (ZCE) over a specific window length.
double calculate_zce_window(const std::vector<int>& streamA, const std::vector<int>& streamB, std::size_t sample_length);
// Forward declare the enum class so the compiler knows it exists right here
enum class StreamLength;

// Calculates actual probability from a bitstream
double calculate_probability(const std::vector<bool>& stream);

// Generates a valid bitstream satisfying the stringent 0.05% error bounds
std::vector<bool> generate_valid_stream(double target_prob, StreamLength length_mode, uint16_t max_cycles);

// Quick helper to stringify the vector of bits for the console display
std::string format_bit_vector(const std::vector<bool>& stream, size_t max_chars = 20);

// Quick helper to compute raw number of set bits (ones)
size_t count_ones(const std::vector<bool>& stream);

}  // namespace stochastic