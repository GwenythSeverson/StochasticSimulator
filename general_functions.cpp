#pragma once

#include "general_functions.hpp"
#include "bsg/sng.hpp"
#include <cmath>
#include <algorithm>
#include <random>
#include "bsg/lfsr.hpp"
#include <string>
#include <sstream>

namespace StochasticSimulator {


bool can_terminate_early(double confidence) {
    return false;
}



double calculate_zce_window(const std::vector<int>& streamA, const std::vector<int>& streamB, std::size_t sample_length) {
    return false;
}

double calculate_probability(const std::vector<bool>& stream) {
    if (stream.empty()) return 0.0;
    size_t ones = std::count(stream.begin(), stream.end(), true);
    return static_cast<double>(ones) / stream.size();
}

std::vector<bool> generate_valid_stream(double target_prob, StochasticSimulator::StreamLength length_mode, uint16_t max_cycles) {
    BitstreamGenerator bsg;
    const double ALLOWED_ERROR = 0.0005; // 0.05% error margin bounds

    while (true) {
        // Seed with 0 for timed hardware random-seeding path
        FlexibleLFSR lfsr(length_mode, 0);
        std::vector<bool> candidate_stream;
        candidate_stream.reserve(max_cycles);

        // Step through the precise sequence length
        for (uint16_t i = 0; i < max_cycles; ++i) {
            uint16_t rand_val = lfsr.next();
            bool bit = bsg.generate_bit(target_prob, rand_val, max_cycles);
            candidate_stream.push_back(bit);
        }

        double actual_prob = calculate_probability(candidate_stream);
        double error = std::abs(target_prob - actual_prob);

        // Keep the stream only if it lands safely within the 0.05% target accuracy
        if (error <= ALLOWED_ERROR) {
            return candidate_stream;
        }
    }
}

std::string format_bit_vector(const std::vector<bool>& stream, size_t max_chars) {
    std::stringstream ss;
    size_t display_len = std::min(stream.size(), max_chars);
    for (size_t i = 0; i < display_len; ++i) {
        ss << (stream[i] ? "1" : "0");
    }
    if (stream.size() > max_chars) {
        ss << "...";
    }
    return ss.str();
}

size_t count_ones(const std::vector<bool>& stream) {
    return std::count(stream.begin(), stream.end(), true);
}

}  // namespace stochastic