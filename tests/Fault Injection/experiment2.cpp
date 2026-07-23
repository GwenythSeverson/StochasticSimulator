#include <gtest/gtest.h>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>

namespace {

constexpr size_t STREAM_LEN = 32;
constexpr int COUNTER_MAX = 8; // Saturation boundary for the internal feedback counter loop

// Cycle-by-cycle hardware simulation of your Up/Down Counter Divider
double simulate_hardware_ud_divider(const std::vector<bool>& stream_X, const std::vector<bool>& stream_Y) {
    int counter = COUNTER_MAX / 2; // Initialize to mid-scale
    size_t output_ones = 0;

    for (size_t i = 0; i < STREAM_LEN; ++i) {
        // Output bit is determined by the upper half of the counter range
        bool z_out = (counter >= COUNTER_MAX / 2);
        if (z_out) {
            output_ones++;
        }

        // Feedback gating path: feedback stream ANDed with Denominator (Y)
        bool feedback_gate = z_out && stream_Y[i];

        // Up/Down Counter state update logic
        if (stream_X[i] && !feedback_gate) {
            if (counter < COUNTER_MAX) counter++; // Increment
        } else if (!stream_X[i] && feedback_gate) {
            if (counter > 0) counter--;           // Decrement
        }
    }

    return static_cast<double>(output_ones) / static_cast<double>(STREAM_LEN);
}

// Generates a stable, canonical sequence for an exact probability count 
std::vector<bool> generate_canonical_stream(size_t active_bits) {
    std::vector<bool> stream(STREAM_LEN, false);
    for (size_t i = 0; i < active_bits; ++i) {
        stream[i] = true;
    }
    return stream;
}

} // namespace stochasticSimulator

TEST(StochasticDividerExhaustiveTest, Sweep32BitExhaustiveDividerFault) {
    const std::string filename = "divider_exhaustive_data.csv";

    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to open output CSV path: " << filename;

    csv << "CountA,CountB,Flips,CleanOut,FaultedOut,Error\n";

    std::cout << "[INFO] Commencing exhaustive cycle-accurate hardware divider sweep..." << std::endl;

    // Sweep all 0-32 value levels for numerator and denominator
    for (size_t countA = 0; countA <= STREAM_LEN; ++countA) {
       // std::cout << "[PROGRESS] Simulating Numerator Level (CountA): " << countA << "/" << STREAM_LEN << std::endl;
        
        std::vector<bool> stream_X_clean = generate_canonical_stream(countA);

        for (size_t countB = 0; countB <= STREAM_LEN; ++countB) {
            std::vector<bool> stream_Y = generate_canonical_stream(countB);

            // Calculate base clean hardware operation output
            double clean_out = simulate_hardware_ud_divider(stream_X_clean, stream_Y);

            // Inject sequential deterministic bitwise corruption depths
            for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {
                
                std::vector<bool> stream_X_faulted = stream_X_clean;
                
                // Invert exactly 'flips' positions to test state tracking degradation
                for (size_t f = 0; f < flips; ++f) {
                    stream_X_faulted[f] = !stream_X_faulted[f];
                }

                // Route directly through your feedback division logic
                double faulted_out = simulate_hardware_ud_divider(stream_X_faulted, stream_Y);
                double error = faulted_out - clean_out;

                csv << countA << "," << countB << "," << flips << "," 
                    << clean_out << "," << faulted_out << "," << error << "\n";
            }
        }
    }

    csv.close();
    std::cout << "[SUCCESS] Exhaustive architecture behavior saved to " << filename << std::endl;
}