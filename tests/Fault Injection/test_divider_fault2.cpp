/*
 * =========================================================================================
 * test_divider_fault2.cpp
 * =========================================================================================
 * EXHAUSTIVE 32-BIT STOCHASTIC DIVIDER (UP/DOWN COUNTER) FAULT-INJECTION SIMULATOR
 *
 * Mathematical & Architectural Overview:
 * 1. Unary representation in Stochastic Computing:
 *    A value in the interval [0, 1) is represented as a bitstream of length N.
 *    For a 32-bit stream, there are 32 possible unique fractional numbers represented by
 *    ones counts ranging from 0/32 to 31/32 inclusive.
 *
 * 2. Stochastic Division via Up/Down Counter:
 *    The stochastic divider uses an up/down counter feedback loop (Gaines division).
 *    It computes Z = X / Y in the stochastic domain. For accurate unipolar division,
 *    CountX should be <= CountY so that the result doesn't saturate at 1.0.
 *    The divider is STATEFUL and processes entire streams at once, returning a
 *    decoded output probability (double).
 *
 * 3. Exhaustive Input Coverage:
 *    We sweep numerator CountX in [0, 31] and denominator CountY in [1, 31].
 *    CountY starts at 1 to avoid division by zero. We only test pairs where
 *    CountX <= CountY to stay within the valid unipolar division range.
 *    This covers all valid (X, Y) pairs within the 32-bit representation.
 *
 * 4. Fault Injection Sweep:
 *    For each organization, we inject fault intensities from 1 to 32 bits flipped.
 *    Only stream X (numerator) is faulted; stream Y (denominator) remains clean.
 *
 * 5. Single-Open File Output (High Performance):
 *    The CSV file is opened ONCE at the start, all trial records are streamed
 *    directly to the file buffer, and it is closed ONCE at the end.
 *
 * 6. Error Components:
 *    - Ideal Quotient (real): CountX / CountY (as a ratio of ones counts)
 *    - Clean Output Probability: result of ud_counter_division on clean streams
 *    - Faulted Output Probability: result of ud_counter_division on faulted streams
 *    - Precision Error: CleanOutputProb - IdealQuotient
 *    - Bit Flip Error: FaultedOutputProb - CleanOutputProb
 *    - Total Error: FaultedOutputProb - IdealQuotient [Total = Precision + Bit Flip]
 * =========================================================================================
 */

#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>

#include "../../modules/divider.hpp"

using namespace StochasticSimulator;

namespace {

// Definition of constants for the 32-bit simulator sweep.
constexpr size_t STREAM_LEN = 32;            // Length of the stochastic bitstreams
constexpr size_t ORGANIZATIONS_PER_PAIR = 5; // Unique bit permutations per input value pair

// Builds a stream with exactly 'count' ones, shuffled randomly.
std::vector<bool> build_stream(size_t count, std::mt19937& rng) {
    // Step 1: Initialize a stream of all false (zeros)
    std::vector<bool> stream(STREAM_LEN, false);
    // Step 2: Set the first 'count' positions to true (ones)
    for (size_t i = 0; i < count; ++i) {
        stream[i] = true;
    }
    // Step 3: Shuffle to distribute the ones randomly across the stream
    std::shuffle(stream.begin(), stream.end(), rng);
    return stream;
}

// Converts a boolean bitstream vector to a CSV-friendly string: "[1, 0, 1, ...]"
std::string vector_to_csv_string(const std::vector<bool>& vec) {
    std::string s = "\"[";
    for (size_t i = 0; i < vec.size(); ++i) {
        s += (vec[i] ? "1" : "0");
        if (i < vec.size() - 1) s += ", ";
    }
    s += "]\"";
    return s;
}

// Computes the number of ones (true values) inside a boolean stream
size_t count_ones_local(const std::vector<bool>& vec) {
    return static_cast<size_t>(std::count(vec.begin(), vec.end(), true));
}

// Flips exactly num_flips distinct random bit positions in-place
void flip_random_bits(std::vector<bool>& stream, size_t num_flips, std::mt19937& rng) {
    // If flips request exceeds stream size, simply invert the entire stream
    if (num_flips >= stream.size()) {
        stream.flip();
        return;
    }
    // Create bit indices list
    std::vector<size_t> indices(stream.size());
    std::iota(indices.begin(), indices.end(), 0);
    // Shuffle the index list to select which bit positions to flip
    std::shuffle(indices.begin(), indices.end(), rng);
    // Flip the bits at the first num_flips indices
    for (size_t i = 0; i < num_flips; ++i) {
        stream[indices[i]] = !stream[indices[i]];
    }
}

} // namespace

// Setup Google Test fixture class for the exhaustive divider fault injection sweep
class DividerExhaustiveInjectionTest : public ::testing::Test {
protected:
    std::mt19937 rng;      // Pseudo-random number generator engine

    void SetUp() override {
        rng.seed(1337);    // Fixed seed to guarantee identical generation across CI/CD environments
    }
};

// Google Test targeting exhaustive divider fault analysis
TEST_F(DividerExhaustiveInjectionTest, Sweep32BitExhaustiveDividerFaultInjection) {
    // Step 1: Define output path
    const std::string filename = "32bit_exhaustive_divider_trials.csv";

    // Step 2: Open the CSV file ONCE at the start of the test.
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to open output CSV path: " << filename;

    // Step 3: Write the CSV column header row
    csv << "Trial,CountX,CountY,CleanOutputProb,Organization,BitsFlipped,"
           "Vector_X_Original,Vector_Y_Original,Vector_X_Faulted,"
           "FaultedOutputProb,TotalError,PrecisionError,BitFlipError\n";

    size_t trial = 0;             // Running tally of total trials simulated
    size_t qualifying_pairs = 0;  // Count of X/Y ones count combinations evaluated

    // Step 4: Loop through all possible ones counts for numerator Stream X (0 to 31)
    for (size_t countX = 0; countX < STREAM_LEN; ++countX) {

        // Step 5: Loop through denominator Stream Y (1 to 31, skip 0 to avoid division by zero)
        // Only test where countX <= countY (valid unipolar division constraint)
        for (size_t countY = 1; countY < STREAM_LEN; ++countY) {

            // Step 6: Skip invalid division cases where numerator > denominator
            if (countX > countY) continue;

            ++qualifying_pairs;

            // Step 7: Compute the ideal quotient (unquantized mathematical result)
            double ideal_quotient = static_cast<double>(countX) / static_cast<double>(countY);

            // Step 8: Sweep multiple spatial organizations for this combination
            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {

                // Step 9: Build clean numerator and denominator streams
                std::vector<bool> vecX_orig = build_stream(countX, rng);
                std::vector<bool> vecY_orig = build_stream(countY, rng);

                // Step 10: Compute clean output probability by running divider without faults
                double clean_output_prob = ud_counter_division(vecX_orig, vecY_orig);

                // Step 11: Verify stream generation counts (safety assertions)
                ASSERT_EQ(count_ones_local(vecX_orig), countX)
                    << "Stream X generation mismatch for CountX=" << countX;
                ASSERT_EQ(count_ones_local(vecY_orig), countY)
                    << "Stream Y generation mismatch for CountY=" << countY;

                // Step 12: Sweep all possible fault intensities (1 to 32 bit flips)
                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {

                    // Copy original Stream X to inject faults
                    std::vector<bool> vecX_faulted = vecX_orig;

                    // Step 13: Flip the requested number of bits in stream X
                    flip_random_bits(vecX_faulted, flips, rng);

                    // Step 14: Run the divider simulation on faulted streams
                    double faulted_output_prob = ud_counter_division(vecX_faulted, vecY_orig);

                    // Step 15: Compute different error metrics (all in probability domain)
                    // Total Error: difference between faulted output probability and ideal quotient
                    double total_error = faulted_output_prob - ideal_quotient;

                    // Precision Error: difference between clean output probability and ideal quotient
                    double precision_error = clean_output_prob - ideal_quotient;

                    // Bit Flip Error: difference between faulted output and clean output
                    double bit_flip_error = faulted_output_prob - clean_output_prob;

                    // Step 16: Log results directly to the open CSV stream buffer
                    ++trial;
                    csv << trial << "," << countX << "," << countY << "," << clean_output_prob << ","
                        << org << "," << flips << ","
                        << vector_to_csv_string(vecX_orig) << ","
                        << vector_to_csv_string(vecY_orig) << ","
                        << vector_to_csv_string(vecX_faulted) << ","
                        << faulted_output_prob << ","
                        << total_error << ","
                        << precision_error << ","
                        << bit_flip_error << "\n";
                }
            }
        }
    }

    // Step 17: Close the CSV file handle ONCE at the end of the test.
    csv.close();

    // Print summary stats to the console
    std::cout << "[INFO] Evaluated " << qualifying_pairs << " valid (X <= Y) input count pairs." << std::endl;
    std::cout << "[SUCCESS] Exported " << trial << " divider fault-injection trials to " << filename << std::endl;
}
