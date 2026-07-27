/*
 * =========================================================================================
 * test_multiplier_fault2.cpp
 * =========================================================================================
 * EXHAUSTIVE 32-BIT STOCHASTIC MULTIPLIER FAULT-INJECTION SIMULATOR
 *
 * Mathematical & Architectural Overview:
 * 1. Unary representation in Stochastic Computing:
 *    A value in the interval [0, 1) is represented as a bitstream of length N.
 *    For a 32-bit stream, there are 32 possible unique fractional numbers represented by
 *    ones counts ranging from 0/32 to 31/32 inclusive.
 * 
 * 2. Exhaustive Input Coverage:
 *    To ensure "every possible generation possible" is simulated, we sweep all possible
 *    ones counts for both inputs: CountA in [0, 31] and CountB in [0, 31].
 *    This covers all 32 x 32 = 1024 input probability pairs.
 * 
 * 3. Spatial Organizations and Trials:
 *    For each input pair (CountA, CountB), we generate 5 unique spatial organizations
 *    (permutations) of the bit positions using a shuffled index map. This checks different
 *    bit alignments to test correlation sensitivity.
 * 
 * 4. Fault Injection Sweep:
 *    For each organization, we inject fault intensities from 1 to 32 bits flipped.
 *    Only stream A is faulted; stream B remains clean.
 * 
 * 5. Single-Open File Output (High Performance):
 *    Opening and closing a file for append in a loop creates severe filesystem overhead.
 *    To ensure a runtime of seconds rather than hours, we open the file stream ONCE at the
 *    very beginning of the test, write all 163,840 trial records directly to the file buffer,
 *    and close it ONCE at the end.
 *
 * 6. Error Components:
 *    - Ideal Product (real): (CountA / 32.0) * (CountB / 32.0)
 *    - Clean Output Ones Count (quantized): ExpectedCleanOnes = round(CountA * CountB / 32.0)
 *    - Faulted Output Ones Count: FaultedOnesCount
 *    - Precision Error: ExpectedCleanOnes - (CountA * CountB / 32.0)
 *    - Bit Flip Error: FaultedOnesCount - ExpectedCleanOnes
 *    - Total Error: FaultedOnesCount - (CountA * CountB / 32.0) [Total = Precision + Bit Flip]
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
#include <set>
#include <utility>
#include <iostream>

#include "../../modules/multiplier.hpp"

using namespace StochasticSimulator;

namespace {

// Definition of constants for the 32-bit simulator sweep.
constexpr size_t STREAM_LEN = 32;          // Length of the stochastic bitstreams
constexpr size_t ORGANIZATIONS_PER_PAIR = 5; // Unique bit permutations per input value pair

// Builds a pair of 32-bit streams with an EXACT or near-exact zero-correlation-error
// overlap for the given ones-counts, using a shared shuffled index map.
std::pair<std::vector<bool>, std::vector<bool>> build_zce_pair(
    size_t countA, size_t countB, size_t expected_overlap, std::mt19937& rng) {

    // Step 1: Initialize clean boolean arrays of length 32 with false (all zeros)
    std::vector<bool> vecA(STREAM_LEN, false);
    std::vector<bool> vecB(STREAM_LEN, false);

    // Step 2: Create a sequence of indices [0, 1, 2, ..., 31]
    std::vector<size_t> indices(STREAM_LEN);
    std::iota(indices.begin(), indices.end(), 0);

    // Step 3: Randomly shuffle indices to ensure bits are spatially distributed stochastically
    std::shuffle(indices.begin(), indices.end(), rng);

    // Step 4: Map overlapping ones, disjoint ones, and zeros onto the shuffled indices
    for (size_t i = 0; i < STREAM_LEN; ++i) {
        size_t idx = indices[i];
        if (i < expected_overlap) {
            // Overlapping region: both A and B are true (product ones)
            vecA[idx] = true;
            vecB[idx] = true;
        } else if (i < countA) {
            // A-only region: A is true, B is false
            vecA[idx] = true;
            vecB[idx] = false;
        } else if (i < countA + (countB - expected_overlap)) {
            // B-only region: A is false, B is true
            vecA[idx] = false;
            vecB[idx] = true;
        } else {
            // Zero-overlap region: both A and B are false
            vecA[idx] = false;
            vecB[idx] = false;
        }
    }
    return {vecA, vecB};
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

// Setup Google Test fixture class for the exhaustive fault injection sweep
class MultiplierExhaustiveInjectionTest : public ::testing::Test {
protected:
    Multiplier multiplier; // The hardware multiplier under test
    std::mt19937 rng;      // Pseudo-random number generator engine

    void SetUp() override {
        rng.seed(1337);    // Fixed seed to guarantee identical generation across CI/CD environments
    }
};

// Google Test targeting exhaustive fault analysis
TEST_F(MultiplierExhaustiveInjectionTest, Sweep32BitExhaustiveFaultInjection) {
    // Step 1: Define output path
    const std::string filename = "32bit_exhaustive_multiplier_trials.csv";
    
    // Step 2: Open the CSV file ONCE at the start of the test.
    // Keeping this handle open during the loop prevents heavy disk seek/write latency.
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to open output CSV path: " << filename;

    // Step 3: Write the CSV column header row
    csv << "Trial,CountA,CountB,ExpectedCleanOnes,Organization,BitsFlipped,"
           "Vector_A_Original,Vector_B_Original,Vector_A_Faulted,Vector_Output,"
           "FaultedOnesCount,TotalError,PrecisionError,BitFlipError\n";

    size_t trial = 0;             // Running tally of total trials simulated
    size_t qualifying_pairs = 0;  // Count of A/B ones count combinations evaluated

    // Step 4: Loop exhaustively through all possible ones counts for Stream A (0 to 31)
    for (size_t countA = 0; countA < STREAM_LEN; ++countA) {
        
        // Step 5: Loop exhaustively through all possible ones counts for Stream B (0 to 31)
        for (size_t countB = 0; countB < STREAM_LEN; ++countB) {
            
            // Step 6: Determine expected clean overlap (ideal rounded ones count)
            // This is the optimal ones count assuming zero correlation error.
            size_t expected_overlap = static_cast<size_t>(
                std::round((static_cast<double>(countA) * static_cast<double>(countB)) / STREAM_LEN)
            );
            ++qualifying_pairs;

            // Step 7: Sweep multiple spatial organizations for this combination
            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {
                
                // Step 8: Build the clean stream pair (vecA_orig and vecB_orig) using ZCE mapping
                auto zce_pair = build_zce_pair(countA, countB, expected_overlap, rng);
                const std::vector<bool>& vecA_orig = zce_pair.first;
                const std::vector<bool>& vecB_orig = zce_pair.second;

                // Step 9: Verify mathematical counts of generated streams (safety assertions)
                ASSERT_EQ(count_ones_local(vecA_orig), countA) 
                    << "Stream A generation mismatch for CountA=" << countA;
                ASSERT_EQ(count_ones_local(vecB_orig), countB) 
                    << "Stream B generation mismatch for CountB=" << countB;

                // Step 10: Sweep all possible fault intensities (1 to 32 bit flips)
                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {
                    
                    // Copy original Stream A to inject faults
                    std::vector<bool> vecA_faulted = vecA_orig;
                    
                    // Step 11: Flip the requested number of bits in stream A
                    flip_random_bits(vecA_faulted, flips, rng);

                    // Step 12: Run bit-by-bit hardware multiplication simulation
                    std::vector<bool> vecOut(STREAM_LEN);
                    size_t faulted_ones = 0;
                    for (size_t i = 0; i < STREAM_LEN; ++i) {
                        vecOut[i] = multiplier.multiply(vecA_faulted[i], vecB_orig[i]);
                        if (vecOut[i]) {
                            faulted_ones++;
                        }
                    }

                    // Step 13: Compute different error metrics
                    // Ideal product (unquantized) = (countA * countB) / 32.0
                    double ideal_val = (static_cast<double>(countA) * static_cast<double>(countB)) / static_cast<double>(STREAM_LEN);
                    
                    // Total Error is the difference between the faulted count and the ideal mathematical product
                    double total_error = static_cast<double>(faulted_ones) - ideal_val;
                    
                    // Precision Error is the difference between the clean quantized count and the ideal mathematical product
                    double precision_error = static_cast<double>(expected_overlap) - ideal_val;
                    
                    // Bit Flip Error is the difference between the faulted count and the clean quantized count
                    double bit_flip_error = static_cast<double>(faulted_ones) - static_cast<double>(expected_overlap);

                    // Step 14: Log results directly to the open CSV stream buffer
                    ++trial;
                    csv << trial << "," << countA << "," << countB << "," << expected_overlap << ","
                        << org << "," << flips << ","
                        << vector_to_csv_string(vecA_orig) << ","
                        << vector_to_csv_string(vecB_orig) << ","
                        << vector_to_csv_string(vecA_faulted) << ","
                        << vector_to_csv_string(vecOut) << ","
                        << faulted_ones << "," 
                        << total_error << "," 
                        << precision_error << "," 
                        << bit_flip_error << "\n";
                }
            }
        }
    }

    // Step 15: Close the CSV file handle ONCE at the end of the test.
    csv.close();

    // Print summary stats to the console
    std::cout << "[INFO] Evaluated " << qualifying_pairs << " input count pairs exhaustively." << std::endl;
    std::cout << "[SUCCESS] Exported " << trial << " fault-injection trials to " << filename << std::endl;
}
