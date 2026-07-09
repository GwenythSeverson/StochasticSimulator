/*
 * =========================================================================================
 * test_adder_fault2.cpp
 * =========================================================================================
 * EXHAUSTIVE 32-BIT STOCHASTIC ADDER (MUX) FAULT-INJECTION SIMULATOR
 *
 * Mathematical & Architectural Overview:
 * 1. Unary representation in Stochastic Computing:
 *    A value in the interval [0, 1) is represented as a bitstream of length N.
 *    For a 32-bit stream, there are 32 possible unique fractional numbers represented by
 *    ones counts ranging from 0/32 to 31/32 inclusive.
 *
 * 2. Stochastic Addition via MUX:
 *    The stochastic adder uses a Multiplexer (MUX) gate. A third "select" stream S
 *    chooses between inputs A and B on each clock cycle:
 *      Output[i] = S[i] ? A[i] : B[i]
 *    When S has a 50% density (16 ones out of 32), the output represents (A + B) / 2.
 *    This is scaled addition: the result is always halved.
 *
 * 3. Exhaustive Input Coverage:
 *    We sweep all possible ones counts for both inputs:
 *      CountA in [0, 31] and CountB in [0, 31].
 *    This covers all 32 x 32 = 1024 input probability pairs.
 *
 * 4. Select Stream:
 *    The select stream S is generated with exactly 16 ones (50% density) for each
 *    organization, then shuffled randomly. This ensures unbiased scaled addition.
 *
 * 5. Fault Injection Sweep:
 *    For each organization, we inject fault intensities from 1 to 32 bits flipped.
 *    Only stream A is faulted; stream B and S remain clean.
 *
 * 6. Single-Open File Output (High Performance):
 *    The CSV file is opened ONCE at the start, all 163,840 trial records are streamed
 *    directly to the file buffer, and it is closed ONCE at the end.
 *
 * 7. Error Components:
 *    - Ideal Sum (real): (CountA + CountB) / (2 * 32.0)   [scaled by MUX halving]
 *    - Ideal Ones Count: IdealVal = (CountA + CountB) / 2.0
 *    - Clean Output Ones Count (quantized): ExpectedCleanOnes = counted from clean sim
 *    - Faulted Output Ones Count: FaultedOnesCount
 *    - Precision Error: ExpectedCleanOnes - IdealVal
 *    - Bit Flip Error: FaultedOnesCount - ExpectedCleanOnes
 *    - Total Error: FaultedOnesCount - IdealVal [Total = Precision + Bit Flip]
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

#include "../../modules/adder.hpp"

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

// Setup Google Test fixture class for the exhaustive adder fault injection sweep
class AdderExhaustiveInjectionTest : public ::testing::Test {
protected:
    Adder adder;           // The hardware adder (MUX) under test
    std::mt19937 rng;      // Pseudo-random number generator engine

    void SetUp() override {
        rng.seed(1337);    // Fixed seed to guarantee identical generation across CI/CD environments
    }
};

// Google Test targeting exhaustive adder fault analysis
TEST_F(AdderExhaustiveInjectionTest, Sweep32BitExhaustiveAdderFaultInjection) {
    // Step 1: Define output path
    const std::string filename = "32bit_exhaustive_adder_trials.csv";

    // Step 2: Open the CSV file ONCE at the start of the test.
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

            ++qualifying_pairs;

            // Step 6: Sweep multiple spatial organizations for this combination
            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {

                // Step 7: Build clean streams A, B, and select stream S (50% density = 16 ones)
                std::vector<bool> vecA_orig = build_stream(countA, rng);
                std::vector<bool> vecB_orig = build_stream(countB, rng);
                std::vector<bool> vecS = build_stream(STREAM_LEN / 2, rng); // 16 ones for 50% select

                // Step 8: Compute clean output ones count by running the MUX adder without faults
                size_t clean_ones = 0;
                for (size_t i = 0; i < STREAM_LEN; ++i) {
                    if (adder.add_scaled_or_weighted(vecA_orig[i], vecB_orig[i], vecS[i])) {
                        clean_ones++;
                    }
                }
                size_t expected_clean_ones = clean_ones;

                // Step 9: Compute the ideal (unquantized) ones count for scaled addition
                // MUX addition: output = (A + B) / 2, so ideal ones = (countA + countB) / 2.0
                double ideal_val = (static_cast<double>(countA) + static_cast<double>(countB)) / 2.0;

                // Step 10: Verify stream generation counts (safety assertions)
                ASSERT_EQ(count_ones_local(vecA_orig), countA)
                    << "Stream A generation mismatch for CountA=" << countA;
                ASSERT_EQ(count_ones_local(vecB_orig), countB)
                    << "Stream B generation mismatch for CountB=" << countB;

                // Step 11: Sweep all possible fault intensities (1 to 32 bit flips)
                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {

                    // Copy original Stream A to inject faults
                    std::vector<bool> vecA_faulted = vecA_orig;

                    // Step 12: Flip the requested number of bits in stream A
                    flip_random_bits(vecA_faulted, flips, rng);

                    // Step 13: Run bit-by-bit hardware MUX addition simulation
                    std::vector<bool> vecOut(STREAM_LEN);
                    size_t faulted_ones = 0;
                    for (size_t i = 0; i < STREAM_LEN; ++i) {
                        vecOut[i] = adder.add_scaled_or_weighted(vecA_faulted[i], vecB_orig[i], vecS[i]);
                        if (vecOut[i]) {
                            faulted_ones++;
                        }
                    }

                    // Step 14: Compute different error metrics
                    // Total Error: difference between faulted output and ideal mathematical sum
                    double total_error = static_cast<double>(faulted_ones) - ideal_val;

                    // Precision Error: difference between clean quantized output and ideal
                    double precision_error = static_cast<double>(expected_clean_ones) - ideal_val;

                    // Bit Flip Error: difference between faulted output and clean quantized output
                    double bit_flip_error = static_cast<double>(faulted_ones) - static_cast<double>(expected_clean_ones);

                    // Step 15: Log results directly to the open CSV stream buffer
                    ++trial;
                    csv << trial << "," << countA << "," << countB << "," << expected_clean_ones << ","
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

    // Step 16: Close the CSV file handle ONCE at the end of the test.
    csv.close();

    // Print summary stats to the console
    std::cout << "[INFO] Evaluated " << qualifying_pairs << " input count pairs exhaustively." << std::endl;
    std::cout << "[SUCCESS] Exported " << trial << " adder fault-injection trials to " << filename << std::endl;
}
