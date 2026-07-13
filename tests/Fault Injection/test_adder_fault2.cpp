/*
 * =========================================================================================
 * test_adder_fault2.cpp
 * =========================================================================================
 * EXHAUSTIVE 32-BIT STOCHASTIC ADDER (MUX) FAULT-INJECTION SIMULATOR
 *
 * Mathematical & Architectural Overview:
 * 1. Unary representation in Stochastic Computing:
 *    A value in the interval [0, 1) is represented as a bitstream of length N.
 *    For a 32-bit stream, there are 32 possible unique fractional numbers represented
 *    by ones counts ranging from 0/32 to 31/32 inclusive.
 *
 * 2. Stochastic Addition via MUX:
 *    The stochastic adder uses a Multiplexer (MUX) gate. A third "select" stream S
 *    chooses between inputs A and B on each clock cycle:
 *      Output[i] = S[i] ? A[i] : B[i]
 *    When S has a 50% density (16 ones out of 32), the output represents (A + B) / 2.
 *    This is scaled addition — the result is always halved.
 *
 * 3. ZCE Co-Construction (Key Difference from Previous Version):
 *    Streams A, B, and S are built TOGETHER from a single shared shuffled index map
 *    (identical strategy to the multiplier's build_zce_pair). This guarantees the
 *    clean output has exactly the expected ideal ones count — zero correlation error
 *    in the baseline. Without this, PrecisionError mixes quantization noise AND
 *    spatial correlation noise, making the two error sources indistinguishable.
 *
 *    Construction layout (shared index, first 16 = S=1 positions, last 16 = S=0):
 *      S=1 positions: A[i]=1 for the first a_in_s slots; B[i]=1 for first b_not_in_ns slots
 *      S=0 positions: A[i]=1 for the first a_not_in_s slots; B[i]=1 for first b_in_ns slots
 *      ExpectedCleanOnes = a_in_s + b_in_ns (exact, no rounding variance in baseline)
 *
 * 4. Exhaustive Input Coverage:
 *    We sweep all ones counts for both inputs:
 *      CountA in [0, 31] and CountB in [0, 31].
 *    This covers all 32 x 32 = 1024 input probability pairs.
 *
 * 5. Fault Injection Sweep:
 *    For each organization, we inject fault intensities from 1 to 32 bits flipped.
 *    Only stream A is faulted; streams B and S remain clean.
 *
 * 6. Single-Open File Output (High Performance):
 *    The CSV file is opened ONCE at the start, all trial records are streamed
 *    directly to the file buffer, and it is closed ONCE at the end.
 *
 * 7. Error Components:
 *    - Ideal Sum (real): (CountA + CountB) / 2.0   [MUX halving]
 *    - Expected Clean Ones Count: a_in_s + b_in_ns  [ZCE-locked, no variance]
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
#include <tuple>
#include <iostream>

#include "../../modules/adder.hpp"

using namespace StochasticSimulator;

namespace {

// Definition of constants for the 32-bit simulator sweep.
constexpr size_t STREAM_LEN = 32;             // Length of the stochastic bitstreams
constexpr size_t S_ONES = STREAM_LEN / 2;     // Select stream density = 50% = 16 ones
constexpr size_t ORGANIZATIONS_PER_PAIR = 100;  // Unique bit permutations per input pair

/*
 * build_mux_zce_triple()
 *
 * Constructs streams A, B, and S together from a single shared shuffled index so
 * that the MUX output has EXACTLY (a_in_s + b_in_ns) ones — zero correlation error
 * in the baseline. This mirrors what build_zce_pair does for the multiplier.
 *
 * Layout of the shared shuffled indices:
 *   indices[0 .. S_ONES-1]     → S=1 positions (MUX selects A)
 *   indices[S_ONES .. LEN-1]   → S=0 positions (MUX selects B)
 *
 * Within S=1 positions:
 *   first a_in_s      → A=1  (these contribute 1s to output via MUX selecting A)
 *   first b_not_in_ns → B=1  (B has 1s here but MUX ignores them; keeps countB correct)
 *
 * Within S=0 positions:
 *   first a_not_in_s  → A=1  (A has 1s here but MUX ignores them; keeps countA correct)
 *   first b_in_ns     → B=1  (these contribute 1s to output via MUX selecting B)
 *
 * Returns: {vecA, vecB, vecS, expected_clean_ones}
 */
std::tuple<std::vector<bool>, std::vector<bool>, std::vector<bool>, size_t>
build_mux_zce_triple(size_t countA, size_t countB, std::mt19937& rng) {

    std::vector<bool> vecA(STREAM_LEN, false);
    std::vector<bool> vecB(STREAM_LEN, false);
    std::vector<bool> vecS(STREAM_LEN, false);

    // Shared index shuffle — all three streams are laid out from this single map
    std::vector<size_t> indices(STREAM_LEN);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    // How many of A's ones go in S=1 positions (these flow through to the output)
    size_t a_in_s = static_cast<size_t>(
        std::round(static_cast<double>(countA) * S_ONES / STREAM_LEN));
    // Clamp so we never ask for more ones than slots exist in each region
    a_in_s = std::min(a_in_s, std::min(countA, S_ONES));
    size_t a_not_in_s = countA - a_in_s;
    if (a_not_in_s > S_ONES) { a_not_in_s = S_ONES; a_in_s = countA - a_not_in_s; }

    // How many of B's ones go in S=0 positions (these flow through to the output)
    size_t b_in_ns = static_cast<size_t>(
        std::round(static_cast<double>(countB) * S_ONES / STREAM_LEN));
    b_in_ns = std::min(b_in_ns, std::min(countB, S_ONES));
    size_t b_not_in_ns = countB - b_in_ns;
    if (b_not_in_ns > S_ONES) { b_not_in_ns = S_ONES; b_in_ns = countB - b_not_in_ns; }

    // --- Lay out S=1 positions (indices[0 .. S_ONES-1]) ---
    for (size_t i = 0; i < S_ONES; ++i) {
        size_t idx = indices[i];
        vecS[idx] = true;
        vecA[idx] = (i < a_in_s);        // contributing A ones
        vecB[idx] = (i < b_not_in_ns);   // non-contributing B ones (MUX ignores)
    }

    // --- Lay out S=0 positions (indices[S_ONES .. STREAM_LEN-1]) ---
    for (size_t i = S_ONES; i < STREAM_LEN; ++i) {
        size_t j   = i - S_ONES;
        size_t idx = indices[i];
        vecS[idx] = false;
        vecA[idx] = (j < a_not_in_s);    // non-contributing A ones (MUX ignores)
        vecB[idx] = (j < b_in_ns);       // contributing B ones
    }

    // The clean output ones count is exact — no rounding variance
    size_t expected_clean_ones = a_in_s + b_in_ns;

    return {vecA, vecB, vecS, expected_clean_ones};
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
    if (num_flips >= stream.size()) {
        stream.flip();
        return;
    }
    std::vector<size_t> indices(stream.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);
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
        rng.seed(1337);    // Fixed seed for reproducible results
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
           "Vector_A_Original,Vector_B_Original,Vector_S,Vector_A_Faulted,Vector_Output,"
           "FaultedOnesCount,TotalError,PrecisionError,BitFlipError\n";

    size_t trial = 0;
    size_t qualifying_pairs = 0;

    // Step 4: Loop exhaustively through all possible ones counts for Stream A (0 to 31)
    for (size_t countA = 0; countA < STREAM_LEN; ++countA) {

        // Step 5: Loop exhaustively through all possible ones counts for Stream B (0 to 31)
        for (size_t countB = 0; countB < STREAM_LEN; ++countB) {

            ++qualifying_pairs;

            // Step 6: Sweep multiple spatial organizations for this combination
            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {

                // Step 7: Build the ZCE-locked (A, B, S) triple from a shared index map.
                // expected_clean_ones is exact — no correlation noise in the baseline.
                auto [vecA_orig, vecB_orig, vecS, expected_clean_ones] =
                    build_mux_zce_triple(countA, countB, rng);

                // Step 8: Mathematical ideal for a perfect 50-50 MUX adder
                double ideal_val = (static_cast<double>(countA) + static_cast<double>(countB)) / 2.0;

                // Step 9: Verify stream generation counts (safety assertions)
                ASSERT_EQ(count_ones_local(vecA_orig), countA)
                    << "Stream A mismatch: CountA=" << countA;
                ASSERT_EQ(count_ones_local(vecB_orig), countB)
                    << "Stream B mismatch: CountB=" << countB;
                ASSERT_EQ(count_ones_local(vecS), S_ONES)
                    << "Select stream S mismatch";

                // Step 10: Sweep all possible fault intensities (1 to 32 bit flips)
                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {

                    // Copy original Stream A to inject faults; B and S stay clean
                    std::vector<bool> vecA_faulted = vecA_orig;
                    flip_random_bits(vecA_faulted, flips, rng);

                    // Step 11: Run bit-by-bit hardware MUX addition simulation
                    std::vector<bool> vecOut(STREAM_LEN);
                    size_t faulted_ones = 0;
                    for (size_t i = 0; i < STREAM_LEN; ++i) {
                        vecOut[i] = adder.add_scaled_or_weighted(
                            vecA_faulted[i], vecB_orig[i], vecS[i]);
                        if (vecOut[i]) faulted_ones++;
                    }

                    // Step 12: Compute error metrics
                    // Total Error: faulted output vs. mathematical ideal (A+B)/2
                    double total_error    = static_cast<double>(faulted_ones) - ideal_val;
                    // Precision Error: clean quantized output vs. mathematical ideal
                    // (now zero-correlation, so this is pure quantization rounding only)
                    double precision_error = static_cast<double>(expected_clean_ones) - ideal_val;
                    // Bit Flip Error: faulted output vs. clean quantized output
                    double bit_flip_error  = static_cast<double>(faulted_ones)
                                           - static_cast<double>(expected_clean_ones);

                    // Step 13: Log results
                    ++trial;
                    csv << trial << "," << countA << "," << countB << ","
                        << expected_clean_ones << "," << org << "," << flips << ","
                        << vector_to_csv_string(vecA_orig) << ","
                        << vector_to_csv_string(vecB_orig) << ","
                        << vector_to_csv_string(vecS) << ","
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

    csv.close();

    std::cout << "[INFO] Evaluated " << qualifying_pairs
              << " input count pairs exhaustively." << std::endl;
    std::cout << "[SUCCESS] Exported " << trial
              << " adder fault-injection trials to " << filename << std::endl;
}
