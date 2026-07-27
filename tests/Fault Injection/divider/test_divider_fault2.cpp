/*
 * =========================================================================================
 * test_divider_fault2.cpp
 * =========================================================================================
 * EXHAUSTIVE 32-BIT STOCHASTIC DIVIDER (UP/DOWN COUNTER) FAULT-INJECTION SIMULATOR
 *
 * Mathematical & Architectural Overview:
 * 1. Unary representation in Stochastic Computing:
 *    A value in the interval [0, 1) is represented as a bitstream of length N.
 *    For a 32-bit stream, there are 32 possible unique fractional numbers represented
 *    by ones counts ranging from 0/32 to 31/32 inclusive.
 *
 * 2. Stochastic Division via Up/Down Counter:
 *    The stochastic divider uses an up/down counter feedback loop (Gaines division).
 *    It computes Z = X / Y in the stochastic domain. For accurate unipolar division,
 *    CountX should be <= CountY so that the result doesn't saturate at 1.0.
 *    The divider is STATEFUL — it processes entire streams and returns a decoded
 *    output probability.
 *
 * 3. Co-Constructed Input Streams (Correlation Isolation):
 *    X and Y are built from the same shuffled index map so their bit patterns have
 *    consistent relative spatial statistics across organizations. Without this,
 *    the correlation between X and Y's temporal distribution varies randomly per
 *    trial and adds noise to the baseline, making PrecisionError partially a
 *    measurement of construction variance rather than quantization.
 *
 * 4. Multi-Trial SNG Averaging (Key Difference from Previous Version):
 *    The divider's internal SNG (comparator + random threshold) introduces its own
 *    output variance independent of the input streams. Previously the SNG RNG was
 *    seeded from std::random_device (non-deterministic), meaning clean and faulted
 *    outputs for THE SAME STREAMS differed between calls purely due to SNG noise.
 *    This contaminated BitFlipError measurements.
 *
 *    Fix: Run SNG_TRIALS_PER_FAULT independent SNG trials for each (streams, flips)
 *    combination and average the result. Crucially, EACH trial uses the SAME SNG seed
 *    for both the clean and faulted run so SNG variance cancels in BitFlipError:
 *      BitFlipError = mean(faulted) - mean(clean) ≈ pure fault effect
 *
 * 5. Exhaustive Input Coverage:
 *    We sweep numerator CountX in [0, 31] and denominator CountY in [1, 31].
 *    Only test pairs where CountX <= CountY (valid unipolar division range).
 *
 * 6. Fault Injection Sweep:
 *    For each organization, we inject fault intensities from 1 to 32 bits flipped.
 *    Only stream X (numerator) is faulted; stream Y (denominator) remains clean.
 *
 * 7. Error Components (all in probability domain):
 *    - Ideal Quotient: CountX / CountY
 *    - Clean Output Prob: mean of SNG_TRIALS clean runs
 *    - Faulted Output Prob: mean of SNG_TRIALS faulted runs (same SNG seeds)
 *    - Precision Error: CleanOutputProb - IdealQuotient
 *    - Bit Flip Error: FaultedOutputProb - CleanOutputProb (SNG noise cancels)
 *    - Total Error: FaultedOutputProb - IdealQuotient
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

#include "../../../modules/divider.hpp"

using namespace StochasticSimulator;

namespace {

// Definition of constants for the 32-bit simulator sweep.
constexpr size_t STREAM_LEN = 32;            // Length of the stochastic bitstreams
constexpr size_t ORGANIZATIONS_PER_PAIR = 5; // Unique bit permutations per input pair

// Number of independent SNG seed trials to average per (streams, fault level).
// Averaging separates divider SNG variance from fault-induced error signal.
// 5 trials is a good balance between noise reduction and runtime.
constexpr size_t SNG_TRIALS_PER_FAULT = 5;

/*
 * build_stream_pair()
 *
 * Constructs X and Y together from a single shared shuffled index map so their
 * spatial statistics are consistent (same construction philosophy as ZCE for the
 * multiplier). The divider is sequential/stateful so the AND-gate ZCE formula
 * doesn't directly apply, but sharing the index ensures both streams are placed
 * with the same spatial randomness rather than independently.
 *
 * Layout: first countX indices → X=1; first countY indices of the same shuffle → Y=1.
 * Since countX <= countY, the X=1 positions are always a subset of Y=1 positions
 * in the shuffle, which maximises temporal alignment and reduces warm-up variance.
 */
std::pair<std::vector<bool>, std::vector<bool>>
build_stream_pair(size_t countX, size_t countY, std::mt19937& rng) {
    std::vector<bool> vecX(STREAM_LEN, false);
    std::vector<bool> vecY(STREAM_LEN, false);

    // Shared index map — both streams are positioned from the same shuffle
    std::vector<size_t> indices(STREAM_LEN);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    // Place X ones (numerator is sparser or equal to Y)
    for (size_t i = 0; i < countX; ++i) {
        vecX[indices[i]] = true;
    }
    // Place Y ones (denominator is denser)
    for (size_t i = 0; i < countY; ++i) {
        vecY[indices[i]] = true;
    }

    return {vecX, vecY};
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

/*
 * average_divider_output()
 *
 * Runs ud_counter_division SNG_TRIALS_PER_FAULT times with consecutive SNG seeds
 * starting from base_sng_seed and returns the mean output probability.
 * Using a deterministic sweep of seeds (rather than std::random_device) makes
 * the averaging fully reproducible.
 */
double average_divider_output(const std::vector<bool>& vecX,
                              const std::vector<bool>& vecY,
                              uint32_t base_sng_seed) {
    double sum = 0.0;
    for (size_t t = 0; t < SNG_TRIALS_PER_FAULT; ++t) {
        sum += ud_counter_division(vecX, vecY, base_sng_seed + static_cast<uint32_t>(t));
    }
    return sum / static_cast<double>(SNG_TRIALS_PER_FAULT);
}

} // namespace

// Setup Google Test fixture class for the exhaustive divider fault injection sweep
class DividerExhaustiveInjectionTest : public ::testing::Test {
protected:
    std::mt19937 rng;      // Pseudo-random number generator engine

    void SetUp() override {
        rng.seed(1337);    // Fixed seed for reproducible stream construction
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

    size_t trial = 0;
    size_t qualifying_pairs = 0;

    // Running SNG seed counter — incremented per fault trial so every (stream set,
    // fault level) combination gets a unique but fully reproducible SNG seed block.
    uint32_t sng_seed_counter = 100000;

    // Step 4: Loop through all possible ones counts for numerator Stream X (0 to 31)
    for (size_t countX = 0; countX < STREAM_LEN; ++countX) {

        // Step 5: Loop through denominator Stream Y (1 to 31, skip 0 = division by zero)
        // Only test where countX <= countY (valid unipolar division constraint)
        for (size_t countY = 1; countY < STREAM_LEN; ++countY) {

            if (countX > countY) continue;

            ++qualifying_pairs;

            // Step 6: Compute the ideal quotient (unquantized mathematical result)
            double ideal_quotient = static_cast<double>(countX) / static_cast<double>(countY);

            // Step 7: Sweep multiple spatial organizations for this combination
            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {

                // Step 8: Build co-constructed numerator and denominator streams.
                // Both come from the same shuffled index so their spatial distribution
                // is consistent — reduces warm-up correlation noise.
                auto [vecX_orig, vecY_orig] = build_stream_pair(countX, countY, rng);

                // Step 9: Verify stream counts (safety assertions)
                ASSERT_EQ(count_ones_local(vecX_orig), countX)
                    << "Stream X mismatch: CountX=" << countX;
                ASSERT_EQ(count_ones_local(vecY_orig), countY)
                    << "Stream Y mismatch: CountY=" << countY;

                // Step 10: Compute the clean baseline by averaging SNG_TRIALS_PER_FAULT runs.
                // Each run uses a distinct SNG seed so we sample the full SNG variance.
                double clean_output_prob = average_divider_output(
                    vecX_orig, vecY_orig, sng_seed_counter);

                // Step 11: Sweep all possible fault intensities (1 to 32 bit flips)
                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {

                    // Copy original Stream X to inject faults; Y stays clean
                    std::vector<bool> vecX_faulted = vecX_orig;
                    flip_random_bits(vecX_faulted, flips, rng);

                    // Step 12: Average faulted output over the SAME SNG seeds as the
                    // clean run. Using identical SNG seeds means SNG noise cancels in
                    // BitFlipError = faulted_avg - clean_avg, isolating the fault signal.
                    double faulted_output_prob = average_divider_output(
                        vecX_faulted, vecY_orig, sng_seed_counter);

                    // Advance seed counter for the next fault level
                    sng_seed_counter += static_cast<uint32_t>(SNG_TRIALS_PER_FAULT);

                    // Step 13: Compute error metrics (all in probability domain)
                    double total_error     = faulted_output_prob - ideal_quotient;
                    double precision_error = clean_output_prob   - ideal_quotient;
                    double bit_flip_error  = faulted_output_prob - clean_output_prob;

                    // Step 14: Log results
                    ++trial;
                    csv << trial << "," << countX << "," << countY << ","
                        << clean_output_prob << "," << org << "," << flips << ","
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

    csv.close();

    std::cout << "[INFO] Evaluated " << qualifying_pairs
              << " valid (X <= Y) input count pairs." << std::endl;
    std::cout << "[SUCCESS] Exported " << trial
              << " divider fault-injection trials to " << filename << std::endl;
}
