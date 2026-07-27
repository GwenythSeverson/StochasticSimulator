/*#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <string>
#include <cmath>
#include <set>
#include <utility>

#include "../../../modules/multiplier.hpp"

/*
 * 32-bit exhaustive-fraction zero-correlation-error (ZCE) fault injection sweep
 * for the AND-gate multiplier.
 *
 * OUTPUT: 32bit_ZCE_multiplier_trials.csv
 *   One row per (fraction pair x organization x bits-flipped) trial. Faults are
 *   injected into stream A only; stream B remains clean. MATLAB can then:
 *     - group by BitsFlipped and average Error -> error vs. #bits-flipped curve
 *     - group by (CountA, CountB) and average Error across organizations
 *     - confirm coverage of all 33 possible 32-bit fractions
 */

/*using namespace StochasticSimulator;

namespace {

constexpr size_t STREAM_LEN = 32;
constexpr size_t ORGANIZATIONS_PER_PAIR = 5;

// Builds a pair of 32-bit streams with an EXACT zero-correlation-error overlap
// for the given ones-counts, using a shared shuffled index map. Deterministic
// construction guarantees ZCE=0 rather than relying on random search.
std::pair<std::vector<bool>, std::vector<bool>> build_zce_pair(
    size_t countA, size_t countB, size_t expected_overlap, std::mt19937& rng) {

    std::vector<bool> vecA(STREAM_LEN, false);
    std::vector<bool> vecB(STREAM_LEN, false);

    std::vector<size_t> indices(STREAM_LEN);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    for (size_t i = 0; i < STREAM_LEN; ++i) {
        size_t idx = indices[i];
        if (i < expected_overlap) {
            vecA[idx] = true;
            vecB[idx] = true;
        } else if (i < countA) {
            vecA[idx] = true;
            vecB[idx] = false;
        } else if (i < countA + (countB - expected_overlap)) {
            vecA[idx] = false;
            vecB[idx] = true;
        } else {
            vecA[idx] = false;
            vecB[idx] = false;
        }
    }
    return {vecA, vecB};
}

std::string vector_to_csv_string(const std::vector<bool>& vec) {
    std::string s = "\"[";
    for (size_t i = 0; i < vec.size(); ++i) {
        s += (vec[i] ? "1" : "0");
        if (i < vec.size() - 1) s += ", ";
    }
    s += "]\"";
    return s;
}

size_t count_ones_local(const std::vector<bool>& vec) {
    return static_cast<size_t>(std::count(vec.begin(), vec.end(), true));
}

// Flips exactly num_flips distinct random bit positions in-place.
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

class MultiplierInjectionTest : public ::testing::Test {
protected:
    Multiplier multiplier;
    std::mt19937 rng;

    void SetUp() override {
        rng.seed(1337); // deterministic across CI runs
    }
};

TEST_F(MultiplierInjectionTest, Sweep32BitZeroCorrelationFaultInjection) {
    const std::string filename = "32bit_ZCE_multiplier_trials.csv";
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to open output CSV path.";

    csv << "Trial,CountA,CountB,ExpectedCleanOnes,Organization,BitsFlipped,"
           "Vector_A_Original,Vector_B_Original,Vector_A_Faulted,Vector_Output,"
           "FaultedOnesCount,Error\n";

    std::set<size_t> countA_values_seen;
    size_t trial = 0;
    size_t qualifying_pairs = 0;

    for (size_t countA = 0; countA <= STREAM_LEN; ++countA) {
        for (size_t countB = 0; countB <= STREAM_LEN; ++countB) {
            size_t product = countA * countB;
            if (product % STREAM_LEN != 0) {
                continue; // no exact zero-correlation-error version exists at 32 bits
            }
            size_t expected_overlap = product / STREAM_LEN;
            ++qualifying_pairs;

            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {
                auto zce_pair = build_zce_pair(countA, countB, expected_overlap, rng);
                const std::vector<bool>& vecA_orig = zce_pair.first;
                const std::vector<bool>& vecB_orig = zce_pair.second;

                // Sanity: verify this construction really is zero-correlation-error
                size_t clean_ones = 0;
                for (size_t i = 0; i < STREAM_LEN; ++i) {
                    if (multiplier.multiply(vecA_orig[i], vecB_orig[i])) clean_ones++;
                }
                ASSERT_EQ(clean_ones, expected_overlap)
                    << "ZCE construction failed for countA=" << countA
                    << ", countB=" << countB << ", org=" << org;
                // Also confirm the ones counts themselves are exactly as requested.
                ASSERT_EQ(count_ones_local(vecA_orig), countA);
                ASSERT_EQ(count_ones_local(vecB_orig), countB);

                countA_values_seen.insert(countA);

                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {
                    std::vector<bool> vecA_faulted = vecA_orig;
                    flip_random_bits(vecA_faulted, flips, rng);

                    std::vector<bool> vecOut(STREAM_LEN);
                    size_t faulted_ones = 0;
                    for (size_t i = 0; i < STREAM_LEN; ++i) {
                        vecOut[i] = multiplier.multiply(vecA_faulted[i], vecB_orig[i]);
                        if (vecOut[i]) faulted_ones++;
                    }

                    double error = static_cast<double>(faulted_ones) -
                                    static_cast<double>(expected_overlap);

                    ++trial;
                    csv << trial << "," << countA << "," << countB << "," << expected_overlap << ","
                        << org << "," << flips << ","
                        << vector_to_csv_string(vecA_orig) << ","
                        << vector_to_csv_string(vecB_orig) << ","
                        << vector_to_csv_string(vecA_faulted) << ","
                        << vector_to_csv_string(vecOut) << ","
                        << faulted_ones << "," << error << "\n";
                }
            }
        }
    }

    csv.close();

    // Confirm every possible 32-bit probability (0/32 ... 32/32) appears as countA
    // in at least one valid zero-correlation trial ("does every representable
    // fraction actually exist in the data").
    for (size_t k = 0; k <= STREAM_LEN; ++k) {
        EXPECT_TRUE(countA_values_seen.count(k))
            << "Fraction " << k << "/32 never appeared as a valid countA - coverage gap!";
    }

    std::cout << "[INFO] Qualifying (countA,countB) pairs at 32 bits: " << qualifying_pairs << std::endl;
    std::cout << "[SUCCESS] Exported " << trial << " fault-injection trials to " << filename << std::endl;
}*/

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

#include "../../../modules/multiplier.hpp"

/*
 * 32-bit exhaustive-fraction near zero-correlation-error (ZCE) fault injection sweep
 * for the AND-gate multiplier with perfectly balanced trial counts (~1,000 per group).
 *
 * OUTPUT: 32bit_ZCE_multiplier_trials.csv
 * Balanced layout: 7 (B values) * 5 (orgs) * 32 (flips) = 1120 trials per CountA step.
 */

using namespace StochasticSimulator;

namespace {

constexpr size_t STREAM_LEN = 32;
constexpr size_t UNIQUE_B_SAMPLES_PER_A = 7;
constexpr size_t ORGANIZATIONS_PER_PAIR = 5;

// Builds a pair of 32-bit streams with an EXACT or near-exact zero-correlation-error
// overlap for the given ones-counts, using a shared shuffled index map.
std::pair<std::vector<bool>, std::vector<bool>> build_zce_pair(
    size_t countA, size_t countB, size_t expected_overlap, std::mt19937& rng) {

    std::vector<bool> vecA(STREAM_LEN, false);
    std::vector<bool> vecB(STREAM_LEN, false);

    std::vector<size_t> indices(STREAM_LEN);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    for (size_t i = 0; i < STREAM_LEN; ++i) {
        size_t idx = indices[i];
        if (i < expected_overlap) {
            vecA[idx] = true;
            vecB[idx] = true;
        } else if (i < countA) {
            vecA[idx] = true;
            vecB[idx] = false;
        } else if (i < countA + (countB - expected_overlap)) {
            vecA[idx] = false;
            vecB[idx] = true;
        } else {
            vecA[idx] = false;
            vecB[idx] = false;
        }
    }
    return {vecA, vecB};
}

std::string vector_to_csv_string(const std::vector<bool>& vec) {
    std::string s = "\"[";
    for (size_t i = 0; i < vec.size(); ++i) {
        s += (vec[i] ? "1" : "0");
        if (i < vec.size() - 1) s += ", ";
    }
    s += "]\"";
    return s;
}

size_t count_ones_local(const std::vector<bool>& vec) {
    return static_cast<size_t>(std::count(vec.begin(), vec.end(), true));
}

// Flips exactly num_flips distinct random bit positions in-place.
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

class MultiplierInjectionTest : public ::testing::Test {
protected:
    Multiplier multiplier;
    std::mt19937 rng;

    void SetUp() override {
        rng.seed(1337); // deterministic seed across runs
    }
};

TEST_F(MultiplierInjectionTest, Sweep32BitZeroCorrelationFaultInjection) {
    const std::string filename = "32bit_ZCE_multiplier_trials.csv";
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to open output CSV path.";

    csv << "Trial,CountA,CountB,ExpectedCleanOnes,Organization,BitsFlipped,"
           "Vector_A_Original,Vector_B_Original,Vector_A_Faulted,Vector_Output,"
           "FaultedOnesCount,Error\n";

    std::set<size_t> countA_values_seen;
    size_t trial = 0;
    size_t qualifying_pairs = 0;

    // Scan mid-range active probabilities uniformly
    for (size_t countA = 1; countA < STREAM_LEN; ++countA) {
        
        // 1. Generate available active CountB choices (1 to 31)
        std::vector<size_t> possible_B_values(STREAM_LEN - 1);
        std::iota(possible_B_values.begin(), possible_B_values.end(), 1);
        
        // 2. Shuffle to select an un-biased mix of cross-stream densities
        std::shuffle(possible_B_values.begin(), possible_B_values.end(), rng);

        for (size_t b_idx = 0; b_idx < UNIQUE_B_SAMPLES_PER_A; ++b_idx) {
            size_t countB = possible_B_values[b_idx];
            
            // Calculate best mathematical ZCE overlap baseline using rounding
            size_t expected_overlap = static_cast<size_t>(
                std::round((static_cast<double>(countA) * static_cast<double>(countB)) / STREAM_LEN)
            );
            ++qualifying_pairs;

            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {
                auto zce_pair = build_zce_pair(countA, countB, expected_overlap, rng);
                const std::vector<bool>& vecA_orig = zce_pair.first;
                const std::vector<bool>& vecB_orig = zce_pair.second;

                // Validate original stream bit properties
                ASSERT_EQ(count_ones_local(vecA_orig), countA);
                ASSERT_EQ(count_ones_local(vecB_orig), countB);
                countA_values_seen.insert(countA);

                // Run uniform fault injection sweep
                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {
                    std::vector<bool> vecA_faulted = vecA_orig;
                    flip_random_bits(vecA_faulted, flips, rng);

                    std::vector<bool> vecOut(STREAM_LEN);
                    size_t faulted_ones = 0;
                    for (size_t i = 0; i < STREAM_LEN; ++i) {
                        vecOut[i] = multiplier.multiply(vecA_faulted[i], vecB_orig[i]);
                        if (vecOut[i]) faulted_ones++;
                    }

                    double error = static_cast<double>(faulted_ones) -
                                   static_cast<double>(expected_overlap);

                    ++trial;
                    csv << trial << "," << countA << "," << countB << "," << expected_overlap << ","
                        << org << "," << flips << ","
                        << vector_to_csv_string(vecA_orig) << ","
                        << vector_to_csv_string(vecB_orig) << ","
                        << vector_to_csv_string(vecA_faulted) << ","
                        << vector_to_csv_string(vecOut) << ","
                        << faulted_ones << "," << error << "\n";
                }
            }
        }
    }

    csv.close();

    // Verify mathematical coverage properties for active domain (1 to 31)
    for (size_t k = 1; k < STREAM_LEN; ++k) {
        EXPECT_TRUE(countA_values_seen.count(k))
            << "Fraction " << k << "/32 missing from balanced execution run!";
    }

    std::cout << "[INFO] Uniformly logged " << UNIQUE_B_SAMPLES_PER_A 
              << " random cross-configurations per probability step." << std::endl;
    std::cout << "[INFO] Qualifying (countA,countB) stream test pairs: " << qualifying_pairs << std::endl;
    std::cout << "[SUCCESS] Exported " << trial << " fault-injection trials to " << filename << std::endl;
}
