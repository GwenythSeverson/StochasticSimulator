/*
 * =========================================================================================
 * test_adder_fault2.cpp
 * =========================================================================================
 * EXHAUSTIVE 32-BIT STOCHASTIC ADDER (MUX) FAULT-INJECTION SIMULATOR (DEDUPLICATED WITH PROGRESS)
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
#include <set>
#include <iomanip> // For nice progress formatting

#include "../../modules/adder.hpp"

using namespace StochasticSimulator;

namespace {

constexpr size_t STREAM_LEN = 32;             
constexpr size_t S_ONES = STREAM_LEN / 2;     
constexpr size_t ORGANIZATIONS_PER_PAIR = 1000;  

uint32_t pack_stream_to_uint32(const std::vector<bool>& vec) {
    uint32_t val = 0;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i]) {
            val |= (1U << i);
        }
    }
    return val;
}

std::tuple<std::vector<bool>, std::vector<bool>, std::vector<bool>, size_t>
build_mux_zce_triple(size_t countA, size_t countB, std::mt19937& rng) {
    std::vector<bool> vecA(STREAM_LEN, false);
    std::vector<bool> vecB(STREAM_LEN, false);
    std::vector<bool> vecS(STREAM_LEN, false);

    std::vector<size_t> indices(STREAM_LEN);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), rng);

    std::uniform_int_distribution<int> coin_flip(0, 1);

    // --- Unbiased Boundary-Aware Rounding for Stream A ---
    double exact_a = static_cast<double>(countA) * S_ONES / STREAM_LEN;
    size_t a_in_s;
    if (std::fmod(exact_a, 1.0) == 0.5) {
        size_t down_val = static_cast<size_t>(std::floor(exact_a));
        size_t up_val   = static_cast<size_t>(std::ceil(exact_a));
        
        bool down_legal = (down_val <= countA) && (down_val <= S_ONES) && ((countA - down_val) <= S_ONES);
        bool up_legal   = (up_val <= countA) && (up_val <= S_ONES) && ((countA - up_val) <= S_ONES);
        
        if (down_legal && up_legal) {
            a_in_s = (coin_flip(rng) == 1) ? up_val : down_val;
        } else if (up_legal) {
            a_in_s = up_val;
        } else {
            a_in_s = down_val;
        }
    } else {
        a_in_s = static_cast<size_t>(std::round(exact_a));
    }
    a_in_s = std::min(a_in_s, std::min(countA, S_ONES));
    size_t a_not_in_s = countA - a_in_s;
    if (a_not_in_s > S_ONES) { a_not_in_s = S_ONES; a_in_s = countA - a_not_in_s; }

    // --- Unbiased Boundary-Aware Rounding for Stream B ---
    double exact_b = static_cast<double>(countB) * S_ONES / STREAM_LEN;
    size_t b_in_ns;
    if (std::fmod(exact_b, 1.0) == 0.5) {
        size_t down_val = static_cast<size_t>(std::floor(exact_b));
        size_t up_val   = static_cast<size_t>(std::ceil(exact_b));
        
        bool down_legal = (down_val <= countB) && (down_val <= S_ONES) && ((countB - down_val) <= S_ONES);
        bool up_legal   = (up_val <= countB) && (up_val <= S_ONES) && ((countB - up_val) <= S_ONES);
        
        if (down_legal && up_legal) {
            b_in_ns = (coin_flip(rng) == 1) ? up_val : down_val;
        } else if (up_legal) {
            b_in_ns = up_val;
        } else {
            b_in_ns = down_val;
        }
    } else {
        b_in_ns = static_cast<size_t>(std::round(exact_b));
    }
    b_in_ns = std::min(b_in_ns, std::min(countB, S_ONES));
    size_t b_not_in_ns = countB - b_in_ns;
    if (b_not_in_ns > S_ONES) { b_not_in_ns = S_ONES; b_in_ns = countB - b_not_in_ns; }

    // --- Lay out S=1 positions ---
    for (size_t i = 0; i < S_ONES; ++i) {
        size_t idx = indices[i];
        vecS[idx] = true;
        vecA[idx] = (i < a_in_s);
        vecB[idx] = (i < b_not_in_ns);
    }

    // --- Lay out S=0 positions ---
    for (size_t i = S_ONES; i < STREAM_LEN; ++i) {
        size_t j   = i - S_ONES;
        size_t idx = indices[i];
        vecS[idx] = false;
        vecA[idx] = (j < a_not_in_s);
        vecB[idx] = (j < b_in_ns);
    }

    size_t expected_clean_ones = a_in_s + b_in_ns;
    return {vecA, vecB, vecS, expected_clean_ones};
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

class AdderExhaustiveInjectionTest : public ::testing::Test {
protected:
    Adder adder;
    std::mt19937 rng;

    void SetUp() override {
        rng.seed(1337);
    }
};

TEST_F(AdderExhaustiveInjectionTest, Sweep32BitExhaustiveAdderFaultInjection) {
    const std::string filename = "32bit_exhaustive_adder_trials.csv";

    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to open output CSV path: " << filename;

    csv << "Trial,CountA,CountB,ExpectedCleanOnes,Organization,BitsFlipped,"
           "Vector_A_Original,Vector_B_Original,Vector_S,Vector_A_Faulted,Vector_Output,"
           "FaultedOnesCount,TotalError,PrecisionError,BitFlipError\n";

    size_t trial = 0;
    
    // Progress calculation variables
    constexpr size_t TOTAL_PAIRS = (STREAM_LEN + 1) * (STREAM_LEN + 1); // 33 * 33 = 1089 pairs
    size_t completed_pairs = 0;

    std::cout << "[START] Starting sweep of " << TOTAL_PAIRS << " input pairs..." << std::endl;

    for (size_t countA = 0; countA <= STREAM_LEN; ++countA) {
        for (size_t countB = 0; countB <= STREAM_LEN; ++countB) {

            std::set<std::tuple<uint32_t, uint32_t, uint32_t>> seen_configurations;
            size_t consecutive_duplicate_failures = 0;
            constexpr size_t MAX_DUPLICATE_ATTEMPTS = 200; 

            for (size_t org = 0; org < ORGANIZATIONS_PER_PAIR; ++org) {
                
                std::vector<bool> vecA_orig, vecB_orig, vecS;
                size_t expected_clean_ones = 0;
                bool unique_found = false;

                while (consecutive_duplicate_failures < MAX_DUPLICATE_ATTEMPTS) {
                    auto [temp_A, temp_B, temp_S, clean_ones] = build_mux_zce_triple(countA, countB, rng);
                    
                    uint32_t packedA = pack_stream_to_uint32(temp_A);
                    uint32_t packedB = pack_stream_to_uint32(temp_B);
                    uint32_t packedS = pack_stream_to_uint32(temp_S);
                    
                    auto config_key = std::make_tuple(packedA, packedB, packedS);
                    
                    if (seen_configurations.insert(config_key).second) {
                        vecA_orig = std::move(temp_A);
                        vecB_orig = std::move(temp_B);
                        vecS = std::move(temp_S);
                        expected_clean_ones = clean_ones;
                        unique_found = true;
                        consecutive_duplicate_failures = 0; 
                        break;
                    } else {
                        consecutive_duplicate_failures++;
                    }
                }

                if (!unique_found) {
                    break;
                }

                double ideal_val = (static_cast<double>(countA) + static_cast<double>(countB)) / 2.0;

                ASSERT_EQ(count_ones_local(vecA_orig), countA) << "Stream A mismatch: CountA=" << countA;
                ASSERT_EQ(count_ones_local(vecB_orig), countB) << "Stream B mismatch: CountB=" << countB;
                ASSERT_EQ(count_ones_local(vecS), S_ONES) << "Select stream S mismatch";

                for (size_t flips = 1; flips <= STREAM_LEN; ++flips) {
                    std::vector<bool> vecA_faulted = vecA_orig;
                    flip_random_bits(vecA_faulted, flips, rng);

                    std::vector<bool> vecOut(STREAM_LEN);
                    size_t faulted_ones = 0;
                    for (size_t i = 0; i < STREAM_LEN; ++i) {
                        vecOut[i] = adder.add_scaled_or_weighted(
                            vecA_faulted[i], vecB_orig[i], vecS[i]);
                        if (vecOut[i]) faulted_ones++;
                    }

                    double total_error     = static_cast<double>(faulted_ones) - ideal_val;
                    double precision_error = static_cast<double>(expected_clean_ones) - ideal_val;
                    double bit_flip_error  = static_cast<double>(faulted_ones) - static_cast<double>(expected_clean_ones);

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

            // --- PROGRESS UPDATER ---
            ++completed_pairs;
            if (completed_pairs % 25 == 0 || completed_pairs == TOTAL_PAIRS) {
                double percent = (static_cast<double>(completed_pairs) / TOTAL_PAIRS) * 100.0;
                std::cout << "[PROGRESS] " << std::setw(6) << std::fixed << std::setprecision(2) << percent << "% | "
                          << "Processed " << completed_pairs << "/" << TOTAL_PAIRS << " pairs | "
                          << "Active: (A=" << countA << ", B=" << countB << ") | "
                          << "Total Trials: " << trial << std::endl;
            }
        }
    }

    csv.close();

    std::cout << "[SUCCESS] Completed execution!" << std::endl;
    std::cout << "[SUCCESS] Exported " << trial
              << " unique trials to " << filename << std::endl;
}