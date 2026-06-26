#include <gtest/gtest.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <random>

// Include overall project utility definitions
#include "../../general_functions.hpp"

namespace StochasticSimulator {

TEST(DividerAccuracyTest, GenerateMatlabCorrelationMatrixData) {
    // =========================================================================
    // 1. CONFIGURATION BLOCK
    // =========================================================================
    double target_pa = 0.2500;
    double target_pb = 0.5000;
    const int NUM_TRIALS = 1000;

    StreamLength selected_length = StreamLength::Length_1024;
    uint16_t max_cycles = 1024;

    // Target perfect count: (0.25 / 0.50) * 1024 = 512 ones
    double expected_pout = target_pa / target_pb;
    size_t expected_ones = static_cast<size_t>(std::round(expected_pout * max_cycles));

    std::vector<std::vector<bool>> streams_A;
    std::vector<std::vector<bool>> streams_B;
    std::vector<std::vector<bool>> results;

    streams_A.reserve(NUM_TRIALS);
    streams_B.reserve(NUM_TRIALS);
    results.reserve(NUM_TRIALS);

    std::cout << "[INFO] Initiating C++ divider accuracy sweep pipeline (1000 perfect trials)..." << std::endl;

    // =========================================================================
    // 2. GENERATION LOOP
    // =========================================================================

    std::random_device rd;
    std::mt19937 gen(rd());
    const int COUNTER_MAX = 32;
    std::uniform_int_distribution<int> r_num_gen(0, COUNTER_MAX - 1);

    int successful_trials = 0;
    while (successful_trials < NUM_TRIALS) {
        // Generate streams ONCE per keeper
        std::vector<bool> trial_A = generate_valid_stream(target_pa, selected_length, max_cycles);
        std::vector<bool> trial_B = generate_valid_stream(target_pb, selected_length, max_cycles);

        // Re-run the stochastic simulation on the SAME streams until it hits
        bool found = false;
        for (int attempt = 0; attempt < 200 && !found; ++attempt) {
            int counter = 0;
            std::vector<bool> trial_out;
            trial_out.reserve(max_cycles);
            bool prev_z_bit = false;

            for (uint16_t i = 0; i < max_cycles; ++i) {
                bool x_bit   = trial_A[i];
                bool y_bit   = trial_B[i];
                bool mul_bit = prev_z_bit && y_bit;

                if      (x_bit && !mul_bit && counter < COUNTER_MAX) counter++;
                else if (!x_bit && mul_bit && counter > 0)           counter--;

                int rn     = r_num_gen(gen);
                bool z_bit = (counter > rn);
                trial_out.push_back(z_bit);
                prev_z_bit = z_bit;
            }

            if (count_ones(trial_out) == expected_ones) {
                streams_A.push_back(trial_A);
                streams_B.push_back(trial_B);
                results.push_back(std::move(trial_out));
                successful_trials++;
                found = true;
            }
        }
        // If no hit in 200 attempts, outer loop regenerates fresh streams
    }

    // =========================================================================
    // 3. SINGLE-PASS CSV GENERATION (High Efficiency Batch Write)
    // =========================================================================
    std::string filename = "divider_accuracy_trials.csv";
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to write CSV file path.";

    csv << "Trial,StreamA_Ones_Count_Length,Vector_A,Vector_B,Vector_Output,Result_Ones_Count_Length\n";

    for (int t = 0; t < NUM_TRIALS; ++t) {
        size_t count_a   = count_ones(streams_A[t]);
        size_t count_out = count_ones(results[t]);

        std::stringstream mat_vec_a, mat_vec_b, mat_vec_out;
        mat_vec_a   << "\"[";
        mat_vec_b   << "\"[";
        mat_vec_out << "\"[";

        for (uint16_t i = 0; i < max_cycles; ++i) {
            std::string comma_space = (i < max_cycles - 1) ? ", " : "";
            mat_vec_a   << streams_A[t][i] << comma_space;
            mat_vec_b   << streams_B[t][i] << comma_space;
            mat_vec_out << results[t][i]   << comma_space;
        }

        mat_vec_a   << "]\"";
        mat_vec_b   << "]\"";
        mat_vec_out << "]\"";

        csv << (t + 1) << ","
            << count_a << "/" << max_cycles << ","
            << mat_vec_a.str() << ","
            << mat_vec_b.str() << ","
            << mat_vec_out.str() << ","
            << count_out << "/" << max_cycles << "\n";
    }

    csv.close();
    std::cout << "[SUCCESS] Exported " << NUM_TRIALS << " perfect divider traces to: " << filename << std::endl;
}

} // namespace StochasticSimulator