#include <gtest/gtest.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

// Include overall project utility definitions
#include "../../general_functions.hpp"
#include "../../modules/multiplier.hpp"

namespace StochasticSimulator {

TEST(MultiplierAccuracyTest, GenerateMatlabCorrelationMatrixData) {
    // =========================================================================
    // 1. CONFIGURATION BLOCK
    // =========================================================================
    double target_pa = 0.5000;            
    double target_pb = 0.5000;            
    const int NUM_TRIALS = 10;            
    
    StreamLength selected_length = StreamLength::Length_1024;
    uint16_t max_cycles = 1024; 

    Multiplier multiplier;

    std::vector<std::vector<bool>> streams_A(NUM_TRIALS);
    std::vector<std::vector<bool>> streams_B(NUM_TRIALS);
    std::vector<std::vector<bool>> results(NUM_TRIALS);

    std::cout << "[INFO] Initiating C++ accuracy sweep pipeline..." << std::endl;

    // =========================================================================
    // 2. GENERATION LOOP
    // =========================================================================
    for (int t = 0; t < NUM_TRIALS; ++t) {
        streams_A[t] = generate_valid_stream(target_pa, selected_length, max_cycles);
        streams_B[t] = generate_valid_stream(target_pb, selected_length, max_cycles);
        results[t].reserve(max_cycles);

        for (uint16_t i = 0; i < max_cycles; ++i) {
            bool out_bit = multiplier.multiply(streams_A[t][i], streams_B[t][i]);
            results[t].push_back(out_bit);
        }
    }

    // =========================================================================
    // 3. ACCURACY TRIALS CSV GENERATION (One Line Per Trial)
    // =========================================================================
    std::string filename = "multiplier_accuracy_trials.csv";
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to write CSV file path.";

    // Header row outlining the exact requested layout
    csv << "Trial,StreamA_Ones_Count_Length,Vector_A,Vector_B,Vector_Output,Result_Ones_Count_Length\n";

    for (int t = 0; t < NUM_TRIALS; ++t) {
        size_t count_a = count_ones(streams_A[t]);
        size_t count_out = count_ones(results[t]);

        // Build fully-expanded bracketed arrays: "[1, 0, 1, ...]"
        std::stringstream mat_vec_a, mat_vec_b, mat_vec_out;
        mat_vec_a << "\"[";
        mat_vec_b << "\"[";
        mat_vec_out << "\"[";
        
        for (uint16_t i = 0; i < max_cycles; ++i) {
            std::string comma_space = (i < max_cycles - 1) ? ", " : "";
            mat_vec_a << streams_A[t][i] << comma_space;
            mat_vec_b << streams_B[t][i] << comma_space;
            mat_vec_out << results[t][i] << comma_space;
        }
        
        mat_vec_a << "]\"";
        mat_vec_b << "]\"";
        mat_vec_out << "]\"";

        // Write row elements mapped strictly to your requested format
        csv << (t + 1) << "," 
            << count_a << "/" << max_cycles << "," 
            << mat_vec_a.str() << "," 
            << mat_vec_b.str() << "," 
            << mat_vec_out.str() << "," 
            << count_out << "/" << max_cycles << "\n";
    }
    csv.close();

    // =========================================================================
    // 4. CONSOLE DISPLAY REPORT (Preserved for human review)
    // =========================================================================
    std::cout << "\n============= SIMULATION EXPERIMENT SUMMARY =============" << std::endl;
    std::cout << "Trial |   Count A   |   Count B   | Count Out | Vector A (First 20b) | Vector B (First 20b)" << std::endl;
    std::cout << "-----------------------------------------------------------------------------------" << std::endl;

    for (int t = 0; t < NUM_TRIALS; ++t) {
        size_t count_a = count_ones(streams_A[t]);
        size_t count_b = count_ones(streams_B[t]);
        size_t count_out = count_ones(results[t]);

        std::string vec_a_short = format_bit_vector(streams_A[t], 20);
        std::string vec_b_short = format_bit_vector(streams_B[t], 20);

        std::cout << std::setw(3) << (t + 1) << "   | " 
                  << std::setw(4) << count_a << "/" << max_cycles << " | " 
                  << std::setw(4) << count_b << "/" << max_cycles << " | " 
                  << std::setw(4) << count_out << "/" << max_cycles << " | "
                  << std::setw(20) << std::left << vec_a_short << " | " 
                  << vec_b_short << std::right << std::endl;
    }

    std::cout << "===================================================================================" << std::endl;
    std::cout << "[SUCCESS] Exported single-row accuracy traces to: " << filename << std::endl;
}

} // namespace StochasticSimulator