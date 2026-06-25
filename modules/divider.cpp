#include "general_functions.hpp"
#include "bsg/sng.hpp"
#include "bsg/lfsr.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include "general_functions.hpp"
#include <random>


// https://jsm.ece.wisc.edu/docs/wu-ieeedt2021.pdf divider from figure 2a
// "Gaines division (GDIV) [1] for unipolar and bipolar SC is shown in Figure 2a and b."
namespace StochasticSimulator {

    /*
    * up down counter division
    */
double ud_counter_division(const std::vector<bool>& stream_X, const std::vector<bool>& stream_Y) {
    if (stream_X.empty() || stream_X.size() != stream_Y.size()) {
        std::cerr << "Error: Input streams must be non-empty and of identical length.\n";
        return 0.0;
    }

    size_t max_cycles = stream_X.size();
    
    // Using a 0-indexed counter to match hardware SNG behavior cleanly
    int counter = 0; 
    const int COUNTER_MAX = 32; 
    
    std::vector<bool> stream_Z;
    stream_Z.reserve(max_cycles);

    // Setup random number generation to act as the hardware RN block (e.g., an LFSR)
    std::random_device rd;
    std::mt19937 gen(rd());
    // Generates a uniform random threshold integer between 0 and COUNTER_MAX - 1
    std::uniform_int_distribution<int> r_num_gen(0, COUNTER_MAX - 1);

    // Track z_bit from the *previous* cycle for the feedback loop
    bool prev_z_bit = false; 

    for (size_t i = 0; i < max_cycles; ++i) {
        bool x_bit = stream_X[i];
        bool y_bit = stream_Y[i];

        // Feedback uses the registered state from the previous clock cycle
        bool mul_bit = prev_z_bit && y_bit;

        // 1. Update the counter state (Saturating bounds: 0 to COUNTER_MAX)
        if (x_bit && !mul_bit) {
            if (counter < COUNTER_MAX) counter++;
        } else if (!x_bit && mul_bit) {
            if (counter > 0) counter--;
        }

        // 2. Stochastic Number Generator (SNG) Block Replacement
        // Generate a fresh random number (RN) for this clock cycle
        int rn = r_num_gen(gen);

        // Comparator logic: If the integrated count is greater than the random threshold,
        // we output a 1. This spreads the 1s out randomly instead of clumping.
        bool current_z_bit = (counter >= rn);
        stream_Z.push_back(current_z_bit);

        // Save for the next cycle's feedback loop
        prev_z_bit = current_z_bit;
    }

    return calculate_probability(stream_Z);
}
}  // namespace StochasticSimulator
