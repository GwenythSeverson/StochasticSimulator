#include "general_functions.hpp"
#include "bsg/sng.hpp"
#include "bsg/lfsr.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include "general_functions.hpp"



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
    int counter = 0; 
    const int COUNTER_MAX = 32; 
    
    std::vector<bool> stream_Z;
    stream_Z.reserve(max_cycles);

    // Track z_bit from the *previous* cycle for the feedback loop
    bool prev_z_bit = false; 

    for (size_t i = 0; i < max_cycles; ++i) {
        bool x_bit = stream_X[i];
        bool y_bit = stream_Y[i];

        // Feedback uses the registered state from the previous clock cycle
        bool mul_bit = prev_z_bit && y_bit;

        // 1. Update the counter state
        if (x_bit && !mul_bit) {
            if (counter < COUNTER_MAX) counter++;
        } else if (!x_bit && mul_bit) {
            if (counter > 0) counter--;
        }

        // 2. Combinational Output: Determine current bit based on the updated state
        bool current_z_bit = (counter > 0);
        stream_Z.push_back(current_z_bit);

        // Save for the next cycle's feedback loop
        prev_z_bit = current_z_bit;
    }

    return calculate_probability(stream_Z);
}
}  // namespace StochasticSimulator
