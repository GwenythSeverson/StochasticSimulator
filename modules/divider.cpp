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
double ud_counter_division(const std::vector<bool>& stream_X, const std::vector<bool>& stream_Y, uint32_t sng_seed) {

    if (stream_X.empty() || stream_X.size() != stream_Y.size()) {
        std::cerr << "Error: Input streams must be non-empty and of identical length.\n";
        return 0.0;
    }

    size_t max_cycles = stream_X.size();
    

// extra note here to know this is where you would probably input new code for any pre loading like this 
// counter max/2 section that causes the warm up tim eot be cut in half in this situation so the proability starts around 1/2 
// and can get its proabbility close and closer to the actual needed values faster at the start. but for now 
// its set to counter starts at zero with a larger warmup time for higher proability values. but it allows for an asnwer to be always found and not stuck at zero for a long time.
// zeros in the end. 
    const int COUNTER_MAX = 32; 
    int counter = 0;
    
    std::vector<bool> stream_Z;
    stream_Z.reserve(max_cycles);

  // We use a 16-bit LFSR (Length_65536) to ensure the random sequence does not loop for 
    // at least 65,536 cycles, which is critical for long 1024+ bit streams.
    // We ensure the seed is non-zero (LFSRs hate 0).
    uint16_t safe_seed = (sng_seed % 65535) + 1; 
    FlexibleLFSR lfsr(StreamLength::Length_65536, safe_seed);
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

            // Hardware LFSR draw. We need a 5-bit value (0-31) for the COUNTER_MAX=32 comparator.
        // We take the lower 5 bits of the 16-bit LFSR output (equivalent to tapping 5 wires).
        int rn = lfsr.next() % 32;

        // Comparator logic: If the integrated count is greater than the random threshold,
        // we output a 1. This spreads the 1s out randomly instead of clumping.
        bool current_z_bit = (counter > rn);
        stream_Z.push_back(current_z_bit);

        // Save for the next cycle's feedback loop
        prev_z_bit = current_z_bit;
    }

    return calculate_probability(stream_Z);
}
}  // namespace StochasticSimulator
