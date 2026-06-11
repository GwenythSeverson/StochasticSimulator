// Stochastic number generator (SNG) implementation using a linear feedback shift register (LFSR)
#include "sng.hpp"
#include <stdexcept>

namespace StochasticSimulator {

bool BitstreamGenerator::generate_bit(double probability, uint16_t random_val, uint16_t max_scale) const {
    if (probability < 0.0 || probability > 1.0) {
        throw std::invalid_argument("Probability must be between 0.0 and 1.0");
    }

    // Dynamic scale mapping (maps 0.0-1.0 to the current LFSR bit-width bounds)
    uint16_t threshold = static_cast<uint16_t>(probability * static_cast<double>(max_scale));
    
    // Hardware digital comparator logic
    return threshold > random_val;
}

} // namespace StochasticSimulator