#pragma once
#include <cstdint>

namespace StochasticSimulator {

class BitstreamGenerator {
public:
    BitstreamGenerator() = default;

    /**
     * @brief Generates a single stochastic bit based on a probability threshold.
     * @param probability A value between 0.0 and 1.0.
     * @param random_val A pseudo-random integer supplied by the LFSR.
     * @param max_scale The maximum possible value the LFSR can output (e.g., 128, 1024).
     * @return true (1) if the threshold > random_val, false (0) otherwise.
     */
    bool generate_bit(double probability, uint16_t random_val, uint16_t max_scale) const;
};

} // namespace StochasticSimulator