#include "lfsr.hpp"
#include <stdexcept>

namespace StochasticSimulator {

FlexibleLFSR::FlexibleLFSR(StreamLength lengthMode, uint16_t seed) {
    if (seed == 0) {
        throw std::invalid_argument("Seed cannot be 0. :(");
    }

    initial_seed = seed;
    state = seed;
    cycle_count = 0;

    switch (lengthMode) {
        case StreamLength::Length_128:
            max_cycles = 128;
            polynomial_mask = 0x0041; // 7-bit
            break;
        case StreamLength::Length_256:
            max_cycles = 256;
            polynomial_mask = 0x008C; // 8-bit
            break;
        case StreamLength::Length_512:
            max_cycles = 512;
            polynomial_mask = 0x0120; // 9-bit
            break;
        case StreamLength::Length_1024:
            max_cycles = 1024;
            polynomial_mask = 0x0204; // 10-bit
            break;
        case StreamLength::Length_4096:
            max_cycles = 4096;
            polynomial_mask = 0x0829; // 12-bit
            break;
        case StreamLength::Length_16384:
            max_cycles = 16384;
            polynomial_mask = 0x202C; // 14-bit
            break;
        case StreamLength::Length_65536:
            max_cycles = 65536;
            polynomial_mask = 0xB400; // 16-bit prolly slow
            break;
    }
    
    if (seed >= max_cycles) {
        throw std::invalid_argument("Seed is out of bounds for the selected stream length.");
    }
}

uint16_t FlexibleLFSR::next() {
    // Artificial 0 injection to cleanly complete the power-of-2 sequence length
    if (cycle_count == (max_cycles - 1)) {
        cycle_count = 0;
        return 0;
    }

    bool output_bit = state & 1;
    state >>= 1;
    if (output_bit) {
        state ^= polynomial_mask;
    }

    cycle_count++;
    return state;
}

uint16_t FlexibleLFSR::get_max_cycles() const {
    return max_cycles;
}

void FlexibleLFSR::reset() {
    state = initial_seed;
    cycle_count = 0;
}

}  // namespace StochasticSimulator