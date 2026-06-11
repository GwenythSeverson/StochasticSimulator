
/**
 * Linear feedback shift register (LFSR) implementation for stochastic number generation. 
 *Between fibbonacci and Galois LFSR implementations, Ill start with Galois Configuration
 * because the the XOR gates are in between the shoft registers rather than rippling through
 * multiple external ones in the fibb lfsr design. 
 * to be expanded upon. 
 */

 #pragma once
#include <cstdint>

namespace StochasticSimulator {

enum class StreamLength {
    Length_128,    // 7-bit LFSR
    Length_256,    // 8-bit LFSR
    Length_512,    // 9-bit LFSR
    Length_1024,   // 10-bit LFSR
    Length_4096,   // 12-bit LFSR
    Length_16384,  // 14-bit LFSR
    Length_65536   // 16-bit LFSR
};

class FlexibleLFSR {
public:
    explicit FlexibleLFSR(StreamLength lengthMode = StreamLength::Length_1024, uint16_t seed = 0x01);

    uint16_t next();
    uint16_t get_max_cycles() const;
    void reset();

private:
    uint16_t state;
    uint16_t initial_seed;
    uint16_t cycle_count;
    uint16_t max_cycles;
    uint16_t polynomial_mask;
};

}  // namespace StochasticSimulator