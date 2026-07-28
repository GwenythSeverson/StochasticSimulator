
/**
 * Linear feedback shift register (LFSR) implementation for stochastic number generation. 
 *Between fibbonacci and Galois LFSR implementations, Ill start with Galois Configuration
 * because the the XOR gates are in between the shoft registers rather than rippling through
 * multiple external ones in the fibb lfsr design. 
 * to be expanded upon. 
 */

 #pragma once
#include <cstdint>
// LFSR MODE SELECTION:
// 1. Auto-Random Mode (For multi-run accuracy tests / avoiding correlation bias):
//    -> FlexibleLFSR stream(StreamLength::Length_256); // Omit seed to auto-randomize
//
// 2. Deterministic Mode (For debugging, fault injection, & early termination tests):
//    -> FlexibleLFSR stream(StreamLength::Length_256, 42); // Pass explicit hardware seed
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
private:
    uint16_t initial_seed;
    uint16_t state;
    uint16_t polynomial_mask;
    // MUST be 32-bit: Length_65536 sets this to 65536, which does not fit in a uint16_t and
    // silently truncated to 0. That made every seed fail the bounds check below, so every
    // Length_65536 construction with an explicit seed threw.
    uint32_t max_cycles;
    uint32_t cycle_count;

public:
    // Setting default seed to 0 triggers automatic time-randomization
    FlexibleLFSR(StreamLength lengthMode, uint16_t seed = 0);

    uint16_t next();
    uint32_t get_max_cycles() const;
    void reset();
};

} // namespace StochasticSimulator