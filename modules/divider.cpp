// Gaines division (GDIV), unipolar -- see divider.hpp for the loop diagram and the equilibrium
// derivation. https://jsm.ece.wisc.edu/docs/wu-ieeedt2021.pdf figure 2a.
//
// "Gaines division (GDIV) [1] for unipolar and bipolar SC is shown in Figure 2a and b."

#include "divider.hpp"

#include "bsg/lfsr.hpp"
#include "general_functions.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace StochasticSimulator {
namespace {

/**
 * The comparator RNG runs on Length_16384, deliberately NOT Length_65536.
 *
 * FlexibleLFSR stores max_cycles in a uint16_t, so Length_65536 assigns 65536 into 16 bits and
 * silently truncates it to 0. That makes the constructor's "seed >= max_cycles" bounds check true
 * for every seed, so the old Length_65536 line here threw std::invalid_argument on every single
 * call -- the divider could not run at all. A 14-bit LFSR has a 16384-cycle period, which is
 * ample for the comparator, and sidesteps the truncation entirely.
 */
constexpr StreamLength RNG_MODE = StreamLength::Length_16384;
constexpr uint32_t RNG_SCALE = 16384;  // LFSR outputs land in [0, RNG_SCALE)

bool is_power_of_two(unsigned v) {
    return v != 0 && (v & (v - 1)) == 0;
}

// Shared front-end check. Kept as a plain bool so both entry points report the same way.
bool inputs_are_valid(const std::vector<bool>& stream_X, const std::vector<bool>& stream_Y) {
    if (stream_X.empty() || stream_X.size() != stream_Y.size()) {
        std::cerr << "Error: Input streams must be non-empty and of identical length.\n";
        return false;
    }
    return true;
}

}  // namespace

unsigned gaines_counter_depth_for_length(std::size_t stream_length) {
    if (stream_length == 0) {
        throw std::invalid_argument("Gaines divider stream length must be at least 1.");
    }

    // Smallest power of two d with 2*d*d >= stream_length, i.e. d ~ sqrt(N/2). Measured against
    // an exhaustive depth sweep this picked the best depth at N = 512, 1024, 4096 and 16384, and
    // tied with the best at N = 50000.
    unsigned depth = 1;
    while (static_cast<std::size_t>(depth) * depth * 2 < stream_length && depth < 1024) {
        depth *= 2;
    }

    if (depth < 16) depth = 16;        // below this the 1/D quantization dominates everything
    if (depth > 1024) depth = 1024;    // above this nothing realistic settles in time
    return depth;
}

std::vector<bool> gaines_division_stream(const std::vector<bool>& stream_X,
                                         const std::vector<bool>& stream_Y,
                                         uint32_t sng_seed,
                                         unsigned counter_depth) {
    if (!inputs_are_valid(stream_X, stream_Y)) {
        return {};
    }

    const std::size_t max_cycles = stream_X.size();

    if (counter_depth == GAINES_AUTO_DEPTH) {
        counter_depth = gaines_counter_depth_for_length(max_cycles);
    }
    // A power of two that divides RNG_SCALE keeps the comparator scaling exact -- no rounding
    // sneaks into the one place where an off-by-one would bias every output bit.
    if (!is_power_of_two(counter_depth) || counter_depth > RNG_SCALE) {
        throw std::invalid_argument(
            "Gaines counter depth must be a power of two no greater than 16384.");
    }

    const int depth = static_cast<int>(counter_depth);
    const uint32_t step = RNG_SCALE / counter_depth;  // scales the counter up to the RNG's range

    // Start at mid-scale rather than 0. The counter has to walk to its equilibrium either way,
    // and starting halfway means the worst-case walk is D/2 instead of D. Measured at N = 1024,
    // depth 32: RMSE 0.028 starting at zero against 0.021 starting at mid.
    int counter = depth / 2;

    // Seed must stay inside the LFSR's own bounds, and an LFSR seeded with 0 is stuck forever.
    uint16_t safe_seed = static_cast<uint16_t>((sng_seed % (RNG_SCALE - 1)) + 1);
    FlexibleLFSR lfsr(RNG_MODE, safe_seed);

    std::vector<bool> stream_Z;
    stream_Z.reserve(max_cycles);

    // The feedback path is registered: this cycle's AND uses LAST cycle's Z, which is what the
    // hardware does and what keeps the loop from being combinationally circular.
    bool prev_z_bit = false;

    for (std::size_t i = 0; i < max_cycles; ++i) {
        bool x_bit = stream_X[i];
        bool y_bit = stream_Y[i];

        bool mul_bit = prev_z_bit && y_bit;

        // Saturating up/down counter. Holding when x_bit == mul_bit is the whole point: the
        // counter integrates the DIFFERENCE between X and Z*Y, and parks where they agree.
        if (x_bit && !mul_bit) {
            if (counter < depth) ++counter;
        } else if (!x_bit && mul_bit) {
            if (counter > 0) --counter;
        }

        // SNG comparator over the LFSR's full width. counter == depth scales to exactly
        // RNG_SCALE, which beats every possible draw, so a saturated counter pins Z to 1.
        uint32_t random_value = lfsr.next();
        bool current_z_bit = static_cast<uint32_t>(counter) * step > random_value;

        stream_Z.push_back(current_z_bit);
        prev_z_bit = current_z_bit;
    }

    return stream_Z;
}

double ud_counter_division(const std::vector<bool>& stream_X,
                           const std::vector<bool>& stream_Y,
                           uint32_t sng_seed,
                           unsigned counter_depth) {
    std::vector<bool> stream_Z =
        gaines_division_stream(stream_X, stream_Y, sng_seed, counter_depth);

    // Empty means the front-end check rejected the inputs; see the note in divider.hpp about why
    // this collapses to the same 0.0 a genuine zero quotient would produce.
    if (stream_Z.empty()) {
        return 0.0;
    }
    return calculate_probability(stream_Z);
}

}  // namespace StochasticSimulator
