// uMUL: unipolar unary multiplier with conditional bit stream generation.
// uGEMM Figure 3(a) -- see uMUL.hpp for the datapath sketch and the rationale.
// https://jsm.ece.wisc.edu/docs/wu-toppicks2021.pdf 
//
// "only when input 0 is logic one, the RNG inside the bit stream generator for input 1 will
//  update, and the generator will eventually generate a new bit, i.e., input 0 is an enable
//  signal to the bit stream generator. As such, we thoroughly eliminate the correlation problem
//  and achieve high accuracy with merely an extra enable signal."

#include "uMUL.hpp"

#include "general_functions.hpp"

#include <stdexcept>
#include <string>

namespace StochasticSimulator {
namespace {

unsigned validated_width(unsigned width_bits) {
    if (width_bits == 0 || width_bits > UnaryMultiplier::MAX_WIDTH) {
        throw std::invalid_argument(
            "uMUL width must be between 1 and " + std::to_string(UnaryMultiplier::MAX_WIDTH) +
            "; the counter is 2^width flip-flops, so wider than that is not hardware any more.");
    }
    return width_bits;
}

}  // namespace

// Members initialize in declaration order, so `entry` and `window` can lean on `width`.
UnaryMultiplier::UnaryMultiplier(unsigned width_bits, unsigned sobol_dimension)
    : width(validated_width(width_bits)),
      entry(1ull << width),
      window(static_cast<std::size_t>(entry)),
      oldest(0),
      ones_count(0),
      rng(width, sobol_dimension),
      rng_index(0),
      enabled_cycles(0) {
    reset();
}

// Delegates through SobolRNG so the StreamLength -> width mapping lives in exactly one place.
UnaryMultiplier::UnaryMultiplier(StreamLength lengthMode, unsigned sobol_dimension)
    : UnaryMultiplier(SobolRNG(lengthMode, sobol_dimension).get_width(), sobol_dimension) {}

/**
 * Two error terms pull in opposite directions (see SIZING in the header): warm-up costs roughly
 * 2^width / N, counter quantization costs roughly 1 / 2^width. Their sum is smallest where
 * 2^width ~ sqrt(N); measuring the actual RMSE across N = 256..4096 puts the optimum a factor of
 * two above that, so the window is sized to 2 * sqrt(N):
 *     width = 1 + floor(log2(N) / 2)
 * That reproduced the measured best width at N = 256, 1024 and 4096, and came within one step
 * (0.007 vs 0.006 RMSE) at N = 512.
 */
unsigned UnaryMultiplier::width_for_stream_length(std::size_t stream_length) {
    if (stream_length == 0) {
        throw std::invalid_argument("uMUL stream length must be at least 1.");
    }

    unsigned log2_length = 0;
    while ((static_cast<std::size_t>(1) << (log2_length + 1)) <= stream_length) {
        ++log2_length;
    }

    unsigned width_bits = 1 + log2_length / 2;
    return (width_bits > MAX_WIDTH) ? MAX_WIDTH : width_bits;
}

UnaryMultiplier UnaryMultiplier::for_stream_length(std::size_t stream_length,
                                                   unsigned sobol_dimension) {
    return UnaryMultiplier(width_for_stream_length(stream_length), sobol_dimension);
}

void UnaryMultiplier::reset() {
    // Power-on contents of the shift register are 0,1,0,1,... -- exactly half ones, so the unit
    // starts out assuming p1 = 0.5. The window cannot be right before it has seen 2^width bits of
    // in_1, and a neutral guess costs far less warm-up error than starting the tally at zero
    // would (which would hold the output at 0 for the whole fill).
    for (std::size_t i = 0; i < window.size(); ++i) {
        window[i] = (i % 2 == 1);
    }
    oldest = 0;  // slot 0 is the bit that shifts out first
    ones_count = static_cast<uint32_t>(entry / 2);

    rng.reset();
    rng_index = 0;
    enabled_cycles = 0;
}

bool UnaryMultiplier::multiply(bool in_0, bool in_1) {
    // ---- C: counter -------------------------------------------------------------------------
    // Shifts on every cycle whether or not the enable fires -- C is tracking in_1's value, it is
    // not producing output. The oldest bit falls out, in_1 falls in, and the running tally
    // follows that one swap, which is what keeps the read O(1) instead of a popcount.
    bool departing = window[oldest];
    window[oldest] = in_1;
    oldest = (oldest + 1) % window.size();

    if (in_1 && !departing) {
        ++ones_count;
    } else if (!in_1 && departing) {
        --ones_count;
    }

    // ---- G: bit stream generator ------------------------------------------------------------
    // ones_count is BINARY here (the thick line in Figure 3(a)) -- every other wire in this unit
    // carries a single unary bit, this one carries an integer in [0, 2^width]. The comparator is
    // the same one BitstreamGenerator and SobolRNG::next_bit use: truncating threshold, strict
    // ">". ones_count == entry therefore pins the regenerated stream to all ones.
    uint32_t random_value = rng.value_at(rng_index);
    bool regenerated_bit = ones_count > random_value;

    // ---- en: the enable ---------------------------------------------------------------------
    // The whole trick, in one branch. in_0 clocks G's RNG index forward; on an in_0 == 0 cycle the
    // RNG holds its value, so no random number is burned on a cycle the AND gate was going to
    // zero anyway. G therefore walks its sequence in order across exactly the cycles that can
    // carry a 1, and lands p1 of them -- regardless of how in_0's and in_1's ones line up. That
    // is why uMUL has no correlation requirement and the plain AND-gate Multiplier does.
    // Read first, then advance: this cycle's output used the value the register already held.
    if (in_0) {
        rng_index = (rng_index + 1) % entry;
        ++enabled_cycles;
    }

    // ---- ANDMUL -----------------------------------------------------------------------------
    return in_0 && regenerated_bit;
}

uint32_t UnaryMultiplier::get_counter_value() const {
    return ones_count;
}

double UnaryMultiplier::get_tracked_probability() const {
    return static_cast<double>(ones_count) / static_cast<double>(entry);
}

uint64_t UnaryMultiplier::get_enabled_cycles() const {
    return enabled_cycles;
}

unsigned UnaryMultiplier::get_width() const {
    return width;
}

uint64_t UnaryMultiplier::get_entry() const {
    return entry;
}

std::vector<bool> umul_stream(const std::vector<bool>& stream_0,
                              const std::vector<bool>& stream_1,
                              unsigned width,
                              unsigned sobol_dimension) {
    if (stream_0.empty() || stream_0.size() != stream_1.size()) {
        throw std::invalid_argument("uMUL input streams must be non-empty and of identical length.");
    }

    if (width == UnaryMultiplier::AUTO_WIDTH) {
        width = UnaryMultiplier::width_for_stream_length(stream_0.size());
    }
    UnaryMultiplier umul(width, sobol_dimension);

    std::vector<bool> stream_out;
    stream_out.reserve(stream_0.size());
    for (std::size_t i = 0; i < stream_0.size(); ++i) {
        stream_out.push_back(umul.multiply(stream_0[i], stream_1[i]));
    }
    return stream_out;
}

double umul_probability(const std::vector<bool>& stream_0,
                        const std::vector<bool>& stream_1,
                        unsigned width,
                        unsigned sobol_dimension) {
    return calculate_probability(umul_stream(stream_0, stream_1, width, sobol_dimension));
}

}  // namespace StochasticSimulator
