// SuMUL: unipolar unary multiplier with conditional bit stream generation.
// uGEMM Figure 3(a) -- see SuMUL.hpp for the datapath sketch and the rationale.
// https://jsm.ece.wisc.edu/docs/wu-toppicks2021.pdf 
//
// "only when input 0 is logic one, the RNG inside the bit stream generator for input 1 will
//  update, and the generator will eventually generate a new bit, i.e., input 0 is an enable
//  signal to the bit stream generator. As such, we thoroughly eliminate the correlation problem
//  and achieve high accuracy with merely an extra enable signal."

#include "SuMUL.hpp"

#include "general_functions.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
#endif

namespace StochasticSimulator {
namespace {

// Population count of one 64-bit word -- the leaf of the parallel-counter tree. Uses the CPU
// instruction where available and falls back to the standard SWAR reduction otherwise, so the
// counter's arithmetic is identical on every target.
inline uint32_t popcount64(uint64_t v) {
#if defined(_MSC_VER) && defined(_M_X64)
    return static_cast<uint32_t>(__popcnt64(v));
#elif defined(__GNUC__) || defined(__clang__)
    return static_cast<uint32_t>(__builtin_popcountll(v));
#else
    v = v - ((v >> 1) & 0x5555555555555555ull);
    v = (v & 0x3333333333333333ull) + ((v >> 2) & 0x3333333333333333ull);
    v = (v + (v >> 4)) & 0x0F0F0F0F0F0F0F0Full;
    return static_cast<uint32_t>((v * 0x0101010101010101ull) >> 56);
#endif
}

inline bool read_bit(const std::vector<uint64_t>& words, std::size_t bit) {
    return ((words[bit >> 6] >> (bit & 63)) & 1ull) != 0ull;
}

inline void write_bit(std::vector<uint64_t>& words, std::size_t bit, bool value) {
    uint64_t mask = 1ull << (bit & 63);
    if (value) {
        words[bit >> 6] |= mask;
    } else {
        words[bit >> 6] &= ~mask;
    }
}

std::size_t words_for(uint64_t bits) {
    return static_cast<std::size_t>((bits + 63) / 64);
}

unsigned validated_width(unsigned width_bits) {
    if (width_bits == 0 || width_bits > SUnaryMultiplier::MAX_WIDTH) {
        throw std::invalid_argument(
            "SuMUL width must be between 1 and " + std::to_string(SUnaryMultiplier::MAX_WIDTH) +
            "; the counter is 2^width flip-flops, so wider than that is not hardware any more.");
    }
    return width_bits;
}

}  // namespace

// Members initialize in declaration order, so `entry` and `window` can lean on `width`.
SUnaryMultiplier::SUnaryMultiplier(unsigned width_bits, unsigned sobol_dimension)
    : width(validated_width(width_bits)),
      entry(1ull << width),
      window(words_for(entry), 0ull),
      oldest(0),
      rng(width, sobol_dimension),
      rng_index(0),
      enabled_cycles(0) {
    reset();
}

// Delegates through SobolRNG so the StreamLength -> width mapping lives in exactly one place.
SUnaryMultiplier::SUnaryMultiplier(StreamLength lengthMode, unsigned sobol_dimension)
    : SUnaryMultiplier(SobolRNG(lengthMode, sobol_dimension).get_width(), sobol_dimension) {}

/**
 * Two error terms pull in opposite directions (see SIZING in the header): warm-up costs roughly
 * 2^width / N, counter quantization costs roughly 1 / 2^width. Their sum is smallest where
 * 2^width ~ sqrt(N); measuring the actual RMSE across N = 256..4096 puts the optimum a factor of
 * two above that, so the window is sized to 2 * sqrt(N):
 *     width = 1 + floor(log2(N) / 2)
 * That reproduced the measured best width at N = 256, 1024 and 4096, and came within one step
 * (0.007 vs 0.006 RMSE) at N = 512.
 */
unsigned SUnaryMultiplier::width_for_stream_length(std::size_t stream_length) {
    if (stream_length == 0) {
        throw std::invalid_argument("SuMUL stream length must be at least 1.");
    }

    unsigned log2_length = 0;
    while ((static_cast<std::size_t>(1) << (log2_length + 1)) <= stream_length) {
        ++log2_length;
    }

    unsigned width_bits = 1 + log2_length / 2;
    return (width_bits > MAX_WIDTH) ? MAX_WIDTH : width_bits;
}

SUnaryMultiplier SUnaryMultiplier::for_stream_length(std::size_t stream_length,
                                                   unsigned sobol_dimension) {
    return SUnaryMultiplier(width_for_stream_length(stream_length), sobol_dimension);
}

void SUnaryMultiplier::reset() {
    // Power-on contents of the shift register are 0,1,0,1,... -- exactly half ones, so the unit
    // starts out assuming p1 = 0.5. The window cannot be right before it has seen 2^width bits of
    // in_1, and a neutral guess costs far less warm-up error than starting the tally at zero
    // would (which would hold the output at 0 for the whole fill).
    // Clear the padding bits in the last word too, so the popcount only ever sees real cells.
    std::fill(window.begin(), window.end(), 0ull);
    for (std::size_t i = 0; i < static_cast<std::size_t>(entry); ++i) {
        write_bit(window, i, (i % 2 == 1));
    }
    oldest = 0;  // slot 0 is the bit that shifts out first

    rng.reset();
    rng_index = 0;
    enabled_cycles = 0;
}

bool SUnaryMultiplier::multiply(bool in_0, bool in_1) {
    // ---- C: counter -------------------------------------------------------------------------
    // Shifts on every cycle whether or not the enable fires -- C is tracking in_1's value, it is
    // not producing output. Writing in_1 over the oldest cell and advancing the head IS the shift:
    // the departing bit is simply overwritten.
    write_bit(window, oldest, in_1);
    oldest = (oldest + 1) % static_cast<std::size_t>(entry);

    // The count is derived from the cells, never carried forward -- see the note on `window` in
    // the header. This is the parallel counter, and it is what makes C self-scrubbing.
    uint32_t ones_count = get_counter_value();

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
    // is why SuMUL has no correlation requirement and the plain AND-gate Multiplier does.
    // Read first, then advance: this cycle's output used the value the register already held.
    if (in_0) {
        rng_index = (rng_index + 1) % entry;
        ++enabled_cycles;
    }

    // ---- ANDMUL -----------------------------------------------------------------------------
    return in_0 && regenerated_bit;
}

uint32_t SUnaryMultiplier::get_counter_value() const {
    // Parallel counter over the shift register. Padding bits above `entry` are held at zero by
    // reset() and never written, so summing whole words is exact.
    uint32_t total = 0;
    for (std::size_t w = 0; w < window.size(); ++w) {
        total += popcount64(window[w]);
    }
    return total;
}

std::size_t SUnaryMultiplier::physical_index(std::size_t cell) const {
    if (cell >= static_cast<std::size_t>(entry)) {
        throw std::out_of_range("SuMUL window cell index is past the end of the shift register.");
    }
    // Cell 0 is the bit that shifts out next, which lives at `oldest`.
    return (oldest + cell) % static_cast<std::size_t>(entry);
}

std::size_t SUnaryMultiplier::window_cells() const {
    return static_cast<std::size_t>(entry);
}

bool SUnaryMultiplier::get_window_cell(std::size_t cell) const {
    return read_bit(window, physical_index(cell));
}

void SUnaryMultiplier::set_window_cell(std::size_t cell, bool value) {
    write_bit(window, physical_index(cell), value);
}

void SUnaryMultiplier::flip_window_cell(std::size_t cell) {
    std::size_t bit = physical_index(cell);
    write_bit(window, bit, !read_bit(window, bit));
}

uint64_t SUnaryMultiplier::get_rng_index() const {
    return rng_index;
}

void SUnaryMultiplier::set_rng_index(uint64_t index) {
    rng_index = index % entry;
}

double SUnaryMultiplier::get_tracked_probability() const {
    return static_cast<double>(get_counter_value()) / static_cast<double>(entry);
}

uint64_t SUnaryMultiplier::get_enabled_cycles() const {
    return enabled_cycles;
}

unsigned SUnaryMultiplier::get_width() const {
    return width;
}

uint64_t SUnaryMultiplier::get_entry() const {
    return entry;
}

std::vector<bool> sumul_stream(const std::vector<bool>& stream_0,
                              const std::vector<bool>& stream_1,
                              unsigned width,
                              unsigned sobol_dimension) {
    if (stream_0.empty() || stream_0.size() != stream_1.size()) {
        throw std::invalid_argument("SuMUL input streams must be non-empty and of identical length.");
    }

    if (width == SUnaryMultiplier::AUTO_WIDTH) {
        width = SUnaryMultiplier::width_for_stream_length(stream_0.size());
    }
    SUnaryMultiplier sumul(width, sobol_dimension);

    std::vector<bool> stream_out;
    stream_out.reserve(stream_0.size());
    for (std::size_t i = 0; i < stream_0.size(); ++i) {
        stream_out.push_back(sumul.multiply(stream_0[i], stream_1[i]));
    }
    return stream_out;
}

double sumul_probability(const std::vector<bool>& stream_0,
                        const std::vector<bool>& stream_1,
                        unsigned width,
                        unsigned sobol_dimension) {
    return calculate_probability(sumul_stream(stream_0, stream_1, width, sobol_dimension));
}

}  // namespace StochasticSimulator
