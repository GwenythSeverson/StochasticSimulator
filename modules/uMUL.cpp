// uMUL: unipolar unary multiplier with conditional bit stream generation.
// uGEMM Figure 3(a) -- see uMUL.hpp for the datapath sketch and the rationale.
// https://jsm.ece.wisc.edu/docs/wu-toppicks2021.pdf
//
// "only when input 0 is logic one, the RNG inside the bit stream generator for input 1 will
//  update, and the generator will eventually generate a new bit, i.e., input 0 is an enable
//  signal to the bit stream generator. As such, we thoroughly eliminate the correlation problem
//  and achieve high accuracy with merely an extra enable signal."
//
// in_0 is a unary bitstream; in_1 is a binary value loaded into a register, matching the
// synthesized uMUL_uni.sv (`iA`, `iB`, `loadB`, `oC <= iA & (iB_buf > sobolSeq)`).

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
            "; that is the Sobol generator's own ceiling.");
    }
    return width_bits;
}

}  // namespace

// Members initialize in declaration order, so `entry` can lean on `width`.
UnaryMultiplier::UnaryMultiplier(unsigned width_bits, unsigned sobol_dimension)
    : width(validated_width(width_bits)),
      entry(1ull << width),
      value(0),
      rng(width, sobol_dimension),
      rng_index(0),
      enabled_cycles(0) {
    reset();
}

// Delegates through SobolRNG so the StreamLength -> width mapping lives in exactly one place.
UnaryMultiplier::UnaryMultiplier(StreamLength lengthMode, unsigned sobol_dimension)
    : UnaryMultiplier(SobolRNG(lengthMode, sobol_dimension).get_width(), sobol_dimension) {}

/**
 * ceil(log2(stream_length)). There is no tradeoff to balance here -- unlike SuMUL's sliding
 * window, width costs `width` flip-flops rather than 2^width, so the only thing to get right is
 * that G's RNG must not wrap before the run ends. Overshooting is free and buys finer value
 * quantization; undershooting makes the Sobol sequence repeat mid-stream.
 */
unsigned UnaryMultiplier::width_for_stream_length(std::size_t stream_length) {
    if (stream_length == 0) {
        throw std::invalid_argument("uMUL stream length must be at least 1.");
    }

    unsigned width_bits = 1;
    while (width_bits < MAX_WIDTH &&
           (static_cast<std::size_t>(1) << width_bits) < stream_length) {
        ++width_bits;
    }
    return width_bits;
}

UnaryMultiplier UnaryMultiplier::for_stream_length(std::size_t stream_length,
                                                  unsigned sobol_dimension) {
    return UnaryMultiplier(width_for_stream_length(stream_length), sobol_dimension);
}

void UnaryMultiplier::load_value(uint64_t new_value) {
    // 2^width is a legal value, not an off-by-one: the comparator is a strict ">", so 2^width is
    // the only setting that beats every random number and therefore the only way to say p1 = 1.
    if (new_value > entry) {
        throw std::invalid_argument("uMUL value must be in [0, 2^width].");
    }
    value = new_value;
}

void UnaryMultiplier::load_probability(double probability) {
    if (!(probability >= 0.0 && probability <= 1.0)) {
        throw std::invalid_argument("uMUL probability must be between 0.0 and 1.0.");
    }
    // Truncating threshold, same as BitstreamGenerator and SobolRNG::next_bit. p = 1.0 lands
    // exactly on 2^width, which is what pins the generated stream to all ones.
    value = static_cast<uint64_t>(probability * static_cast<double>(entry));
    if (value > entry) {
        value = entry;  // guard the float edge case
    }
}

void UnaryMultiplier::load_from_stream(const std::vector<bool>& stream) {
    if (stream.empty()) {
        throw std::invalid_argument("uMUL cannot load an operand from an empty stream.");
    }
    uint64_t length = static_cast<uint64_t>(stream.size());
    if (width < 64 && length > (1ull << (64 - width))) {
        throw std::invalid_argument("uMUL operand stream is too long to scale without overflow.");
    }

    uint64_t ones = 0;
    for (bool bit : stream) {
        if (bit) ++ones;
    }
    // Integer throughout: the counter of Figure 3(a) run to completion, then scaled once. All
    // ones gives exactly `entry`, no ones gives exactly 0, and the only rounding anywhere is this
    // single truncating divide.
    value = (ones << width) / length;
}

bool UnaryMultiplier::multiply(bool in_0) {
    return multiply(in_0, CycleUpset{});
}

bool UnaryMultiplier::multiply(bool in_0, const CycleUpset& upset) {
    // ---- en: the enable wire ------------------------------------------------------------------
    // Resolved FIRST, because there is physically one wire and both the AND gate and the index
    // register hang off it. A glitch here is therefore not a one-cycle error: it also shifts G's
    // position in the Sobol sequence for every cycle that follows.
    bool enable = (in_0 != upset.in_0_flip);

    // ---- G: bit stream generator ------------------------------------------------------------
    // C's output is just the loaded register -- in_1 arrived binary, so there is nothing to
    // count. The comparator is the SNG's: truncating threshold, strict ">", the same convention
    // BitstreamGenerator and SobolRNG::next_bit use, which is what makes the rails exact.
    // The bus mask is clipped to `width` wires; a w-bit bus has no more than that to hit.
    uint32_t bus_mask = static_cast<uint32_t>(upset.rng_bus_flip & (entry - 1));
    uint32_t random_value = rng.value_at(rng_index) ^ bus_mask;
    bool generated_bit = value > static_cast<uint64_t>(random_value);

    // in_0 clocks G's RNG index forward; on an enable-low cycle the RNG holds its value, so no
    // random number is burned on a cycle the AND gate was going to zero anyway. G therefore walks
    // its sequence in order across exactly the cycles that can carry a 1, and lands p1 of them --
    // regardless of what in_0's stream looks like. That is why uMUL has no correlation
    // requirement and the plain AND-gate Multiplier does.
    // Read first, then advance: this cycle's output used the value the register already held.
    if (enable) {
        rng_index = (rng_index + 1) % entry;
        ++enabled_cycles;
    }

    // ---- ANDMUL, then the output wire ---------------------------------------------------------
    bool out = enable && generated_bit;
    return out != upset.out_flip;
}

void UnaryMultiplier::reset() {
    // Run state only. The loaded value survives, because it is an operand and not part of the
    // run -- one weight against many activation streams is the whole point of the loaded form.
    rng.reset();
    rng_index = 0;
    enabled_cycles = 0;
}

uint64_t UnaryMultiplier::get_value() const {
    return value;
}

double UnaryMultiplier::get_probability() const {
    return static_cast<double>(value) / static_cast<double>(entry);
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

void UnaryMultiplier::flip_value_bit(unsigned bit) {
    if (bit >= width) {
        throw std::out_of_range("uMUL value bit index is past the width of the register.");
    }
    // Deliberately NOT re-clamped to [0, 2^width]: a flipped register is allowed to hold whatever
    // the upset made it hold. Clamping here would be the model repairing a fault it is supposed
    // to be measuring.
    value ^= (1ull << bit);
}

void UnaryMultiplier::flip_rng_index_bit(unsigned bit) {
    if (bit >= width) {
        throw std::out_of_range("uMUL rng_index bit index is past the width of the register.");
    }
    // Stays inside [0, entry) for any bit < width, so no wrap is needed.
    rng_index ^= (1ull << bit);
}

uint64_t UnaryMultiplier::get_rng_index() const {
    return rng_index;
}

void UnaryMultiplier::set_rng_index(uint64_t index) {
    rng_index = index % entry;
}

std::vector<bool> umul_stream(const std::vector<bool>& stream_0,
                              double probability_1,
                              unsigned width,
                              unsigned sobol_dimension) {
    if (stream_0.empty()) {
        throw std::invalid_argument("uMUL input stream must be non-empty.");
    }

    if (width == UnaryMultiplier::AUTO_WIDTH) {
        width = UnaryMultiplier::width_for_stream_length(stream_0.size());
    }
    UnaryMultiplier umul(width, sobol_dimension);
    umul.load_probability(probability_1);

    std::vector<bool> stream_out;
    stream_out.reserve(stream_0.size());
    for (std::size_t i = 0; i < stream_0.size(); ++i) {
        stream_out.push_back(umul.multiply(stream_0[i]));
    }
    return stream_out;
}

std::vector<bool> umul_stream(const std::vector<bool>& stream_0,
                              const std::vector<bool>& stream_1,
                              unsigned width,
                              unsigned sobol_dimension) {
    if (stream_0.empty() || stream_1.empty()) {
        throw std::invalid_argument("uMUL input streams must be non-empty.");
    }

    if (width == UnaryMultiplier::AUTO_WIDTH) {
        width = UnaryMultiplier::width_for_stream_length(stream_0.size());
    }
    UnaryMultiplier umul(width, sobol_dimension);
    umul.load_from_stream(stream_1);

    std::vector<bool> stream_out;
    stream_out.reserve(stream_0.size());
    for (std::size_t i = 0; i < stream_0.size(); ++i) {
        stream_out.push_back(umul.multiply(stream_0[i]));
    }
    return stream_out;
}

double umul_probability(const std::vector<bool>& stream_0,
                        double probability_1,
                        unsigned width,
                        unsigned sobol_dimension) {
    return calculate_probability(umul_stream(stream_0, probability_1, width, sobol_dimension));
}

double umul_probability(const std::vector<bool>& stream_0,
                        const std::vector<bool>& stream_1,
                        unsigned width,
                        unsigned sobol_dimension) {
    return calculate_probability(umul_stream(stream_0, stream_1, width, sobol_dimension));
}

}  // namespace StochasticSimulator
