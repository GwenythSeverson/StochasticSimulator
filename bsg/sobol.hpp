/**
 * Sobol low-discrepancy sequence generator for stochastic number generation.
 *
 * Where the LFSR gives you a *pseudo-random* permutation of [0, 2^n), Sobol gives you a
 * *progressively stratified* one. That is what makes it work "for any length" --
 * truncating a Sobol-driven bitstream early still leaves you with a near-optimal estimate of p,
 * instead of the random-walk error an LFSR leaves behind.
 * 
 * * (best fast ex[planation) s an evenly-spread Sobol sequence of any requested length, 
 * and thresholds it into a unipolar stochastic bitstream — exact at power-of-two lengths,
 *  within about one bit otherwise.
 *
 * SIZING (this is the whole "any length" mechanism):
 *   The generator has a period of 2^width. You do not make a fixed-width generator arbitrarily
 *   long -- you size the width to the precision you need:
 *      - Pick your precision  -> that is 2^width.
 *      - Build SobolRNG(width) (or SobolRNG::for_length(N), which does the ceil(log2 N) for you).
 *      - Run up to 2^width generator-advancing steps for an exact result; any shorter run is a
 *        progressively-refined approximation that is already evenly spread.
 *      - Need more precision later? Bump width. Nothing else changes.
 *   Past 2^width steps the sequence simply repeats (the stream is periodic, same as the LFSR).
 *
 * CORRELATION:
 *   Two operands driven by the SAME dimension are perfectly correlated and will wreck any
 *   multiplier / ZCE measurement. Give every independent operand its own dimension:
 *      SobolRNG rng_a(10, 1);   // dimension 1
 *      SobolRNG rng_b(10, 2);   // dimension 2  -> low cross-correlation
 *
 *   Narrow generators cannot tell every dimension apart -- below a certain width the polynomial
 *   recurrence never runs and higher rows collapse onto lower ones. Asking for such a dimension
 *   throws rather than silently returning a duplicate stream; see minimum_width_for_dimension().
 *   Every dimension is distinct from width 7 up, so all StreamLength modes are safe.
 *
 * USAGE:
 *   1. Drop-in for FlexibleLFSR (same next() / get_max_cycles() / reset() shape):
 *        SobolRNG rng(StreamLength::Length_1024, 1);
 *        uint16_t r = static_cast<uint16_t>(rng.next());
 *   2. Self-contained SNG (no separate BitstreamGenerator needed, any width):
 *        SobolRNG rng(12, 1);
 *        bool bit = rng.next_bit(0.375);
 *   3. One-shot stream of arbitrary length:
 *        std::vector<bool> s = generate_sobol_stream(0.375, 1000, 1);
 */

#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include "lfsr.hpp"  // StreamLength, for interop with the existing LFSR-based pipeline

namespace StochasticSimulator {

class SobolRNG {
public:
    // 32-bit direction vectors; period 2^32 is already far past any practical stream length.
    static constexpr unsigned MAX_WIDTH = 32;
    // Number of tabulated dimensions. Extend at runtime with set_direction_numbers().
    static constexpr unsigned MAX_DIMENSION = 19;

    /**
     * @param width      Bit precision. Period (and exact-result stream length) is 2^width.
     * @param dimension  1-based Sobol dimension. Use a distinct dimension per independent operand.
     */
    explicit SobolRNG(unsigned width, unsigned dimension = 1);

    // Interop constructor: matches the widths the existing FlexibleLFSR stream modes use.
    explicit SobolRNG(StreamLength lengthMode, unsigned dimension = 1);

    // Sizes the generator to hold `stream_length` distinct points (rounds up to a power of two).
    static SobolRNG for_length(uint64_t stream_length, unsigned dimension = 1);
    static unsigned width_for_length(uint64_t stream_length);

    /**
     * Narrow generators cannot tell every dimension apart: below a certain width the polynomial
     * recurrence never runs, so higher tabulated rows collapse onto lower ones and would hand you
     * two perfectly correlated streams. Constructing such a pair throws; this reports the smallest
     * width at which `dimension` is genuinely distinct. Dimension 1 is always valid (returns 1).
     * At width >= 7 every tabulated dimension is distinct, so the StreamLength modes are all safe.
     */
    static unsigned minimum_width_for_dimension(unsigned dimension);

    // Next point of the sequence, an integer in [0, 2^width). Wraps to the start after 2^width calls.
    uint32_t next();

    // Same point mapped to [0.0, 1.0).
    double next_uniform();

    /**
     * Full SNG in one call: advances the sequence and applies the hardware comparator.
     * Uses the same truncating threshold + strict ">" comparison as BitstreamGenerator so
     * Sobol and LFSR streams stay directly comparable.
     */
    bool next_bit(double probability);

    // Random access without stepping (Gray-code evaluation) -- handy for fault injection replay.
    uint32_t value_at(uint64_t index) const;
    void jump_to(uint64_t index);

    uint64_t get_max_cycles() const;  // 2^width -- the period, i.e. the exact-result length
    unsigned get_width() const;
    unsigned get_dimension() const;
    uint64_t get_position() const;    // how many points have been consumed since reset
    void reset();

    /**
     * Escape hatch for a hand-supplied direction-number row (e.g. a line from the official
     * Joe-Kuo table) when you need a dimension past MAX_DIMENSION or a specific tuned row.
     * @param s       Degree of the primitive polynomial.
     * @param a       Its inner coefficients a_1..a_{s-1} packed with a_1 as the MSB.
     * @param m_init  s initial direction numbers; each must be odd with m_i < 2^i.
     */
    void set_direction_numbers(unsigned s, uint32_t a, const std::vector<uint32_t>& m_init);

private:
    unsigned width;
    unsigned dimension;
    uint64_t period;
    uint64_t counter;
    uint32_t state;
    std::vector<uint32_t> direction;  // direction[i] = v_{i+1}

    void build_direction_vectors(unsigned s, uint32_t a, const std::vector<uint32_t>& m_init);
    void load_tabulated_dimension(unsigned dim);
};

/**
 * One-shot bitstream of any length, width auto-sized to the length.
 * Give each operand its own `dimension` so the streams are not correlated.
 */
std::vector<bool> generate_sobol_stream(double probability, std::size_t length, unsigned dimension = 1);

}  // namespace StochasticSimulator
