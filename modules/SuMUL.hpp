/**
 * uMUL -- the unary multiplier from uGEMM, Figure 3(a) (unipolar).
 *   Wu, Li, Yin, Hsiao, Kim, San Miguel. "uGEMM: Unary Computing for GEMM Applications."
 *   https://jsm.ece.wisc.edu/docs/wu-toppicks2021.pdf
 * Unipolar only -- the bipolar XNOR decomposition in Figure 3(b) is deliberately not built here.
 *
 * WHAT IT FIXES:
 *   The plain AND-gate Multiplier next door is correct only if its two input streams are
 *   uncorrelated. Feed it two streams that happen to line their 1s up and it reports garbage --
 *   AND(s, s) = p, not p^2. The AND gate is not really the problem, though. Look at what it does:
 *   the output can only be 1 on the cycles where in_0 is 1. On every cycle where in_0 is 0 the
 *   output is forced to 0 no matter what in_1 was, so in_1's bit that cycle was generated and
 *   thrown away. Correlation error is precisely the accounting error left by those wasted bits.
 *
 *   uMUL stops wasting them. It does not take in_1's bitstream straight into the gate; it
 *   regenerates in_1 locally, and only lets that generator step forward on the cycles the AND
 *   gate can actually use -- so in_0 becomes an ENABLE signal [(1) in Figure 3(a)]. The generator
 *   then walks its RNG sequence in order across exactly the in_0-is-1 cycles, and lays down p1 of
 *   them as ones. Nothing is discarded, so nothing is miscounted, and the result is p0 * p1 no
 *   matter how the two input streams are aligned.
 *
 * THE DATAPATH (one clock cycle):
 *
 *      in_0 ---------------------------+-----------------,
 *                                      |                  )-- out      (2) ANDMUL
 *                                   (enable)          ,--'
 *                                      v              |
 *      in_1 --> [ C: counter ] ==bin==> [ G: bit stream gen ]           (1) conditional bs gen
 *                                        comparator vs RNG
 *
 *   C   Counter. A 2^width-deep shift register holding the most recent 2^width bits of in_1;
 *       its output is the NUMBER OF ONES in that window. That is a plain binary integer in
 *       [0, 2^width] -- the thick line in Figure 3(a) -- not a unary bit. It is the only
 *       binary signal in the whole unit; in_0, in_1 and out are all unary bits.
 *   G   Bit stream generator. The usual SNG comparator, "binary value > random value", using the
 *       same truncating-threshold, strict-">" convention as BitstreamGenerator and
 *       SobolRNG::next_bit so its streams stay comparable with the rest of the pipeline.
 *   en  in_0. The RNG index inside G advances by 1 on a cycle only if in_0 is 1. C keeps
 *       shifting every cycle regardless -- it is tracking in_1's value, not producing output.
 *
 * WHY SOBOL:
 *   G consumes only as many random numbers as in_0 has ones, so with p0 = 0.1 a 1024-cycle run
 *   gives G about 100 draws. A Sobol prefix of any length is already evenly spread, so those 100
 *   draws still land at p1; an LFSR prefix is a random walk and would not. uGEMM's own default is
 *   Sobol dimension 1, which is what this defaults to.
 *
 * SIZING (`width` is the one knob, and it is a real tradeoff -- measured, not guessed):
 *   The counter window is 2^width deep, and that single number sets both error terms:
 *     - Warm-up.   The window boots half full (p1 = 0.5) and cannot be right until it has seen
 *                  2^width bits of in_1. That contaminates the first 2^width of N cycles, so the
 *                  error it contributes grows with 2^width / N. Set width = log2(N) and the
 *                  window never fills at all -- measured RMSE 0.073, i.e. the unit is useless.
 *     - Precision. C can only express p1 in steps of 1 / 2^width, so a narrow window quantizes.
 *                  This error shrinks with 2^width.
 *   Balancing the two puts the best window near 2 * sqrt(N), which is what
 *   width_for_stream_length() returns, and what the free functions below use unless told
 *   otherwise. What the rule picks against what an exhaustive width sweep found to be best:
 *        N = 256  -> 5 (best 5)      N = 1024 -> 6 (best 6)
 *        N = 512  -> 5 (best 6)      N = 4096 -> 7 (best 7)
 *   Only N = 512 is off, and by one step -- 0.0073 RMSE against the best 0.0061.
 *   Rule of thumb if you are picking by hand: keep N / 2^width somewhere around 8 to 32.
 *
 * MEASURED, N = 1024, auto width, sweeping p0,p1 over 0.05..0.95:
 *     operands sharing an RNG (SCC = 1)   uMUL RMSE 0.0044 max 0.012
 *                                          AND RMSE 0.1113 max 0.250   <- collapses to min(p0,p1)
 *     independent operands                uMUL RMSE 0.0061 max 0.016
 *                                          AND RMSE 0.0010 max 0.003
 *   That is the trade in one table. Against correlated inputs uMUL is ~25x better and the AND
 *   gate is simply wrong; against genuinely independent inputs the AND gate wins, because it
 *   reads in_1 directly and pays no warm-up while uMUL has to learn p1 first.
 *
 * EXACT 0 AND EXACT 1 are the warm-up's worst case, for the same reason: the window boots at 0.5,
 *   so a few cycles' worth of the wrong bit gets out before it fills. At N = 1024 that is
 *   1.0 * 0.0 -> 0.015 and 1.0 * 1.0 -> 0.984 rather than a clean 0 and 1. The plain AND gate is
 *   exact at both. Widen N or narrow width to shrink it; it does not go to zero.
 *
 * WHAT IT DOES NOT FIX:
 *   C is a SLIDING window, so it tracks in_1's *local* density. That is exactly right for the
 *   well-mixed streams an SNG produces, and wrong for a non-stationary stream -- feed it in_1
 *   with all its ones bunched at the front and C will faithfully report ~1.0 during the bunch.
 *   uMUL removes the correlation requirement, not the requirement that a stream be well mixed.
 *
 * USAGE:
 *   1. Whole stream in one call -- width sized to the stream automatically:
 *        std::vector<bool> out = umul_stream(stream_a, stream_b);
 *        double z = umul_probability(stream_a, stream_b);
 *   2. Per-cycle, mirroring the AND-gate Multiplier's loop shape (but stateful -- see below):
 *        UnaryMultiplier umul = UnaryMultiplier::for_stream_length(N);
 *        for (size_t i = 0; i < N; ++i) out.push_back(umul.multiply(stream_a[i], stream_b[i]));
 *   3. Explicit width when you want to control the tradeoff yourself:
 *        UnaryMultiplier umul(6);
 *
 * NOT STATELESS -- unlike Multiplier, and like Adder is not either in spirit. It carries the
 * counter window and the RNG index across cycles, so one instance belongs to one stream pair.
 * Reuse it for a second pair only after reset().
 *
 * ASYMMETRY: multiply(a, b) and multiply(b, a) are both correct but are not the same circuit.
 * in_0 is the enable and passes through untouched; in_1 gets regenerated. Put the operand you
 * trust more, or the sparser one, on in_0.
 */

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>

#include "bsg/lfsr.hpp"   // StreamLength
#include "bsg/sobol.hpp"  // the RNG inside G

namespace StochasticSimulator {

class UnaryMultiplier {
public:
    // The counter is 2^width physical flip-flops, so this stops being hardware past ~16.
    static constexpr unsigned MAX_WIDTH = 16;
    // uGEMM's own default configuration is an 8-bit width on Sobol dimension 1. Good when the
    // stream is comfortably longer than the 256-deep window; otherwise prefer the auto-sizing.
    static constexpr unsigned DEFAULT_WIDTH = 8;
    // Sentinel for the free functions: size the width from the stream length instead.
    static constexpr unsigned AUTO_WIDTH = 0;

    /**
     * @param width            Precision of the unit. The counter window is 2^width bits deep and
     *                         the RNG inside G has period 2^width. See SIZING above.
     * @param sobol_dimension  Dimension for G's RNG. Give each uMUL in a design its own dimension
     *                         so their regenerated streams are not correlated with each other.
     */
    explicit UnaryMultiplier(unsigned width = DEFAULT_WIDTH, unsigned sobol_dimension = 1);

    // Interop constructor: sizes the unit to one of the existing stream-length modes. Note this
    // sets width = log2(length), which is the degenerate case SIZING warns about -- it is here
    // for symmetry with SobolRNG, not because it is the right width for a stream that long.
    explicit UnaryMultiplier(StreamLength lengthMode, unsigned sobol_dimension = 1);

    // Balances warm-up against counter precision: the window lands near 2 * sqrt(stream_length).
    static unsigned width_for_stream_length(std::size_t stream_length);

    // Same thing as a factory, mirroring SobolRNG::for_length.
    static UnaryMultiplier for_stream_length(std::size_t stream_length, unsigned sobol_dimension = 1);

    /**
     * @brief One clock cycle of Figure 3(a).
     * @param in_0  The enable-side operand. Passes into the AND gate and gates G's RNG.
     * @param in_1  The regenerated-side operand. Shifts into the counter C.
     * @return The output bit, in_0 AND G(counter).
     */
    bool multiply(bool in_0, bool in_1);

    // Returns the unit to its power-on state: counter window half full, RNG index at 0.
    void reset();

    // C's binary output this cycle: ones currently in the window, in [0, 2^width].
    uint32_t get_counter_value() const;

    // The same thing as a probability -- what the unit currently believes p1 to be.
    double get_tracked_probability() const;

    // How many times the enable has fired, i.e. how many random numbers G has consumed.
    uint64_t get_enabled_cycles() const;

    unsigned get_width() const;
    uint64_t get_entry() const;  // 2^width: window depth and RNG period

private:
    unsigned width;
    uint64_t entry;

    // Counter C. A circular buffer standing in for the shift register: `oldest` is the slot that
    // falls off next, and `ones_count` is maintained incrementally so the read stays O(1).
    std::vector<bool> window;
    std::size_t oldest;
    uint32_t ones_count;

    // Bit stream generator G. `rng_index` is the register the enable signal clocks; the RNG's
    // value this cycle is Sobol point number rng_index. Addressing the sequence by index rather
    // than stepping it keeps "the RNG only moves when enabled" literal in the code, and keeps the
    // unit replayable for fault injection.
    SobolRNG rng;
    uint64_t rng_index;
    uint64_t enabled_cycles;
};

/**
 * Runs a whole stream pair through one uMUL and returns the output stream.
 * Leaving `width` at AUTO_WIDTH sizes the counter to the stream via width_for_stream_length().
 * Throws if the streams are empty or of unequal length.
 */
std::vector<bool> umul_stream(const std::vector<bool>& stream_0,
                              const std::vector<bool>& stream_1,
                              unsigned width = UnaryMultiplier::AUTO_WIDTH,
                              unsigned sobol_dimension = 1);

/** The same run, decoded straight to a probability. Expect p0 * p1. */
double umul_probability(const std::vector<bool>& stream_0,
                        const std::vector<bool>& stream_1,
                        unsigned width = UnaryMultiplier::AUTO_WIDTH,
                        unsigned sobol_dimension = 1);

}  // namespace StochasticSimulator
