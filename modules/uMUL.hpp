/**
 * uMUL -- the unary multiplier from uGEMM, Figure 3(a) (unipolar).
 *   Wu, Li, Yin, Hsiao, Kim, San Miguel. "uGEMM: Unary Computing for GEMM Applications."
 *   https://jsm.ece.wisc.edu/docs/wu-toppicks2021.pdf
 * Unipolar only -- the decomposed-XNORMUL bipolar form in Figure 3(b) is not built here.
 *
 * ONE STREAM AND ONE BINARY NUMBER. in_0 arrives as a unary bitstream; in_1 is a binary value
 * LOADED into a register before the run. This is the form uGEMM actually synthesized --
 * hw/kernel/.../uMUL_uni.sv takes `iA` (one bit), `iB` (8 bits) and `loadB`, and computes
 * `oC <= iA & (iB_buf > sobolSeq)` -- and it is the form a GEMM wants, because the weight operand
 * is known before the layer runs.
 *
 * WHAT IT FIXES:
 *   The plain AND-gate Multiplier next door is correct only if its two input streams are
 *   uncorrelated. Feed it two streams that line their 1s up and it reports garbage --
 *   AND(s, s) = p, not p^2. The AND gate is not really the problem, though. Look at what it does:
 *   the output can only be 1 on the cycles where in_0 is 1. On every cycle where in_0 is 0 the
 *   output is forced to 0 no matter what in_1 was, so in_1's bit that cycle was generated and
 *   thrown away. Correlation error is precisely the accounting error left by those wasted bits.
 *
 *   uMUL stops wasting them. It generates in_1's stream locally from the loaded value, and only
 *   lets that generator step forward on the cycles the AND gate can actually use -- so in_0
 *   becomes an ENABLE signal [(1) in Figure 3(a)]. The generator then walks its RNG sequence in
 *   order across exactly the in_0-is-1 cycles and lands p1 of them as ones. Nothing is discarded,
 *   so nothing is miscounted, and the result is p0 * p1 no matter what in_0's stream looks like.
 *
 * THE DATAPATH (one clock cycle):
 *
 *      in_0 ---------------------------+-----------------,
 *                                      |                  )-- out      (2) ANDMUL
 *                                   (enable)          ,--'
 *                                      v              |
 *      in_1 --> [ C: value register ] ==bin==> [ G: bit stream gen ]   (1) conditional bs gen
 *               (loaded once, w bits)           comparator vs RNG
 *
 *   C   The counter position in Figure 3(a), held as a plain w-bit register. Because in_1 is
 *       supplied as a binary number there is nothing to count at run time -- the value is loaded
 *       and then just sits there. It is the thick line in Figure 3(a), and the only binary signal
 *       in the unit; in_0 and out are unary bits.
 *   G   Bit stream generator. The usual SNG comparator, "value > random value", using the same
 *       truncating-threshold, strict-">" convention as BitstreamGenerator and SobolRNG::next_bit
 *       so its streams stay comparable with the rest of the pipeline.
 *   en  in_0. The RNG index inside G advances by 1 on a cycle only if in_0 is 1.
 *
 * WHY THE LOADED FORM IS THE ACCURATE ONE:
 *   in_1 is known exactly, so there is no estimator anywhere in the unit and therefore no
 *   estimator error. Every failure mode that comes from watching a stream and guessing its
 *   density is simply absent:
 *     - no warm-up. The value is right on cycle 1.
 *     - the zero rail is exact. value = 0 never fires the strict ">", so p1 = 0 comes out as a
 *       bit-perfect zero.
 *     - no dependence on in_1's stream quality, because there is no in_1 stream.
 *     - no dependence on in_1 being stationary, well-mixed, or long enough to measure.
 *
 * THE ONE RAIL IS NOT REACHABLE, and that is the RTL's behaviour rather than an oversight.
 *   The value register is exactly `width` bits, holding [0, 2^width - 1], the same as `iB` in
 *   uMUL_uni.sv. The comparator is `value > random`, and random also spans [0, 2^width - 1], so
 *   the largest expressible density is
 *       (2^width - 1) / 2^width      -- 1023/1024 at width 10, 255/256 at the RTL's width 8
 *   and p1 = 1.0 simply does not exist as an operand. load_probability(1.0) saturates to
 *   2^width - 1 rather than pretending otherwise.
 *
 *   Representing p1 = 1.0 exactly would take a (width + 1)-bit register so the value 2^width
 *   could be stored -- one extra flip-flop, purely to buy one extra operand at the top of the
 *   range. uGEMM does not spend it, so neither does this. The consequence is worth stating
 *   plainly: multiplying by "one" gives you 1023/1024, not a pass-through.
 *   The only quantization left anywhere is the 1/2^width step of the value register, and since
 *   that costs `width` flip-flops rather than 2^width you can simply make it big. Compare
 *   SuMUL next door, whose whole SIZING section is a two-sided tradeoff forced on it by having to
 *   estimate p1 from a sliding window.
 *
 * MEASURED, N = 1024, width 10, 361-point sweep of p0,p1 over 0.05..0.95
 * (regenerate with --gtest_filter=*CharacterizeAccuracy*):
 *
 *                                        RMSE      max
 *     uMUL, in_0 from an LFSR          0.0007   0.0021
 *     uMUL, in_0 from Sobol            0.0007   0.0021
 *     AND gate, independent operands   0.0011   0.0030
 *     AND gate, SCC = 1                0.1110   0.2500   <- collapses to min(p0, p1)
 *
 *   Two things to read off that. First, uMUL beats the plain AND gate even on the AND gate's
 *   BEST case -- independent operands -- which SuMUL never managed, because SuMUL had to pay
 *   warm-up for the privilege of estimating p1. Remove the estimator and the debt goes with it.
 *   Second, the two uMUL rows are identical: in_0's generator makes no difference whatsoever,
 *   because the enable never reasons about in_0's structure, only about which cycles it passes.
 *   There is no correlated-operand row because there is no second operand stream to correlate.
 *
 * WHY SOBOL:
 *   G consumes only as many random numbers as in_0 has ones, so with p0 = 0.1 a 1024-cycle run
 *   gives G about 100 draws. A Sobol prefix of any length is already evenly spread, so those 100
 *   draws still land at p1; an LFSR prefix is a random walk and would not. uGEMM's own default is
 *   Sobol dimension 1, which is what this defaults to.
 *
 * SIZING -- `width` is not a tradeoff here:
 *   It sets two things, and both want it large: the value register's precision (steps of
 *   1/2^width) and G's RNG period 2^width, which must cover the enabled cycles or the Sobol
 *   sequence repeats mid-run. Cost is linear in width, so overshooting is free.
 *       width_for_stream_length(N) = ceil(log2(N))
 *   is the floor set by the RNG period; go wider if you want finer value quantization. MAX_WIDTH
 *   is 32 -- the Sobol generator's own ceiling -- rather than SuMUL's 16, precisely because there
 *   is no 2^width shift register to pay for.
 *
 * IF YOU ONLY HAVE in_1 AS A STREAM: load_from_stream() counts its ones once, up front, and
 *   loads the exact ratio. That is the counter of Figure 3(a) run to completion before the
 *   multiply starts, instead of incrementally during it, and it is strictly better -- the value
 *   is exact from cycle 1 rather than converging. umul_stream()'s two-stream overload does this
 *   for you. It costs one pass over in_1 and no per-cycle hardware.
 *
 * FAULT BEHAVIOUR -- every register here is PERSISTENT. There is no self-scrubbing state at all,
 *   which is the sharpest possible contrast with SuMUL, whose window cells shift an upset off the
 *   end within 2^width cycles and heal with no intervention. Here:
 *     value      -- an upset moves the output density by 2^bit / 2^width for the ENTIRE run and
 *                   is never repaired. A flip in the MSB halves or doubles the operand.
 *     rng_index  -- an upset permanently offsets G's position in the Sobol sequence.
 *   That is the price of having no redundancy: the unit is exact when healthy and has nothing to
 *   fall back on when it is not. SuMUL trades accuracy for a datapath that repairs itself. This
 *   asymmetry is the point of the bit-level accessors below.
 *
 * USAGE:
 *   1. Whole stream in one call, in_1 as a probability -- width sized to the stream:
 *        std::vector<bool> out = umul_stream(stream_a, 0.375);
 *        double z = umul_probability(stream_a, 0.375);
 *   2. in_1 as a stream you want counted exactly first:
 *        double z = umul_probability(stream_a, stream_b);
 *   3. Per-cycle, loading the operand once and reusing it -- the GEMM shape, where one weight
 *      multiplies many activations:
 *        UnaryMultiplier umul = UnaryMultiplier::for_stream_length(N);
 *        umul.load_probability(0.375);
 *        for (size_t i = 0; i < N; ++i) out.push_back(umul.multiply(stream_a[i]));
 *        umul.reset();                      // clears the run, KEEPS the loaded value
 *        for (size_t i = 0; i < N; ++i) out2.push_back(umul.multiply(other[i]));
 *
 * NOT STATELESS -- it carries the loaded value and the RNG index across cycles. reset() clears
 * the run state (RNG index, enabled-cycle count) and deliberately PRESERVES the loaded value, so
 * a loaded operand can be reused across streams the way a weight is.
 *
 * ASYMMETRY is now structural rather than a convention: in_0 is the stream and the enable, in_1
 * is the loaded number. They are not interchangeable and there is no overload that pretends
 * otherwise.
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
    // Only the value register and the RNG scale with width -- `width` flip-flops, not 2^width --
    // so this runs all the way to the Sobol generator's own ceiling.
    static constexpr unsigned MAX_WIDTH = SobolRNG::MAX_WIDTH;
    // uGEMM's own default configuration is an 8-bit width on Sobol dimension 1, matching the
    // 8-bit `iB` port of uMUL_uni.sv. Fine for streams up to 256 cycles; past that the RNG wraps.
    static constexpr unsigned DEFAULT_WIDTH = 8;
    // Sentinel for the free functions: size the width from the stream length instead.
    static constexpr unsigned AUTO_WIDTH = 0;

    /**
     * @param width            Precision of the value register and of G's comparator. The RNG
     *                         inside G has period 2^width. Bigger is simply better here.
     * @param sobol_dimension  Dimension for G's RNG. Give each uMUL in a design its own dimension
     *                         so their generated streams are not correlated with each other.
     */
    explicit UnaryMultiplier(unsigned width = DEFAULT_WIDTH, unsigned sobol_dimension = 1);

    // Interop constructor: sizes the unit to one of the existing stream-length modes.
    explicit UnaryMultiplier(StreamLength lengthMode, unsigned sobol_dimension = 1);

    // ceil(log2(stream_length)) -- just enough RNG period to cover the run without wrapping.
    static unsigned width_for_stream_length(std::size_t stream_length);

    // Same thing as a factory, mirroring SobolRNG::for_length.
    static UnaryMultiplier for_stream_length(std::size_t stream_length, unsigned sobol_dimension = 1);

    // ---- Loading in_1 -- the `loadB` port ----------------------------------------------------

    /**
     * Loads the raw register value. The register is `width` bits, so the range is
     * [0, 2^width - 1] -- exactly `iB` in uMUL_uni.sv. Throws above that; 2^width is NOT a legal
     * operand and p1 = 1.0 is therefore not representable. See the note in the file header.
     */
    void load_value(uint64_t value);

    /**
     * Same, as a probability in [0, 1]. Truncating threshold, matching BitstreamGenerator.
     * SATURATES: probability 1.0 loads 2^width - 1, the largest the register holds, which
     * generates a density of (2^width - 1)/2^width rather than a stream of solid ones.
     */
    void load_probability(double probability);

    /**
     * Counts the ones in `stream` and loads the ratio -- Figure 3(a)'s counter run to completion
     * before the multiply begins. Integer throughout, so the only rounding is the register's own
     * 1/2^width truncation. An all-ones stream saturates to 2^width - 1, same as
     * load_probability(1.0). Throws on an empty stream.
     */
    void load_from_stream(const std::vector<bool>& stream);

    // ---- Running ------------------------------------------------------------------------------

    /**
     * @brief One clock cycle of Figure 3(a).
     * @param in_0  The enable-side operand: the unary bit. Passes into the AND gate and gates G.
     * @return The output bit, in_0 AND G(value).
     */
    bool multiply(bool in_0);

    /**
     * Single-event TRANSIENTS on this cycle's combinational paths -- the wires between the
     * registers. The registers themselves are hit with flip_value_bit() / flip_rng_index_bit(),
     * whose effects persist; everything in here lasts exactly one cycle.
     *
     * Together those two mechanisms cover the unit's entire radiation cross-section:
     *
     *      SITE                       COUNT      LIFETIME     HOW TO INJECT
     *      value register             width      permanent    flip_value_bit(b)
     *      rng_index register         width      permanent    flip_rng_index_bit(b)
     *      RNG output bus             width      one cycle    CycleUpset::rng_bus_flip
     *      in_0 / enable wire         1          one cycle    CycleUpset::in_0_flip
     *      output wire                1          one cycle    CycleUpset::out_flip
     *      ----------------------------------------------------------------------
     *      total sites              3*width + 2
     *
     * The Sobol direction vectors are deliberately NOT a site: they are constants, and in
     * hardware they are ROM or tied-off wiring rather than flip-flops. Add them if your fault
     * model includes configuration memory.
     */
    struct CycleUpset {
        uint32_t rng_bus_flip = 0;      // XOR mask on G's RNG output, masked to `width` wires
        bool     in_0_flip    = false;  // glitch on the enable wire
        bool     out_flip     = false;  // glitch on the output wire
    };

    /**
     * @brief The same clock cycle, with transients applied.
     *
     * Note that in_0_flip is seen by BOTH the AND gate and the index register, because there is
     * only one enable wire. A glitch there does not merely corrupt one output bit -- it can also
     * desynchronise G's position in the Sobol sequence for the rest of the run, which is why a
     * one-cycle transient on this unit is not necessarily a one-cycle error.
     */
    bool multiply(bool in_0, const CycleUpset& upset);

    /**
     * Clears the RNG index and the enabled-cycle count. PRESERVES the loaded value, so one loaded
     * operand can be run against several in_0 streams -- the GEMM case. Use load_value(0) if you
     * actually want the operand cleared too.
     */
    void reset();

    // ---- Inspection ---------------------------------------------------------------------------

    uint64_t get_value() const;         // C's binary output: the loaded register, in [0, 2^width]
    double get_probability() const;     // the same thing as p1
    uint64_t get_enabled_cycles() const;  // how many random numbers G has consumed
    unsigned get_width() const;
    uint64_t get_entry() const;         // 2^width: the RNG period and the value register's full scale

    // ---- Bit-level state access, for fault-injection campaigns -------------------------------
    // Both registers are PERSISTENT -- nothing here self-scrubs, so no upset is ever repaired.
    // See the note in the file header; this is the opposite of SuMUL's behaviour and the reason
    // both units are worth keeping side by side.
    void flip_value_bit(unsigned bit);       // inject an upset into the operand register
    void flip_rng_index_bit(unsigned bit);   // inject an upset into the index register
    uint64_t get_rng_index() const;
    void set_rng_index(uint64_t index);

private:
    unsigned width;
    uint64_t entry;

    // Counter C, degenerate: in_1 arrives already binary, so the register is written once by
    // load_*() and read every cycle. No accumulation, no estimation, no per-cycle state change.
    uint64_t value;

    // Bit stream generator G. `rng_index` is the register the enable signal clocks; the RNG's
    // value this cycle is Sobol point number rng_index. Addressing the sequence by index rather
    // than stepping it keeps "the RNG only moves when enabled" literal in the code, and keeps the
    // unit replayable for fault injection.
    SobolRNG rng;
    uint64_t rng_index;
    uint64_t enabled_cycles;
};

/**
 * Runs a whole stream through one uMUL against a known probability, and returns the output stream.
 * Leaving `width` at AUTO_WIDTH sizes the RNG to the stream via width_for_stream_length().
 * Throws if the stream is empty or the probability is outside [0, 1].
 */
std::vector<bool> umul_stream(const std::vector<bool>& stream_0,
                              double probability_1,
                              unsigned width = UnaryMultiplier::AUTO_WIDTH,
                              unsigned sobol_dimension = 1);

/**
 * Convenience overload for when in_1 is only available as a stream: counts its ones once, loads
 * the exact ratio, then runs. Equivalent to load_from_stream() -- NOT a per-cycle estimator.
 * The two streams need not be the same length.
 */
std::vector<bool> umul_stream(const std::vector<bool>& stream_0,
                              const std::vector<bool>& stream_1,
                              unsigned width = UnaryMultiplier::AUTO_WIDTH,
                              unsigned sobol_dimension = 1);

/** The same runs, decoded straight to a probability. Expect p0 * p1. */
double umul_probability(const std::vector<bool>& stream_0,
                        double probability_1,
                        unsigned width = UnaryMultiplier::AUTO_WIDTH,
                        unsigned sobol_dimension = 1);

double umul_probability(const std::vector<bool>& stream_0,
                        const std::vector<bool>& stream_1,
                        unsigned width = UnaryMultiplier::AUTO_WIDTH,
                        unsigned sobol_dimension = 1);

}  // namespace StochasticSimulator
