/**
 * uSADD -- the unary SCALED ADDER from uGEMM.
 *   Wu, Li, Yin, Hsiao, Kim, San Miguel. "uGEMM: Unary Computing for GEMM Applications."
 *   https://jsm.ece.wisc.edu/docs/wu-toppicks2021.pdf
 * Unipolar only. Computes the SCALED sum of n unary streams:
 *
 *      out  =  (p_0 + p_1 + ... + p_{n-1}) / n
 *
 * -----------------------------------------------------------------------------------------
 * WHAT IT FIXES -- THE MUX ADDER THROWS AWAY n-1 BITS EVERY CYCLE
 *
 *   The standard stochastic adder is a MUX: pick one input at random each cycle and pass it
 *   through. With a 0.5 select stream, out = (a + b)/2 in expectation, and that "in expectation"
 *   is the whole problem. On every cycle the MUX routes ONE input and DISCARDS the other. The
 *   discarded bit carried real information about its operand, and throwing it away is what turns
 *   an exact quantity into a sampled one:
 *
 *       MUX      out = sel ? a : b        1 of 2 input bits used per cycle, 1 discarded
 *       uSADD    out = f(a + b)           2 of 2 used
 *
 *   The MUX's error therefore has TWO independent sources -- how the operands happen to line up
 *   with each other, and how the select stream happens to fall -- and neither shrinks by making
 *   the operands better. uSADD removes both by never discarding a bit.
 *
 * -----------------------------------------------------------------------------------------
 * THE DATAPATH (one clock cycle)
 *
 *      in_0 --,                                    (1) ACCUMULATION      (2) AVERAGE BY N
 *      in_1 --+--> [ PC: parallel counter ] --s--> (+) --> [ A ] --carry--> out
 *       ...  -'         sums the n bits             ^        |
 *      in_n-1                                       |        v
 *                                                   '----- residue, [0, n)
 *
 *   Matching uGEMM figure 3(c): block (1) is the parallel counter feeding the accumulator,
 *   block (2) is the carry-out that performs the divide-by-N.
 *
 *   PC  Parallel counter. Combinational popcount of the n input bits, so s = number of ones
 *       this cycle, in [0, n]. Needs ceil(log2(n+1)) bits. For n = 2 it is a half adder.
 *   A   Accumulator. A += s every cycle, and THE OUTPUT IS THE CARRY OUT of that addition:
 *       reaching n emits a 1 and leaves the residue behind. AT MOST ONE BIT IS EMITTED PER
 *       CYCLE, because the output is a unary stream at one bit per cycle.
 *
 *       uGEMM states this as: "the output is set to the carry bit of the accumulator... only
 *       when there are N ones in the input, a logic one at the output will be generated,
 *       exactly the output scaling."
 *
 *       THE REGISTER IS NARROWER THAN THE SUM, and that is the point of the carry formulation.
 *       A + s reaches 2n-1 at most, but that value lives only inside the adder; what is clocked
 *       back is the residue, always in [0, n). So A is ceil(log2(n)) bits -- for n = 2, ONE BIT
 *       holding {0, 1} -- not the ceil(log2(2n)) an earlier revision of this file used. That
 *       sizing error did not change any healthy result, since compare-and-subtract is the same
 *       arithmetic as a carry, but it doubled the accumulator's apparent fault cross-section.
 *
 * WHY THAT IS EXACTLY THE SCALED SUM. A is a running credit counter, so over N cycles
 *
 *       output ones  =  floor( (total input ones) / n )
 *
 * and dividing by N gives (SUM p_i)/n up to the final residue, which is < n ones out of N. The
 * only error the healthy unit can make is that leftover residue: at most n-1 output ones, i.e.
 * (n-1)/N, and it is a DEFICIT never a surplus. There is no estimator and no RNG anywhere.
 *
 * -----------------------------------------------------------------------------------------
 * THE COLLAPSE -- STRONGER THAN uMUL'S, AND IT IS WHY THIS UNIT IS EXHAUSTIVELY TESTABLE
 *
 *   floor(total ones / n) does not mention WHERE any one sits. Two consequences:
 *
 *     1. The output ones-count depends ONLY on the total ones across all inputs. Every one of
 *        the C(256,128)^2 = 3.3e151 arrangements of two 0.5 streams gives the SAME answer.
 *        There is no correlation requirement, so unlike the MUX there is nothing to filter for.
 *     2. A fault anywhere in the input streams matters only through its NET effect on that
 *        total: -1 per flipped one, +1 per flipped zero. The whole 768-bit fault space of a
 *        two-input adder collapses onto a single integer.
 *
 *   Compare the MUX, whose output is |{i : sel_i=1, a_i=1}| + |{i : sel_i=0, b_i=1}| -- a
 *   position-dependent quantity that needs all eight (a,b,sel) region sizes to describe.
 *
 * -----------------------------------------------------------------------------------------
 * FAULT BEHAVIOUR -- THE ACCUMULATOR IS A RESIDUE, NOT AN OPERAND
 *
 *   This is the sharpest contrast with uMUL, and the reason both are worth measuring.
 *
 *     uMUL's value register HOLDS AN OPERAND. An upset there changes the number being
 *       multiplied for the entire run: unbounded, permanent, and a flip of the MSB takes the
 *       answer to zero. Nothing repairs it.
 *     uSADD's accumulator HOLDS A RESIDUE -- credit not yet spent, always in [0, n). An upset
 *       adds or removes credit ONCE. The unit keeps integrating correctly afterwards, so the
 *       total output moves by about (perturbation / n) bits and then the error stops growing.
 *       The damage is BOUNDED BY THE SIZE OF THE REGISTER, not by the length of the run.
 *
 *   PC upsets behave like accumulator upsets, because A integrates whatever PC hands it: a
 *   one-cycle bus glitch of +2^b feeds 2^b spurious credit into the running sum, exactly as if
 *   2^b extra input ones had arrived. Also bounded, also one-time.
 *
 *   SITE INVENTORY, per cycle:
 *       SITE                     COUNT              LIFETIME     HOW TO INJECT
 *       input streams            n                  one cycle    flip the caller's bits
 *       PC output bus            ceil(log2(n+1))    one cycle    CycleUpset::pc_flip
 *       A register               acc_width          persistent   flip_accumulator_bit(b)
 *       output wire              1                  one cycle    CycleUpset::out_flip
 *
 * -----------------------------------------------------------------------------------------
 * SIZING. acc_width must hold the RESIDUE, [0, n-1], so ceil(log2(n)) bits; AUTO_WIDTH picks
 * that. The pre-drain sum is wider but never reaches the register -- the carry takes it. A wider
 * accumulator buys nothing when healthy; it only ENLARGES THE FAULT CROSS-SECTION, since
 * arithmetic WRAPS at 2^acc_width rather than saturating -- see the note on add().
 *
 * USAGE:
 *   1. Whole streams in one call:
 *        std::vector<bool> out = usadd_stream({stream_a, stream_b});
 *        double z = usadd_probability({stream_a, stream_b});      // expect (pa + pb)/2
 *   2. Per cycle:
 *        UnaryScaledAdder add(2);
 *        for (size_t i = 0; i < N; ++i) out.push_back(add.add(a[i], b[i]));
 *        add.reset();                                             // clears A and the counters
 *
 * NOT STATELESS -- it carries the accumulator across cycles. That is the entire unit.
 */

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>

namespace StochasticSimulator {

class UnaryScaledAdder {
public:
    // Two is the interesting case and the one the poster measures, but the unit is general.
    static constexpr unsigned MAX_INPUTS = 32;
    static constexpr unsigned DEFAULT_INPUTS = 2;
    // Sentinel: size the accumulator from the input count via width_for_inputs().
    static constexpr unsigned AUTO_WIDTH = 0;

    /**
     * @param inputs     How many unary streams are summed. The scale factor is 1/inputs.
     * @param acc_width  Bits in the accumulator A. AUTO_WIDTH = ceil(log2(inputs)), the
     *                   smallest that holds the residue the carry leaves behind.
     */
    explicit UnaryScaledAdder(unsigned inputs = DEFAULT_INPUTS,
                              unsigned acc_width = AUTO_WIDTH);

    // ceil(log2(n)) -- just enough to hold the residue [0, n-1] that the carry leaves behind.
    static unsigned width_for_inputs(unsigned inputs);

    // ---- Running ------------------------------------------------------------------------------

    /**
     * @brief One clock cycle.
     * @param in  The n input bits this cycle. Size must equal get_inputs().
     * @return The output bit: 1 if the accumulator reached n, else 0.
     */
    bool add(const std::vector<bool>& in);

    // Convenience for the two-input case, which is what a scaled adder usually is.
    bool add(bool in_0, bool in_1);

    /**
     * Single-event TRANSIENTS on this cycle's combinational paths. The accumulator itself is hit
     * with flip_accumulator_bit(), whose effect persists; everything in here lasts one cycle.
     *
     * Note that a pc_flip is NOT a one-cycle error in the output. A integrates whatever PC hands
     * it, so spurious credit stays in the running sum and shifts every later emission. The total
     * output ones moves by about (flip magnitude / n) -- bounded, but not confined to one bit.
     */
    struct CycleUpset {
        uint32_t pc_flip  = 0;      // XOR mask on the parallel counter's output bus
        bool     out_flip = false;  // glitch on the output wire
    };

    bool add(const std::vector<bool>& in, const CycleUpset& upset);

    // Clears A and the run counters. There is no loaded operand to preserve.
    void reset();

    // ---- Inspection ---------------------------------------------------------------------------

    uint32_t get_accumulator() const;     // A's current contents, the unspent credit
    unsigned get_inputs() const;          // n
    unsigned get_acc_width() const;       // bits in A
    unsigned get_pc_width() const;        // bits on the PC output bus, ceil(log2(n+1))
    uint64_t get_cycles() const;          // cycles clocked since reset
    uint64_t get_emitted_ones() const;    // output ones emitted since reset

    // ---- Bit-level state access, for fault-injection campaigns -------------------------------
    // A is the unit's ONLY sequential state. Flipping a bit here changes the unspent credit, so
    // the unit continues integrating correctly from a wrong starting point -- a bounded, one-time
    // error rather than uMUL's permanent operand corruption. That contrast is the point.
    void flip_accumulator_bit(unsigned bit);
    void set_accumulator(uint32_t value);

private:
    unsigned inputs;
    unsigned acc_width;
    uint32_t acc_mask;    // 2^acc_width - 1, applied after every write so A wraps like hardware
    uint32_t acc;
    uint64_t cycles;
    uint64_t emitted_ones;
};

/**
 * Runs whole streams through one uSADD and returns the output stream. All inputs must be the
 * same length and non-empty. Expect a density of (SUM p_i)/n.
 */
std::vector<bool> usadd_stream(const std::vector<std::vector<bool>>& streams,
                               unsigned acc_width = UnaryScaledAdder::AUTO_WIDTH);

/** Two-stream convenience overload. */
std::vector<bool> usadd_stream(const std::vector<bool>& stream_0,
                               const std::vector<bool>& stream_1,
                               unsigned acc_width = UnaryScaledAdder::AUTO_WIDTH);

/** The same runs, decoded straight to a probability. */
double usadd_probability(const std::vector<std::vector<bool>>& streams,
                         unsigned acc_width = UnaryScaledAdder::AUTO_WIDTH);

double usadd_probability(const std::vector<bool>& stream_0,
                         const std::vector<bool>& stream_1,
                         unsigned acc_width = UnaryScaledAdder::AUTO_WIDTH);

}  // namespace StochasticSimulator
