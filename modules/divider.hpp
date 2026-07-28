/**
 * Gaines division (GDIV), unipolar -- Figure 2(a) of
 *   https://jsm.ece.wisc.edu/docs/wu-ieeedt2021.pdf
 * Computes Z = X / Y with an up/down counter in a feedback loop (Gaines' ADDIE).
 *
 * THE LOOP:
 *
 *     X ------------------> [ +/- ]
 *                              |  up/down counter, depth D, saturating at 0 and D
 *                              v
 *                          [ counter ] ==bin==> [ comparator vs RNG ] --+--> Z
 *                              ^                                        |
 *                              |                                        |
 *                          [ AND ] <----- Y                             |
 *                              ^------------- registered Z (1 cycle) ---+
 *
 *   Per cycle: mul = Z(t-1) AND Y(t). The counter counts UP when X is 1 and mul is 0, DOWN when
 *   X is 0 and mul is 1, and holds when they agree. The counter's binary value drives an SNG
 *   comparator that emits Z.
 *
 * WHY IT LANDS ON X/Y:
 *   The counter stops drifting when up-events and down-events balance:
 *       P(X=1, mul=0) = P(X=0, mul=1)
 *       px(1 - pm)    = (1 - px)pm      ->    px = pm
 *   and pm = P(Z=1 AND Y=1) = pz * py when Z and Y are uncorrelated, so pz = px / py.
 *   That "uncorrelated" is a real precondition: Z is regenerated from an independent LFSR, but
 *   do not hand this function a Y that shares an RNG with the Z it will produce.
 *
 * DOMAIN: unipolar, so Z is clamped to [0, 1] by the counter's own saturation. X must be <= Y or
 *   the true quotient exceeds 1 and the counter simply pins at D -- 0.8 / 0.4 returns ~0.995,
 *   not 2.0. That is the design behaving correctly, not an error, so it is not reported as one.
 *
 * COUNTER DEPTH is the one knob, and it trades the same two errors uMUL's width does:
 *     - Settling.  The counter starts at mid-scale and has to walk to its equilibrium. A deeper
 *                  counter walks further, so it needs a longer stream before it arrives.
 *     - Precision. The counter resolves Z only to 1/D, and dithers around the true value by
 *                  about that much. A shallow counter is coarse and noisy.
 *   Measured over a px,py sweep at 3 seeds, the best depth tracks sqrt(N) closely:
 *        N =   512 -> 16      N = 16384 -> 128
 *        N =  1024 -> 32      N = 50000 -> 128..256 (tied)
 *        N =  4096 -> 64
 *   gaines_counter_depth_for_length() implements that, and is the default. The old hard-coded
 *   depth of 32 was right for N ~ 1024 and left ~4x accuracy on the table at N = 50000.
 */

#pragma once

#include "general_functions.hpp"
#include <cstdint>
#include <cstddef>
#include <vector>

namespace StochasticSimulator {

// Pass as `counter_depth` to size the counter from the stream length instead of fixing it.
constexpr unsigned GAINES_AUTO_DEPTH = 0;

/**
 * Counter depth for a given stream length: a power of two near sqrt(N), clamped to [16, 1024].
 */
unsigned gaines_counter_depth_for_length(std::size_t stream_length);

/**
 * @brief The divider's actual output BITSTREAM, one bit per input cycle.
 * Use this when the quotient has to feed another stochastic unit, or when you want to measure
 * the stream itself (ZCE, early termination, fault injection) rather than just its mean.
 *
 * @param stream_X      Numerator bitstream.
 * @param stream_Y      Denominator bitstream. Should satisfy p(X) <= p(Y); see DOMAIN above.
 * @param sng_seed      Seed for the internal comparator RNG. Fixed value -> reproducible output.
 * @param counter_depth Up/down counter depth, or GAINES_AUTO_DEPTH to size it to the stream.
 * @return The Z bitstream, or an empty vector if the inputs are malformed.
 */
std::vector<bool> gaines_division_stream(const std::vector<bool>& stream_X,
                                         const std::vector<bool>& stream_Y,
                                         uint32_t sng_seed = 42,
                                         unsigned counter_depth = GAINES_AUTO_DEPTH);

/**
 * @brief The same division, decoded straight to a probability.
 * @return double The decoded output probability, or 0.0 if the inputs are malformed.
 *
 * NOTE ON THE ERROR CONTRACT: 0.0 is also a perfectly valid quotient (X = 0), so a caller cannot
 * tell a failure from a result. Every other module in this project throws instead. This one keeps
 * returning 0.0 only because the existing unit test pins that behaviour down; prefer
 * gaines_division_stream() and check for an empty vector if you need to distinguish the two.
 */
double ud_counter_division(const std::vector<bool>& stream_X,
                           const std::vector<bool>& stream_Y,
                           uint32_t sng_seed = 42,
                           unsigned counter_depth = GAINES_AUTO_DEPTH);

}  // namespace StochasticSimulator
