// uSADD: unipolar unary scaled adder, parallel counter into an accumulator.
// uGEMM's scaled-addition unit -- see uSADD.hpp for the datapath sketch and the rationale.
// https://jsm.ece.wisc.edu/docs/wu-toppicks2021.pdf
//
//     out = (p_0 + p_1 + ... + p_{n-1}) / n
//
// Every input bit is counted. Nothing is discarded, which is the whole difference from the MUX
// adder next door, and the reason this unit has no correlation requirement.

#include "uSADD.hpp"

#include "general_functions.hpp"

#include <stdexcept>
#include <string>

namespace StochasticSimulator {
namespace {

unsigned validated_inputs(unsigned n) {
    if (n < 1 || n > UnaryScaledAdder::MAX_INPUTS) {
        throw std::invalid_argument(
            "uSADD input count must be between 1 and " +
            std::to_string(UnaryScaledAdder::MAX_INPUTS) + ".");
    }
    return n;
}

}  // namespace

unsigned UnaryScaledAdder::width_for_inputs(unsigned n) {
    validated_inputs(n);
    // THE REGISTER HOLDS THE RESIDUE, [0, n-1], so it needs ceil(log2(n)) bits -- NOT enough to
    // hold the pre-drain sum. uGEMM takes the output from the CARRY BIT of the accumulator, so
    // the wide value (up to 2n-1) exists only inside the adder for one gate delay; what gets
    // clocked into the register is what is left after the carry leaves. For n = 2 that is a
    // ONE-BIT accumulator holding {0, 1}.
    unsigned bits = 1;
    while ((1u << bits) < n) {
        ++bits;
    }
    return bits;
}

// Members initialize in declaration order, so acc_mask can lean on acc_width.
UnaryScaledAdder::UnaryScaledAdder(unsigned input_count, unsigned width_bits)
    : inputs(validated_inputs(input_count)),
      acc_width(width_bits == AUTO_WIDTH ? width_for_inputs(input_count) : width_bits),
      acc_mask(0),
      acc(0),
      cycles(0),
      emitted_ones(0) {
    if (acc_width == 0 || acc_width > 31) {
        throw std::invalid_argument("uSADD accumulator width must be between 1 and 31 bits.");
    }
    // An accumulator too narrow to hold the residue [0, n-1] loses credit on a HEALTHY run,
    // which would be a silently wrong unit rather than an interesting fault. Refuse it.
    if (acc_width < width_for_inputs(inputs)) {
        throw std::invalid_argument(
            "uSADD accumulator needs at least " + std::to_string(width_for_inputs(inputs)) +
            " bits to hold the residue for " + std::to_string(inputs) + " inputs; a narrower "
            "one drops credit even with no faults.");
    }
    acc_mask = (acc_width >= 32) ? 0xFFFFFFFFu : ((1u << acc_width) - 1u);
    reset();
}

unsigned UnaryScaledAdder::get_pc_width() const {
    // The PC output holds 0..n, so it needs ceil(log2(n+1)) bits.
    unsigned bits = 1;
    while ((1u << bits) < inputs + 1u) {
        ++bits;
    }
    return bits;
}

bool UnaryScaledAdder::add(const std::vector<bool>& in) {
    return add(in, CycleUpset{});
}

bool UnaryScaledAdder::add(bool in_0, bool in_1) {
    if (inputs != 2) {
        throw std::invalid_argument(
            "uSADD::add(bool, bool) is the two-input convenience overload, but this unit has " +
            std::to_string(inputs) + " inputs; pass a vector instead.");
    }
    std::vector<bool> pair{in_0, in_1};
    return add(pair, CycleUpset{});
}

bool UnaryScaledAdder::add(const std::vector<bool>& in, const CycleUpset& upset) {
    if (in.size() != inputs) {
        throw std::invalid_argument(
            "uSADD received " + std::to_string(in.size()) + " input bits but has " +
            std::to_string(inputs) + " inputs.");
    }

    // ---- PC: parallel counter -----------------------------------------------------------------
    // Combinational popcount of this cycle's input bits. In hardware this is an adder tree, not
    // sequential state, so an upset here is a one-cycle bus glitch -- but see below for why its
    // CONSEQUENCE is not confined to one cycle.
    uint32_t s = 0;
    for (bool bit : in) {
        if (bit) ++s;
    }
    // The bus is only get_pc_width() wires wide, so a strike cannot flip anything above that.
    const uint32_t pc_mask = (1u << get_pc_width()) - 1u;
    s ^= (upset.pc_flip & pc_mask);

    // ---- A: accumulator, and the CARRY that is the output ---------------------------------
    // uGEMM: "the output is set to the carry bit of the accumulator... only when there are N
    // ones in the input, a logic one at the output will be generated, exactly the output
    // scaling."
    //
    // So the sum is formed at FULL width inside the adder -- up to (n-1) + n = 2n-1 -- and the
    // carry out of the "reached n" boundary becomes the output bit. Only the RESIDUE is clocked
    // back into the register, which is why the register is just ceil(log2(n)) bits wide and not
    // wide enough to hold the sum. Computing the sum in a local rather than in `acc` is what
    // makes that distinction real instead of cosmetic.
    const uint32_t sum = acc + s;

    // The carry. At most ONE bit leaves per cycle, because the output is a unary stream at one
    // bit per cycle; surplus credit stays behind and is spent later. That is exactly why the
    // unit is self-correcting -- a perturbation to A is spent down rather than compounded.
    const bool out = (sum >= inputs);

    // Masked on write, so the register WRAPS at 2^acc_width exactly as hardware would --
    // deliberately not saturated, because clamping here would be the model repairing a fault it
    // is supposed to be measuring. A healthy run never reaches the mask; only an upset can.
    acc = (out ? (sum - inputs) : sum) & acc_mask;

    ++cycles;
    if (out) ++emitted_ones;

    // ---- output wire --------------------------------------------------------------------------
    return out != upset.out_flip;
}

void UnaryScaledAdder::reset() {
    acc = 0;
    cycles = 0;
    emitted_ones = 0;
}

uint32_t UnaryScaledAdder::get_accumulator() const { return acc; }
unsigned UnaryScaledAdder::get_inputs() const { return inputs; }
unsigned UnaryScaledAdder::get_acc_width() const { return acc_width; }
uint64_t UnaryScaledAdder::get_cycles() const { return cycles; }
uint64_t UnaryScaledAdder::get_emitted_ones() const { return emitted_ones; }

void UnaryScaledAdder::flip_accumulator_bit(unsigned bit) {
    if (bit >= acc_width) {
        throw std::out_of_range("uSADD accumulator bit index is past the width of the register.");
    }
    // Deliberately NOT re-clamped into [0, n): a struck register is allowed to hold whatever the
    // upset made it hold. The unit then keeps integrating correctly from a wrong residue, which
    // is the behaviour the fault campaigns exist to measure.
    acc = (acc ^ (1u << bit)) & acc_mask;
}

void UnaryScaledAdder::set_accumulator(uint32_t value) {
    acc = value & acc_mask;
}

// ---- Free functions ---------------------------------------------------------------------------

std::vector<bool> usadd_stream(const std::vector<std::vector<bool>>& streams,
                               unsigned acc_width) {
    if (streams.empty()) {
        throw std::invalid_argument("uSADD needs at least one input stream.");
    }
    const std::size_t length = streams.front().size();
    if (length == 0) {
        throw std::invalid_argument("uSADD input streams must be non-empty.");
    }
    for (const std::vector<bool>& s : streams) {
        if (s.size() != length) {
            throw std::invalid_argument("uSADD input streams must all be the same length.");
        }
    }

    UnaryScaledAdder unit(static_cast<unsigned>(streams.size()), acc_width);

    std::vector<bool> out;
    out.reserve(length);
    std::vector<bool> cycle(streams.size());
    for (std::size_t i = 0; i < length; ++i) {
        for (std::size_t k = 0; k < streams.size(); ++k) {
            cycle[k] = streams[k][i];
        }
        out.push_back(unit.add(cycle));
    }
    return out;
}

std::vector<bool> usadd_stream(const std::vector<bool>& stream_0,
                               const std::vector<bool>& stream_1,
                               unsigned acc_width) {
    return usadd_stream(std::vector<std::vector<bool>>{stream_0, stream_1}, acc_width);
}

double usadd_probability(const std::vector<std::vector<bool>>& streams, unsigned acc_width) {
    return calculate_probability(usadd_stream(streams, acc_width));
}

double usadd_probability(const std::vector<bool>& stream_0,
                         const std::vector<bool>& stream_1,
                         unsigned acc_width) {
    return calculate_probability(usadd_stream(stream_0, stream_1, acc_width));
}

}  // namespace StochasticSimulator
