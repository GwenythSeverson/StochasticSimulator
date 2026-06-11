#include "adder.hpp"

#include "general_functions.hpp"

namespace StochasticSimulator {

// Pure, stateless adder mirroring a physical Multiplexer (MUX)
bool Adder::add_scaled_or_weighted(bool lhs, bool rhs, bool select_bit) const {
    // A standard MUX acts as an adder in stochastic computing.
    // If select_bit has a 50% probability, output value is (lhs + rhs) / 2
    // Basically, weighted addition abilities. 
    return select_bit ? lhs : rhs;
}
}  // namespace StochasticSimulator
