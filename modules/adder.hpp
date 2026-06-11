#pragma once

#include <vector>

namespace StochasticSimulator {


class Adder {
public:
    bool add_scaled_or_weighted(bool lhs, bool rhs, bool select_bit) const; // stateless. select bit cant be assymmetric slect bit for weighted addition
};

}  // namespace StochasticSimulator
