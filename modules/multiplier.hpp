#pragma once

#include <vector>

namespace StochasticSimulator {

  //COMPLETELY STATELESS. NOT LIKE ADDERS AND DIVIDERS.
class Multiplier {
public:
    // Default constructor
    Multiplier() = default;

    /**
     * @brief Performs stochastic multiplication on two single bits.
     *  standard hardware AND gate.
     * * @param lhs The bit from the left-hand side stream.
     * @param rhs The bit from the right-hand side stream.
     * @return true if both bits are 1, false otherwise.
     */
    bool multiply(bool lhs, bool rhs) const;
};

}  // namespace StochasticSimulator
