#pragma once

#include <vector>

namespace StochasticSimulator {

class Multiplier {
public:
  // todo: consider making a struct for return data types like ideal, observed, error, and zce?

    /**
     * Computes stochastic multiplication by passing two bitstreams through
     * an AND gate and returning the estimated output probability.
     */
    double multiply(const std::vector<bool>& lhs, const std::vector<bool>& rhs) const;
};

}  // namespace StochasticSimulator
