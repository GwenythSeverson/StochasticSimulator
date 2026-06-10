#pragma once

#include <vector>

namespace StochasticSimulator {

class Divider {
public:
    // Pure logic execution: route numerator bits through the denominator stream. more dividers to be expanded on
    std::vector<int> divide(const std::vector<int>& numerator, const std::vector<int>& denominator) const;
};

}  // namespace StochasticSimulator
