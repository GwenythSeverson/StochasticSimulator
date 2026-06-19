#pragma once

#include <vector>

namespace StochasticSimulator {

class Divider {
public:
    // Pure logic execution: route numerator bits through the denominator stream. more dividers to be expanded on
    std::vector<bool> divide(const std::vector<bool>& numerator, const std::vector<bool>& denominator) const;
};

}  // namespace StochasticSimulator
