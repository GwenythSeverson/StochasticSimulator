#pragma once

#include <vector>

namespace StochasticSimulator {

class Adder {
public:
// Takes the two data streams AND the select stream (scaling factor). more adder types to be expanded on
    std::vector<int> add(const std::vector<int>& lhs, 
                         const std::vector<int>& rhs, 
                         const std::vector<int>& select) const;
};

}  // namespace StochasticSimulator
