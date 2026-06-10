#include "multiplier.hpp"
#include <stdexcept>
#include <vector>
#include "general_functions.hpp"

namespace StochasticSimulator {

double Multiplier::multiply(const std::vector<bool>& lhs, const std::vector<bool>& rhs) const {
    const size_t len = std::min(lhs.size(), rhs.size());
    if (len == 0){
        throw std::invalid_argument("Input bitstreams must not be empty.");
    } 
    
    int count = 0;
    for (size_t i = 0; i < len; ++i)
        count += (lhs[i] & rhs[i]) ? 1 : 0;
    
    return static_cast<double>(count) / static_cast<double>(len);
}

}  // namespace StochasticSimulator
