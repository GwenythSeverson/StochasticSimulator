#include "multiplier.hpp"
#include <stdexcept>
#include <vector>
#include "general_functions.hpp"

namespace StochasticSimulator {

class Multiplier {
public:
    // A pure, stateless hardware simulation step
    bool multiply(bool lhs, bool rhs) const {
        return lhs && rhs;
    }
};

}  // namespace StochasticSimulator
