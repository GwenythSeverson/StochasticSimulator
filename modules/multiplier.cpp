#include "multiplier.hpp"
#include <stdexcept>
#include <vector>
#include "general_functions.hpp"

namespace StochasticSimulator {

    //be careful when calling, to make sure that each stream is pulled from 
    // the same clock. (dont pull lhs(1) and rhs(2)) 
class Multiplier {
public:
    // A pure, stateless hardware simulation step
    bool multiply(bool lhs, bool rhs) const {
        return lhs && rhs;
    }
};

}  // namespace StochasticSimulator
