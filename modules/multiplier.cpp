#pragma once
#include "multiplier.hpp"


namespace StochasticSimulator {

    //be careful when calling, to make sure that each stream is pulled from 
    // the same clock. (dont pull lhs(1) and rhs(2)) 
    bool Multiplier::multiply(bool lhs, bool rhs) const {
        return lhs && rhs;
    }
}  // namespace StochasticSimulator

