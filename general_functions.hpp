#pragma once

#include <cstddef>
#include <vector>

namespace stochastic {

inline bool terminate_early_flag = false;



// Determines whether an extreme probability can short-circuit the stream.
bool can_terminate_early(double confidence);

// Calculates Julie Hsiao's Zero Correlation Error (ZCE) over a specific window length.
double calculate_zce_window(const std::vector<int>& streamA, const std::vector<int>& streamB, std::size_t sample_length);

}  // namespace stochastic