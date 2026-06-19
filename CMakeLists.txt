cmake_minimum_required(VERSION 3.16)

project(stochastic_computer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Automatically fetch and build Google Test from source
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# 2. Hardware simulator source files (NO main.cpp here!)
add_executable(stochastic_computer
  general_functions.cpp
  bsg/lfsr.cpp
  bsg/sng.cpp
  modules/adder.cpp
  modules/multiplier.cpp
  modules/divider.cpp
  tests/Unit\ Tests/test_adders.cpp
  tests/Unit\ Tests/test_multipliers.cpp
 "tests/Functional Tests/test_multiplier_accuracy.cpp"
)

# 3. Map header file locations
target_include_directories(stochastic_computer PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/bsg
  ${CMAKE_CURRENT_SOURCE_DIR}/modules
  ${gtest_SOURCE_DIR}/include
)

# 4. Link against both gtest and gtest_main (gtest_main provides the auto-main function!)
target_link_libraries(stochastic_computer PRIVATE
    gtest
    gtest_main
)