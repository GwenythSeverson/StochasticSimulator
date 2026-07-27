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
  bsg/sobol.cpp
  modules/adder.cpp
  modules/multiplier.cpp
  modules/divider.cpp
  tests/Unit\ Tests/test_adders.cpp
  tests/Unit\ Tests/test_multipliers.cpp
  "tests/Unit\ Tests/test_dividers.cpp"
 "tests/Functional Tests/multiplier/test_multiplier_accuracy.cpp"
 "tests/Functional Tests/gen_fault_tests.cpp"

 #"tests/Fault Injection/multiplier/multiplier_injections.cpp"
# "tests/Fault Injection/multiplier/test_multiplier_fault2.cpp"
 #"tests/Fault Injection/adder/test_adder_fault2.cpp"
 "tests/Fault Injection/divider/test_divider_fault2.cpp"
 #"tests/Fault Injection/experiment2.cpp"
 )
# cut out ones for speed
# "tests/Functional Tests/divider/test_dividers_accuracy.cpp"


# 3. Map header file locations
target_include_directories(stochastic_computer PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/bsg
  ${CMAKE_CURRENT_SOURCE_DIR}/modules
  ${CMAKE_CURRENT_SOURCE_DIR}/tests
  ${gtest_SOURCE_DIR}/include
)

# 4. Link against both gtest and gtest_main (gtest_main provides the auto-main function!)
target_link_libraries(stochastic_computer PRIVATE
    gtest
    gtest_main
)