cmake_minimum_required(VERSION 3.16)

project(stochastic_computer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# =====================================================================================
# HOW TO USE THIS FILE
#
# Sources are grouped into named lists below. To drop a whole category from the build,
# comment out its set(...) block -- an undefined variable expands to nothing, so the
# add_executable() at the bottom needs no edit. To drop one test, comment out its line.
#
#   set(UNIT_TESTS ...)        <- comment the whole block to skip every unit test
#     # "tests/Unit Tests/x"   <- or comment one line to skip just that file
#
# BEFORE RECOMPILING, though: you usually do not need to. Filter at run time instead --
# it is instant and it does not invalidate the build:
#   stochastic_computer.exe --gtest_list_tests
#   stochastic_computer.exe --gtest_filter=UnaryMultiplierTest.*
#   stochastic_computer.exe --gtest_filter=*uMUL*:*SuMUL*
#   stochastic_computer.exe --gtest_filter=-*Exhaustive*      (everything EXCEPT the slow ones)
# Comment things out here only when you want the compile itself to be faster.
#
# Approximate run times are noted so you know what you are paying for.
# =====================================================================================

# 1. Automatically fetch and build Google Test from source
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# =====================================================================================
# CORE -- the simulator itself. Not optional; the tests will not link without it.
# (NO main.cpp here -- gtest_main provides main.)
# =====================================================================================
set(CORE_SOURCES
  general_functions.cpp

  # Bitstream generation
  bsg/lfsr.cpp
  bsg/sng.cpp
  bsg/sobol.cpp

  # Arithmetic units
  modules/adder.cpp
  modules/multiplier.cpp          # plain AND-gate multiplier
  modules/SuMUL.cpp               # uGEMM Fig 3(a), sliding-window counter
  modules/uMUL.cpp                # uGEMM Fig 3(a), loaded binary operand
  modules/divider.cpp
)

# =====================================================================================
# UNIT TESTS -- fast, structural, run these constantly. Whole section is ~0.2 s.
# =====================================================================================
set(UNIT_TESTS
  "tests/Unit Tests/test_adders.cpp"
  "tests/Unit Tests/test_multipliers.cpp"
  "tests/Unit Tests/test_SuMUL.cpp"                 # ~0.02 s
  "tests/Unit Tests/test_uMUL.cpp"                  # ~0.08 s
  "tests/Unit Tests/test_dividers.cpp"
)

# =====================================================================================
# FUNCTIONAL TESTS -- accuracy sweeps and characterization. Slow, and several write CSVs
# into the working directory. Comment this block out for a quick edit-compile-run loop.
# =====================================================================================
set(FUNCTIONAL_TESTS
  "tests/Functional Tests/multiplier/ANDmul/test_multiplier_accuracy.cpp"
  "tests/Functional Tests/multiplier/SuMUL/test_SuMUL_accuracy.cpp"    # -> sumul_exhaustive_trials.csv
  "tests/Functional Tests/multiplier/uMUL/test_uMUL_behavior.cpp"      # ~41 s, -> umul_behavior_exhaustive.csv (29 MB)

  # Cut for speed -- re-enable when you actually want the divider numbers.
  # "tests/Functional Tests/divider/test_dividers_accuracy.cpp"
)

# =====================================================================================
# FAULT INJECTION -- single-event-upset campaigns. These are the slow, CSV-heavy ones.
# gen_fault_tests.cpp lives under "Functional Tests/" on disk but belongs here by purpose.
# =====================================================================================
set(FAULT_INJECTION_TESTS
  "tests/Functional Tests/gen_fault_tests.cpp"
  "tests/Fault Injection/divider/test_divider_fault2.cpp"              # ~1 s, -> 32bit_exhaustive_divider_trials.csv

  # Currently disabled -- re-enable individually as needed.
  # "tests/Fault Injection/multiplier/multiplier_injections.cpp"
  # "tests/Fault Injection/multiplier/test_multiplier_fault2.cpp"
  # "tests/Fault Injection/adder/test_adder_fault2.cpp"
  # "tests/Fault Injection/experiment2.cpp"
)

# =====================================================================================
# TARGET -- comment a list out above and it simply vanishes from here.
# =====================================================================================
add_executable(stochastic_computer
  ${CORE_SOURCES}
  ${UNIT_TESTS}
  ${FUNCTIONAL_TESTS}
  ${FAULT_INJECTION_TESTS}
)

# Header search paths
target_include_directories(stochastic_computer PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/bsg
  ${CMAKE_CURRENT_SOURCE_DIR}/modules
  ${CMAKE_CURRENT_SOURCE_DIR}/tests
  ${gtest_SOURCE_DIR}/include
)

# gtest_main provides the auto-main function
target_link_libraries(stochastic_computer PRIVATE
    gtest
    gtest_main
)

# Report what actually went into the build, so a commented-out block is visible at
# configure time instead of being discovered when a test mysteriously stops running.
list(LENGTH UNIT_TESTS _n_unit)
list(LENGTH FUNCTIONAL_TESTS _n_func)
list(LENGTH FAULT_INJECTION_TESTS _n_fault)
message(STATUS "stochastic_computer: ${_n_unit} unit, ${_n_func} functional, ${_n_fault} fault-injection test files")
