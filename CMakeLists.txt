cmake_minimum_required(VERSION 3.16)
project(stochastic_computer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(stochastic_computer
  main.cpp
  general_functions.cpp
  modules/adder.cpp
  modules/multiplier.cpp
  modules/divider.cpp
  tests/fault_tests.cpp
  tests/early_termination_tests.cpp
)

target_include_directories(stochastic_computer PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/modules/include
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/include
)
