cmake_minimum_required(VERSION 3.16)

set(CMAKE_TOOLCHAIN_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/../vcpkg/scripts/buildsystems/vcpkg.cmake"
    CACHE STRING "")

project(stochastic_computer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# GUnit
add_subdirectory(libs/GUnit)

# find GTest (installed via vcpkg)
find_package(GTest CONFIG REQUIRED)

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
  ${CMAKE_CURRENT_SOURCE_DIR}/libs/GUnit/include
)

target_link_libraries(stochastic_computer PRIVATE
    GTest::gtest
    GTest::gtest_main
    gunit
)