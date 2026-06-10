# StochasticSimulator
A C++ based stochastic computing simulator. 

Organization in modules is a C++ representation of the architecture of each kind of specified part. If an implimentation of a specific adder is missing, add it in the adder files. testing files live under tests and are named as the main testing point. 

GUnit/Gtests will be implimented, with a current plan to create structs that hold trial records, bit stream records, and maybe a function that records a snapshot at every step to put into a CSV to generate grpahs with in MATLAB. 

## Setup
1. Clone the repo with submodules:
   git clone --recurse-submodules <your-repo-url>

2. Install vcpkg:
   git clone https://github.com/microsoft/vcpkg
   cd vcpkg
   ./bootstrap-vcpkg.sh

3. Install dependencies via vcpkg:
   vcpkg install gtest

4. Build:
   mkdir build && cd build
   cmake ..
   cmake --build .