# StochasticSimulator
A C++ based stochastic computing simulator. 

Organization in modules is a C++ representation of the architecture of each kind of specified part. If an implimentation of a specific adder is missing, add it in the adder files. testing files live under tests and are named as the main testing point. 

GUnit/Gtests will be implimented, with a current plan to create structs that hold trial records, bit stream records, and maybe a function that records a snapshot at every step to put into a CSV to generate grpahs with in MATLAB. 

## Setup
 install CMake https://cmake.org/download/

## Tests
```powershell
# 1. Configure the build environment
cmake -B build -S .

# 2. Compile the simulator and test targets
cmake --build build --config Release

# 3. Run the verification test suite
.\build\Release\stochastic_computer.exe
```

## File explanations/ progress
Bit Stream Generator (bsg) folder- 
    Linear Feedback Shift Register (lfsr)-
        - 2 mode galois lfsr. 
            - mode 0: auto random for multi run accuracy tests and avoided correlation bais
            - mode hardware seed: deterministic mode for debugging, fault injection, and early termination tests
        - Stochastic Number Generator (sng)- 

keepdata-
    Divider Data-
        Alternative_divider_data-
        CNT_starting_point_data-
            did 2 trial runs of .25/.5 wiht teh counter initialized at different spots (0 & 16)
            and decided that the warmup wasnt incredibly effected
    Multiplier Data-
        just basic mult behavior graphs and the csv it was pulled from

Modeules folder- 
    adder- scaled/weighted and gate streamed, counter not yet implemented
    multiplier- and gate, streamed

Tests-
    Unit Tests- 
        - test_adders- 
            - currently tests truth table of mux signals
            - TODO test weighted/scaled addition functionality with preset streams
        - test_multipliers- 
            - currently tests truth table of and gate multiplier
            - TODO add preset stream test with loop demo
        - test dividers-
            - checks bitwise behavior accuracy in saturation counting divider
    Functional Tests-
        - test_multipler_accuracy- **includes correaltion error in results 
            - creates cvv for matlab graphs and later analysis
        - fault_tests- evaluates how hardware faults affect elements by controlling variance before injecting radiation errors.
