# StochasticSimulator
A C++ based stochastic computing simulator. 

Organization in modules is a C++ representation of the architecture of each kind of specified part. If an implimentation of a specific adder is missing, add it in the adder files. testing files live under tests and are named as the main testing point. 

GUnit/Gtests will be implimented, with a current plan to create structs that hold trial records, bit stream records, and maybe a function that records a snapshot at every step to put into a CSV to generate grpahs with in MATLAB. 

## Setup
vcpkg for gtest/gunit- will make setup instructions or make easier for someone to use later

## File explanations/ progress
bsg folder- 
    Bit Stream Generator with galois lfsr and sotchastic number generator

Modeules folder- 
    adder- scaled/weighted and gate streamed, counter not yet implemented
    multiplier- and gate, streamed
Tests-
    no unit or flip/Etermination tests. 
