# StochasticSimulator
A C++ based stochastic computing simulator. 

Organization in modules is a C++ representation of the architecture of each kind of specified part. If an implimentation of a specific adder is missing, add it in the adder files. testing files live under tests and are named as the main testing point. 

GUnit/Gtests will be implimented, with a current plan to create structs that hold trial records, bit stream records, and maybe a function that records a snapshot at every step to put into a CSV to generate grpahs with in MATLAB. 