
#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <memory>
#include <iostream>
#include <cmath>
#include <numeric>
#include "adder.hpp"
#include "divider.hpp"
#include "general_functions.hpp"
#include "multiplier.hpp"


/*
* Please note that while I organized this to the best of my ability, this is still a majority of 
* AI generated code, so there may be soem redundancy or inefficiency. I have tried to comment and explain 
* beyond the AI but just message me if you have any questions.
* overview of the test:
* This evaluates how hardware faults affect 
* stochastic computing elements by controlling
*  variance before injecting radiation errors. It 
* uses a deterministic random seed (1337) to 
* make the pseudo-random bit shuffling reproducible across
*  test executions. To eliminate statistical noise,
*  the generateStream and deterministic alignment loop 
*  force exactly 96 overlaps (pre-calculated via expected 
*  cross-product probability) instead of letting standard
*  random generation chuck around an average. This
* ensures that the baseline math matrix has a 
* zero-correlation mathematical error, meaning 
* any deviation found in gunit_output_256bit.csv 
* is strictly caused by the radiation simulation.
* this logs both the original and faulted bitstreams 
* alongside their absolute population counts (a_ones / length)
*  allowing me to see exactly which bit transitions caused a failure if needed.
* the test runner inherits from a base 
* StochasticFaultTest fixture that manages the state of 
* the random number generator and injector. Because the 
* fault configurations are independent of the operational 
* logic, this design abstracts the component under test;
* you can change this code to test adders or 
* dividers simply by replacing the multiplier.multiply() 
* call with the new hardware while reusing the 
* exact same stream generation, logging pipeline, and radiation burst profiles.
*/
enum class FaultType {
    SEU_RANDOM,
    SEU_EVEN_SPACED,
    MEU_BURST
};

class StochasticFaultTest : public ::testing::Test {
protected:
    StochasticSimulator::Multiplier multiplier;
    std::mt19937 rng;

    void SetUp() override {
        rng.seed(1337); // Consistent seed for deterministic testing
    }

    std::string vectorToString(const std::vector<bool>& vec) {
        std::string s = "\"[";
        for (size_t i = 0; i < vec.size(); ++i) {
            s += (vec[i] ? "1" : "0");
            if (i < vec.size() - 1) s += ", ";
        }
        s += "]\"";
        return s;
    }

    // Generates a stream with an exact quantity of '1' bits to eliminate variance
    std::vector<bool> generateStream(size_t length, double probability) {
        std::vector<bool> stream(length, false);
        size_t ones_count = static_cast<size_t>(std::round(length * probability));
        for (size_t i = 0; i < ones_count; ++i) stream[i] = true;
        std::shuffle(stream.begin(), stream.end(), rng);
        return stream;
    }

    void injectFaults(std::vector<bool>& stream, FaultType type, size_t totalFlips, size_t burstSize = 1) {
        size_t length = stream.size();
        if (totalFlips == 0 || length == 0) return;

        if (type == FaultType::SEU_EVEN_SPACED) {
            size_t interval = length / totalFlips;
            if (interval == 0) interval = 1;
            for (size_t i = 0; i < totalFlips && (i * interval) < length; ++i) {
                stream[i * interval] = !stream[i * interval];
            }
        } 
        else if (type == FaultType::SEU_RANDOM) {
            std::uniform_int_distribution<size_t> dist(0, length - 1);
            for (size_t i = 0; i < totalFlips; ++i) {
                size_t idx = dist(rng);
                stream[idx] = !stream[idx];
            }
        } 
        else if (type == FaultType::MEU_BURST) {
            std::uniform_int_distribution<size_t> start_dist(0, length - 1);
            size_t flips_injected = 0;
            while (flips_injected < totalFlips) {
                size_t start_idx = start_dist(rng);
                for (size_t b = 0; b < burstSize && flips_injected < totalFlips; ++b) {
                    stream[(start_idx + b) % length] = !stream[(start_idx + b) % length];
                    flips_injected++;
                }
            }
        }
    }

    // Updated helper to log both original and faulted vectors
    void logToCsv(const std::string& filename, size_t trial, size_t length, 
                  const std::vector<bool>& origA, const std::vector<bool>& origB,
                  const std::vector<bool>& faultedA, const std::vector<bool>& faultedB, 
                  const std::vector<bool>& vecOut, size_t faulted_ones, bool isFirst) {
        std::ofstream csvFile;
        if (isFirst) {
            csvFile.open(filename, std::ios::out);
            csvFile << "Trial,StreamA_Original_Ones_Count_Length,Vector_A_Original,Vector_B_Original,Vector_A_Faulted,Vector_B_Faulted,Vector_Output,Result_Ones_Count_Length\n";
        } else {
            csvFile.open(filename, std::ios::app);
        }

        size_t a_ones = 0;
        for(bool b : origA) if(b) a_ones++;

        csvFile << trial << "," 
                << a_ones << "/" << length << ","
                << vectorToString(origA) << ","
                << vectorToString(origB) << ","
                << vectorToString(faultedA) << ","
                << vectorToString(faultedB) << ","
                << vectorToString(vecOut) << ","
                << faulted_ones << "/" << length << "\n";
        csvFile.close();
    }
};

// Test 1: Verification Accuracy Check using 8 bits
TEST_F(StochasticFaultTest, VerifyAccuracyWith8Bits) {
    size_t bitLength = 8;
    double pA = 0.5; // 4 ones
    double pB = 0.5; // 4 ones
    
    std::vector<bool> vecA = generateStream(bitLength, pA);
    std::vector<bool> vecB = generateStream(bitLength, pB);
    
    // Save original states before injection
    std::vector<bool> origA = vecA;
    std::vector<bool> origB = vecB;
    
    // --- Pre-Fault Zero Correlation Error Validation ---
    size_t expectedCleanOnes = static_cast<size_t>((bitLength * pA) * (bitLength * pB) / bitLength);
    size_t actualCleanOnes = 0;
    for (size_t i = 0; i < bitLength; ++i) {
        if (multiplier.multiply(vecA[i], vecB[i])) actualCleanOnes++;
    }
    
    // Assert zero-correlation mathematical mismatch before faulting
    ASSERT_EQ(actualCleanOnes, expectedCleanOnes) 
        << "Pre-fault correlation error found in 8-bit baseline!";
    
    // Inject a single SEU bitflip to Vector A
    injectFaults(vecA, FaultType::SEU_RANDOM, 1);
    
    std::vector<bool> vecOut(bitLength);
    size_t faulted_ones = 0;
    for (size_t i = 0; i < bitLength; ++i) {
        vecOut[i] = multiplier.multiply(vecA[i], vecB[i]);
        if (vecOut[i]) faulted_ones++;
    }

    logToCsv("gunit_output_8bit.csv", 1, bitLength, origA, origB, vecA, vecB, vecOut, faulted_ones, true);
}

TEST_F(StochasticFaultTest, GenerateGraphDataWith256Bits) {
    size_t bitLength = 256;
    size_t totalTrials = 50;
    double pA = 0.75; 
    double pB = 0.50; 
    size_t radiationIntensity = 12; 
    size_t burstSize = 4;           

    std::string filename = "gunit_output_256bit.csv";

    // Mathematically perfect targets
    size_t countA = static_cast<size_t>(bitLength * pA); // Exactly 192
    size_t countB = static_cast<size_t>(bitLength * pB); // Exactly 128
    size_t expectedCleanOnes = (countA * countB) / bitLength; // Exactly 96

    for (size_t t = 1; t <= totalTrials; ++t) {
        std::vector<bool> vecA(bitLength, false);
        std::vector<bool> vecB(bitLength, false);

        // 1. Create a shared index map
        std::vector<size_t> indices(bitLength);
        std::iota(indices.begin(), indices.end(), 0);
        // Shuffle the indices collectively to keep things stochastic
        std::shuffle(indices.begin(), indices.end(), rng); 

        // 2. Map the bits deterministically relative to each other
        for (size_t i = 0; i < bitLength; ++i) {
            size_t idx = indices[i];
            
            if (i < expectedCleanOnes) {
                vecA[idx] = true;
                vecB[idx] = true;
            } else if (i < countA) {
                vecA[idx] = true;
                vecB[idx] = false;
            } else if (i < countA + (countB - expectedCleanOnes)) {
                vecA[idx] = false;
                vecB[idx] = true;
            } else {
                vecA[idx] = false;
                vecB[idx] = false;
            }
        }

        // 3. HARD TEST: Verify that the pre-fault calculation has absolutely ZERO error.
        size_t actualCleanOnes = 0;
        for (size_t i = 0; i < bitLength; ++i) {
            if (multiplier.multiply(vecA[i], vecB[i])) actualCleanOnes++;
        }
        
        ASSERT_EQ(actualCleanOnes, expectedCleanOnes)
            << "Sanity check failed! Math matrix calculation error.";

        // --- CAPTURE ORIGINALS HERE ---
        std::vector<bool> origA = vecA;
        std::vector<bool> origB = vecB;

        // 4. Inject Radiation Faults
        injectFaults(vecA, FaultType::MEU_BURST, radiationIntensity, burstSize);
        injectFaults(vecB, FaultType::MEU_BURST, radiationIntensity, burstSize);

        // 5. Compute Faulty Output
        std::vector<bool> vecOut(bitLength);
        size_t faulted_ones = 0;
        for (size_t i = 0; i < bitLength; ++i) {
            vecOut[i] = multiplier.multiply(vecA[i], vecB[i]);
            if (vecOut[i]) faulted_ones++;
        }

        logToCsv(filename, t, bitLength, origA, origB, vecA, vecB, vecOut, faulted_ones, (t == 1));
    }
}

void run_fault_tests() {
    std::cout << "To run GUnit tests, execute your primary test binary directly." << std::endl;
}