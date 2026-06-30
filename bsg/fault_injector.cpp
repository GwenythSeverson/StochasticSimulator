#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <numeric>

// Global or static random engine for performance (prevents re-seeding overhead)
static std::mt19937& get_rng() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

namespace FaultInjector {

    // 1. Flip exactly one specific bit by index
    void flipSpecificBit(std::vector<bool>& stream, size_t index) {
        if (index >= stream.size()) {
            throw std::out_of_range("Index out of stream bounds.");
        }
        stream[index] = !stream[index];
    }

    // 2. Flip absolutely all bits in the stream
    void flipAllBits(std::vector<bool>& stream) {
        stream.flip(); // bitwise inversion call
    } // honestly i dont think this can happen irl but if it does, it should be very bad if p!=.50

    // 3. Flip exactly N random bits across the entire stream (Independent Faults)
    void flipRandomBits(std::vector<bool>& stream, size_t num_flips) {
        if (num_flips > stream.size()) {
            num_flips = stream.size(); // Cap at max stream length
        }
        
        auto& gen = get_rng();
        
        // If num_flips is small relative to size, sample random indices
        if (num_flips < stream.size() / 4) {
            std::uniform_int_distribution<size_t> dist(0, stream.size() - 1);
            size_t flips_done = 0;
            while (flips_done < num_flips) {
                size_t idx = dist(gen);
                // Ensure we don't double-flip the same bit in one pass if strict count is needed
                if (stream[idx] == stream[idx]) { // Trivial check, or use a tracking set
                     stream[idx] = !stream[idx];
                     flips_done++;
                }
            }
        } else {
            // For dense faults, shuffle an index vector (Fisher-Yates) to guarantee distinct flips wihtout reselecting already flipped bit
            std::vector<size_t> indices(stream.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::shuffle(indices.begin(), indices.end(), gen);
            
            for (size_t i = 0; i < num_flips; ++i) {
                stream[indices[i]] = !stream[indices[i]];
            }
        }
    }

    // 4. Clumped Multi-Bit Upset (MEU/MBU) at a specified fixed location
    void injectClumpedFaultFixed(std::vector<bool>& stream, size_t start_index, size_t clump_size) {
        if (start_index >= stream.size()) return;
        size_t end_index = std::min(start_index + clump_size, stream.size());
        
        for (size_t i = start_index; i < end_index; ++i) {
            stream[i] = !stream[i];
        }
    }

    // 5. Clumped Multi-Bit Upset (MEU/MBU) at a completely randomized location
    // simulating unpredictable LEO high-energy particle strikes
    void injectClumpedFaultRandom(std::vector<bool>& stream, size_t clump_size) {
        if (stream.empty() || clump_size == 0) return;
        if (clump_size > stream.size()) clump_size = stream.size();

        auto& gen = get_rng();
        std::uniform_int_distribution<size_t> dist(0, stream.size() - clump_size);
        size_t start_index = dist(gen);

        for (size_t i = start_index; i < start_index + clump_size; ++i) {
            stream[i] = !stream[i];
        }
    }
    // 6. Multiple Clumped MEUs/MBUs with variable sizes and random locations
    // Perfect for simulating intense, multi-strike radiation events (e.g., solar flares)
    void injectMultipleClumpedFaultsRandom(std::vector<bool>& stream, size_t num_clumps, size_t min_clump_size, size_t max_clump_size) {
        if (stream.empty() || num_clumps == 0 || min_clump_size == 0 || min_clump_size > max_clump_size) {
            return;
        }

        auto& gen = get_rng();
        std::uniform_int_distribution<size_t> size_dist(min_clump_size, max_clump_size);

        for (size_t c = 0; c < num_clumps; ++c) {
            // Determine the size of this specific clump strike
            size_t current_clump_size = size_dist(gen);
            if (current_clump_size > stream.size()) {
                current_clump_size = stream.size();
            }

            // Pick a random valid starting index where the whole clump fits
            std::uniform_int_distribution<size_t> start_dist(0, stream.size() - current_clump_size);
            size_t start_index = start_dist(gen);

            // Inject the clump by flipping the bits
            for (size_t i = start_index; i < start_index + current_clump_size; ++i) {
                stream[i] = !stream[i];
            }
        }
    }

        ////////////////////////////////////////
        /* please add in future a scalable implementation of LEO radiation rates of flips where all 
        streams are generated and then flipped accourding to poisson distribution of time exposed to 
        see if there is a point where the longer a stream is the more degraded the data is for datacenter 
        sized models - is it even reasonable to think about this? if 4m operations are being done in a nanosecond,
        ~ 100k bits are being flipped in it. How risky is that?  */
}