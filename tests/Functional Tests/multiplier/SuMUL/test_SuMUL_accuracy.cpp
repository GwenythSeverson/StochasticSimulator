#include <gtest/gtest.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

// Four levels up: SuMUL/ -> multiplier/ -> Functional Tests/ -> tests/ -> repo root
#include "../../../../general_functions.hpp"
#include "../../../../modules/SuMUL.hpp"
#include "../../../../bsg/lfsr.hpp"
#include "../../../../bsg/sng.hpp"

// =============================================================================================
// EXHAUSTIVE SuMUL BEHAVIOUR SWEEP
//
// One target probability pair, and EVERY possible stream arrangement that produces it.
//
// A Length_1024 LFSR accepts seeds 1..1023, so there are exactly 1023 distinct streams for a
// given target -- all of them rotations of one underlying sequence. Crossing operand A's 1023
// seeds with operand B's 1023 gives 1023^2 = 1,046,529 arrangements, and this runs all of them.
// Measured at ~12 us/trial optimised, that is ~13 s (slower in a Debug build, see NOTE below).
//
// WHY NOT ALSO SWEEP EVERY TARGET PROBABILITY:
//   There are 1025 reachable targets (k/1024 for k = 0..1024), so the full cross product is
//   1025^2 * 1023^2 = 1.1e12 trials -- about 0.4 CPU-years, and 13 PB if the bit vectors were
//   written out. The compute is merely awful; the storage is impossible. Sweeping the target
//   axis is a separate, much cheaper experiment (1025^2 pairs at one seed each is ~13 s) and
//   belongs in its own test.
//
// WHY THE CSV HOLDS SUMMARIES, NOT BIT VECTORS:
//   1,046,529 rows x 4 bit-vectors is ~12.6 GB and MATLAB would never load it. The heavy lifting
//   is therefore done here in C++ -- every trial is truncated at each length and folded into
//   running mean/min/max/variance accsumulators -- and the CSV carries one row per truncation
//   length. Same "fast C++, MATLAB graphs it" split, just with the aggregation on the fast side.
//
// WHY 0.5 x 0.875 AND NOT 0.5 x 0.5:
//   SuMUL's counter boots half full, i.e. it starts out ASSUMING p1 = 0.5. Run at p1 = 0.5 and
//   that opening guess is already correct, so the warm-up transient -- the most important thing
//   about this unit's short-stream behaviour -- is invisible, and the early points instead show
//   a deterministic startup artifact. At p1 = 0.875 the unit starts out WRONG and has to climb:
//       counter boots at 32/64 = 0.500 and must reach 56/64 = 0.875
//       so the output starts near 0.5 * 0.500 = 0.250 and rises to 0.5 * 0.875 = 0.4375
//
// NOTE ON PICKING YOUR OWN TARGETS: use exact multiples of 1/1024. Values like 0.2, 0.4, 0.7,
//   0.9 and 0.95 are NOT representable, and generate_valid_stream() spins forever on them (over
//   a full period the LFSR visits every value once, so a stream always has exactly
//   floor(p*1024) ones -- retrying changes the phase, never the count). This test sidesteps that
//   entirely by driving the LFSR directly with explicit seeds, which is also what makes it
//   reproducible run to run.
//
// NOTE ON RUNTIME: this is a Debug-build test by default and will take a few minutes there.
//   Build Release, or raise SEED_STRIDE below to subsample the arrangement space.
// =============================================================================================

namespace StochasticSimulator {

namespace {

// 1 = every arrangement (1023^2). Raise to subsample if a Debug run is too slow; the mean is
// stable after only a few dozen samples, though min/max need the full sweep to be exact.
constexpr uint16_t SEED_STRIDE = 1;

constexpr uint16_t STREAM_BITS = 1024;
constexpr uint16_t MIN_SEED = 1;
constexpr uint16_t MAX_SEED = 1023;  // FlexibleLFSR rejects seed >= max_cycles

// Truncation lengths, matching the MATLAB early-termination convention.
const std::vector<uint16_t>& truncation_lengths() {
    static const std::vector<uint16_t> lens = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    return lens;
}

std::vector<bool> stream_for(uint16_t ones_target, uint16_t seed) {
    BitstreamGenerator bsg;
    FlexibleLFSR lfsr(StreamLength::Length_1024, seed);
    const double p = static_cast<double>(ones_target) / STREAM_BITS;

    std::vector<bool> s;
    s.reserve(STREAM_BITS);
    for (uint16_t i = 0; i < STREAM_BITS; ++i) {
        s.push_back(bsg.generate_bit(p, lfsr.next(), STREAM_BITS));
    }
    return s;
}

struct Accsumulator {
    double sum = 0.0;
    double sum_sq = 0.0;
    double sum_abs_err = 0.0;
    double min_est = 2.0;
    double max_est = -1.0;

    void add(double est, double ideal) {
        sum += est;
        sum_sq += est * est;
        sum_abs_err += std::fabs(est - ideal);
        if (est < min_est) min_est = est;
        if (est > max_est) max_est = est;
    }
    double mean(std::size_t n) const { return sum / n; }
    double stddev(std::size_t n) const {
        double m = mean(n);
        double var = sum_sq / n - m * m;
        return var > 0.0 ? std::sqrt(var) : 0.0;
    }
};

}  // namespace

TEST(SUMULAccuracyTest, ExhaustiveArrangementSweep) {
    // =========================================================================
    // 1. CONFIGURATION BLOCK -- targets must be exact multiples of 1/1024
    // =========================================================================
    const uint16_t ones_A = 512;  // p0 = 0.5000, the enable-side operand
    const uint16_t ones_B = 896;  // p1 = 0.8750, the regenerated operand

    const double target_pa = static_cast<double>(ones_A) / STREAM_BITS;
    const double target_pb = static_cast<double>(ones_B) / STREAM_BITS;
    const double ideal_product = target_pa * target_pb;  // 0.4375

    // The operation string is carried in the CSV so the MATLAB title can state exactly what was
    // computed rather than hard-coding it in two places that then drift apart.
    std::ostringstream op;
    op << std::fixed << std::setprecision(4)
       << "SuMUL  " << target_pa << " x " << target_pb << " = " << ideal_product;
    const std::string operation = op.str();

    const std::vector<uint16_t>& lens = truncation_lengths();

    std::cout << "[INFO] " << operation << std::endl;
    std::cout << "[INFO] sweeping every arrangement: seeds " << MIN_SEED << ".." << MAX_SEED
              << " on each operand, stride " << SEED_STRIDE << std::endl;

    // =========================================================================
    // 2. PRE-BUILD EVERY DISTINCT STREAM (1023 per operand, not 1023^2)
    // =========================================================================
    std::vector<uint16_t> seeds;
    for (uint16_t s = MIN_SEED; s <= MAX_SEED; s += SEED_STRIDE) seeds.push_back(s);

    std::vector<std::vector<bool>> A_streams, B_streams;
    A_streams.reserve(seeds.size());
    B_streams.reserve(seeds.size());
    for (uint16_t s : seeds) {
        A_streams.push_back(stream_for(ones_A, s));
        B_streams.push_back(stream_for(ones_B, s));
    }

    // Every arrangement really does carry the intended probability.
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        ASSERT_EQ(count_ones(A_streams[i]), ones_A) << "A seed " << seeds[i];
        ASSERT_EQ(count_ones(B_streams[i]), ones_B) << "B seed " << seeds[i];
    }

    const std::size_t total_trials = seeds.size() * seeds.size();
    std::cout << "[INFO] " << seeds.size() << " streams per operand -> "
              << total_trials << " arrangements" << std::endl;

    // =========================================================================
    // 3. EXHAUSTIVE SWEEP
    // =========================================================================
    std::vector<Accsumulator> acc(lens.size());
    std::size_t done = 0;
    int next_report = 10;

    for (std::size_t ia = 0; ia < A_streams.size(); ++ia) {
        const std::vector<bool>& A = A_streams[ia];
        for (std::size_t ib = 0; ib < B_streams.size(); ++ib) {
            std::vector<bool> Z = sumul_stream(A, B_streams[ib]);

            // One pass over Z, harvesting every truncation length as we go.
            std::size_t ones = 0;
            std::size_t next_len = 0;
            for (uint16_t k = 0; k < STREAM_BITS; ++k) {
                ones += Z[k] ? 1 : 0;
                if (next_len < lens.size() && (k + 1) == lens[next_len]) {
                    acc[next_len].add(static_cast<double>(ones) / lens[next_len], ideal_product);
                    ++next_len;
                }
            }

            ++done;
            int pct = static_cast<int>(100.0 * done / total_trials);
            if (pct >= next_report) {
                std::cout << "[INFO] " << pct << "% (" << done << "/" << total_trials << ")"
                          << std::endl;
                next_report += 10;
            }
        }
    }

    // =========================================================================
    // 4. CSV GENERATION -- one row per truncation length
    // =========================================================================
    std::string filename = "sumul_exhaustive_trials.csv";
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to write CSV file path.";

    csv << "Operation,Target_PA,Target_PB,Ideal_Product,Total_Arrangements,"
           "N,Mean_Est,Mean_AbsError,Min_Est,Max_Est,Std_Est\n";

    csv << std::fixed << std::setprecision(8);
    for (std::size_t i = 0; i < lens.size(); ++i) {
        csv << "\"" << operation << "\"," << target_pa << "," << target_pb << ","
            << ideal_product << "," << total_trials << ","
            << lens[i] << ","
            << acc[i].mean(total_trials) << ","
            << acc[i].sum_abs_err / total_trials << ","
            << acc[i].min_est << "," << acc[i].max_est << ","
            << acc[i].stddev(total_trials) << "\n";
    }
    csv.close();

    // =========================================================================
    // 5. CONSOLE REPORT
    // =========================================================================
    std::cout << std::fixed << std::setprecision(5)
              << "\n     N   mean est   mean |err|     min       max       std\n";
    for (std::size_t i = 0; i < lens.size(); ++i) {
        std::cout << std::setw(6) << lens[i] << "   "
                  << acc[i].mean(total_trials) << "    "
                  << acc[i].sum_abs_err / total_trials << "   "
                  << acc[i].min_est << "   " << acc[i].max_est << "   "
                  << acc[i].stddev(total_trials) << "\n";
    }
    std::cout << "[SUCCESS] Exported aggregate traces to: " << filename << std::endl;

    // =========================================================================
    // 6. GUARDS
    // =========================================================================
    const std::size_t last = lens.size() - 1;
    EXPECT_NEAR(acc[last].mean(total_trials), ideal_product, 0.02)
        << "full-length SuMUL output drifted off the true product";

    // The warm-up must be visible: at N = 16 the counter is only a quarter filled, so the
    // estimate should still be sitting well below the final answer.
    const std::size_t n16 = 3;  // lens = {2,4,8,16,...}
    EXPECT_LT(acc[n16].mean(total_trials), acc[last].mean(total_trials))
        << "expected the counter warm-up to hold the early estimate low";
}

}  // namespace StochasticSimulator
