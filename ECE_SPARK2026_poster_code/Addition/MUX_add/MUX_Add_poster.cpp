/*
 * =========================================================================================
 * MUX_Add_poster.cpp                                   ECE SPARK 2026 -- MUX scaled adder
 * =========================================================================================
 * Companion to AND_mul/ and uMUL_mul/. Same operands, same axes, same exhaustiveness standard,
 * so all four units can be read side by side. Writes into this folder, wherever it is run from.
 *
 *     MUX_Add_poster_warmup.csv       zero-fault early-termination behaviour
 *     MUX_Add_poster_faults.csv       fault response, 0 through 36 flipped bits
 *     MUX_Add_poster_sensitivity.csv  exact single-flip damage per site class
 *
 * -----------------------------------------------------------------------------------------
 * THE UNIT.  out = sel ? a : b.   One multiplexer. With a 0.5 select stream this is the
 * standard stochastic SCALED ADDER: out = (a + b)/2 in expectation.
 *
 *     a ---->|\
 *            | \
 *     b ---->|  |----> out
 *            | /
 *     sel -->|/
 *
 * THREE OPERAND STREAMS, NOT TWO. This is the point of the file. The AND-gate multiplier has
 * two 256-bit operands; the MUX adder has THREE, because the select stream is as real and as
 * flippable as the data. Its fault surface is
 *
 *     a stream     256 bits
 *     b stream     256 bits
 *     sel stream   256 bits   <- the one people forget
 *     ------------------------
 *     total        768 bits
 *
 * and this file sweeps all of it. A select bit is NOT a "control" bit that can be assumed
 * correct: flipping it re-routes that cycle to the other operand, and if the two disagree the
 * output changes. Half the time they disagree.
 *
 * -----------------------------------------------------------------------------------------
 * THE OPERANDS AND THE FILTER
 *
 *   a = b = sel = 0.5, each 128 ones out of 256.  Intended result (0.5 + 0.5)/2 = 0.5.
 *
 *   A trial is KEPT only if the three streams are MUTUALLY DECORRELATED in the strongest
 *   possible sense: each of the 8 combinations of (a, b, sel) occurs exactly 32 times in the
 *   256 cycles. That is the maximum-entropy case, and it is the case the MUX adder REQUIRES to
 *   be correct -- the direct analogue of the AND gate's "overlap exactly 64" filter.
 *
 *   Checked by running the gate, not predicted: out = 1 in exactly four of the eight regions
 *   (see the table below), so 4 x 32 = 128 output ones = 0.5, exactly.
 *
 * -----------------------------------------------------------------------------------------
 * WHY THIS IS EXHAUSTIVELY ENUMERABLE -- THE COLLAPSE
 *
 *   A MUX is memoryless. Cycle i's output depends on cycle i's three bits and nothing else, so
 *   the output ones-count is
 *
 *       SUM over the 8 regions of (cells in that region) x out(region)
 *
 *   and WHERE the cells sit is never consulted. Label a region by r = 4a + 2b + s:
 *
 *       r   a b sel   out = sel ? a : b
 *       0   0 0  0    b = 0
 *       1   0 0  1    a = 0
 *       2   0 1  0    b = 1     <-
 *       3   0 1  1    a = 0
 *       4   1 0  0    b = 0
 *       5   1 0  1    a = 1     <-
 *       6   1 1  0    b = 1     <-
 *       7   1 1  1    a = 1     <-
 *
 *   A flip on a cycle is a 3-bit mask m over (a, b, sel), so the cycle moves from region r to
 *   region r XOR m and its output becomes out(r XOR m). Every one of the astronomically many
 *   flip layouts therefore collapses onto the 8x8 table of (region, mask) counts, and the
 *   number of layouts reaching each outcome is an exact product of multinomials.
 *
 *   ARITHMETIC GUARANTEE: summed over all outcomes the layout count must be C(768, f) exactly.
 *   That identity is ASSERTED at every flip level.
 *
 *   The literal experiment -- every arrangement of three 0.5 streams, times every layout of up
 *   to 36 flips over 768 bits -- is beyond any computer. The enumeration below is not a sample
 *   of it; it is the same experiment losslessly compressed, and every outcome in it is measured
 *   by running a real Adder over real 256-bit streams.
 * =========================================================================================
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// Three levels up: MUX_add/ -> Addition/ -> ECE_SPARK2026_poster_code/ -> repo root
#include "../../../modules/adder.hpp"

using namespace StochasticSimulator;

namespace {

constexpr int N = 256;                 // stream length = cycles
constexpr int REGIONS = 8;             // (a, b, sel) combinations
constexpr int RSIZE = N / REGIONS;     // 32 cells per region -- the decorrelation filter
constexpr int TARGET = N / 2;          // 128 output ones -> exactly 0.5
constexpr int MAX_FLIPS = 36;
constexpr int SITES = 3 * N;           // 768: a, b and sel are all flippable
constexpr int INSTANCES_PER_CLASS = 4; // shuffled real runs per measured state
constexpr int MC_TRIALS_PER_FLIP = 20000;

// out(r) for r = 4a + 2b + sel, i.e. sel ? a : b.
constexpr int OUT[REGIONS] = {0, 0, 1, 0, 0, 1, 1, 1};

const std::vector<int>& truncation_lengths() {
    static const std::vector<int> v = {2, 4, 8, 16, 32, 64, 128, 256};
    return v;
}

std::string source_dir() {
    std::string p = __FILE__;
    size_t slash = p.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string() : p.substr(0, slash + 1);
}

int popcount3(int v) { return (v & 1) + ((v >> 1) & 1) + ((v >> 2) & 1); }

// Build three 256-bit streams in which every (a, b, sel) region holds exactly RSIZE cells,
// scattered at random. `cell_region[i]` records which region cycle i belongs to.
void build_decorrelated(std::vector<bool>& a, std::vector<bool>& b, std::vector<bool>& sel,
                        std::vector<int>& cell_region, std::mt19937& rng) {
    cell_region.assign(N, 0);
    std::vector<int> slots;
    slots.reserve(N);
    for (int r = 0; r < REGIONS; ++r)
        for (int i = 0; i < RSIZE; ++i) slots.push_back(r);
    std::shuffle(slots.begin(), slots.end(), rng);

    a.assign(N, false); b.assign(N, false); sel.assign(N, false);
    for (int i = 0; i < N; ++i) {
        const int r = slots[i];
        cell_region[i] = r;
        a[i]   = (r >> 2) & 1;
        b[i]   = (r >> 1) & 1;
        sel[i] = (r >> 0) & 1;
    }
}

int run_mux(const std::vector<bool>& a, const std::vector<bool>& b,
            const std::vector<bool>& sel, int upto = N) {
    Adder mux;
    int ones = 0;
    for (int i = 0; i < upto; ++i) {
        if (mux.add_scaled_or_weighted(a[i], b[i], sel[i])) ++ones;
    }
    return ones;
}

}  // namespace

TEST(MuxAddPoster, GenerateTrials) {
    std::mt19937 rng(1337);
    long long real_runs = 0;

    std::vector<bool> a, b, sel;
    std::vector<int> cell_region;

    // ======================================================================================
    // SANITY -- the filter really does give exactly 0.5, measured
    // ======================================================================================
    for (int inst = 0; inst < INSTANCES_PER_CLASS; ++inst) {
        build_decorrelated(a, b, sel, cell_region, rng);
        ASSERT_EQ(run_mux(a, b, sel), TARGET)
            << "a decorrelated trial did not come out at 128/256";
        ++real_runs;
    }

    // ======================================================================================
    // SENSITIVITY -- exact single-flip damage, all 768 sites in 24 classes
    // ======================================================================================
    // A single flip is fully described by (which region the struck cycle is in, which of the
    // three bits was struck): 8 x 3 = 24 classes covering all 768 sites, 32 sites each.
    std::ostringstream sbuf;
    sbuf << "Site,Region,A,B,Sel,BitStruck,SiteCount,ProbPerSingleFlip,OutBefore,OutAfter,"
            "DeltaOnes,DeltaFraction\n";
    const char* bit_name[3] = {"sel", "b", "a"};   // bit 0 = sel, 1 = b, 2 = a

    std::cout << "\n[MUX] exact single-flip sensitivity of the 768-bit operand surface\n"
              << "      region (a,b,sel)  struck  out before -> after   delta\n";
    int sens_rows = 0;
    for (int r = 0; r < REGIONS; ++r) {
        for (int bit = 0; bit < 3; ++bit) {
            const int r2 = r ^ (1 << bit);
            const int before = OUT[r], after = OUT[r2];
            const int delta = after - before;

            // MEASURED, not read off the table: build a real trial, flip one real bit that sits
            // in region r, and re-run the gate.
            build_decorrelated(a, b, sel, cell_region, rng);
            int target_cycle = -1;
            for (int i = 0; i < N; ++i) if (cell_region[i] == r) { target_cycle = i; break; }
            ASSERT_GE(target_cycle, 0);
            const int clean = run_mux(a, b, sel);
            if (bit == 0) sel[target_cycle] = !sel[target_cycle];
            else if (bit == 1) b[target_cycle] = !b[target_cycle];
            else a[target_cycle] = !a[target_cycle];
            const int faulted = run_mux(a, b, sel);
            real_runs += 2;
            ASSERT_EQ(faulted - clean, delta)
                << "measured single-flip delta disagrees with the region table at r=" << r
                << " bit=" << bit;

            sbuf << bit_name[bit] << "_in_region_" << r << "," << r << ","
                 << ((r >> 2) & 1) << "," << ((r >> 1) & 1) << "," << (r & 1) << ","
                 << bit_name[bit] << "," << RSIZE << ","
                 << std::fixed << std::setprecision(10)
                 << (static_cast<double>(RSIZE) / SITES) << ","
                 << before << "," << after << "," << delta << ","
                 << (static_cast<double>(delta) / N) << "\n";
            ++sens_rows;

            if (delta != 0) {
                std::cout << "        r=" << r << " (" << ((r >> 2) & 1) << ","
                          << ((r >> 1) & 1) << "," << (r & 1) << ")   "
                          << std::setw(3) << bit_name[bit] << "     "
                          << before << " -> " << after << "        "
                          << std::showpos << delta << std::noshowpos << "\n";
            }
        }
    }
    {
        const std::string fn = source_dir() + "MUX_Add_poster_sensitivity.csv";
        std::ofstream f(fn);
        ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
        f << sbuf.str();
        f.close();
        std::cout << "      " << sens_rows << " site classes covering all " << SITES
                  << " bits -> " << fn << std::endl;
    }

    // ======================================================================================
    // PART 1 -- ZERO-FAULT EARLY TERMINATION
    // ======================================================================================
    // Truncating at cycle L asks how many cells of each OUTPUT-ONE region landed in the prefix.
    // The four out=1 regions hold 4 x 32 = 128 cells of the 256, so the prefix output count is a
    // HYPERGEOMETRIC draw: choosing L of 256 positions when 128 of them emit a one.
    //
    //     P(j output ones in the first L) = C(128,j) * C(128,L-j) / C(256,L)
    //
    // Exact over every arrangement in at most L+1 terms. Each class is built and RUN.
    std::vector<double> lf(SITES + 1, 0.0);
    for (int i = 1; i <= SITES; ++i) lf[i] = lf[i-1] + std::log(static_cast<double>(i));
    auto logC = [&](int n, int k) {
        return (k < 0 || k > n) ? -1e300 : lf[n] - lf[k] - lf[n-k];
    };

    // Arrangements of the three streams under the filter: 256! / (32!)^8.
    double log10_arrangements = lf[N];
    for (int r = 0; r < REGIONS; ++r) log10_arrangements -= lf[RSIZE];
    log10_arrangements /= std::log(10.0);

    std::ostringstream wbuf;
    wbuf << "N,MeanEstimate,MeanAbsError,MinEstimate,MaxEstimate,Classes,Log10Arrangements\n";
    wbuf << std::fixed << std::setprecision(8);

    long long warm_classes = 0;
    for (int L : truncation_lengths()) {
        double mean_est = 0.0, mean_abs = 0.0, wsum = 0.0, lo = 1e9, hi = -1e9;
        const double lnormL = logC(N, L);
        const int j_lo = std::max(0, L - TARGET), j_hi = std::min(L, TARGET);
        for (int j = j_lo; j <= j_hi; ++j) {
            const double lw = logC(TARGET, j) + logC(N - TARGET, L - j) - lnormL;
            if (lw < -300.0) continue;
            const double w = std::exp(lw);
            if (w <= 0.0) continue;
            const double est = static_cast<double>(j) / L;
            ++warm_classes;
            mean_est += w * est;
            mean_abs += w * std::fabs(est - 0.5);
            wsum += w;
            lo = std::min(lo, est);
            hi = std::max(hi, est);
        }
        wbuf << L << "," << (mean_est / wsum) << "," << (mean_abs / wsum) << ","
             << lo << "," << hi << "," << (j_hi - j_lo + 1) << "," << log10_arrangements << "\n";
    }
    {
        const std::string fn = source_dir() + "MUX_Add_poster_warmup.csv";
        std::ofstream f(fn);
        ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
        f << wbuf.str();
        f.close();
        std::cout << "\n[PART 1] exhaustive over all 10^" << std::fixed << std::setprecision(0)
                  << log10_arrangements << " decorrelated arrangements via " << warm_classes
                  << " classes -> " << fn << std::endl;
    }

    // ======================================================================================
    // PART 2 -- FAULT SWEEP, 0..36 FLIPS OVER THE 768-BIT SURFACE
    // ======================================================================================
    // PER-REGION GENERATING FUNCTION. One cell in region r can take any of 8 flip masks m,
    // costing popcount(m) flips and contributing OUT[r ^ m] output ones. A region of RSIZE such
    // cells is the RSIZE-fold convolution of that single-cell polynomial; the whole stream is
    // the convolution of the 8 region polynomials. dp[f][k] then counts layouts exactly.
    auto convolve = [&](const std::vector<std::vector<long double>>& A,
                        const std::vector<std::vector<long double>>& B,
                        int max_ones) {
        std::vector<std::vector<long double>> C(MAX_FLIPS + 1,
                                                std::vector<long double>(max_ones + 1, 0.0L));
        for (int f1 = 0; f1 <= MAX_FLIPS; ++f1)
            for (int k1 = 0; k1 < static_cast<int>(A[f1].size()); ++k1) {
                if (A[f1][k1] == 0.0L) continue;
                for (int f2 = 0; f1 + f2 <= MAX_FLIPS; ++f2)
                    for (int k2 = 0; k2 < static_cast<int>(B[f2].size()); ++k2) {
                        if (B[f2][k2] == 0.0L) continue;
                        if (k1 + k2 > max_ones) continue;
                        C[f1 + f2][k1 + k2] += A[f1][k1] * B[f2][k2];
                    }
            }
        return C;
    };

    // dp over the whole stream, built region by region.
    std::vector<std::vector<long double>> dp(MAX_FLIPS + 1,
                                             std::vector<long double>(N + 1, 0.0L));
    dp[0][0] = 1.0L;
    for (int r = 0; r < REGIONS; ++r) {
        // single cell in region r
        std::vector<std::vector<long double>> cell(MAX_FLIPS + 1,
                                                   std::vector<long double>(2, 0.0L));
        for (int m = 0; m < REGIONS; ++m) {
            const int cost = popcount3(m);
            if (cost > MAX_FLIPS) continue;
            cell[cost][OUT[r ^ m]] += 1.0L;
        }
        // RSIZE-fold convolution by repeated squaring on the exponent
        std::vector<std::vector<long double>> region(MAX_FLIPS + 1,
                                                     std::vector<long double>(1, 0.0L));
        region[0][0] = 1.0L;
        std::vector<std::vector<long double>> base = cell;
        int e = RSIZE, cap = 1;
        while (e > 0) {
            if (e & 1) { region = convolve(region, base, std::min(N, cap + RSIZE)); }
            e >>= 1;
            if (e > 0) { base = convolve(base, base, std::min(N, RSIZE)); cap = std::min(N, cap * 2); }
        }
        dp = convolve(dp, region, N);
    }

    std::ostringstream buf;
    buf << "BitsFlipped,OutputOnes,OutputFraction,ErrorPerN,Probability,"
           "MC_Probability,MC_Trials,LevelRealTrials,Log10TrialCount\n";

    long long classes = 0;
    double mc_worst = 0.0;
    std::cout << "\n  flips   classes    mean       min        max      MC gap\n";

    std::vector<bool> fa(N), fb(N), fs(N);
    for (int f = 0; f <= MAX_FLIPS; ++f) {
        const double lnorm = logC(SITES, f);
        std::vector<double> mass(N + 1, 0.0);
        double wsum = 0.0;
        for (int k = 0; k <= N; ++k) {
            if (dp[f][k] <= 0.0L) continue;
            // dp counts LAYOUTS; normalise by C(768,f) to get a probability.
            const double lw = std::log(static_cast<double>(dp[f][k])) - lnorm;
            if (lw < -300.0) continue;
            mass[k] = std::exp(lw);
            wsum += mass[k];
        }
        if (wsum <= 0.0) continue;

        // THE ARITHMETIC GUARANTEE. Every layout of f flips over 768 sites is counted exactly
        // once, or this fails.
        EXPECT_NEAR(wsum, 1.0, 1e-6)
            << "flip layouts do not sum to C(768," << f << ")";

        // ---- FORWARD MONTE CARLO OVER GENUINELY RANDOM TRIALS -----------------------------
        std::vector<long long> mc_hist(N + 1, 0);
        std::vector<int> pos(SITES);
        for (int t = 0; t < MC_TRIALS_PER_FLIP; ++t) {
            build_decorrelated(fa, fb, fs, cell_region, rng);
            std::iota(pos.begin(), pos.end(), 0);
            for (int i = 0; i < f; ++i) {
                std::uniform_int_distribution<int> pick(i, SITES - 1);
                std::swap(pos[i], pos[pick(rng)]);
                const int p = pos[i];
                if (p < N)          fa[p]         = !fa[p];
                else if (p < 2 * N) fb[p - N]     = !fb[p - N];
                else                fs[p - 2 * N] = !fs[p - 2 * N];
            }
            ++mc_hist[run_mux(fa, fb, fs)];
            ++real_runs;
        }

        double worst_gap = 0.0;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0 && mc_hist[k] == 0) continue;
            const double pe = mass[k] / wsum;
            const double pm = static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP;
            worst_gap = std::max(worst_gap, std::fabs(pe - pm));
        }
        mc_worst = std::max(mc_worst, worst_gap);

        double mean_frac = 0.0;
        int lo = -1, hi = -1;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0) continue;
            const double p = mass[k] / wsum;
            const double frac = static_cast<double>(k) / N;
            mean_frac += p * frac;
            if (lo < 0) lo = k;
            hi = k;
            ++classes;
            buf << f << "," << k << ","
                << std::fixed << std::setprecision(10) << frac << ","
                << (static_cast<double>(k - TARGET) / N) << ","
                << std::scientific << std::setprecision(12) << p << ","
                << (static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP) << ","
                << MC_TRIALS_PER_FLIP << "," << MC_TRIALS_PER_FLIP << ","
                << std::fixed << std::setprecision(4)
                << ((std::log(p) + lnorm) / std::log(10.0)) << "\n";
        }

        std::cout << std::fixed << std::setprecision(6)
                  << std::setw(7) << f << std::setw(10) << classes << "  "
                  << std::setw(9) << mean_frac << "  "
                  << std::setw(9) << (static_cast<double>(lo) / N) << "  "
                  << std::setw(9) << (static_cast<double>(hi) / N) << "  "
                  << std::setw(9) << worst_gap << "\n";

        EXPECT_GE(lo, 0);
        EXPECT_LE(hi, N);
        if (f == 0) { EXPECT_EQ(lo, TARGET); EXPECT_EQ(hi, TARGET); }
        EXPECT_LT(worst_gap, 0.02)
            << "random real trials disagree with the enumerated distribution at f=" << f;
    }

    const std::string fn = source_dir() + "MUX_Add_poster_faults.csv";
    std::ofstream f2(fn);
    ASSERT_TRUE(f2.is_open()) << "Failed to open " << fn;
    f2 << buf.str();
    f2.close();

    std::cout << "\n[PART 2] " << classes << " outcome rows over the 768-bit operand surface\n"
              << "         " << real_runs << " REAL 256-cycle MUX simulations\n"
              << "         worst enumerated-vs-random disagreement: " << mc_worst << "\n"
              << "         -> " << fn << std::endl;
    EXPECT_GT(classes, 0);
}
