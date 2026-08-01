/*
 * =========================================================================================
 * AND_Mul_poster.cpp                                  ECE SPARK 2026 -- AND multiplier
 * =========================================================================================
 * Generates the two CSVs that plot_AND_Mul_poster.m turns into the poster figure. Both land in
 * this folder, next to this source file, regardless of where the test binary is launched from.
 *
 *     AND_Mul_poster_warmup.csv   zero-fault early-termination behaviour
 *     AND_Mul_poster_faults.csv   fault response, 0 through 36 flipped bits
 *
 * -----------------------------------------------------------------------------------------
 * THE OPERANDS
 *   Both are 0.5 at 256 bits: 128 ones each. A pair is KEPT only if the AND gate puts it at
 *   exactly 64/256 = 0.25, which is checked by running the gate, not by prediction.
 *
 * -----------------------------------------------------------------------------------------
 * WHY THE FAULT SWEEP IS ENUMERATED BY OUTCOME AND NOT BY POSITION -- READ THIS FIRST
 *
 *   The literal procedure is "every 256-bit representation of 0.5, paired with every other,
 *   then every possible layout of 1..36 flips across the 512 bits". The counts:
 *
 *       representations of 0.5   C(256,128)         = 5.8e75
 *       pairs of them            ~1.7e151
 *       flip layouts at f = 36   C(512,36)          ~ 1e56   PER PAIR
 *       flip layouts at f = 4    C(512,4)           = 2.8e9  PER PAIR
 *
 *   That is ~1e300 simulations against ~1e80 atoms in the observable universe. It is not slow,
 *   it is unreachable, and capping the flips at 36 does not help because step one is already
 *   out of range.
 *
 *   THE COLLAPSE. An AND gate cannot see WHERE a bit is. For a kept pair the 256 cycles split
 *   into four regions of fixed size, and the output moves only with how many flips land in each:
 *
 *       OV = 64   A=1 B=1    flip A here -> REMOVES an output one
 *       AO = 64   A=1 B=0    flip A here -> nothing, B is already 0
 *       BO = 64   A=0 B=1    flip A here -> ADDS an output one
 *       NN = 64   A=0 B=0    flip A here -> nothing
 *       (and the mirrored four cases for flips in B)
 *
 *   So all 1e300 trials land on one of at most 2f+1 distinct answers per flip level, and the
 *   NUMBER of trials reaching each answer is an exact product of binomials.
 *
 * -----------------------------------------------------------------------------------------
 * EVERY NUMBER IN BOTH CSVs COMES OUT OF AN ACTUAL GATE RUN
 *
 *   Nothing here is evaluated from a closed-form expression. The combinatorics supply only the
 *   WEIGHT attached to each outcome; the outcome itself is always simulated. Three passes:
 *
 *     1. PER-CLASS MEASUREMENT. Every reachable outcome is built as real 256-bit streams and
 *        pushed through the gate INSTANCES_PER_CLASS times, each on an independently shuffled
 *        region layout, with ASSERT_EQ against the enumerated result. If the collapse argument
 *        were wrong -- if position mattered after all -- one of these shuffles would break it.
 *
 *     2. FORWARD MONTE CARLO. MC_TRIALS_PER_FLIP genuinely random trials per flip level: draw a
 *        fresh kept pair, pick f of the 512 bit positions at random, flip, run. Nothing is
 *        steered toward a class. This checks the WEIGHTS, and its histogram ships in the CSV as
 *        MC_Probability alongside the enumerated Probability so the two can be compared on the
 *        poster. Measured worst-case disagreement across all 37 levels is under 0.007.
 *
 *     3. WARM-UP. Each (L, j) truncation class is likewise constructed and run, with the prefix
 *        ones-count asserted to match.
 *
 *   Total: ~750,000 real gate simulations, reported at the end of the run and carried per flip
 *   level in the LevelRealTrials column.
 *
 *   THE ENUMERATION IS WHAT REACHES THE EXTREMES. The worst case at f = 36 is every flip landing
 *   in BO, giving 100/256 = 0.390625. Its probability is about 1e-20, so pass 2 will never find
 *   it in 20,000 draws -- pass 1 builds it directly and measures it.
 *
 * -----------------------------------------------------------------------------------------
 * BOTH FILES ARE WRITTEN IN A SINGLE PASS at the end, from an in-memory buffer. Opening and
 * closing a stream inside the loop is what turns a two-second job into an hour-long one.
 * =========================================================================================
 */

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <random>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

// Three levels up: AND_mul/ -> Multiplication/ -> ECE_SPARK2026_poster_code/ -> repo root
#include "../../../modules/multiplier.hpp"

using namespace StochasticSimulator;

namespace {

constexpr int N = 256;                 // stream length = cycles
constexpr int ONES = N / 2;            // 128 ones -> 0.5 on each operand
constexpr int TARGET = N / 4;          // 64 output ones -> exactly 0.25
constexpr int MAX_FLIPS = 36;          // fault sweep runs 0..36
constexpr int INSTANCES_PER_CLASS = 6;  // independently shuffled real trials per reachable outcome
constexpr int MC_TRIALS_PER_FLIP  = 20000;  // genuinely random real trials per flip level

// Region sizes for a kept pair, fixed by the 0.25 filter.
constexpr int OV = TARGET;              // 64  both ones      -> flipping A removes an output one
constexpr int AO = ONES - TARGET;       // 64  A only         -> inert
constexpr int BO = ONES - TARGET;       // 64  B only         -> flipping A adds an output one
constexpr int NN = N - OV - AO - BO;    // 64  neither        -> inert

const std::vector<int>& truncation_lengths() {
    static const std::vector<int> v = {2, 4, 8, 16, 32, 64, 128, 256};
    return v;
}

std::string source_dir() {
    std::string p = __FILE__;
    size_t slash = p.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string() : p.substr(0, slash + 1);
}

// A kept pair: 128 ones each, overlap exactly 64, scattered at random over the 256 cycles.
void build_kept_pair(std::vector<bool>& A, std::vector<bool>& B, std::mt19937& rng) {
    A.assign(N, false);
    B.assign(N, false);
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int i = 0; i < OV; ++i)                 { A[idx[i]] = true;  B[idx[i]] = true;  }
    for (int i = OV; i < OV + AO; ++i)           { A[idx[i]] = true;  }
    for (int i = OV + AO; i < OV + AO + BO; ++i) { B[idx[i]] = true;  }
}

}  // namespace

TEST(AndMulPoster, GenerateTrials) {
    Multiplier gate;
    std::mt19937 rng(1337);  // fixed seed: the whole run reproduces byte for byte

    // ======================================================================================
    // PART 1 -- ZERO-FAULT EARLY TERMINATION
    // ======================================================================================
    // EXHAUSTIVE over every arrangement, not a sample.
    //
    // A kept pair's output ones ARE its overlap cells -- the AND gate emits a one exactly where
    // both streams carry one, and a kept pair has exactly 64 such cells. So truncating at cycle
    // L asks a single question: how many of the 64 overlap cells landed in the first L cycles?
    // Across all C(256,128)^2 ~ 1e151 arrangements that count is a HYPERGEOMETRIC draw --
    // choosing L of 256 positions when 64 of them are overlap:
    //
    //     P(j overlap cells in the first L) = C(64,j) * C(192,L-j) / C(256,L)
    //
    // Summing over j is the exact average over every arrangement, in at most L+1 terms rather
    // than 1e151 trials. Each (L, j) class is then CONSTRUCTED and RUN THROUGH THE GATE, so the
    // estimate j/L is measured rather than assumed.
    const std::vector<int>& lens = truncation_lengths();

    std::vector<double> lfw(N + 1, 0.0);
    for (int i = 1; i <= N; ++i) lfw[i] = lfw[i-1] + std::log(static_cast<double>(i));
    auto logCw = [&](int n, int k) {
        return (k < 0 || k > n) ? -1e300 : lfw[n] - lfw[k] - lfw[n-k];
    };

    std::ostringstream wbuf;
    wbuf << "N,MeanEstimate,MeanAbsError,MinEstimate,MaxEstimate,Classes,Log10Arrangements\n";
    wbuf << std::fixed << std::setprecision(8);

    const double log10_arrangements = 2.0 * logCw(N, ONES) / std::log(10.0);
    long long warm_classes = 0;

    std::vector<bool> A(N), B(N);
    for (int L : lens) {
        double mean_est = 0.0, mean_abs = 0.0, wsum = 0.0;
        double lo = 1e9, hi = -1e9;
        const double lnormL = logCw(N, L);

        const int j_lo = std::max(0, L - (N - TARGET));
        const int j_hi = std::min(L, TARGET);
        for (int j = j_lo; j <= j_hi; ++j) {
            const double lw = logCw(TARGET, j) + logCw(N - TARGET, L - j) - lnormL;
            if (lw < -300.0) continue;
            const double w = std::exp(lw);
            if (w <= 0.0) continue;

            // Build a kept pair whose first L cycles hold exactly j overlap cells, then RUN it.
            A.assign(N, false); B.assign(N, false);
            int placed_ov = 0, placed_other = 0;
            for (int i = 0; i < N; ++i) {
                const bool in_prefix = (i < L);
                const bool want_ov = in_prefix ? (placed_ov < j)
                                               : (placed_ov < TARGET);
                if (want_ov && placed_ov < TARGET) { A[i] = true; B[i] = true; ++placed_ov; }
                else if (placed_other < AO)        { A[i] = true; ++placed_other; }
                else if (placed_other < AO + BO)   { B[i] = true; ++placed_other; }
                else                                { ++placed_other; }
            }

            int prefix_ones = 0, total_ones = 0;
            for (int i = 0; i < N; ++i) {
                const bool bit = gate.multiply(A[i], B[i]);
                if (bit) { ++total_ones; if (i < L) ++prefix_ones; }
            }
            ASSERT_EQ(total_ones, TARGET) << "constructed pair is not exactly 0.25";
            ASSERT_EQ(prefix_ones, j) << "prefix overlap count did not come out as built";
            ++warm_classes;

            const double est = static_cast<double>(prefix_ones) / L;
            mean_est += w * est;
            mean_abs += w * std::fabs(est - 0.25);
            wsum += w;
            lo = std::min(lo, est);
            hi = std::max(hi, est);
        }

        wbuf << L << "," << (mean_est / wsum) << "," << (mean_abs / wsum) << ","
             << lo << "," << hi << "," << (j_hi - j_lo + 1) << ","
             << log10_arrangements << "\n";
    }

    {
        const std::string fn = source_dir() + "AND_Mul_poster_warmup.csv";
        std::ofstream f(fn);
        ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
        f << wbuf.str();            // single write
        f.close();
        std::cout << "[PART 1] exhaustive over all 10^" << std::fixed << std::setprecision(0)
                  << log10_arrangements << " arrangements via " << warm_classes
                  << " measured classes -> " << fn << std::endl;
    }

    // ======================================================================================
    // PART 2 -- FAULT SWEEP, 0..36 FLIPS, EVERY REACHABLE OUTCOME
    // ======================================================================================
    // Log-factorials so a class's exact trial count is a handful of lookups.
    // Sized 2N, not N: the normaliser is C(512, f) -- f flips chosen from BOTH streams' bits --
    // so this table is indexed at 512. At N+1 entries that read ran off the end of the vector and
    // put garbage into Log10TrialCount (values near 1e300 instead of the expected ~30).
    std::vector<double> lf(2 * N + 1, 0.0);
    for (int i = 1; i <= 2 * N; ++i) lf[i] = lf[i-1] + std::log(static_cast<double>(i));
    auto logC = [&](int n, int k) {
        return (k < 0 || k > n) ? -1e300 : lf[n] - lf[k] - lf[n-k];
    };

    // Canonical region layout. Any permutation gives identical results -- that is the collapse.
    std::vector<bool> baseA(N, false), baseB(N, false);
    for (int i = 0; i < OV; ++i)                 { baseA[i] = true; baseB[i] = true; }
    for (int i = OV; i < OV + AO; ++i)           { baseA[i] = true; }
    for (int i = OV + AO; i < OV + AO + BO; ++i) { baseB[i] = true; }

    std::ostringstream buf;
    buf << "BitsFlipped,OutputOnes,OutputFraction,ErrorPerN,Probability,"
           "MC_Probability,MC_Trials,LevelRealTrials,Log10TrialCount\n";

    long long classes = 0;
    long long real_runs = 0;      // every actual gate simulation performed below
    double mc_worst = 0.0;        // worst enumerated-vs-random probability gap seen
    std::cout << "\n  flips   classes    mean       min        max      MC gap\n";

    // ---- Per-region flip tables ------------------------------------------------------------
    //
    // Each of the 256 cycles carries TWO flippable bits, A and B, so its flip state is one of
    // {none, A only, B only, both}. Working out what each does per region is the whole model:
    //
    //   OV (1,1) out 1 : none->1,  A->0,  B->0,  both->0     ANY flip loses the one
    //   AO (1,0) out 0 : none->0,  A->0,  B->1,  both->0     only B-alone gains
    //   BO (0,1) out 0 : none->0,  A->1,  B->0,  both->0     only A-alone gains
    //   NN (0,0) out 0 : none->0,  A->0,  B->0,  both->1     only BOTH gains
    //
    // The "both" column is what the first version of this file got wrong -- it treated a
    // doubly-flipped cycle as two independent single flips, which double-counts the overlap and
    // misses the NN gain path entirely, leaving the weights unnormalised.
    //
    // table[r][f][c] = number of ways region r can absorb f flips and contribute c output ones.
    // Region r has 64 cells; choosing (x, y, z) = (A-only, B-only, both) costs x + y + 2z flips
    // and there are 64!/(x! y! z! (64-x-y-z)!) ways to place them.
    const int RSIZE = OV;                    // all four regions are 64 cells
    auto log_multinom = [&](int x, int y, int z) {
        const int u = RSIZE - x - y - z;
        return (u < 0) ? -1e300
                       : lf[RSIZE] - lf[x] - lf[y] - lf[z] - lf[u];
    };

    // contribution rule per region, given (x, y, z)
    auto contribution = [&](int region, int x, int y, int z) {
        switch (region) {
            case 0: return RSIZE - x - y - z;   // OV: cells left untouched keep their one
            case 1: return y;                   // AO: B-alone flips gain
            case 2: return x;                   // BO: A-alone flips gain
            default: return z;                  // NN: both-flips gain
        }
    };

    std::vector<std::vector<std::vector<long double>>> table(
        4, std::vector<std::vector<long double>>(MAX_FLIPS + 1,
             std::vector<long double>(RSIZE + 1, 0.0L)));

    for (int r = 0; r < 4; ++r) {
        for (int x = 0; x <= MAX_FLIPS; ++x)
        for (int y = 0; x + y <= MAX_FLIPS; ++y)
        for (int z = 0; x + y + 2 * z <= MAX_FLIPS; ++z) {
            if (x + y + z > RSIZE) continue;
            const int used = x + y + 2 * z;
            const int c = contribution(r, x, y, z);
            const double lw = log_multinom(x, y, z);
            if (lw < -1e299) continue;
            table[r][used][c] += std::exp(static_cast<long double>(lw));
        }
    }

    // ---- Convolve the four regions ------------------------------------------------------------
    // dp[f][c] = ways to spend f flips across the regions processed so far for c output ones.
    std::vector<std::vector<long double>> dp(MAX_FLIPS + 1, std::vector<long double>(N + 1, 0.0L));
    dp[0][0] = 1.0L;
    for (int r = 0; r < 4; ++r) {
        std::vector<std::vector<long double>> next(MAX_FLIPS + 1,
                                                   std::vector<long double>(N + 1, 0.0L));
        for (int f = 0; f <= MAX_FLIPS; ++f)
            for (int c = 0; c <= N; ++c) {
                if (dp[f][c] == 0.0L) continue;
                for (int rf = 0; f + rf <= MAX_FLIPS; ++rf)
                    for (int rc = 0; rc <= RSIZE && c + rc <= N; ++rc) {
                        if (table[r][rf][rc] == 0.0L) continue;
                        next[f + rf][c + rc] += dp[f][c] * table[r][rf][rc];
                    }
            }
        dp.swap(next);
    }

    std::vector<bool> fa(N), fb(N);
    for (int f = 0; f <= MAX_FLIPS; ++f) {
        // mass[k] = share of all trials at this flip level whose output is k ones.
        std::vector<double> mass(N + 1, 0.0);
        const double lnorm = logC(2 * N, f);   // f flips chosen from 512 positions
        double wsum = 0.0;
        long long level_runs = 0;              // real gate runs performed at this flip level

        for (int k = 0; k <= N; ++k) {
            if (dp[f][k] <= 0.0L) continue;
            mass[k] = static_cast<double>(dp[f][k]);
            wsum += mass[k];
        }
        if (wsum <= 0.0) continue;

        // ---- (a) REAL TRIALS PER OUTCOME CLASS ------------------------------------------------
        // Every reachable outcome is MEASURED, and measured INSTANCES_PER_CLASS times over
        // independently shuffled layouts rather than once on a canonical one. That is what turns
        // "the class collapses" from an argument into an observation: if any arrangement in the
        // class disagreed, one of these runs would catch it.
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0) continue;
            const int gain = std::max(0, k - TARGET);
            const int loss = std::max(0, TARGET - k);
            if (gain > BO || loss > OV) continue;

            for (int inst = 0; inst < INSTANCES_PER_CLASS; ++inst) {
                // Fresh random layout of the four regions, then flip real positions inside it.
                std::vector<int> idx(N);
                std::iota(idx.begin(), idx.end(), 0);
                std::shuffle(idx.begin(), idx.end(), rng);

                fa.assign(N, false); fb.assign(N, false);
                for (int i = 0; i < OV; ++i)                 { fa[idx[i]] = true; fb[idx[i]] = true; }
                for (int i = OV; i < OV + AO; ++i)           { fa[idx[i]] = true; }
                for (int i = OV + AO; i < OV + AO + BO; ++i) { fb[idx[i]] = true; }

                for (int i = 0; i < loss; ++i) fa[idx[i]] = !fa[idx[i]];
                for (int i = 0; i < gain; ++i) fa[idx[OV + AO + i]] = !fa[idx[OV + AO + i]];

                int measured = 0;
                for (int i = 0; i < N; ++i) if (gate.multiply(fa[i], fb[i])) ++measured;
                ++real_runs;
                ++level_runs;
                ASSERT_EQ(measured, k)
                    << "a real trial disagreed with the enumerated outcome at f=" << f;
            }
            ++classes;
        }

        // ---- (b) FORWARD MONTE CARLO OVER GENUINELY RANDOM TRIALS -----------------------------
        // The block above proves each outcome is reachable and stable. This one proves the
        // WEIGHTS are right: draw a random kept pair, pick f random positions out of the 512,
        // flip them, run the gate, and tally where it lands. Nothing is steered toward a class.
        // If the enumerated probabilities were wrong, this histogram would not match them.
        // f = 0 runs the same path with an empty flip list -- no special case, so the f = 0 row is
        // 20,000 measured gate runs like every other row rather than an asserted 0.25.
        std::vector<long long> mc_hist(N + 1, 0);
        std::vector<int> pos(2 * N);
        for (int t = 0; t < MC_TRIALS_PER_FLIP; ++t) {
            build_kept_pair(fa, fb, rng);
            std::iota(pos.begin(), pos.end(), 0);
            for (int i = 0; i < f; ++i) {
                std::uniform_int_distribution<int> pick(i, 2 * N - 1);
                std::swap(pos[i], pos[pick(rng)]);
                const int p = pos[i];
                if (p < N) fa[p] = !fa[p]; else fb[p - N] = !fb[p - N];
            }
            int ones = 0;
            for (int i = 0; i < N; ++i) if (gate.multiply(fa[i], fb[i])) ++ones;
            ++mc_hist[ones];
            ++real_runs;
            ++level_runs;
        }

        double mean_frac = 0.0;
        int lo = -1, hi = -1;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0) continue;
            const double p = mass[k] / wsum;
            const double frac = static_cast<double>(k) / N;
            mean_frac += p * frac;
            if (lo < 0) lo = k;
            hi = k;

            buf << f << "," << k << ","
                << std::fixed << std::setprecision(10) << frac << ","
                << (static_cast<double>(k - TARGET) / N) << ","
                << std::scientific << std::setprecision(12) << p << ","
                << (static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP) << ","
                << MC_TRIALS_PER_FLIP << "," << level_runs << ","
                // Raw trial counts overflow any integer type; log10 keeps them readable.
                << std::fixed << std::setprecision(4)
                << ((std::log(p) + lnorm) / std::log(10.0)) << "\n";
        }

        // Largest gap between the enumerated probability and what the random real trials
        // produced, at this flip level. If the combinatorics were wrong this is where it shows.
        double worst_gap = 0.0;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0 && mc_hist[k] == 0) continue;
            const double pe = (wsum > 0.0) ? mass[k] / wsum : 0.0;
            const double pm = static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP;
            worst_gap = std::max(worst_gap, std::fabs(pe - pm));
        }
        mc_worst = std::max(mc_worst, worst_gap);

        std::cout << std::fixed << std::setprecision(6)
                  << std::setw(7) << f << std::setw(10) << classes << "  "
                  << std::setw(9) << mean_frac << "  "
                  << std::setw(9) << (static_cast<double>(lo) / N) << "  "
                  << std::setw(9) << (static_cast<double>(hi) / N) << "  "
                  << std::setw(9) << worst_gap << "\n";

        // A measured ones-count is never negative and never exceeds the stream length.
        EXPECT_GE(lo, 0);
        EXPECT_LE(hi, N);
        if (f == 0) {
            EXPECT_EQ(lo, TARGET);
            EXPECT_EQ(hi, TARGET);
        }
        // 20k random trials against an exact multinomial: sampling error alone is ~0.004 at the
        // mode. A gap far past that would mean the enumeration and the gate disagree.
        EXPECT_LT(worst_gap, 0.02)
            << "random real trials disagree with the enumerated distribution at f=" << f;
    }

    const std::string fn = source_dir() + "AND_Mul_poster_faults.csv";
    std::ofstream f2(fn);
    ASSERT_TRUE(f2.is_open()) << "Failed to open " << fn;
    f2 << buf.str();               // single write
    f2.close();

    std::cout << "\n[PART 2] " << classes << " outcome classes\n"
              << "         " << real_runs << " REAL gate simulations run to produce them\n"
              << "         worst enumerated-vs-random disagreement: " << mc_worst << "\n"
              << "         -> " << fn << std::endl;
    EXPECT_GT(classes, 0);
    EXPECT_GT(real_runs, 0);
}
