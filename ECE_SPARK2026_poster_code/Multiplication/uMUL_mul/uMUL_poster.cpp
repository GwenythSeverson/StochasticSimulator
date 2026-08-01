/*
 * =========================================================================================
 * uMUL_poster.cpp                                      ECE SPARK 2026 -- uMUL multiplier
 * =========================================================================================
 * Twin of AND_mul/AND_Mul_poster.cpp, same operands, same axes, same exhaustiveness standard,
 * so the two figures can be read side by side. Generates the two CSVs that
 * plot_uMUL_poster.m turns into the poster figures. Both land in this folder, next to this
 * source file, regardless of where the test binary is launched from.
 *
 *     uMUL_poster_warmup.csv   zero-fault early-termination behaviour
 *     uMUL_poster_faults.csv   fault response, 0 through 36 flipped bits
 *
 * See HOW_THIS_DATA_WAS_MADE.txt in this folder for the full reasoning, the trial-count
 * arithmetic behind every decision, and the honest limitations.
 *
 * -----------------------------------------------------------------------------------------
 * THE OPERANDS -- AND THE ASYMMETRY THAT IS THE WHOLE POINT
 *
 *   in_0 = 0.5 as a 256-bit unary STREAM        128 ones out of 256 bits
 *   in_1 = 0.5 as an 8-bit BINARY REGISTER      value 128 out of 256   (uMUL_uni.sv's `iB`)
 *   product 0.25                                64 output ones out of 256
 *
 *   The AND gate next door holds both operands as 256-bit unary streams: 512 bits of operand
 *   storage. uMUL holds one stream and one register: 264 bits. Half the storage -- and that
 *   halving is exactly where its accuracy advantage comes from, because the register is exact
 *   while a unary stream is only a sample. This file measures what that trade costs when the
 *   bits get hit.
 *
 * -----------------------------------------------------------------------------------------
 * NO FILTER IS NEEDED HERE, AND THAT IS ALREADY A RESULT
 *
 *   The AND-gate file has to DISCARD operand pairs: two arbitrary 0.5 streams do not multiply
 *   to 0.25, they multiply to whatever their overlap happens to be, so only the pairs the gate
 *   puts at exactly 64/256 are kept. That is the correlation problem.
 *
 *   uMUL has no such problem and therefore no filter. Every one of the C(256,128) = 5.8e75
 *   arrangements of in_0 gives exactly 64 output ones. The reason is the enable:
 *
 *       out[t] = in_0[t] AND (value > sobol[index])       index++ only when in_0[t] == 1
 *
 *   The RNG index advances once per one in in_0, so the enabled cycles consume Sobol points
 *   0, 1, 2, ... n1-1 in order, where n1 is the NUMBER of ones in in_0. Where those ones sit
 *   is never consulted. So
 *
 *       output ones = G(value, n1) = #{ j < n1 : value > sobol[j] }
 *
 *   and that is a function of two integers, full stop. This is a stronger collapse than the
 *   AND gate's: the AND gate collapses over flip POSITIONS but still depends on the operand
 *   pair, while uMUL collapses the entire 5.8e75-arrangement operand space onto one number.
 *
 * -----------------------------------------------------------------------------------------
 * THE FAULT SURFACE COVERED HERE -- 264 BITS, AND WHY THOSE 264
 *
 *   Flips are injected into the OPERAND STORAGE, which is the direct counterpart of the AND
 *   gate's two operand streams:
 *
 *       in_0 stream        256 bits    flipping one toggles whether that cycle is enabled,
 *                                      so n1 moves by +/-1 and the answer moves by <= 1/256
 *       value register       8 bits    flipping one changes the operand for the WHOLE run
 *       ----------------------------
 *       total              264 bits
 *
 *   Both are static for the run, so both preserve the collapse: the faulted output is still
 *   G(value', n1'), just with a different pair of integers. That is what keeps this sweep
 *   exactly enumerable rather than sampled.
 *
 *   THE REGISTER BITS ARE NOT LIKE THE STREAM BITS. Sensitivity, exact, printed at run time:
 *   flipping stream bits moves the answer by at most 1/256 each, while flipping value bit 7
 *   takes the operand to 0 and the answer to a hard zero -- from ONE upset. 8 of the 264 bits
 *   carry essentially all of the risk. That asymmetry is the figure.
 *
 *   OUT OF SCOPE, deliberately, and documented in the .txt rather than hidden: the rng_index
 *   register, the RNG output bus, and the enable/output wires. Those are generator internals
 *   whose upsets are TIME-dependent -- an rng_index hit at cycle t only disturbs the enabled
 *   cycles after t -- so they do not collapse to a pair of integers and cannot be enumerated
 *   exhaustively at this scale. Including them would force this file to become a sampler, and
 *   the whole value of the AND-gate twin is that it is not one.
 *
 * -----------------------------------------------------------------------------------------
 * EVERY NUMBER IN BOTH CSVs COMES OUT OF AN ACTUAL GATE RUN
 *
 *   G(value, n1) is never evaluated from a formula. It is measured by building a real 256-bit
 *   in_0 stream, loading the real value register, and clocking a real UnaryMultiplier 256
 *   times. Each distinct (value, n1) state is measured INSTANCES_PER_CLASS times on
 *   INDEPENDENTLY SHUFFLED arrangements and the results are asserted equal -- which is what
 *   turns "position does not matter" from an argument into an observation.
 *
 *   The combinatorics supply only the WEIGHT on each measured outcome, and those weights are
 *   themselves checked against MC_TRIALS_PER_FLIP genuinely random trials per flip level.
 *
 * -----------------------------------------------------------------------------------------
 * BOTH FILES ARE WRITTEN IN A SINGLE PASS at the end, from an in-memory buffer.
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

// Three levels up: uMUL_mul/ -> Multiplication/ -> ECE_SPARK2026_poster_code/ -> repo root
#include "../../../modules/uMUL.hpp"

using namespace StochasticSimulator;

namespace {

constexpr int N = 256;                  // stream length = cycles
constexpr int ONES = N / 2;             // 128 ones -> in_0 = 0.5
constexpr unsigned WIDTH = 8;           // value register width; matches uMUL_uni.sv's 8-bit iB
constexpr uint64_t ENTRY = 1ull << WIDTH;   // 256 -- the register's full scale and the RNG period
constexpr uint64_t VALUE = ENTRY / 2;   // 128 -> in_1 = 0.5
constexpr int TARGET = N / 4;           // 64 output ones -> exactly 0.25
constexpr int MAX_FLIPS = 36;           // fault sweep runs 0..36
constexpr int REG_BITS = static_cast<int>(WIDTH);   // 8 flippable operand-register bits
constexpr int SITES = N + REG_BITS;     // 264 -- the operand surface swept below

constexpr int INSTANCES_PER_CLASS = 4;      // independently shuffled real runs per (value, n1) state
constexpr int MC_TRIALS_PER_FLIP  = 20000;  // genuinely random real trials per flip level

const std::vector<int>& truncation_lengths() {
    static const std::vector<int> v = {2, 4, 8, 16, 32, 64, 128, 256};
    return v;
}

std::string source_dir() {
    std::string p = __FILE__;
    size_t slash = p.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string() : p.substr(0, slash + 1);
}

// Hand-rolled rather than __builtin_popcount / std::popcount: the former is GCC-and-Clang only
// and this repo builds under MSVC, the latter needs C++20 and this project is on C++17.
int popcount8(uint32_t v) {
    int n = 0;
    while (v) { n += static_cast<int>(v & 1u); v >>= 1; }
    return n;
}

// A 256-bit in_0 stream carrying exactly `ones` ones, scattered at random.
void build_stream(std::vector<bool>& s, int ones, std::mt19937& rng) {
    s.assign(N, false);
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int i = 0; i < ones; ++i) s[idx[i]] = true;
}

// One real run: load `value`, clock the stream through, count output ones.
int run_gate(const std::vector<bool>& in_0, uint64_t value) {
    UnaryMultiplier umul(WIDTH, 1);
    umul.load_value(value);
    int ones = 0;
    for (int i = 0; i < N; ++i) if (umul.multiply(in_0[i])) ++ones;
    return ones;
}

}  // namespace

TEST(uMULPoster, GenerateTrials) {
    std::mt19937 rng(1337);   // fixed seed: the whole run reproduces byte for byte
    long long real_runs = 0;  // every actual 256-cycle gate simulation performed below

    // ======================================================================================
    // G(value, n1) -- MEASURED, MEMOISED, AND VERIFIED
    // ======================================================================================
    // The one primitive everything else is built from. On first request for a (value, n1) pair
    // it builds INSTANCES_PER_CLASS independently shuffled in_0 streams carrying n1 ones, runs
    // each through a real UnaryMultiplier, and ASSERTS all four agree. That assertion IS the
    // collapse claim of the header: if where the ones sat could change the answer, four random
    // shuffles would catch it. After that the value is cached, because it has been established
    // that there is nothing left to vary.
    std::vector<std::vector<int>> memo(ENTRY, std::vector<int>(N + 1, -1));
    std::vector<bool> scratch(N);
    auto G = [&](uint64_t value, int n1) -> int {
        if (n1 < 0 || n1 > N) return -1;
        int& slot = memo[value][n1];
        if (slot >= 0) return slot;
        int first = -1;
        for (int inst = 0; inst < INSTANCES_PER_CLASS; ++inst) {
            build_stream(scratch, n1, rng);
            const int measured = run_gate(scratch, value);
            ++real_runs;
            if (inst == 0) first = measured;
            else EXPECT_EQ(measured, first)
                << "uMUL output depended on WHERE the ones sat: value=" << value
                << " n1=" << n1;
        }
        slot = first;
        return first;
    };

    // The healthy answer, measured rather than assumed.
    ASSERT_EQ(G(VALUE, ONES), TARGET)
        << "0.5 x 0.5 did not come out at 64/256 with a clean unit";

    // ---- Exact single-bit sensitivity of every operand site ----------------------------------
    // The sharpest numbers on the poster: one upset in bit 7 destroys the answer outright, while
    // one upset anywhere in the 256-bit stream cannot move it past 1/256. Written to its own CSV
    // because figure 3 is built on it -- it is the "where does uMUL fail" half of the story.
    std::ostringstream sbuf;
    sbuf << "Site,SiteCount,ProbPerSingleFlip,FaultedValue,OutputOnes,OutputFraction,Delta\n";
    auto sens_row = [&](const std::string& name, int count, long long vshown, int k) {
        sbuf << name << "," << count << ","
             << std::fixed << std::setprecision(10)
             << (static_cast<double>(count) / SITES) << "," << vshown << "," << k << ","
             << (static_cast<double>(k) / N) << ","
             << (static_cast<double>(k - TARGET) / N) << "\n";
    };

    std::cout << "\n[uMUL] exact single-flip sensitivity of the 264-bit operand surface"
              << " (healthy = 0.250000)\n";
    for (int b = REG_BITS - 1; b >= 0; --b) {
        const uint64_t v = VALUE ^ (1ull << b);
        const int k = G(v, ONES);
        sens_row("value_bit_" + std::to_string(b), 1, static_cast<long long>(v), k);
        std::cout << "        value bit " << b << "  value " << std::setw(4) << v
                  << " -> " << std::setw(4) << k << "/256 = " << std::fixed
                  << std::setprecision(6) << (static_cast<double>(k) / N)
                  << "   delta " << std::showpos << (static_cast<double>(k - TARGET) / N)
                  << std::noshowpos << "\n";
    }
    {
        // A stream flip lands on a cycle that held a one (n1 -> 127) or a zero (n1 -> 129).
        // 128 sites each, and the operand register is untouched, so the value column stays 128.
        const int lo = G(VALUE, ONES - 1), hi = G(VALUE, ONES + 1);
        sens_row("stream_bit_on_a_one",  ONES,     static_cast<long long>(VALUE), lo);
        sens_row("stream_bit_on_a_zero", N - ONES, static_cast<long long>(VALUE), hi);
        std::cout << "        stream bit on a one  (128 -> 127 ones) -> " << lo << "/256\n"
                  << "        stream bit on a zero (128 -> 129 ones) -> " << hi << "/256\n"
                  << "        i.e. any one of the 256 stream bits moves the answer by 1/256"
                  << " = 0.003906\n";
    }
    {
        const std::string fn = source_dir() + "uMUL_poster_sensitivity.csv";
        std::ofstream f(fn);
        ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
        f << sbuf.str();            // single write
        f.close();
        std::cout << "        -> " << fn << std::endl;
    }

    // ======================================================================================
    // PART 1 -- ZERO-FAULT EARLY TERMINATION
    // ======================================================================================
    // EXHAUSTIVE over every arrangement of in_0, not a sample.
    //
    // Truncating at cycle L asks one question: how many ones of in_0 landed in the first L
    // cycles? Call it n. The enabled cycles inside the prefix are the FIRST n enabled cycles of
    // the run, so they consume Sobol points 0..n-1 and the prefix output is G(value, n) -- the
    // same primitive as everywhere else. Across all C(256,128) arrangements n is a
    // HYPERGEOMETRIC draw, choosing L of 256 positions when 128 of them are ones:
    //
    //     P(n ones in the first L) = C(128,n) * C(128,L-n) / C(256,L)
    //
    // Summing over n is the exact average over all 5.8e75 arrangements in at most L+1 terms.
    //
    // NOTE THE CONTRAST WITH THE AND GATE. There the prefix count was the OVERLAP of two
    // streams, drawn from only 64 overlap cells; here it is the ones of a single stream, drawn
    // from 128. And G halves it: G(128, n) rises by one every two ones. So uMUL's truncation
    // error is roughly half the AND gate's at the same L, which is what the two figure 1s show.
    const std::vector<int>& lens = truncation_lengths();

    std::vector<double> lf(2 * N + 1, 0.0);
    for (int i = 1; i <= 2 * N; ++i) lf[i] = lf[i-1] + std::log(static_cast<double>(i));
    auto logC = [&](int n, int k) {
        return (k < 0 || k > n) ? -1e300 : lf[n] - lf[k] - lf[n-k];
    };

    // One stream's worth of arrangements this time, not a pair: in_1 is a register, not a stream.
    const double log10_arrangements = logC(N, ONES) / std::log(10.0);

    std::ostringstream wbuf;
    wbuf << "N,MeanEstimate,MeanAbsError,MinEstimate,MaxEstimate,Classes,Log10Arrangements\n";
    wbuf << std::fixed << std::setprecision(8);

    long long warm_classes = 0;
    for (int L : lens) {
        double mean_est = 0.0, mean_abs = 0.0, wsum = 0.0;
        double lo = 1e9, hi = -1e9;
        const double lnormL = logC(N, L);

        const int n_lo = std::max(0, L - (N - ONES));
        const int n_hi = std::min(L, ONES);
        for (int n = n_lo; n <= n_hi; ++n) {
            const double lw = logC(ONES, n) + logC(N - ONES, L - n) - lnormL;
            if (lw < -300.0) continue;
            const double w = std::exp(lw);
            if (w <= 0.0) continue;

            // Prefix output ones = G(value, n): the first n enabled cycles, measured.
            const double est = static_cast<double>(G(VALUE, n)) / L;
            ++warm_classes;

            mean_est += w * est;
            mean_abs += w * std::fabs(est - 0.25);
            wsum += w;
            lo = std::min(lo, est);
            hi = std::max(hi, est);
        }

        wbuf << L << "," << (mean_est / wsum) << "," << (mean_abs / wsum) << ","
             << lo << "," << hi << "," << (n_hi - n_lo + 1) << ","
             << log10_arrangements << "\n";
    }

    {
        const std::string fn = source_dir() + "uMUL_poster_warmup.csv";
        std::ofstream f(fn);
        ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
        f << wbuf.str();            // single write
        f.close();
        std::cout << "\n[PART 1] exhaustive over all 10^" << std::fixed << std::setprecision(0)
                  << log10_arrangements << " in_0 arrangements via " << warm_classes
                  << " measured classes -> " << fn << std::endl;
    }

    // ======================================================================================
    // PART 2 -- FAULT SWEEP, 0..36 FLIPS OVER THE 264-BIT OPERAND SURFACE
    // ======================================================================================
    // THE ENUMERATION. f flips are placed on 264 sites. Split them by where they land:
    //
    //     r  flips on the 8 register bits, choosing a specific SUBSET S -- C(8,r) subsets,
    //        each giving its own faulted operand  value' = 128 XOR mask(S)
    //     a  flips on the 128 stream positions holding a one   -> that cycle stops enabling
    //     b  flips on the 128 stream positions holding a zero  -> that cycle starts enabling
    //        with r + a + b = f, so n1' = 128 - a + b
    //
    // The faulted output is G(value', n1'), measured. The number of layouts in the class is
    //
    //     ways(r-subset, a, b) = C(128,a) * C(128,b)          (the subset itself is one choice)
    //
    // and summing checks out exactly by Vandermonde, which is the arithmetic guarantee that
    // nothing is double counted and nothing is missed:
    //
    //     SUM_r C(8,r) * SUM_a C(128,a)*C(128,f-r-a) = SUM_r C(8,r)*C(256,f-r) = C(264,f)
    //
    // That identity is ASSERTED at every flip level below, in log space.
    std::ostringstream buf;
    buf << "BitsFlipped,OutputOnes,OutputFraction,ErrorPerN,Probability,"
           "MC_Probability,MC_Trials,LevelRealTrials,ProbRegisterClean,ProbRegisterHit,"
           "ProbMSBIntact,ProbMSBLost,Log10TrialCount\n";

    long long classes = 0;
    double mc_worst = 0.0;
    std::cout << "\n  flips   classes    mean       min        max     P(reg hit)   MC gap\n";

    std::vector<bool> fs(N);
    for (int f = 0; f <= MAX_FLIPS; ++f) {
        // mass[k] = share of all layouts at this flip level whose output is k ones, carried twice
        // over in two different decompositions. Both ship in the CSV because between them they
        // are the entire explanation of the figure.
        //
        //   clean / hit   was ANY of the 8 register bits struck?
        //   MSB intact /  was BIT 7 specifically struck? This is the split that matters, and it
        //   MSB lost      is not the same question. value = 128 has exactly ONE bit set, so:
        //                     bit 7 survives -> value' = 128 + mask(S) in [128,255] -> out >= 0.25
        //                     bit 7 dies     -> value' = mask(S\{7})   in [0,127]   -> out <  0.25
        //                 Killing the MSB leaves the operand equal to whatever the surviving low
        //                 bits say, which is a small number, so those trials collapse toward zero.
        //                 That population IS the low shelf of figure 3, and its mass is exactly
        //                 f/264 -- one specific site out of 264 -- NOT the register-hit mass.
        //                 Only about 13-20% of register hits land below 0.25; the rest go up.
        std::vector<double> mass(N + 1, 0.0), mass_clean(N + 1, 0.0), mass_hit(N + 1, 0.0),
                            mass_msb_ok(N + 1, 0.0), mass_msb_lost(N + 1, 0.0);
        const double lnorm = logC(SITES, f);     // f flips chosen from 264 sites
        double wsum = 0.0, wsum_hit = 0.0;
        long long level_runs = 0;
        const long long runs_before = real_runs;
        long long level_classes = 0;

        for (int r = 0; r <= std::min(REG_BITS, f); ++r) {
            // Every subset of exactly r register bits, as a bitmask over the 8 bits.
            for (uint32_t S = 0; S < (1u << REG_BITS); ++S) {
                if (popcount8(S) != r) continue;
                const uint64_t vfault = VALUE ^ static_cast<uint64_t>(S);

                const int rem = f - r;
                for (int a = 0; a <= std::min(ONES, rem); ++a) {
                    const int b = rem - a;
                    if (b > N - ONES) continue;
                    const int n1 = ONES - a + b;

                    const double lw = logC(ONES, a) + logC(N - ONES, b) - lnorm;
                    if (lw < -300.0) continue;
                    const double w = std::exp(lw);
                    if (w <= 0.0) continue;

                    const int k = G(vfault, n1);   // <- the real gate run (memoised)
                    ++level_classes;

                    mass[k] += w;
                    wsum += w;
                    if (r == 0) mass_clean[k] += w;
                    else      { mass_hit[k]  += w; wsum_hit += w; }
                    // The MSB is the top register bit, so its mask bit is 1 << (REG_BITS - 1).
                    if (S & (1u << (REG_BITS - 1))) mass_msb_lost[k] += w;
                    else                            mass_msb_ok[k]   += w;
                }
            }
        }
        if (wsum <= 0.0) continue;

        // The Vandermonde check. Every layout of f flips over 264 sites is accounted for exactly
        // once, or this fails. Tolerance is generous only because the sum is built in log space.
        EXPECT_NEAR(wsum, 1.0, 1e-9)
            << "flip layouts do not sum to C(264," << f << ") at f=" << f;

        // ---- FORWARD MONTE CARLO OVER GENUINELY RANDOM TRIALS --------------------------------
        // The enumeration above proves each outcome is reachable and stable. This proves the
        // WEIGHTS are right: draw a fresh random in_0 with 128 ones, pick f of the 264 sites
        // uniformly at random, apply them for real -- toggling stream bits, XORing register bits
        // -- and clock the whole thing through. Nothing is steered toward a class.
        std::vector<long long> mc_hist(N + 1, 0);
        std::vector<int> pos(SITES);
        for (int t = 0; t < MC_TRIALS_PER_FLIP; ++t) {
            build_stream(fs, ONES, rng);
            uint64_t v = VALUE;
            std::iota(pos.begin(), pos.end(), 0);
            for (int i = 0; i < f; ++i) {
                std::uniform_int_distribution<int> pick(i, SITES - 1);
                std::swap(pos[i], pos[pick(rng)]);
                const int p = pos[i];
                if (p < N) fs[p] = !fs[p];
                else       v ^= (1ull << (p - N));
            }
            ++mc_hist[run_gate(fs, v)];
            ++real_runs;
            ++level_runs;
        }

        // Largest gap between the enumerated probability and what the random trials produced.
        double worst_gap = 0.0;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0 && mc_hist[k] == 0) continue;
            const double pe = mass[k] / wsum;
            const double pm = static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP;
            worst_gap = std::max(worst_gap, std::fabs(pe - pm));
        }
        mc_worst = std::max(mc_worst, worst_gap);

        // Real runs at this level = the MC trials above plus any fresh G() measurements the
        // enumeration triggered. States already measured at a lower f cost nothing again, which
        // is why this falls toward MC_TRIALS_PER_FLIP as the sweep proceeds.
        level_runs = real_runs - runs_before;

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
                << MC_TRIALS_PER_FLIP << "," << level_runs << ","
                << (mass_clean[k] / wsum) << "," << (mass_hit[k] / wsum) << ","
                << (mass_msb_ok[k] / wsum) << "," << (mass_msb_lost[k] / wsum) << ","
                // Raw layout counts overflow any integer type; log10 keeps them readable.
                << std::fixed << std::setprecision(4)
                << ((std::log(p) + lnorm) / std::log(10.0)) << "\n";
        }

        std::cout << std::fixed << std::setprecision(6)
                  << std::setw(7) << f << std::setw(10) << level_classes << "  "
                  << std::setw(9) << mean_frac << "  "
                  << std::setw(9) << (static_cast<double>(lo) / N) << "  "
                  << std::setw(9) << (static_cast<double>(hi) / N) << "  "
                  << std::setw(9) << (wsum_hit / wsum) << "  "
                  << std::setw(9) << worst_gap << "\n";

        EXPECT_GE(lo, 0);
        EXPECT_LE(hi, N);
        if (f == 0) {
            EXPECT_EQ(lo, TARGET);
            EXPECT_EQ(hi, TARGET);
        }
        // 20k random trials against an exact distribution: sampling error alone is of order
        // 0.003. A gap far past that would mean the enumeration and the gate disagree.
        EXPECT_LT(worst_gap, 0.02)
            << "random real trials disagree with the enumerated distribution at f=" << f;
    }

    const std::string fn = source_dir() + "uMUL_poster_faults.csv";
    std::ofstream f2(fn);
    ASSERT_TRUE(f2.is_open()) << "Failed to open " << fn;
    f2 << buf.str();               // single write
    f2.close();

    std::cout << "\n[PART 2] " << classes << " outcome rows over the 264-bit operand surface\n"
              << "         " << real_runs << " REAL 256-cycle gate simulations run to produce them\n"
              << "         worst enumerated-vs-random disagreement: " << mc_worst << "\n"
              << "         -> " << fn << std::endl;
    EXPECT_GT(classes, 0);
    EXPECT_GT(real_runs, 0);
}
