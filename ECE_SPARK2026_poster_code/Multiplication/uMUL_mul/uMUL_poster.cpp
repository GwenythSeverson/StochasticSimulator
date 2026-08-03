/*
 * =========================================================================================
 * uMUL_poster.cpp                                      ECE SPARK 2026 -- uMUL multiplier
 * =========================================================================================
 * Twin of AND_mul/AND_Mul_poster.cpp, same operands, same axes, same exhaustiveness standard,
 * so the two figures can be read side by side. Generates the two CSVs that
 * plot_uMUL_poster.m turns into the poster figures. Both land in this folder, next to this
 * source file, regardless of where the test binary is launched from.
 *
 *     uMUL_poster_warmup.csv           zero-fault early-termination behaviour
 *     uMUL_poster_sensitivity.csv      exact single-flip delta of every operand site
 *     uMUL_poster_faults_stream.csv    fault response, REGISTER HARDENED, f = 0..36 over 256 sites
 *     uMUL_poster_faults_all.csv       fault response, NOTHING HARDENED,  f = 0..36 over 264 sites
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
 * THE FAULT SURFACE COVERED HERE -- 264 BITS, SWEPT AS TWO SEPARATE CAMPAIGNS
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
 *   PART 2 SWEEPS THE SURFACE TWICE, AS A BEFORE-AND-AFTER ON HARDENING THE REGISTER:
 *
 *       2A  register hardened   256 sites, stream only.   Answer bounded to 0.25 +/- f/256.
 *       2B  nothing hardened    264 sites, flips spread UNIFORMLY over stream and register.
 *                               Answer reaches a hard 0.000 from ONE upset, at every flip count.
 *
 *   Same operand, same flip counts, same x axis -- the ONLY difference between the two panels is
 *   whether the 8 register bits are in the target set. That is what makes them a fair pair, and
 *   it is the poster's argument for selective hardening stated as an experiment.
 *
 *   Each is a real independent campaign with its own Monte Carlo and its own trial count; they
 *   share no trials. An earlier revision instead reported ONE joint sweep split into "register
 *   survived" / "register hit" slices, so both figures quoted the same trial count and the
 *   damaged slice held only 8/264 of the mass at f=1.
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
    // PART 2 -- TWO SEPARATE FAULT CAMPAIGNS, NOT ONE CAMPAIGN SPLIT IN TWO
    // ======================================================================================
    // THIS REPLACED A SINGLE JOINT SWEEP, and the reason is a reporting problem, not a physics
    // problem. The old version placed f flips anywhere on the combined 264-bit surface and then
    // sorted the resulting outcomes into "no register bit was hit" and "at least one was",
    // shipping both as separate figures. Both figures were therefore CONDITIONAL SLICES of one
    // population, and both quoted that population's full trial count -- 806,376 real gate runs
    // on each panel, for two panels that between them used those same 806,376 runs once. Worse,
    // the register-hit slice thinned out badly at low f: at f=1 only 8 of 264 layouts land on
    // the register, so a panel drawn from ~3% of the mass was being presented on equal footing
    // with one drawn from the other 97%.
    //
    // The fix is to stop conditioning and start running two independent experiments, each with
    // its own surface, its own flip range, its own Monte Carlo, and its own honest trial count:
    //
    //   2A  STREAM ONLY     surface = the 256 in_0 bits.   The value register is never touched.
    //                       f = 0..36. Answers the question "what does uMUL do when the part it
    //                       shares with the AND gate gets hit?", so it is the panel that is
    //                       directly comparable to the AND gate's fault figure.
    //
    //   2B  REGISTER ONLY   surface = the 8 value-register bits. The stream is never touched.
    //                       f = 0..8 AND NO FURTHER -- there are only 8 bits, so a 9th flip does
    //                       not exist. This sweep is LITERALLY exhaustive: all 2^8 = 256 register
    //                       states are visited, C(8,f) of them at each f, no weighting needed
    //                       beyond uniform-over-subsets.
    //
    // Neither sweep is a subset of the other and neither borrows the other's trials.
    //
    // WHAT IS GIVEN UP BY SPLITTING. The joint sweep answered "if f upsets land somewhere in the
    // operand storage, what happens?" -- which is the right question for a part sitting in a
    // radiation field, where the attacker does not respect module boundaries. These two sweeps
    // answer "what happens if the f upsets land HERE", which is the right question for deciding
    // WHAT TO HARDEN. The joint result is still recoverable from these two by weighting them
    // 256/264 and 8/264 per flip; it is just no longer the thing plotted.
    std::vector<bool> fs(N);

    // --------------------------------------------------------------------------------------
    // 2A -- STREAM-ONLY SWEEP, 0..36 FLIPS OVER 256 SITES
    // --------------------------------------------------------------------------------------
    // f flips land on the stream and nowhere else, so value' = value = 128 throughout and the
    // only thing that moves is the ones-count:
    //
    //     a  flips on the 128 positions holding a one   -> that cycle stops enabling
    //     b  flips on the 128 positions holding a zero  -> that cycle starts enabling
    //        with a + b = f, so n1' = 128 - a + b
    //
    // Layouts in the class = C(128,a) * C(128,b), and Vandermonde gives the guarantee that the
    // decomposition is complete and non-overlapping:
    //
    //     SUM_a C(128,a) * C(128,f-a) = C(256,f)
    //
    // asserted at every flip level, in log space. This is the AND gate's decomposition with one
    // stream instead of two, which is exactly why the two figures can be read against each other.
    std::ostringstream sbuf2;
    sbuf2 << "BitsFlipped,OutputOnes,OutputFraction,ErrorPerN,Probability,"
             "MC_Probability,MC_Trials,LevelRealTrials,Log10TrialCount\n";

    long long stream_classes = 0;
    double stream_mc_worst = 0.0;
    const long long stream_runs_start = real_runs;
    std::cout << "\n[2A] STREAM ONLY -- 256 sites, register never touched\n"
              << "  flips   classes    mean       min        max      MC gap\n";

    for (int f = 0; f <= MAX_FLIPS; ++f) {
        std::vector<double> mass(N + 1, 0.0);
        const double lnorm = logC(N, f);       // f flips chosen from the 256 stream sites
        double wsum = 0.0;
        const long long runs_before = real_runs;
        long long level_classes = 0;

        for (int a = 0; a <= std::min(ONES, f); ++a) {
            const int b = f - a;
            if (b > N - ONES) continue;
            const int n1 = ONES - a + b;

            const double lw = logC(ONES, a) + logC(N - ONES, b) - lnorm;
            if (lw < -300.0) continue;
            const double w = std::exp(lw);
            if (w <= 0.0) continue;

            const int k = G(VALUE, n1);        // <- the real gate run (memoised)
            ++level_classes;
            mass[k] += w;
            wsum += w;
        }
        if (wsum <= 0.0) continue;

        EXPECT_NEAR(wsum, 1.0, 1e-9)
            << "stream flip layouts do not sum to C(256," << f << ") at f=" << f;

        // ---- FORWARD MONTE CARLO, STREAM SITES ONLY ---------------------------------------
        // Fresh random in_0 with 128 ones, f of the 256 STREAM positions picked uniformly, each
        // toggled for real, clocked through a real unit. The register is left alone, which is
        // the entire difference from sweep 2B below.
        std::vector<long long> mc_hist(N + 1, 0);
        std::vector<int> pos(N);
        for (int t = 0; t < MC_TRIALS_PER_FLIP; ++t) {
            build_stream(fs, ONES, rng);
            std::iota(pos.begin(), pos.end(), 0);
            for (int i = 0; i < f; ++i) {
                std::uniform_int_distribution<int> pick(i, N - 1);
                std::swap(pos[i], pos[pick(rng)]);
                fs[pos[i]] = !fs[pos[i]];
            }
            ++mc_hist[run_gate(fs, VALUE)];
            ++real_runs;
        }

        double worst_gap = 0.0;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0 && mc_hist[k] == 0) continue;
            const double pe = mass[k] / wsum;
            const double pm = static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP;
            worst_gap = std::max(worst_gap, std::fabs(pe - pm));
        }
        stream_mc_worst = std::max(stream_mc_worst, worst_gap);

        const long long level_runs = real_runs - runs_before;

        double mean_frac = 0.0;
        int lo = -1, hi = -1;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0) continue;
            const double p = mass[k] / wsum;
            const double frac = static_cast<double>(k) / N;
            mean_frac += p * frac;
            if (lo < 0) lo = k;
            hi = k;
            ++stream_classes;

            sbuf2 << f << "," << k << ","
                  << std::fixed << std::setprecision(10) << frac << ","
                  << (static_cast<double>(k - TARGET) / N) << ","
                  << std::scientific << std::setprecision(12) << p << ","
                  << (static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP) << ","
                  << MC_TRIALS_PER_FLIP << "," << level_runs << ","
                  << std::fixed << std::setprecision(4)
                  << ((std::log(p) + lnorm) / std::log(10.0)) << "\n";
        }

        std::cout << std::fixed << std::setprecision(6)
                  << std::setw(7) << f << std::setw(10) << level_classes << "  "
                  << std::setw(9) << mean_frac << "  "
                  << std::setw(9) << (static_cast<double>(lo) / N) << "  "
                  << std::setw(9) << (static_cast<double>(hi) / N) << "  "
                  << std::setw(9) << worst_gap << "\n";

        EXPECT_GE(lo, 0);
        EXPECT_LE(hi, N);
        if (f == 0) { EXPECT_EQ(lo, TARGET); EXPECT_EQ(hi, TARGET); }
        // With the register protected the answer CANNOT leave 0.25 +/- f/256. Asserting that is
        // the point of the panel: it is the bound the AND gate also obeys.
        EXPECT_GE(lo, TARGET - f) << "stream-only fault fell further than f/256 at f=" << f;
        EXPECT_LE(hi, TARGET + f) << "stream-only fault rose further than f/256 at f=" << f;
        EXPECT_LT(worst_gap, 0.02)
            << "random real trials disagree with the enumerated distribution at f=" << f;
    }

    const long long stream_real_runs = real_runs - stream_runs_start;
    {
        const std::string fn = source_dir() + "uMUL_poster_faults_stream.csv";
        std::ofstream f2(fn);
        ASSERT_TRUE(f2.is_open()) << "Failed to open " << fn;
        f2 << sbuf2.str();          // single write
        f2.close();
        std::cout << "     " << stream_classes << " outcome rows, "
                  << stream_real_runs << " REAL gate simulations, worst MC gap "
                  << stream_mc_worst << "\n     -> " << fn << std::endl;
    }

    // --------------------------------------------------------------------------------------
    // 2B -- WHOLE-SURFACE SWEEP, 0..36 FLIPS SPREAD UNIFORMLY OVER ALL 264 SITES
    // --------------------------------------------------------------------------------------
    // THE PAIR 2A/2B IS A BEFORE-AND-AFTER ON ONE DESIGN DECISION: is the value register
    // protected or not?
    //
    //     2A  register hardened   -> 256-bit surface, answer bounded to 0.25 +/- f/256
    //     2B  nothing hardened    -> 264-bit surface, answer can hit 0.000 from ONE upset
    //
    // Both are real independent campaigns with their own Monte Carlo and their own trial count.
    // 2B is the one that describes a part actually sitting in a radiation field, because an upset
    // does not respect module boundaries -- it lands wherever it lands, and 8 of the 264 places it
    // can land are catastrophic. 2A is the counterfactual you get by spending area on hardening.
    //
    //     (An earlier revision made 2B a register-ONLY sweep, f = 0..8. That answered "what if
    //      you deliberately aim f upsets at the register", which is a fault-injection experiment
    //      rather than a radiation one, and it could not show the thing this panel exists to
    //      show: that spreading damage EVENLY across the whole unit still reaches zero. It also
    //      had the odd property that damage FELL as flips rose, because value' = 0 needs the
    //      subset {bit 7} alone. Uniform spreading removes that artifact -- at any f there are
    //      plenty of layouts that strike bit 7 and put the rest on the stream.)
    //
    // THE ENUMERATION. f flips are placed on 264 sites. Split them by where they land:
    //
    //     r  flips on the 8 register bits, choosing a specific SUBSET S -- C(8,r) subsets,
    //        each giving its own faulted operand  value' = 128 XOR mask(S)
    //     a  flips on the 128 stream positions holding a one   -> that cycle stops enabling
    //     b  flips on the 128 stream positions holding a zero  -> that cycle starts enabling
    //        with r + a + b = f, so n1' = 128 - a + b
    //
    // The faulted output is G(value', n1'), measured. Layouts in the class = C(128,a)*C(128,b),
    // and the sum checks out exactly by Vandermonde -- the arithmetic guarantee that nothing is
    // double counted and nothing is missed:
    //
    //     SUM_r C(8,r) * SUM_a C(128,a)*C(128,f-r-a) = SUM_r C(8,r)*C(256,f-r) = C(264,f)
    //
    // ASSERTED at every flip level below, in log space.
    //
    // TWO DECOMPOSITIONS SHIP IN THE CSV because between them they are the whole explanation:
    //     clean / hit   was ANY of the 8 register bits struck?  P(hit) = 1 - C(256,f)/C(264,f)
    //     MSB intact /  was BIT 7 specifically struck? Not the same question, and this is the one
    //     MSB lost      that matters. value = 128 has exactly ONE bit set, so
    //                       bit 7 survives -> value' in [128,255] -> answer >= 0.25
    //                       bit 7 dies     -> value' in [0,127]   -> answer <= 0.25
    //                   and P(bit 7 struck) = f/264 exactly -- ONE site out of 264, NOT the
    //                   register-hit mass. Only about 13-20% of register hits go DOWN; the rest
    //                   push the answer up.
    std::ostringstream rbuf;
    rbuf << "BitsFlipped,OutputOnes,OutputFraction,ErrorPerN,Probability,"
            "MC_Probability,MC_Trials,LevelRealTrials,ProbRegisterClean,ProbRegisterHit,"
            "ProbMSBIntact,ProbMSBLost,Log10TrialCount\n";

    long long all_classes = 0;
    double all_mc_worst = 0.0;
    const long long all_runs_start = real_runs;
    std::cout << "\n[2B] WHOLE SURFACE -- 264 sites, flips spread evenly over stream + register\n"
              << "  flips   classes    mean       min        max    P(reg hit)  P(MSB hit)  MC gap\n";

    for (int f = 0; f <= MAX_FLIPS; ++f) {
        std::vector<double> mass(N + 1, 0.0), mass_clean(N + 1, 0.0), mass_hit(N + 1, 0.0),
                            mass_msb_ok(N + 1, 0.0), mass_msb_lost(N + 1, 0.0);
        const double lnorm = logC(SITES, f);     // f flips chosen from 264 sites
        double wsum = 0.0, wsum_hit = 0.0, wsum_msb = 0.0;
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
                    if (S & (1u << (REG_BITS - 1))) { mass_msb_lost[k] += w; wsum_msb += w; }
                    else                              mass_msb_ok[k]   += w;
                }
            }
        }
        if (wsum <= 0.0) continue;

        // The Vandermonde check. Every layout of f flips over 264 sites is accounted for exactly
        // once, or this fails. Tolerance is generous only because the sum is built in log space.
        EXPECT_NEAR(wsum, 1.0, 1e-9)
            << "flip layouts do not sum to C(264," << f << ") at f=" << f;
        // P(bit 7 struck) = f/264 exactly -- one site out of 264, asserted rather than asserted
        // in a comment. This is the mass of the low mode, and it is NOT the register-hit mass.
        EXPECT_NEAR(wsum_msb, static_cast<double>(f) / SITES, 1e-9)
            << "P(bit 7 struck) is not f/264 at f=" << f;

        // ---- FORWARD MONTE CARLO OVER GENUINELY RANDOM TRIALS ------------------------------
        // Draw a fresh random in_0 with 128 ones, pick f of the 264 sites uniformly at random,
        // apply them for real -- toggling stream bits, XORing register bits -- and clock the
        // whole thing through. UNIFORM over all 264: nothing steers flips toward the register,
        // which is exactly what makes this the "spread evenly" campaign.
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
        }

        double worst_gap = 0.0;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0 && mc_hist[k] == 0) continue;
            const double pe = mass[k] / wsum;
            const double pm = static_cast<double>(mc_hist[k]) / MC_TRIALS_PER_FLIP;
            worst_gap = std::max(worst_gap, std::fabs(pe - pm));
        }
        all_mc_worst = std::max(all_mc_worst, worst_gap);

        const long long level_runs = real_runs - runs_before;

        double mean_frac = 0.0;
        int lo = -1, hi = -1;
        for (int k = 0; k <= N; ++k) {
            if (mass[k] <= 0.0) continue;
            const double p = mass[k] / wsum;
            const double frac = static_cast<double>(k) / N;
            mean_frac += p * frac;
            if (lo < 0) lo = k;
            hi = k;
            ++all_classes;

            rbuf << f << "," << k << ","
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
                  << std::setw(10) << (wsum_msb / wsum) << "  "
                  << std::setw(9) << worst_gap << "\n";

        EXPECT_GE(lo, 0);
        EXPECT_LE(hi, N);
        if (f == 0) { EXPECT_EQ(lo, TARGET); EXPECT_EQ(hi, TARGET); }

        // THE CLAIM THE PANEL EXISTS TO MAKE, AND IT HOLDS AT EVERY FLIP COUNT FROM ONE.
        // Spreading damage evenly over the whole unit still reaches a hard zero, because some
        // layout at every f puts a flip on bit 7 and the rest anywhere. Contrast 2A, where the
        // register is protected and the floor is asserted to stay above 0.25 - f/256.
        if (f >= 1) {
            EXPECT_EQ(lo, 0)
                << "a whole-surface fault at f=" << f << " could not zero the unit, but the "
                << "layout striking bit 7 exists at every flip count";
        }
        EXPECT_LT(worst_gap, 0.02)
            << "random real trials disagree with the enumerated distribution at f=" << f;
    }

    const long long all_real_runs = real_runs - all_runs_start;
    {
        const std::string fn = source_dir() + "uMUL_poster_faults_all.csv";
        std::ofstream f2(fn);
        ASSERT_TRUE(f2.is_open()) << "Failed to open " << fn;
        f2 << rbuf.str();           // single write
        f2.close();
        std::cout << "     " << all_classes << " outcome rows, "
                  << all_real_runs << " REAL gate simulations, worst MC gap "
                  << all_mc_worst << "\n     -> " << fn << std::endl;
    }

    std::cout << "\n[PART 2] two INDEPENDENT campaigns, no shared trials:\n"
              << "         register hardened   256 sites, f=0..36, " << stream_real_runs
              << " real gate simulations   -> figure 2\n"
              << "         nothing hardened    264 sites, f=0..36, " << all_real_runs
              << " real gate simulations   -> figure 3\n"
              << "         grand total for this file: " << real_runs << "\n";
    EXPECT_GT(stream_classes, 0);
    EXPECT_GT(all_classes, 0);
    EXPECT_GT(real_runs, 0);
}
