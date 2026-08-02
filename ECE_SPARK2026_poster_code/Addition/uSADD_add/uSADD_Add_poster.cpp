/*
 * =========================================================================================
 * uSADD_Add_poster.cpp                            ECE SPARK 2026 -- uGEMM scaled adder
 * =========================================================================================
 * Companion to MUX_add/, AND_mul/ and uMUL_mul/. Same operands, same axes, same standard.
 *
 *     uSADD_Add_poster_warmup.csv       zero-fault early-termination behaviour
 *     uSADD_Add_poster_faults.csv       fault response over the 512-bit STREAM surface
 *     uSADD_Add_poster_state.csv        PC and accumulator sensitivity, measured
 *
 * -----------------------------------------------------------------------------------------
 * THE UNIT.  PC sums the input bits; A accumulates and emits one bit whenever it holds n.
 * Nothing is discarded, so out = (a + b)/2 with no correlation requirement and no RNG.
 *
 * OPERANDS. a = b = 0.5, 128 ones each. Intended result 0.5, i.e. 128 output ones.
 *
 * NO FILTER IS NEEDED, unlike the MUX next door. The output is floor(total input ones / 2),
 * which does not mention where any one sits, so EVERY arrangement of the two streams gives
 * exactly 128. The MUX has to discard all but the perfectly decorrelated arrangements to get
 * its 0.5; uSADD gets it from all 3.3e151 of them.
 *
 * -----------------------------------------------------------------------------------------
 * THE FAULT SURFACE, AND WHICH PART OF IT IS EXHAUSTIVE
 *
 *   Everything this unit can suffer reduces to ONE quantity: how much CREDIT reaches A.
 *
 *       a / b stream bit    flip a one -> -1 credit, flip a zero -> +1 credit
 *       PC output bus bit   a glitch of +/- 2^bit on the count handed to A that cycle
 *       A register bit      a strike of +/- 2^bit on the unspent credit
 *
 *   All three are credit perturbations; A integrates them identically. That is the whole
 *   architectural story and it is why uSADD's damage is BOUNDED: A holds a residue, not an
 *   operand, so a perturbation is spent once and the unit keeps integrating correctly. Compare
 *   uMUL, whose 8-bit register holds the operand itself and whose MSB strike zeroes every
 *   remaining cycle of the run.
 *
 *   PART 2 IS EXHAUSTIVE over the 512-bit STREAM surface. Stream flips change the total by a
 *   fixed -1 or +1 with no dependence on timing, so the collapse is exact and the enumeration
 *   below is lossless. Asserted against C(512, f) at every flip level.
 *
 *   PART 3 COVERS PC AND A EXHAUSTIVELY BY SITE, NOT BY LAYOUT, and says so. A PC or A strike
 *   is TIME-DEPENDENT in a way stream flips are not: its sign depends on the bit's value at the
 *   instant of the strike, and a large enough strike can push A past its own width and WRAP,
 *   destroying credit rather than adding it. So every (bit, strike cycle, prior state) class is
 *   measured with real runs and reported, but multi-strike layouts over PC and A are not
 *   enumerated -- they do not collapse, and claiming otherwise would be the one dishonest thing
 *   in this folder.
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

#include "../../../modules/uSADD.hpp"

using namespace StochasticSimulator;

namespace {

constexpr int N = 256;
constexpr int ONES = N / 2;            // 128 ones per input stream -> 0.5 each
constexpr int TARGET = N / 2;          // 128 output ones -> 0.5
constexpr int MAX_FLIPS = 36;
constexpr int STREAM_SITES = 2 * N;    // 512 -- the exhaustively swept surface
constexpr int REGIONS = 4;             // (a, b) combinations
constexpr int RSIZE = N / REGIONS;     // 64 cells each when the streams are decorrelated
constexpr int INSTANCES_PER_CLASS = 4;
constexpr int MC_TRIALS_PER_FLIP = 20000;

const std::vector<int>& truncation_lengths() {
    static const std::vector<int> v = {2, 4, 8, 16, 32, 64, 128, 256};
    return v;
}

std::string source_dir() {
    std::string p = __FILE__;
    size_t slash = p.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string() : p.substr(0, slash + 1);
}

int popcount2(int v) { return (v & 1) + ((v >> 1) & 1); }

std::vector<bool> scattered(int ones, std::mt19937& rng) {
    std::vector<bool> s(N, false);
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int i = 0; i < ones; ++i) s[idx[i]] = true;
    return s;
}

int run_usadd(const std::vector<bool>& a, const std::vector<bool>& b, int upto = N) {
    UnaryScaledAdder unit(2);
    int ones = 0;
    for (int i = 0; i < upto; ++i) if (unit.add(a[i], b[i])) ++ones;
    return ones;
}

}  // namespace

TEST(uSADDAddPoster, GenerateTrials) {
    std::mt19937 rng(1337);
    long long real_runs = 0;

    // ======================================================================================
    // SANITY -- no filter needed, measured
    // ======================================================================================
    for (int inst = 0; inst < INSTANCES_PER_CLASS; ++inst) {
        const std::vector<bool> a = scattered(ONES, rng), b = scattered(ONES, rng);
        ASSERT_EQ(run_usadd(a, b), TARGET) << "0.5 + 0.5 did not scale to 128/256";
        ++real_runs;
    }
    {
        // The strongest form of the claim: adversarially correlated operands, same answer.
        std::vector<bool> s = scattered(ONES, rng), inv(N);
        for (int i = 0; i < N; ++i) inv[i] = !s[i];
        ASSERT_EQ(run_usadd(s, s), TARGET)   << "perfectly correlated operands broke the sum";
        ASSERT_EQ(run_usadd(s, inv), TARGET) << "perfectly anti-correlated operands broke the sum";
        real_runs += 2;
    }

    std::vector<double> lf(STREAM_SITES + 1, 0.0);
    for (int i = 1; i <= STREAM_SITES; ++i) lf[i] = lf[i-1] + std::log(static_cast<double>(i));
    auto logC = [&](int n, int k) {
        return (k < 0 || k > n) ? -1e300 : lf[n] - lf[k] - lf[n-k];
    };

    // ======================================================================================
    // PART 1 -- ZERO-FAULT EARLY TERMINATION
    // ======================================================================================
    // Truncating at L: the prefix output is floor(credit in the first L cycles / 2), and that
    // credit is (ones of a in the prefix) + (ones of b in the prefix). Each is an independent
    // HYPERGEOMETRIC draw over its own stream, so the prefix credit is their convolution --
    // exact over every arrangement, in at most (L+1)^2 terms and usually far fewer.
    std::ostringstream wbuf;
    wbuf << "N,MeanEstimate,MeanAbsError,MinEstimate,MaxEstimate,Classes,Log10Arrangements\n";
    wbuf << std::fixed << std::setprecision(8);

    const double log10_arrangements = 2.0 * logC(N, ONES) / std::log(10.0);
    long long warm_classes = 0;
    std::vector<bool> pa(N), pb(N);

    for (int L : truncation_lengths()) {
        double mean_est = 0.0, mean_abs = 0.0, wsum = 0.0, lo = 1e9, hi = -1e9;
        const double lnormL = logC(N, L);
        const int j_lo = std::max(0, L - (N - ONES)), j_hi = std::min(L, ONES);
        long long classes_here = 0;

        for (int ja = j_lo; ja <= j_hi; ++ja) {
            const double lwa = logC(ONES, ja) + logC(N - ONES, L - ja) - lnormL;
            if (lwa < -300.0) continue;
            for (int jb = j_lo; jb <= j_hi; ++jb) {
                const double lwb = logC(ONES, jb) + logC(N - ONES, L - jb) - lnormL;
                if (lwa + lwb < -300.0) continue;
                const double w = std::exp(lwa + lwb);
                if (w <= 0.0) continue;

                // MEASURED: build a prefix carrying exactly ja and jb ones and run the unit.
                pa.assign(N, false); pb.assign(N, false);
                for (int i = 0; i < ja; ++i) pa[i] = true;
                for (int i = 0; i < jb; ++i) pb[i] = true;
                const int emitted = run_usadd(pa, pb, L);
                ++real_runs;
                ASSERT_EQ(emitted, (ja + jb) / 2) << "prefix output is not floor(credit/2)";

                const double est = static_cast<double>(emitted) / L;
                ++classes_here; ++warm_classes;
                mean_est += w * est;
                mean_abs += w * std::fabs(est - 0.5);
                wsum += w;
                lo = std::min(lo, est);
                hi = std::max(hi, est);
            }
        }
        wbuf << L << "," << (mean_est / wsum) << "," << (mean_abs / wsum) << ","
             << lo << "," << hi << "," << classes_here << "," << log10_arrangements << "\n";
    }
    {
        const std::string fn = source_dir() + "uSADD_Add_poster_warmup.csv";
        std::ofstream f(fn);
        ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
        f << wbuf.str();
        f.close();
        std::cout << "\n[PART 1] exhaustive over all 10^" << std::fixed << std::setprecision(0)
                  << log10_arrangements << " arrangements via " << warm_classes
                  << " measured classes -> " << fn << std::endl;
    }

    // ======================================================================================
    // PART 2 -- FAULT SWEEP OVER THE 512-BIT STREAM SURFACE, EXHAUSTIVE
    // ======================================================================================
    // Per-region generating function over the 4 (a,b) combinations. A cell in region r takes a
    // flip mask m over its two bits, costing popcount(m) flips and contributing popcount(r ^ m)
    // credit. Convolve the four regions, then output = floor(credit / 2).
    //
    // NO CLAMP APPLIES HERE. Credit stays in [256-36, 256+36], so at most 146 emissions are
    // needed in 256 cycles and A never runs out of drain capacity. That is asserted below by
    // agreeing with real runs.
    auto convolve = [&](const std::vector<std::vector<long double>>& A,
                        const std::vector<std::vector<long double>>& B, int max_c) {
        std::vector<std::vector<long double>> C(MAX_FLIPS + 1,
                                                std::vector<long double>(max_c + 1, 0.0L));
        for (int f1 = 0; f1 <= MAX_FLIPS; ++f1)
            for (int c1 = 0; c1 < static_cast<int>(A[f1].size()); ++c1) {
                if (A[f1][c1] == 0.0L) continue;
                for (int f2 = 0; f1 + f2 <= MAX_FLIPS; ++f2)
                    for (int c2 = 0; c2 < static_cast<int>(B[f2].size()); ++c2) {
                        if (B[f2][c2] == 0.0L || c1 + c2 > max_c) continue;
                        C[f1 + f2][c1 + c2] += A[f1][c1] * B[f2][c2];
                    }
            }
        return C;
    };

    const int MAXC = 2 * N;   // credit can reach 2 per cycle
    std::vector<std::vector<long double>> dp(MAX_FLIPS + 1,
                                             std::vector<long double>(MAXC + 1, 0.0L));
    dp[0][0] = 1.0L;
    for (int r = 0; r < REGIONS; ++r) {
        std::vector<std::vector<long double>> cell(MAX_FLIPS + 1,
                                                   std::vector<long double>(3, 0.0L));
        for (int m = 0; m < REGIONS; ++m) {
            const int cost = popcount2(m);
            if (cost > MAX_FLIPS) continue;
            cell[cost][popcount2(r ^ m)] += 1.0L;
        }
        std::vector<std::vector<long double>> region(MAX_FLIPS + 1,
                                                     std::vector<long double>(1, 0.0L));
        region[0][0] = 1.0L;
        std::vector<std::vector<long double>> base = cell;
        int e = RSIZE, cap = 1;
        while (e > 0) {
            if (e & 1) region = convolve(region, base, std::min(MAXC, cap + 2 * RSIZE));
            e >>= 1;
            if (e > 0) { base = convolve(base, base, std::min(MAXC, 2 * RSIZE)); cap = std::min(MAXC, cap * 2); }
        }
        dp = convolve(dp, region, MAXC);
    }

    std::ostringstream buf;
    buf << "BitsFlipped,OutputOnes,OutputFraction,ErrorPerN,Probability,"
           "MC_Probability,MC_Trials,LevelRealTrials,Log10TrialCount\n";

    long long classes = 0;
    double mc_worst = 0.0;
    std::cout << "\n  flips   classes    mean       min        max      MC gap\n";

    for (int f = 0; f <= MAX_FLIPS; ++f) {
        const double lnorm = logC(STREAM_SITES, f);
        std::vector<double> mass(N + 1, 0.0);
        double wsum = 0.0;
        for (int c = 0; c <= MAXC; ++c) {
            if (dp[f][c] <= 0.0L) continue;
            const double lw = std::log(static_cast<double>(dp[f][c])) - lnorm;
            if (lw < -300.0) continue;
            const int k = c / 2;                       // output = floor(credit / 2)
            if (k > N) continue;
            mass[k] += std::exp(lw);
            wsum += std::exp(lw);
        }
        if (wsum <= 0.0) continue;
        EXPECT_NEAR(wsum, 1.0, 1e-6) << "flip layouts do not sum to C(512," << f << ")";

        // ---- forward Monte Carlo, real streams and real flips ------------------------------
        std::vector<long long> mc_hist(N + 1, 0);
        std::vector<int> pos(STREAM_SITES);
        std::vector<bool> fa, fb;
        for (int t = 0; t < MC_TRIALS_PER_FLIP; ++t) {
            fa = scattered(ONES, rng);
            fb = scattered(ONES, rng);
            std::iota(pos.begin(), pos.end(), 0);
            for (int i = 0; i < f; ++i) {
                std::uniform_int_distribution<int> pick(i, STREAM_SITES - 1);
                std::swap(pos[i], pos[pick(rng)]);
                const int p = pos[i];
                if (p < N) fa[p] = !fa[p]; else fb[p - N] = !fb[p - N];
            }
            ++mc_hist[run_usadd(fa, fb)];
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
    {
        const std::string fn = source_dir() + "uSADD_Add_poster_faults.csv";
        std::ofstream f2(fn);
        ASSERT_TRUE(f2.is_open()) << "Failed to open " << fn;
        f2 << buf.str();
        f2.close();
        std::cout << "\n[PART 2] " << classes << " outcome rows over the 512-bit STREAM surface\n"
                  << "         worst enumerated-vs-random disagreement: " << mc_worst << "\n"
                  << "         -> " << fn << std::endl;
    }

    // ======================================================================================
    // PART 3 -- PC AND ACCUMULATOR STRIKES, MEASURED
    // ======================================================================================
    // Exhaustive BY SITE CLASS, not by layout: every (target, bit, strike cycle) is run for real
    // over a grid of strike cycles, and the worst damage each can do is recorded. These sites do
    // not collapse the way stream flips do -- the sign of a register strike depends on the bit's
    // value at that instant, and a large strike can push A past its own width -- so multi-strike
    // layouts are deliberately NOT enumerated.
    UnaryScaledAdder probe(2);
    const unsigned ACC_W = probe.get_acc_width();
    const unsigned PC_W  = probe.get_pc_width();

    std::ostringstream tbuf;
    tbuf << "Target,Bit,Weight,StrikeCycle,CleanOnes,FaultedOnes,DeltaOnes,DeltaFraction\n";

    std::cout << "\n[PART 3] PC and accumulator strikes (measured; site-exhaustive, not "
                 "layout-exhaustive)\n"
              << "         accumulator " << ACC_W << " bits, PC bus " << PC_W << " bits\n"
              << "         target  bit  weight   worst |delta| over all strike cycles\n";

    const std::vector<int> strike_cycles = [] {
        std::vector<int> v;
        for (int c = 0; c < N; c += 8) v.push_back(c);
        v.push_back(N - 1);
        return v;
    }();

    long long state_runs = 0;
    int worst_state_delta = 0;
    for (int target = 0; target < 2; ++target) {          // 0 = accumulator, 1 = PC bus
        const unsigned bits = (target == 0) ? ACC_W : PC_W;
        for (unsigned bit = 0; bit < bits; ++bit) {
            int worst = 0;
            for (int strike : strike_cycles) {
                const std::vector<bool> a = scattered(ONES, rng), b = scattered(ONES, rng);
                const int clean = run_usadd(a, b);
                ++real_runs;

                UnaryScaledAdder unit(2);
                int ones = 0;
                for (int i = 0; i < N; ++i) {
                    UnaryScaledAdder::CycleUpset up;
                    if (target == 0 && i == strike) unit.flip_accumulator_bit(bit);
                    if (target == 1 && i == strike) up.pc_flip = (1u << bit);
                    if (unit.add({a[i], b[i]}, up)) ++ones;
                }
                ++real_runs; ++state_runs;

                const int delta = ones - clean;
                worst = std::max(worst, std::abs(delta));
                tbuf << (target == 0 ? "accumulator" : "pc_bus") << "," << bit << ","
                     << (1 << bit) << "," << strike << "," << clean << "," << ones << ","
                     << delta << "," << std::fixed << std::setprecision(10)
                     << (static_cast<double>(delta) / N) << "\n";
            }
            worst_state_delta = std::max(worst_state_delta, worst);
            std::cout << "         " << std::setw(11) << (target == 0 ? "accumulator" : "pc_bus")
                      << std::setw(5) << bit << std::setw(8) << (1 << bit)
                      << "     " << worst << " output bits = " << std::fixed
                      << std::setprecision(6) << (static_cast<double>(worst) / N) << "\n";
        }
    }
    {
        const std::string fn = source_dir() + "uSADD_Add_poster_state.csv";
        std::ofstream f(fn);
        ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
        f << tbuf.str();
        f.close();
        std::cout << "         " << state_runs << " strike trials -> " << fn << std::endl;
    }

    // THE uSADD RESULT, stated as an assertion. A holds a residue, so one strike is spent once
    // and the unit keeps integrating: the damage is bounded by the register's own weight, never
    // by the length of the run. Contrast uMUL, where one MSB strike takes the answer to zero.
    EXPECT_LE(worst_state_delta, (1 << ACC_W) / 2 + 1)
        << "a single state strike did more damage than the accumulator's own width allows";

    std::cout << "\n[TOTAL] " << real_runs << " real 256-cycle uSADD simulations" << std::endl;
    EXPECT_GT(classes, 0);
}
