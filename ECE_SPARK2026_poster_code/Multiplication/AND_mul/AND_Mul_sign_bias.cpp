/*
 * =========================================================================================
 * AND_Mul_sign_bias.cpp                              ECE SPARK 2026 -- AND multiplier
 * =========================================================================================
 * Generates AND_Mul_sign_bias.csv, which plot_AND_Mul_fault_sign_bias.m turns into the
 * fault-sign-bias figure. Lands in this folder regardless of where the binary is launched.
 *
 * -----------------------------------------------------------------------------------------
 * WHY THIS FILE EXISTS -- THE OLD PANEL WAS SAMPLED, AND IT LOOKED IT
 *
 *   The original "Graph 2: Fault Sign Bias" was drawn from
 *   32bit_exhaustive_multiplier_trials.csv, which sweeps only ORGANIZATIONS_PER_PAIR sampled
 *   arrangements per operand pair and only part of the flip range. The result wobbled around
 *   the trend by a few tenths of a bit -- and that wobble was sampling noise, not physics.
 *
 *   Done EXHAUSTIVELY the answer is provably a straight line, and this file computes it that
 *   way: every operand pair, every flip count from 0 to 32, and every combination of flips
 *   within each count, each carrying its exact multiplicity.
 *
 * -----------------------------------------------------------------------------------------
 * WHAT IS AVERAGED OVER
 *
 *     Count A      the ones in the faulted stream          0 .. 32   (the x axis)
 *     Count B      the ones in the clean stream            0 .. 32   averaged over
 *     f            how many of the 32 bits are flipped     0 .. 32   averaged over
 *     which bits   every C(32,f) combination at that f               averaged over, exactly
 *
 *   Faults hit stream A only, matching the original campaign.
 *
 * -----------------------------------------------------------------------------------------
 * THE COLLAPSE. An AND gate cannot see WHERE a bit is, so for a given operand pair the 32
 * cycles split into four regions and only two of them can move the output:
 *
 *     OV  A=1 B=1   ov = round(cA*cB/32)   flipping A REMOVES an output one
 *     AO  A=1 B=0   cA - ov                inert: B is already 0
 *     BO  A=0 B=1   cB - ov                flipping A ADDS an output one
 *     NN  A=0 B=0   the rest               inert
 *
 *   So with x flips landing in OV and z landing in BO,
 *
 *       output ones = ov - x + z        and        ones delta = z - x
 *
 *   That primitive is MEASURED, not assumed: it is checked against a real Multiplier over a
 *   grid of constructed cases before any averaging happens.
 *
 *   x and z are then hypergeometric over the 32 positions, and this file sums their
 *   distributions term by term rather than using a closed form -- every reachable (x, z) is
 *   enumerated with its exact C(ov,x)C(32-ov,f-x)/C(32,f) weight.
 *
 * -----------------------------------------------------------------------------------------
 * THE RESULT IS EXACTLY LINEAR, and it is worth knowing why before looking at the figure.
 * A flip at a B=1 position adds a one if A was 0 there and removes one if A was 1, so the
 * expected delta per flipped position is proportional to (1 - 2*cA/32). Averaged over f and
 * over Count B that gives
 *
 *       mean ones delta  =  8 - cA/2
 *
 * running from +8 at cA = 0 to -8 at cA = 32 and crossing zero at cA = 16 -- exactly half
 * density, where a flip is equally likely to add or remove. The CSV is checked against that
 * line at the end, and the residual is reported.
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

#include "../../../modules/multiplier.hpp"

using namespace StochasticSimulator;

namespace {

constexpr int N = 32;              // stream length of the original 32-bit campaign

std::string source_dir() {
    std::string p = __FILE__;
    size_t slash = p.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string() : p.substr(0, slash + 1);
}

// The zero-correlation-error overlap the original campaign builds: the ones-count the output
// would have if the two streams were perfectly decorrelated.
int zce_overlap(int cA, int cB) {
    int ov = static_cast<int>(std::lround(static_cast<double>(cA) * cB / N));
    ov = std::max(ov, std::max(0, cA + cB - N));   // cannot be smaller than forced overlap
    ov = std::min(ov, std::min(cA, cB));           // nor larger than either operand
    return ov;
}

}  // namespace

TEST(AndMulSignBias, ExhaustiveOverEveryFaultCombination) {
    Multiplier gate;
    std::mt19937 rng(1337);
    long long real_runs = 0;

    // log-factorials for exact hypergeometric weights
    std::vector<double> lf(N + 1, 0.0);
    for (int i = 1; i <= N; ++i) lf[i] = lf[i-1] + std::log(static_cast<double>(i));
    auto logC = [&](int n, int k) {
        return (k < 0 || k > n) ? -1e300 : lf[n] - lf[k] - lf[n-k];
    };

    // =====================================================================================
    // STEP 1 -- MEASURE THE PRIMITIVE
    // =====================================================================================
    // Verify on real streams that a pair with regions (ov, AO, BO, NN), with x flips placed in
    // OV and z in BO, comes out at ov - x + z. Everything below rests on this, so it is checked
    // before it is used rather than after.
    std::cout << "\n[STEP 1] verifying `ones delta = z - x` against a real AND gate\n";
    int checks = 0;
    for (int cA : {0, 5, 11, 16, 21, 27, 32}) {
        for (int cB : {0, 7, 16, 23, 32}) {
            const int ov = zce_overlap(cA, cB);
            const int AO = cA - ov, BO = cB - ov, NN = N - ov - AO - BO;
            ASSERT_GE(NN, 0) << "region sizes do not partition the stream at cA=" << cA
                             << " cB=" << cB;

            for (int x = 0; x <= std::min(ov, 3); ++x) {
                for (int z = 0; z <= std::min(BO, 3); ++z) {
                    // Build the pair in canonical region order, then shuffle the positions so
                    // the check is not accidentally passing on one convenient layout.
                    std::vector<int> slot(N);
                    std::iota(slot.begin(), slot.end(), 0);
                    std::shuffle(slot.begin(), slot.end(), rng);

                    std::vector<bool> A(N, false), B(N, false);
                    int p = 0;
                    std::vector<int> ovPos, boPos;
                    for (int i = 0; i < ov; ++i, ++p) { A[slot[p]] = true; B[slot[p]] = true;
                                                        ovPos.push_back(slot[p]); }
                    for (int i = 0; i < AO; ++i, ++p) { A[slot[p]] = true; }
                    for (int i = 0; i < BO; ++i, ++p) { B[slot[p]] = true;
                                                        boPos.push_back(slot[p]); }

                    int clean = 0;
                    for (int i = 0; i < N; ++i) if (gate.multiply(A[i], B[i])) ++clean;
                    ASSERT_EQ(clean, ov) << "clean overlap is not ov at cA=" << cA << " cB=" << cB;

                    // Apply the flips to stream A only, as the original campaign does.
                    std::vector<bool> Af = A;
                    for (int i = 0; i < x; ++i) Af[ovPos[i]] = !Af[ovPos[i]];
                    for (int i = 0; i < z; ++i) Af[boPos[i]] = !Af[boPos[i]];

                    int faulted = 0;
                    for (int i = 0; i < N; ++i) if (gate.multiply(Af[i], B[i])) ++faulted;
                    real_runs += 2;
                    ++checks;

                    ASSERT_EQ(faulted - clean, z - x)
                        << "measured ones delta disagrees with (z - x) at cA=" << cA
                        << " cB=" << cB << " x=" << x << " z=" << z;
                }
            }
        }
    }
    std::cout << "         " << checks << " region/flip configurations, all agreeing ("
              << real_runs << " real gate runs)\n";

    // =====================================================================================
    // STEP 2 -- EXHAUSTIVE ENUMERATION, INCLUDING OVER ARRANGEMENTS
    // =====================================================================================
    // AVERAGE OVER EVERY ARRANGEMENT, NOT ONE CANONICAL OVERLAP. An earlier version of this
    // file used the original campaign's zce_overlap() -- round(cA*cB/32) -- as THE overlap for
    // each operand pair. That is one specific layout, and forcing it to a whole number of ones
    // puts a systematic tilt into the answer: lround sends every half-integer up, so every odd
    // Count B is biased by -1, and the average inherits
    //         16 odd values out of 33, times E[f]/32  =  16/(2*33)  =  0.2424
    // of spurious negative bias at cA = 16. That is exactly the deviation the check below
    // caught, and it is also what bent the line.
    //
    // The real quantity is the average over ALL arrangements of the two streams. The overlap is
    // then not a chosen number but a HYPERGEOMETRIC random variable in its own right, and no
    // rounding happens anywhere. Enumerating it is both more exhaustive and exactly linear.
    //
    // Ehit[R][f] = expected flips landing in a region of size R when f flips are spread over
    // the 32 positions, summed term by term over its hypergeometric distribution. It depends
    // only on (R, f), so 33 x 33 entries cover every case the sweep needs.
    std::vector<std::vector<double>> Ehit(N + 1, std::vector<double>(N + 1, 0.0));
    long long classes = 0;
    for (int R = 0; R <= N; ++R) {
        for (int f = 0; f <= N; ++f) {
            const double lnorm = logC(N, f);
            double e = 0.0;
            for (int k = std::max(0, f - (N - R)); k <= std::min(f, R); ++k) {
                const double lw = logC(R, k) + logC(N - R, f - k) - lnorm;
                if (lw < -300.0) continue;
                e += k * std::exp(lw);
                ++classes;
            }
            Ehit[R][f] = e;
        }
    }

    std::ostringstream buf;
    buf << "CountA,BitsFlipped,MeanOnesDelta,Log10Combinations\n";
    std::vector<double> perCountA(N + 1, 0.0);

    for (int cA = 0; cA <= N; ++cA) {
        double acc_over_f = 0.0;
        for (int f = 0; f <= N; ++f) {
            double acc_over_cB = 0.0;

            for (int cB = 0; cB <= N; ++cB) {
                // Enumerate the overlap over every arrangement: choosing which cB of the 32
                // positions carry B's ones, ov of them landing on A's ones.
                const double lnormB = logC(N, cB);
                double delta = 0.0;
                for (int ov = std::max(0, cA + cB - N); ov <= std::min(cA, cB); ++ov) {
                    const double lw = logC(cA, ov) + logC(N - cA, cB - ov) - lnormB;
                    if (lw < -300.0) continue;
                    const double w = std::exp(lw);
                    ++classes;
                    // flips into BO add a one, flips into OV remove one
                    delta += w * (Ehit[cB - ov][f] - Ehit[ov][f]);
                }
                acc_over_cB += delta;
            }

            const double mean_f = acc_over_cB / (N + 1);   // average over Count B
            acc_over_f += mean_f;

            buf << cA << "," << f << ","
                << std::fixed << std::setprecision(10) << mean_f << ","
                << std::setprecision(4) << (logC(N, f) / std::log(10.0)) << "\n";
        }
        perCountA[cA] = acc_over_f / (N + 1);              // average over flip count
    }

    // =====================================================================================
    // STEP 3 -- CHECK IT AGAINST THE LINE IT SHOULD BE
    // =====================================================================================
    // Expected delta per flipped B=1 position is proportional to (1 - 2*cA/32); averaged over f
    // and Count B that is 8 - cA/2. If the enumeration is right this is exact, not a fit.
    double worst = 0.0;
    for (int cA = 0; cA <= N; ++cA) {
        worst = std::max(worst, std::fabs(perCountA[cA] - (8.0 - cA / 2.0)));
    }

    std::cout << "\n[STEP 2] exhaustive: every Count A and Count B (0.." << N << "), every flip\n"
              << "         count 0.." << N << ", every combination within each count\n"
              << "         " << classes << " enumerated outcome terms\n"
              << "\n[STEP 3] mean ones delta vs the exact line 8 - cA/2\n"
              << "         cA      measured        8 - cA/2\n";
    for (int cA : {0, 8, 15, 16, 17, 24, 32}) {
        std::cout << "         " << std::setw(3) << cA << std::fixed << std::setprecision(6)
                  << std::setw(14) << perCountA[cA]
                  << std::setw(14) << (8.0 - cA / 2.0) << "\n";
    }
    std::cout << "         worst deviation from the line: " << std::scientific
              << std::setprecision(3) << worst << std::endl;

    // Rounding the ZCE overlap to a whole number of ones perturbs this very slightly; anything
    // beyond that would mean the enumeration is wrong.
    EXPECT_LT(worst, 0.05) << "the exhaustive mean is not the line it must be";

    const std::string fn = source_dir() + "AND_Mul_sign_bias.csv";
    std::ofstream f(fn);
    ASSERT_TRUE(f.is_open()) << "Failed to open " << fn;
    f << buf.str();
    f.close();
    std::cout << "         -> " << fn << std::endl;
}
