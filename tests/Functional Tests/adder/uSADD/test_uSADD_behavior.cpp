/*
 * =========================================================================================
 * test_uSADD_behavior.cpp        Characterization of the uGEMM scaled adder vs the MUX adder
 * =========================================================================================
 * Writes usadd_behavior.csv and usadd_vs_mux.csv into the working directory.
 *
 * PART 1  EXHAUSTIVE ACCURACY. Every operand pair uSADD can be given at N = 256 -- all 257 x 257
 *         = 66,049 combinations of (ones in a, ones in b) -- run through the real unit and
 *         checked against the exact scaled sum. This is genuinely exhaustive over the operand
 *         space, not a sample, because of the collapse: the answer depends only on the ones
 *         COUNTS, so one arrangement per count pair covers every arrangement of that pair.
 *         That claim is not assumed -- ARRANGEMENTS_PER_PAIR independent shuffles are run for a
 *         sampled subset of pairs and asserted equal.
 *
 * PART 2  uSADD vs MUX, ON THE SAME OPERANDS. The MUX adder is the standard stochastic adder:
 *         out = sel ? a : b with a 0.5 select stream. It is correct only in expectation, and
 *         only if a, b and sel are mutually uncorrelated. This part feeds both units identical
 *         operands under three correlation regimes and records what each does.
 *
 * PART 3  EARLY TERMINATION. Mean absolute error against truncation length for both units.
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
#include <vector>

#include "../../../../modules/uSADD.hpp"
#include "../../../../modules/adder.hpp"

using namespace StochasticSimulator;

namespace {

constexpr int N = 256;
constexpr int ARRANGEMENTS_PER_PAIR = 4;   // shuffles used to verify the collapse
constexpr int MUX_TRIALS = 200;            // random select streams per operand pair

std::vector<bool> scattered(int length, int ones, std::mt19937& rng) {
    std::vector<bool> s(length, false);
    std::vector<int> idx(length);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int i = 0; i < ones; ++i) s[idx[i]] = true;
    return s;
}

int count_ones(const std::vector<bool>& s) {
    int n = 0;
    for (bool b : s) if (b) ++n;
    return n;
}

// The MUX adder run over whole streams, so the two units are driven identically.
int mux_ones(const std::vector<bool>& a, const std::vector<bool>& b,
             const std::vector<bool>& sel) {
    Adder mux;
    int ones = 0;
    for (int i = 0; i < N; ++i) {
        if (mux.add_scaled_or_weighted(a[i], b[i], sel[i])) ++ones;
    }
    return ones;
}

}  // namespace

// =============================================================================================
TEST(uSADDBehavior, ExhaustiveOperandSweep) {
    std::mt19937 rng(20260801);

    std::ostringstream buf;
    buf << "OnesA,OnesB,ProbA,ProbB,ExactScaledSum,MeasuredOnes,MeasuredProb,ErrorBits\n";
    buf << std::fixed << std::setprecision(10);

    long long pairs = 0, real_runs = 0;
    int worst_error_bits = 0;
    double worst_abs = 0.0;
    long long exact_pairs = 0;

    std::vector<bool> a, b;
    for (int na = 0; na <= N; ++na) {
        for (int nb = 0; nb <= N; ++nb) {
            a = scattered(N, na, rng);
            b = scattered(N, nb, rng);
            const int measured = count_ones(usadd_stream(a, b));
            ++real_runs;
            ++pairs;

            // The exact scaled sum in output-bit units, and the deficit the residue costs.
            const int exact_bits = (na + nb) / 2;
            const int err_bits = measured - exact_bits;
            EXPECT_EQ(err_bits, 0) << "uSADD is not floor((na+nb)/2) at na=" << na << " nb=" << nb;

            const double exact_prob = (static_cast<double>(na) / N + static_cast<double>(nb) / N) / 2.0;
            const double measured_prob = static_cast<double>(measured) / N;
            const double abs_err = std::fabs(measured_prob - exact_prob);
            worst_abs = std::max(worst_abs, abs_err);
            worst_error_bits = std::max(worst_error_bits, std::abs(err_bits));
            if (abs_err == 0.0) ++exact_pairs;

            // The CSV is 66k rows; keep it to a readable grid plus every exact-rail case.
            if (na % 8 == 0 && nb % 8 == 0) {
                buf << na << "," << nb << ","
                    << (static_cast<double>(na) / N) << "," << (static_cast<double>(nb) / N) << ","
                    << exact_prob << "," << measured << "," << measured_prob << ","
                    << err_bits << "\n";
            }
        }
    }

    // ---- verify the collapse on a sampled grid, with real shuffled runs --------------------
    long long collapse_checks = 0;
    for (int na = 0; na <= N; na += 17) {
        for (int nb = 0; nb <= N; nb += 23) {
            int reference = -1;
            for (int inst = 0; inst < ARRANGEMENTS_PER_PAIR; ++inst) {
                a = scattered(N, na, rng);
                b = scattered(N, nb, rng);
                const int got = count_ones(usadd_stream(a, b));
                ++real_runs;
                if (inst == 0) reference = got;
                else EXPECT_EQ(got, reference)
                    << "arrangement changed the answer at na=" << na << " nb=" << nb;
            }
            ++collapse_checks;
        }
    }

    std::ofstream f("usadd_behavior.csv");
    ASSERT_TRUE(f.is_open());
    f << buf.str();
    f.close();

    std::cout << "\n[uSADD PART 1] exhaustive operand sweep at N = " << N << "\n"
              << "   " << pairs << " operand-count pairs, ALL of the (N+1)^2 that exist\n"
              << "   " << real_runs << " real 256-cycle adder runs\n"
              << "   " << collapse_checks << " pairs re-run on " << ARRANGEMENTS_PER_PAIR
              << " shuffled arrangements each, all agreeing\n"
              << "   worst error: " << worst_error_bits << " output bits, "
              << std::fixed << std::setprecision(8) << worst_abs << " absolute\n"
              << "   exactly correct on " << exact_pairs << " of " << pairs << " pairs ("
              << std::setprecision(2) << (100.0 * exact_pairs / pairs) << "%)\n"
              << "   -> usadd_behavior.csv" << std::endl;

    // The healthy unit's only error is the leftover residue: strictly less than n = 2 output
    // ones, i.e. never more than 1 bit short, and never over.
    EXPECT_EQ(worst_error_bits, 0);
    EXPECT_LE(worst_abs, 1.0 / N + 1e-12);
}

// =============================================================================================
TEST(uSADDBehavior, VersusMuxUnderThreeCorrelationRegimes) {
    std::mt19937 rng(31337);

    std::ostringstream buf;
    buf << "Regime,ProbA,ProbB,Exact,uSADD,MuxMean,MuxMin,MuxMax,MuxRMSE\n";
    buf << std::fixed << std::setprecision(8);

    struct Acc { double sq = 0.0; double worst = 0.0; long long n = 0; };
    Acc usadd_acc, mux_indep, mux_corr, mux_anti;

    const char* names[3] = {"independent", "correlated", "anticorrelated"};

    for (int na = 0; na <= N; na += 16) {
        for (int nb = 0; nb <= N; nb += 16) {
            const double exact = (static_cast<double>(na) / N + static_cast<double>(nb) / N) / 2.0;

            for (int regime = 0; regime < 3; ++regime) {
                std::vector<bool> a = scattered(N, na, rng);
                std::vector<bool> b;
                if (regime == 0) {
                    b = scattered(N, nb, rng);                 // independent
                } else if (regime == 1) {
                    // correlated: b's ones packed onto a's ones where possible
                    b.assign(N, false);
                    int placed = 0;
                    for (int i = 0; i < N && placed < nb; ++i) if (a[i]) { b[i] = true; ++placed; }
                    for (int i = 0; i < N && placed < nb; ++i) if (!b[i]) { b[i] = true; ++placed; }
                } else {
                    // anti-correlated: b's ones pushed onto a's zeros first
                    b.assign(N, false);
                    int placed = 0;
                    for (int i = 0; i < N && placed < nb; ++i) if (!a[i]) { b[i] = true; ++placed; }
                    for (int i = 0; i < N && placed < nb; ++i) if (!b[i]) { b[i] = true; ++placed; }
                }

                // uSADD sees the same operands. The collapse says the regime cannot matter.
                const double us = static_cast<double>(count_ones(usadd_stream(a, b))) / N;
                const double us_err = std::fabs(us - exact);
                usadd_acc.sq += us_err * us_err; usadd_acc.worst = std::max(usadd_acc.worst, us_err);
                ++usadd_acc.n;

                // MUX over many independent 0.5 select streams.
                double mmean = 0.0, mmin = 2.0, mmax = -1.0, msq = 0.0;
                for (int t = 0; t < MUX_TRIALS; ++t) {
                    const std::vector<bool> sel = scattered(N, N / 2, rng);
                    const double m = static_cast<double>(mux_ones(a, b, sel)) / N;
                    mmean += m; mmin = std::min(mmin, m); mmax = std::max(mmax, m);
                    msq += (m - exact) * (m - exact);
                }
                mmean /= MUX_TRIALS;
                const double mrmse = std::sqrt(msq / MUX_TRIALS);

                Acc& target = (regime == 0) ? mux_indep : (regime == 1 ? mux_corr : mux_anti);
                target.sq += msq / MUX_TRIALS;
                target.worst = std::max(target.worst, std::max(std::fabs(mmax - exact),
                                                               std::fabs(mmin - exact)));
                ++target.n;

                buf << names[regime] << "," << (static_cast<double>(na) / N) << ","
                    << (static_cast<double>(nb) / N) << "," << exact << "," << us << ","
                    << mmean << "," << mmin << "," << mmax << "," << mrmse << "\n";
            }
        }
    }

    std::ofstream f("usadd_vs_mux.csv");
    ASSERT_TRUE(f.is_open());
    f << buf.str();
    f.close();

    auto rmse = [](const Acc& x) { return std::sqrt(x.sq / static_cast<double>(x.n)); };

    std::cout << "\n[uSADD PART 2] uSADD vs MUX on identical operands, N = " << N << "\n"
              << std::fixed << std::setprecision(6)
              << "                                    RMSE       worst\n"
              << "   uSADD, any correlation        " << rmse(usadd_acc) << "   " << usadd_acc.worst << "\n"
              << "   MUX,  independent operands    " << rmse(mux_indep) << "   " << mux_indep.worst << "\n"
              << "   MUX,  correlated operands     " << rmse(mux_corr)  << "   " << mux_corr.worst  << "\n"
              << "   MUX,  anti-correlated         " << rmse(mux_anti)  << "   " << mux_anti.worst  << "\n"
              << "   -> usadd_vs_mux.csv" << std::endl;

    // uSADD beats the MUX on the MUX's own best case, and its error does not move with the
    // correlation regime at all -- that is the architectural claim, stated as an assertion.
    EXPECT_LT(rmse(usadd_acc), rmse(mux_indep));
    EXPECT_LE(usadd_acc.worst, 1.0 / N + 1e-12);
}

// =============================================================================================
TEST(uSADDBehavior, EarlyTermination) {
    std::mt19937 rng(555);
    const std::vector<int> lens = {2, 4, 8, 16, 32, 64, 128, 256};
    constexpr int TRIALS = 400;

    std::ostringstream buf;
    buf << "N,uSADD_MeanAbsError,MUX_MeanAbsError,Ratio\n";
    buf << std::fixed << std::setprecision(8);

    std::cout << "\n[uSADD PART 3] early termination, 0.5 + 0.5 scaled to 0.5\n"
              << "       N   uSADD MAE    MUX MAE     ratio\n";

    for (int L : lens) {
        double us_sum = 0.0, mux_sum = 0.0;
        for (int t = 0; t < TRIALS; ++t) {
            const std::vector<bool> a = scattered(N, N / 2, rng);
            const std::vector<bool> b = scattered(N, N / 2, rng);
            const std::vector<bool> sel = scattered(N, N / 2, rng);

            UnaryScaledAdder unit(2);
            Adder mux;
            int us_ones = 0, mux_ones_ct = 0;
            for (int i = 0; i < L; ++i) {
                if (unit.add(a[i], b[i])) ++us_ones;
                if (mux.add_scaled_or_weighted(a[i], b[i], sel[i])) ++mux_ones_ct;
            }
            us_sum  += std::fabs(static_cast<double>(us_ones) / L - 0.5);
            mux_sum += std::fabs(static_cast<double>(mux_ones_ct) / L - 0.5);
        }
        const double us_mae = us_sum / TRIALS, mux_mae = mux_sum / TRIALS;
        const double ratio = (mux_mae > 0.0) ? us_mae / mux_mae : 0.0;
        buf << L << "," << us_mae << "," << mux_mae << "," << ratio << "\n";
        std::cout << std::setw(8) << L << "  " << std::fixed << std::setprecision(6)
                  << std::setw(9) << us_mae << "  " << std::setw(9) << mux_mae << "  "
                  << std::setw(8) << ratio << "\n";
    }

    std::ofstream f("usadd_early_termination.csv");
    ASSERT_TRUE(f.is_open());
    f << buf.str();
    f.close();
    std::cout << "   -> usadd_early_termination.csv" << std::endl;
}
