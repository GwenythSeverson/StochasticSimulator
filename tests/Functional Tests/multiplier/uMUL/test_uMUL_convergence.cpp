#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// Four levels up: uMUL/ -> multiplier/ -> Functional Tests/ -> tests/ -> repo root
#include "../../../../modules/uMUL.hpp"
#include "../../../../general_functions.hpp"
#include "../../../../bsg/lfsr.hpp"
#include "../../../../bsg/sng.hpp"

// =============================================================================================
// uMUL CONVERGENCE / EARLY-TERMINATION SWEEPS -- THREE OPERATING POINTS
//
// The companion to test_uMUL_behavior.cpp. That one asks "how close is the FINAL answer, for
// every operand pair"; these ask "how does ONE answer get closer as the stream runs", which is
// the plot shape the AND-gate and mux-adder scripts already produce.
//
// THREE IDENTICAL EXPERIMENTS at three output levels, so the shape can be compared near the
// extrema rather than only in the comfortable middle:
//
//     LOW    0.5000 x 0.2002 = 0.1001      umul_convergence_low.csv
//     MID    0.5000 x 0.8750 = 0.4375      umul_convergence_mid.csv
//     HIGH   0.9375 x 0.9600 = 0.9000      umul_convergence_high.csv
//
// Same method, same truncation lengths, same CSV schema -- only the operands move. Each CSV's
// Operation field names its own level so a plot cannot be mistaken for one of its siblings.
// MID uses the same target as test_SuMUL_accuracy.cpp on purpose, so those two CSVs drop into the
// same axes and the difference between the units is the difference between the curves.
//
// WHY THE EXTREMA LOOK DIFFERENT, and it is worth predicting before you plot:
//   Short truncations do not measure the multiplier, they measure how well G's Sobol PREFIX can
//   express p1 in the handful of enabled cycles that have happened. The first eight points of
//   Sobol dimension 1 at width 10 are 0, 512, 768, 256, 384, 896, 640, 128, so
//     LOW  (threshold 205): only j = 0 and j = 7 pass -> 2/8 = 0.25 against a target of 0.20.
//                           The short-N curve OVERSHOOTS.
//     MID  (threshold 896): j = 5 is the only miss -> 7/8 = 0.875, exactly the target, from the
//                           moment eight enabled cycles exist.
//     HIGH (threshold 983): every one of the first eight passes -> 8/8 = 1.0 against 0.96. The
//                           short-N curve reads p0 rather than p0 * p1, then settles.
//   None of that is warm-up in SuMUL's sense; there is no estimator converging on anything. It is
//   quantization -- you cannot express 7/8, or 0.2, or 0.96 in four samples.
//
// WHAT DIFFERS FROM THE SuMUL SWEEP, and it is structural rather than a choice:
//   SuMUL takes two STREAMS, so its arrangement space is 1023 x 1023 = 1,046,529 pairs.
//   uMUL takes one stream and one NUMBER, so only in_0 has arrangements: 1023 of them. There is
//   no second stream to arrange. The mean is just as trustworthy; the min/max envelope is drawn
//   from a thousand samples rather than a million, so it is a fair envelope over a smaller space.
// =============================================================================================

namespace StochasticSimulator {
namespace {

constexpr uint16_t STREAM_BITS = 1024;
constexpr uint16_t MIN_SEED = 1;
constexpr uint16_t MAX_SEED = 1023;   // FlexibleLFSR rejects seed >= max_cycles
constexpr uint16_t SEED_STRIDE = 1;   // 1 = every LFSR seed
constexpr unsigned WIDTH = 10;        // log2(1024): RNG period exactly covers the run
// Every tabulated Sobol dimension is distinct at width >= 7, so at width 10 all 19 are usable.
// Averaging the exact result over all of them is what separates "how uMUL behaves" from "what
// dimension 1's prefix happens to do" -- and the prefix structure is the ONLY thing left that
// makes the error curve kink, now that arrangement sampling has been removed entirely.
constexpr unsigned DIM_COUNT = 19;

// Monte Carlo over uniformly random arrangements. This is NOT how the published numbers are
// produced -- the hypergeometric sum below is exact and needs no sampling. It exists to VALIDATE
// that closed form: if randomly shuffled streams do not land on the exact answer, one of the two
// is wrong.
//
// WHAT RAISING THIS DOES, AND WHAT IT DOES NOT. There are two unrelated random generators in this
// experiment and only one of them responds to trial count:
//   * the SHUFFLE RNG that picks arrangements -- pure sampling noise, falls as 1/sqrt(trials).
//     More trials tightens the validation. That is this knob.
//   * the SOBOL RNG INSIDE G -- the circuit's own generator. Its prefix has clustered pass
//     positions, and that is what puts the kinks in the error curve. It is not noise and no
//     number of trials touches it, because every trial runs the SAME generator. The lever for
//     that is DIM_COUNT above: average over independent dimensions.
// So raising this will NOT flatten the bumps. It will only make the agreement column tighter.
//
// At 50k the standard error is ~0.0005 at the noisiest truncation, well under the structure being
// measured. Costs roughly 50 ms per level.
constexpr int MC_TRIALS = 50000;

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

struct Accumulator {
    double sum = 0.0;
    double sum_sq = 0.0;
    double sum_abs_err = 0.0;
    double min_est = 1e9;
    double max_est = -1e9;

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

// ---- Exhaustive averaging over EVERY arrangement of in_0 ------------------------------------
//
// There are C(1024, 512) ~ 10^306 ways to arrange in_0's ones, so enumerating them is not an
// option -- but it is also not necessary. The output ones count after k enabled cycles is
//
//     out(k) = #{ j < k : sobol[j] < value }
//
// which depends ONLY on k, never on WHERE in_0's ones sit. So truncating at N, the single thing
// that varies across arrangements is how many of in_0's ones landed in the first N cycles, and
// that integer has an exact distribution: drawing N positions out of 1024 from a stream carrying
// ones_A ones is a HYPERGEOMETRIC draw. Weighting out(k)/N by that pmf and summing over
// k = 0..N gives the true expectation over all 10^306 arrangements in at most 1025 terms, with
// no sampling and no seed bias anywhere.
//
// This is what "exhaustive" means here, and it is strictly stronger than enumerating arrangements
// would be -- it is the same answer, computed exactly rather than estimated.

double log_choose(double n, double k) {
    if (k < 0.0 || k > n) return -std::numeric_limits<double>::infinity();
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

// P(exactly k of the first N cycles are ones | the whole stream carries K ones out of n)
double hypergeometric_pmf(long n, long K, long N, long k) {
    if (k < 0 || k > K || (N - k) < 0 || (N - k) > (n - K)) return 0.0;
    double lp = log_choose(K, k) + log_choose(n - K, N - k) - log_choose(n, N);
    return std::exp(lp);
}

struct ExactStats {
    double mean_est = 0.0;
    double mean_abs_err = 0.0;
    double stddev = 0.0;
    double min_est = 0.0;   // true extremes over the whole support
    double max_est = 0.0;
    double p01_est = 0.0;   // 1st / 99th percentile -- what is worth shading on a plot, since the
    double p99_est = 0.0;   // support extremes are reachable but astronomically improbable
};

/**
 * Exact statistics over every arrangement, given the generator's pass-prefix.
 *
 * @param pass_prefix  pass_prefix[k] = output ones after k enabled cycles. Built by running the
 *                     real unit with the enable held high, so this is the circuit's own behaviour
 *                     rather than a reimplementation of the comparator.
 */
ExactStats exact_stats_for(const std::vector<uint32_t>& pass_prefix,
                           uint16_t ones_A, uint16_t trunc_N, double ideal) {
    ExactStats s;
    const long n = STREAM_BITS, K = ones_A, N = trunc_N;
    const long k_lo = std::max(0L, N - (n - K));
    const long k_hi = std::min<long>(N, K);

    s.min_est = 1e9;
    s.max_est = -1e9;

    double cumulative = 0.0;
    bool have_p01 = false, have_p99 = false;
    double sum_w = 0.0, sum_sq = 0.0;

    for (long k = k_lo; k <= k_hi; ++k) {
        double w = hypergeometric_pmf(n, K, N, k);
        if (w <= 0.0) continue;
        double est = static_cast<double>(pass_prefix[static_cast<std::size_t>(k)]) /
                     static_cast<double>(N);

        s.mean_est += w * est;
        s.mean_abs_err += w * std::fabs(est - ideal);
        sum_sq += w * est * est;
        sum_w += w;

        if (est < s.min_est) s.min_est = est;
        if (est > s.max_est) s.max_est = est;

        // est is non-decreasing in k, since pass_prefix is, so the first crossing of each
        // cumulative threshold is the quantile. Taking the last point BELOW 0.99 instead would
        // exclude the top probability mass entirely -- at HIGH, N = 2, that is 88% of it.
        cumulative += w;
        if (!have_p01 && cumulative >= 0.01) { s.p01_est = est; have_p01 = true; }
        if (!have_p99 && cumulative >= 0.99) { s.p99_est = est; have_p99 = true; }
    }

    // Renormalise against the summed weight rather than assuming it is exactly 1: the pmf is
    // evaluated through lgamma, so it carries a little floating-point slack.
    if (sum_w > 0.0) {
        s.mean_est /= sum_w;
        s.mean_abs_err /= sum_w;
        double var = sum_sq / sum_w - s.mean_est * s.mean_est;
        s.stddev = var > 0.0 ? std::sqrt(var) : 0.0;
    }
    if (!have_p01) s.p01_est = s.min_est;
    if (s.p99_est < s.p01_est) s.p99_est = s.max_est;
    return s;
}

/**
 * pass_prefix[k] for one Sobol dimension, measured by running the real UnaryMultiplier with the
 * enable held high. Element k is the output ones count after k enabled cycles.
 */
std::vector<uint32_t> pass_prefix_for(uint16_t value_B, unsigned dimension) {
    UnaryMultiplier umul(WIDTH, dimension);
    umul.load_value(value_B);
    umul.reset();

    std::vector<uint32_t> prefix(STREAM_BITS + 1, 0);
    for (std::size_t t = 0; t < STREAM_BITS; ++t) {
        prefix[t + 1] = prefix[t] + (umul.multiply(true) ? 1u : 0u);
    }
    return prefix;
}

/**
 * One complete sweep at one operating point. Everything that differs between LOW, MID and HIGH is
 * a parameter here, so the three runs are genuinely the same experiment rather than three
 * hand-copied ones that quietly drift apart.
 *
 * @param level     "LOW" / "MID" / "HIGH" -- goes into the CSV so a plot names its own level.
 * @param ones_A    in_0's ones count. in_0 is always a STREAM.
 * @param value_B   in_1's register value. in_1 is always a NUMBER.
 * @param filename  CSV to write.
 */
void run_convergence_sweep(const char* level, uint16_t ones_A, uint16_t value_B,
                           const std::string& filename) {
    const double target_pa = static_cast<double>(ones_A) / STREAM_BITS;
    const double target_pb = static_cast<double>(value_B) / STREAM_BITS;
    const double ideal_product = target_pa * target_pb;

    // Carried in the CSV so the MATLAB titles state exactly what was computed, and which of the
    // three siblings this is, rather than hard-coding it in places that then drift apart.
    std::ostringstream op;
    op << std::fixed << std::setprecision(4)
       << "uMUL " << level << " output   " << target_pa << " x " << target_pb
       << " = " << ideal_product << "   (1 of 3: LOW / MID / HIGH)";
    const std::string operation = op.str();

    const std::vector<uint16_t>& lens = truncation_lengths();

    std::cout << "\n[INFO] " << operation << std::endl;
    std::cout << "[INFO] in_1 is LOADED as the binary value " << value_B
              << "/1024 -- no stream, so no arrangements on that operand" << std::endl;

    // ---- every distinct in_0 arrangement ----
    std::vector<uint16_t> seeds;
    for (uint16_t s = MIN_SEED; s <= MAX_SEED; s += SEED_STRIDE) seeds.push_back(s);

    std::vector<std::vector<bool>> A_streams;
    A_streams.reserve(seeds.size());
    for (uint16_t s : seeds) A_streams.push_back(stream_for(ones_A, s));
    for (std::size_t i = 0; i < seeds.size(); ++i) {
        ASSERT_EQ(count_ones(A_streams[i]), ones_A) << "A seed " << seeds[i];
    }
    const std::size_t total_trials = seeds.size();

    // ---- sweep ----
    std::vector<Accumulator> acc(lens.size());
    UnaryMultiplier umul(WIDTH, 1);
    umul.load_value(value_B);
    ASSERT_EQ(umul.get_value(), value_B);

    for (std::size_t ia = 0; ia < A_streams.size(); ++ia) {
        const std::vector<bool>& A = A_streams[ia];
        umul.reset();  // clears the run, keeps the loaded operand

        std::size_t ones = 0;
        std::size_t next_len = 0;
        for (uint16_t k = 0; k < STREAM_BITS; ++k) {
            if (umul.multiply(A[k])) ++ones;
            if (next_len < lens.size() && (k + 1) == lens[next_len]) {
                acc[next_len].add(static_cast<double>(ones) / lens[next_len], ideal_product);
                ++next_len;
            }
        }
    }

    // ---- EXHAUSTIVE pass: exact over every arrangement, no sampling ----
    std::vector<uint32_t> prefix_d1 = pass_prefix_for(value_B, 1);
    std::vector<ExactStats> exact(lens.size());
    for (std::size_t i = 0; i < lens.size(); ++i) {
        exact[i] = exact_stats_for(prefix_d1, ones_A, lens[i], ideal_product);
    }

    // ---- and again averaged over independent generators ----
    // One Sobol dimension is one particular circuit, and its prefix has structure of its own --
    // clustered pass positions that make the error curve kink. Averaging the exact result over
    // every usable dimension is what a whole ARRAY of these units looks like, and it is the only
    // remaining way to smooth the curve, because arrangement sampling has already been removed.
    // A dimension that collides with a lower one at this width is skipped rather than silently
    // double-counted; SobolRNG throws on those, which is exactly the behaviour we want.
    std::vector<double> dim_mean(lens.size(), 0.0), dim_abs(lens.size(), 0.0),
                        dim_std(lens.size(), 0.0);
    unsigned dims_used = 0;
    for (unsigned d = 1; d <= DIM_COUNT; ++d) {
        std::vector<uint32_t> prefix_d;
        try {
            prefix_d = pass_prefix_for(value_B, d);
        } catch (const std::invalid_argument&) {
            continue;  // dimension not distinct at this width
        }
        ++dims_used;
        for (std::size_t i = 0; i < lens.size(); ++i) {
            ExactStats e = exact_stats_for(prefix_d, ones_A, lens[i], ideal_product);
            dim_mean[i] += e.mean_est;
            dim_abs[i] += e.mean_abs_err;
            dim_std[i] += e.stddev;
        }
    }
    ASSERT_GT(dims_used, 0u);
    for (std::size_t i = 0; i < lens.size(); ++i) {
        dim_mean[i] /= dims_used;
        dim_abs[i] /= dims_used;
        dim_std[i] /= dims_used;
    }

    // ---- MONTE CARLO VALIDATION over uniformly random arrangements ----
    // The exact pass above claims to average over all ~10^307 arrangements without visiting any
    // of them. This shuffles MC_TRIALS genuinely random arrangements of the same ones_A ones,
    // runs each through the real unit, and checks the two agree. It is a CHECK, not a result --
    // if it ever disagrees, the hypergeometric weighting is wrong and everything above is void.
    std::vector<Accumulator> mc(lens.size());
    {
        std::mt19937_64 gen(20260730ull);  // fixed seed: the validation must be reproducible
        std::vector<bool> shuffled(STREAM_BITS, false);
        std::vector<std::size_t> positions(STREAM_BITS);
        for (std::size_t i = 0; i < STREAM_BITS; ++i) positions[i] = i;

        UnaryMultiplier mc_umul(WIDTH, 1);
        mc_umul.load_value(value_B);

        for (int trial = 0; trial < MC_TRIALS; ++trial) {
            std::shuffle(positions.begin(), positions.end(), gen);
            std::fill(shuffled.begin(), shuffled.end(), false);
            for (uint16_t i = 0; i < ones_A; ++i) shuffled[positions[i]] = true;

            mc_umul.reset();
            std::size_t ones = 0, next_len = 0;
            for (uint16_t k = 0; k < STREAM_BITS; ++k) {
                if (mc_umul.multiply(shuffled[k])) ++ones;
                if (next_len < lens.size() && (k + 1) == lens[next_len]) {
                    mc[next_len].add(static_cast<double>(ones) / lens[next_len], ideal_product);
                    ++next_len;
                }
            }
        }
    }

    // ---- CSV ----
    // Mean_Est / Mean_AbsError / Min_Est / Max_Est / Std_Est now carry the EXHAUSTIVE result, so
    // anything reading those column names gets the unbiased numbers. The seed-sampled values are
    // kept alongside under Seed_* so the sampling bias can be seen rather than assumed away.
    std::ofstream csv(filename);
    ASSERT_TRUE(csv.is_open()) << "CRITICAL: Failed to write " << filename;

    const double log10_arrangements =
        log_choose(static_cast<double>(STREAM_BITS), static_cast<double>(ones_A)) / std::log(10.0);

    csv << "Operation,Target_PA,Target_PB,Ideal_Product,Log10_Arrangements,Seed_Samples,Dims,"
           "MC_Trials,N,"
           "Mean_Est,Mean_AbsError,Min_Est,Max_Est,Std_Est,P01_Est,P99_Est,"
           "Seed_Mean_Est,Seed_Mean_AbsError,Seed_Min_Est,Seed_Max_Est,Seed_Std_Est,"
           "DimAvg_Mean_Est,DimAvg_Mean_AbsError,DimAvg_Std_Est,"
           "MC_Mean_Est,MC_Mean_AbsError,MC_Std_Est\n";
    csv << std::fixed << std::setprecision(8);
    for (std::size_t i = 0; i < lens.size(); ++i) {
        csv << "\"" << operation << "\"," << target_pa << "," << target_pb << ","
            << ideal_product << "," << log10_arrangements << "," << total_trials << ","
            << dims_used << "," << MC_TRIALS << "," << lens[i] << ","
            << exact[i].mean_est << "," << exact[i].mean_abs_err << ","
            << exact[i].min_est << "," << exact[i].max_est << "," << exact[i].stddev << ","
            << exact[i].p01_est << "," << exact[i].p99_est << ","
            << acc[i].mean(total_trials) << "," << acc[i].sum_abs_err / total_trials << ","
            << acc[i].min_est << "," << acc[i].max_est << "," << acc[i].stddev(total_trials) << ","
            << dim_mean[i] << "," << dim_abs[i] << "," << dim_std[i] << ","
            << mc[i].mean(MC_TRIALS) << "," << mc[i].sum_abs_err / MC_TRIALS << ","
            << mc[i].stddev(MC_TRIALS) << "\n";
    }
    csv.close();

    // ---- console ----
    std::cout << std::fixed << std::setprecision(0)
              << "        EXACT over all 10^" << log10_arrangements << " arrangements"
              << "   |  " << dims_used << "-DIM AVG  |  MC " << MC_TRIALS
              << " shuffles (validation)\n" << std::setprecision(5)
              << "     N   mean est  mean |err|     p01       p99    "
              << "|  mean |err|  |  mean est   delta vs exact\n";
    double worst_mc_gap = 0.0;
    for (std::size_t i = 0; i < lens.size(); ++i) {
        double gap = std::fabs(mc[i].mean(MC_TRIALS) - exact[i].mean_est);
        worst_mc_gap = std::max(worst_mc_gap, gap);
        std::cout << std::setw(6) << lens[i] << "   "
                  << exact[i].mean_est << "    " << exact[i].mean_abs_err << "   "
                  << exact[i].p01_est << "   " << exact[i].p99_est << "  |   "
                  << dim_abs[i] << "   |  " << mc[i].mean(MC_TRIALS) << "    "
                  << std::scientific << std::setprecision(2) << gap
                  << std::fixed << std::setprecision(5) << "\n";
    }
    std::cout << "  worst MC-vs-exact gap: " << std::scientific << std::setprecision(2)
              << worst_mc_gap << std::fixed << "\n"
              << "[SUCCESS] Exported to: " << filename << std::endl;

    // ---- guards, identical in form at all three levels, read off the EXHAUSTIVE numbers ----
    const std::size_t last = lens.size() - 1;
    const std::size_t n64 = 5;
    ASSERT_EQ(lens[n64], 64);

    // At full length the answer must be right to well under one bit of the 1024.
    EXPECT_NEAR(exact[last].mean_est, ideal_product, 0.005)
        << level << ": full-length output drifted off the true product";

    // No estimator anywhere in the unit, so by N = 64 the mean must already be on target. SuMUL's
    // companion test asserts the OPPOSITE at short N -- its window is still filling. If this ever
    // fails at any of the three levels, an estimator has crept back in.
    EXPECT_NEAR(exact[n64].mean_est, ideal_product, 0.03)
        << level << ": should be unbiased well before full length -- there is no warm-up to pay";

    // Every arrangement carries exactly ones_A ones, so at full length they all enable the same
    // number of cycles and G returns the same count. The support must therefore collapse to a
    // single point: the output depends on HOW MANY ones in_0 has, never on WHERE they are. This
    // is now a statement about ALL 10^306 arrangements, not about 1023 sampled ones.
    EXPECT_DOUBLE_EQ(exact[last].min_est, exact[last].max_est)
        << level << ": arrangement must not matter once the whole stream has run";
    EXPECT_DOUBLE_EQ(exact[last].stddev, 0.0)
        << level << ": zero spread over every arrangement at full length";
    EXPECT_GT(exact[n64].max_est - exact[n64].min_est, 0.0)
        << level << ": the support should still be open at N = 64";

    // The 1023-seed sweep is a SAMPLE of that exact distribution, so it must agree with it -- but
    // only loosely, and the gap is the sampling bias the exhaustive pass exists to remove.
    EXPECT_NEAR(acc[last].mean(total_trials), exact[last].mean_est, 0.005)
        << level << ": seed sample and exact result disagree at full length";

    // THE VALIDATION THAT MAKES THE CLOSED FORM TRUSTWORTHY. MC_TRIALS genuinely random
    // arrangements, run through the real unit, must land on the hypergeometric answer at every
    // truncation. The tolerance is a few standard errors of the MC estimate itself, so this gets
    // tighter -- not looser -- as MC_TRIALS rises. If it ever fails, the weighting is wrong and
    // every "exhaustive" number in the CSV is void.
    for (std::size_t i = 0; i < lens.size(); ++i) {
        double se = exact[i].stddev / std::sqrt(static_cast<double>(MC_TRIALS));
        double tol = std::max(4.0 * se, 1e-6);
        EXPECT_NEAR(mc[i].mean(MC_TRIALS), exact[i].mean_est, tol)
            << level << ": " << MC_TRIALS << " random arrangements disagree with the exact "
            << "hypergeometric result at N = " << lens[i];
    }
}

}  // namespace

// =============================================================================================
// EVERY OPERAND AGAINST EVERY OPERAND
//
// The three named levels above are three points in a space of 1025 x 1024 = 1,049,600 operand
// pairs. This sweeps all of them, and at each pair computes the SAME exact hypergeometric
// expectation over all ~10^307 in_0 arrangements. So the result is the whole unit's convergence,
// not one operating point's:
//
//     1,049,600 operand pairs  x  ~10^307 arrangements each  x  10 truncation lengths
//
// and none of it is sampled. Reported per truncation length:
//     mean    -- the unit's average convergence error, over every operand pair
//     max     -- the worst pair at that length
//     min     -- the best pair, which is exactly 0 and always will be: k1 = 0 and k0 = 0 emit
//                nothing and are error-free by construction, so the true minimum is degenerate
//     p05/p95 -- the band an arbitrary operand pair actually falls in, which is the useful bound
//
// COST: the naive form is ~1e9 inner iterations. Two things keep it tractable -- the pass-prefix
// table is built ONCE per operand value by running the real unit (1024 runs, not 1e6), and the
// hypergeometric support is truncated where the weight falls below 1e-13, which shrinks the inner
// loop from O(N) to O(sqrt(N)) around the mode.
// =============================================================================================
TEST(UMULConvergenceTest, ExhaustiveOverEveryOperandPair) {
    // EVERY termination point, not just the ten powers of two. Early termination means stopping
    // at some cycle, and every cycle is a legal stopping point, so the curve is only complete if
    // all 1024 of them are computed.
    std::vector<uint16_t> lens;
    lens.reserve(STREAM_BITS);
    for (uint16_t n = 1; n <= STREAM_BITS; ++n) lens.push_back(n);

    constexpr std::size_t STREAM_STEPS = STREAM_BITS + 1;  // k0 = 0..1024, a tally over 1024 bits
    constexpr std::size_t VALUE_STEPS = STREAM_BITS;       // k1 = 0..1023, a 10-bit register

    std::cout << "\n[INFO] exhaustive over every operand pair: " << STREAM_STEPS << " x "
              << VALUE_STEPS << " = " << (STREAM_STEPS * VALUE_STEPS) << " pairs,\n"
              << "       each averaged exactly over all ~10^307 in_0 arrangements,\n"
              << "       at every one of " << lens.size() << " termination points" << std::endl;

    // ---- pass-prefix table, built by running the REAL unit once per operand value ----
    // prefix[k1][k] = output ones after k enabled cycles with the register holding k1.
    std::vector<std::vector<uint16_t>> prefix(VALUE_STEPS);
    for (std::size_t k1 = 0; k1 < VALUE_STEPS; ++k1) {
        std::vector<uint32_t> p = pass_prefix_for(static_cast<uint16_t>(k1), 1);
        prefix[k1].assign(p.begin(), p.end());
    }

    // ---- log-factorial table, so a hypergeometric weight is table lookups rather than lgamma ----
    std::vector<double> logfact(STREAM_BITS + 1, 0.0);
    for (std::size_t i = 1; i <= STREAM_BITS; ++i) {
        logfact[i] = logfact[i - 1] + std::log(static_cast<double>(i));
    }
    auto logC = [&](long n, long k) {
        return logfact[n] - logfact[k] - logfact[n - k];
    };

    std::ofstream csv("umul_convergence_allpairs.csv");
    ASSERT_TRUE(csv.is_open());
    csv << "Pairs,Arrangements_Log10,N,Mean_AbsError,Min_AbsError,Max_AbsError,"
           "P05_AbsError,P95_AbsError,Median_AbsError\n";
    csv << std::fixed << std::setprecision(8);

    std::vector<double> per_pair;
    per_pair.reserve(STREAM_STEPS * VALUE_STEPS);
    std::vector<long> ks;
    std::vector<double> ws;

    std::cout << std::fixed << std::setprecision(6)
              << "\n     N        mean        p05      median        p95         max\n";

    for (std::size_t li = 0; li < lens.size(); ++li) {
        const long N = lens[li];
        const double invN = 1.0 / static_cast<double>(N);
        per_pair.clear();

        for (long k0 = 0; k0 < static_cast<long>(STREAM_STEPS); ++k0) {
            // Compact hypergeometric support for this (N, k0): drop weights below 1e-13, which
            // costs nothing measurable and turns an O(N) inner loop into O(sqrt(N)).
            ks.clear();
            ws.clear();
            const long k_lo = std::max(0L, N - (static_cast<long>(STREAM_BITS) - k0));
            const long k_hi = std::min(N, k0);
            const double lnorm = logC(STREAM_BITS, N);
            double wsum = 0.0;
            for (long k = k_lo; k <= k_hi; ++k) {
                double lw = logC(k0, k) + logC(STREAM_BITS - k0, N - k) - lnorm;
                if (lw < -30.0) continue;  // exp(-30) ~ 1e-13
                double w = std::exp(lw);
                ks.push_back(k);
                ws.push_back(w);
                wsum += w;
            }
            if (ws.empty()) continue;
            for (double& w : ws) w /= wsum;

            for (long k1 = 0; k1 < static_cast<long>(VALUE_STEPS); ++k1) {
                const double ideal = (static_cast<double>(k0) / STREAM_BITS) *
                                     (static_cast<double>(k1) / STREAM_BITS);
                const std::vector<uint16_t>& pk = prefix[k1];
                double err = 0.0;
                for (std::size_t t = 0; t < ks.size(); ++t) {
                    err += ws[t] * std::fabs(static_cast<double>(pk[ks[t]]) * invN - ideal);
                }
                per_pair.push_back(err);
            }
        }

        const std::size_t n = per_pair.size();
        double mean = 0.0, lo = per_pair[0], hi = per_pair[0];
        for (double v : per_pair) {
            mean += v;
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
        mean /= static_cast<double>(n);

        // nth_element rather than a full sort: three O(n) partial selections instead of one
        // O(n log n) sort, done 1024 times. Over the whole sweep that is the difference between
        // a couple of minutes and most of an hour.
        auto quantile = [&](double q) {
            std::size_t idx = static_cast<std::size_t>(q * (n - 1));
            std::nth_element(per_pair.begin(), per_pair.begin() + idx, per_pair.end());
            return per_pair[idx];
        };
        const double p05 = quantile(0.05), p50 = quantile(0.50), p95 = quantile(0.95);

        csv << n << ",307," << N << "," << mean << "," << lo << "," << hi << ","
            << p05 << "," << p95 << "," << p50 << "\n";

        // Print only the powers of two -- 1024 console rows would bury everything else.
        if ((N & (N - 1)) == 0) {
            std::cout << std::setw(6) << N << "  " << std::setw(10) << mean << "  "
                      << std::setw(10) << p05 << "  " << std::setw(10) << p50 << "  "
                      << std::setw(10) << p95 << "  " << std::setw(10) << hi << "\n";
        }

        // The best pair is exactly zero at every length, and that is structural rather than a
        // rounding artifact: with k1 = 0 the strict ">" never fires, so the unit emits nothing
        // and matches an ideal product of zero perfectly.
        EXPECT_DOUBLE_EQ(lo, 0.0) << "the k1 = 0 rail should be error-free at N = " << N;
    }
    csv.close();
    std::cout << "[SUCCESS] Exported " << lens.size()
              << " termination points to: umul_convergence_allpairs.csv" << std::endl;
}

// LOW output -- 0.5000 x 0.2002 = 0.1001
TEST(UMULConvergenceTest, LowOutputProbability) {
    run_convergence_sweep("LOW", 512, 205, "umul_convergence_low.csv");
}

// MID output -- 0.5000 x 0.8750 = 0.4375, matching test_SuMUL_accuracy.cpp
TEST(UMULConvergenceTest, MidOutputProbability) {
    run_convergence_sweep("MID", 512, 896, "umul_convergence_mid.csv");
}

// HIGH output -- 0.9375 x 0.9600 = 0.9000
TEST(UMULConvergenceTest, HighOutputProbability) {
    run_convergence_sweep("HIGH", 960, 983, "umul_convergence_high.csv");
}

}  // namespace StochasticSimulator
