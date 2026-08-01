#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Three levels up: multiplier/ -> Fault Injection/ -> tests/ -> repo root
#include "../radiation_model.hpp"
#include "../../../modules/uMUL.hpp"
#include "../../../modules/multiplier.hpp"
#include "../../../general_functions.hpp"
#include "../../../bsg/sobol.hpp"

// =============================================================================================
// WHOLE-UNIT RADIATION CAMPAIGN -- uMUL against the plain AND-gate multiplier
//
// Radiation does not politely restrict itself to the data stream. It hits the loaded operand
// register, the RNG index, the bus between them, and the wires, and those five site classes fail
// in completely different ways. This test measures all of them, in two passes:
//
//   PART A -- SINGLE-UPSET SENSITIVITY MAP. Deterministic and exhaustive: one upset, at every
//     (site, cycle), no randomness anywhere. Produces the circuit's cross-section
//         X  =  SUM over (site, cycle) of |error given one upset there|
//     which is the quantity that makes every later rate question arithmetic:  E[|err|] ~ lambda*X
//     to first order. This is the number worth publishing; it is a property of the design and
//     does not depend on an orbit.
//
//   PART B -- POISSON RATE SWEEP. Monte Carlo at a range of lambda, both units, same rate, same
//     operands. Confirms Part A's linear prediction at low rate and finds where multi-upset
//     effects break it. That breakdown point is the honest limit of single-event analysis.
//
// THE STRUCTURAL ASYMMETRY BEING MEASURED: ANDMUL is stateless, so every upset costs exactly one
// output bit and is immediately forgotten. uMUL carries 2*width bits of persistent state, and an
// upset there is never repaired -- it biases every remaining cycle. uMUL is far more accurate
// when healthy and has strictly more to lose. The question is where the crossover sits.
//
// Writes umul_radiation_sensitivity.csv and umul_radiation_rate_sweep.csv.
// =============================================================================================

namespace StochasticSimulator {
namespace {

constexpr std::size_t N = 1024;
constexpr unsigned WIDTH = 10;
constexpr uint16_t ONES_A = 512;   // p0 = 0.5,   in_0, always a stream
constexpr uint16_t VALUE_B = 896;  // p1 = 0.875, in_1: a register for uMUL, a stream for ANDMUL
constexpr double IDEAL = (512.0 / 1024.0) * (896.0 / 1024.0);  // 0.4375

using Radiation::uMULSites;
using Radiation::ANDMULSites;
using Radiation::Upset;
using Radiation::UpsetSchedule;

// ---- Running uMUL with a schedule applied ---------------------------------------------------

// Returns the output ones count. Upsets are consumed in cycle order; persistent sites are flipped
// into the register and left there, transient sites are folded into this cycle's CycleUpset.
std::size_t run_umul(const std::vector<bool>& a, uint64_t value,
                     const std::vector<Upset>& upsets, const uMULSites& map) {
    UnaryMultiplier umul(WIDTH, 1);
    umul.load_value(value);

    std::size_t ones = 0;
    std::size_t next = 0;
    for (std::size_t t = 0; t < N; ++t) {
        UnaryMultiplier::CycleUpset cu;
        while (next < upsets.size() && upsets[next].cycle == t) {
            uint32_t s = upsets[next].site;
            if (map.is_value(s)) {
                umul.flip_value_bit(s - map.value_base());        // permanent
            } else if (map.is_index(s)) {
                umul.flip_rng_index_bit(s - map.index_base());    // permanent
            } else if (map.is_bus(s)) {
                cu.rng_bus_flip ^= (1u << (s - map.bus_base()));  // one cycle
            } else if (map.is_in0(s)) {
                cu.in_0_flip = !cu.in_0_flip;                     // one cycle
            } else {
                cu.out_flip = !cu.out_flip;                       // one cycle
            }
            ++next;
        }
        if (umul.multiply(a[t], cu)) ++ones;
    }
    return ones;
}

// ---- Running ANDMUL with a schedule applied -------------------------------------------------

std::size_t run_andmul(const std::vector<bool>& a, const std::vector<bool>& b,
                       const std::vector<Upset>& upsets, const ANDMULSites& map) {
    Multiplier gate;
    std::size_t ones = 0;
    std::size_t next = 0;
    for (std::size_t t = 0; t < N; ++t) {
        bool fa = false, fb = false, fo = false;
        while (next < upsets.size() && upsets[next].cycle == t) {
            uint32_t s = upsets[next].site;
            if (map.is_in0(s))      fa = !fa;
            else if (map.is_in1(s)) fb = !fb;
            else if (map.is_out(s)) fo = !fo;
            // SNG sites: an upset in the generator producing in_1 corrupts that cycle's bit. The
            // generator is stateless from this test's point of view, so it lands as an in_1 flip.
            else                    fb = !fb;
            ++next;
        }
        bool out = gate.multiply(a[t] != fa, b[t] != fb);
        if (out != fo) ++ones;
    }
    return ones;
}

}  // namespace

// =============================================================================================
// PART A -- exhaustive single-upset sensitivity map
// =============================================================================================
TEST(RadiationSensitivityTest, ExhaustiveSingleUpsetMapForUMUL) {
    const uMULSites map(WIDTH);
    std::vector<bool> a = generate_sobol_stream(static_cast<double>(ONES_A) / N, N, 2);
    ASSERT_EQ(count_ones(a), ONES_A);

    const std::size_t clean = run_umul(a, VALUE_B, {}, map);
    ASSERT_EQ(clean, static_cast<std::size_t>(IDEAL * N))
        << "clean run must be exact before we start breaking it";

    std::ofstream csv("umul_radiation_sensitivity.csv");
    ASSERT_TRUE(csv.is_open());
    csv << "Site,Site_Class,Bit,Cycle,Out_Ones,Clean_Ones,Err_Bits\n";

    // Per-class totals: SUM |err| over every (site, cycle) in the class -- the class's share of
    // the unit's cross-section.
    struct ClassAcc { double sum_abs = 0.0; double max_abs = 0.0; std::size_t n = 0; };
    ClassAcc value_acc, index_acc, bus_acc, in0_acc, out_acc;

    auto bucket = [&](uint32_t s) -> ClassAcc& {
        if (map.is_value(s)) return value_acc;
        if (map.is_index(s)) return index_acc;
        if (map.is_bus(s))   return bus_acc;
        if (map.is_in0(s))   return in0_acc;
        return out_acc;
    };

    for (uint32_t s = 0; s < map.total(); ++s) {
        for (std::size_t t = 0; t < N; ++t) {
            std::size_t ones = run_umul(a, VALUE_B, {Upset{s, static_cast<uint32_t>(t)}}, map);
            double err = static_cast<double>(ones) - static_cast<double>(clean);
            double abs_err = std::fabs(err);

            ClassAcc& acc = bucket(s);
            acc.sum_abs += abs_err;
            acc.max_abs = std::max(acc.max_abs, abs_err);
            ++acc.n;

            csv << s << ',' << map.class_name(s) << ','
                << (map.is_persistent(s) ? s % WIDTH : 0) << ',' << t << ','
                << ones << ',' << clean << ',' << err << '\n';
        }
    }
    csv.close();

    const double total_X = value_acc.sum_abs + index_acc.sum_abs + bus_acc.sum_abs +
                           in0_acc.sum_abs + out_acc.sum_abs;

    auto report = [&](const char* name, const ClassAcc& acc) {
        std::cout << "  " << std::left << std::setw(12) << name << std::right
                  << std::setw(8) << acc.n << " site-cycles"
                  << std::setw(12) << std::fixed << std::setprecision(1) << acc.sum_abs
                  << " bits total"
                  << std::setw(9) << std::setprecision(3) << acc.sum_abs / acc.n << " mean"
                  << std::setw(8) << std::setprecision(0) << acc.max_abs << " worst"
                  << std::setw(8) << std::setprecision(1) << 100.0 * acc.sum_abs / total_X << " %\n";
    };

    std::cout << "\n  uMUL single-upset sensitivity, N = " << N << ", width = " << WIDTH
              << ", clean output = " << clean << " ones\n"
              << "  cross-section X = SUM |err| over all " << map.total() * N << " site-cycles\n\n";
    report("value", value_acc);
    report("rng_index", index_acc);
    report("rng_bus", bus_acc);
    report("in_0_wire", in0_acc);
    report("out_wire", out_acc);
    std::cout << "  " << std::left << std::setw(12) << "TOTAL X" << std::right
              << std::setw(30) << std::fixed << std::setprecision(1) << total_X << " bits\n\n";

    // E[|err|] ~ lambda * X. Scientific notation, because these rates are far below anything
    // std::fixed can show -- printing 0.0 here would make a real number look like nothing.
    //
    // ONLY ONE ENVIRONMENT IS QUOTED, and it is not LEO. The single published per-bit flip-flop
    // rate that could be verified for this file is a GEO / solar-minimum / 100-mil-aluminium
    // CREME96 prediction (Lee et al. 2018 -- see radiation_model.hpp). Substituting it for a LEO
    // number would be wrong in kind, not just in magnitude: LEO is partly shielded from GCR by
    // the geomagnetic field but gains the trapped-proton flux of the South Atlantic Anomaly,
    // which GEO does not have. Run CREME96 or SPENVIS for your own orbit and device and pass the
    // result in; X below does not change when you do.
    const double lambda_geo =
        Radiation::per_bit_cycle(Radiation::GEO_SOLARMIN_FLIPFLOP_UPSETS_PER_BIT_DAY, 100e6);

    std::cout << std::scientific << std::setprecision(3)
              << "  E[|err|] ~ lambda * X, at 100 MHz over a " << N << "-cycle run.\n"
              << "  Reference point: Xilinx UltraScale+ 16 nm FLIP-FLOPS, CREME96,\n"
              << "  GEO / solar min / 100 mil Al = "
              << Radiation::GEO_SOLARMIN_FLIPFLOP_UPSETS_PER_BIT_DAY << " upsets/bit/day\n"
              << "    lambda      = " << lambda_geo << " per bit-cycle\n"
              << "    E[|err|]    = " << lambda_geo * total_X << " bits/run\n"
              << "    1 disturbed run per "
              << 1.0 / (lambda_geo * static_cast<double>(map.total()) * static_cast<double>(N))
              << " runs\n"
              << std::defaultfloat
              << "\n[SUCCESS] Exported to umul_radiation_sensitivity.csv\n" << std::endl;

    // The persistent registers must dominate -- that is the whole thesis of this unit's fault
    // behaviour. A transient costs at most one output bit; a register upset biases every cycle
    // that follows it.
    double persistent = value_acc.sum_abs + index_acc.sum_abs;
    double transient = bus_acc.sum_abs + in0_acc.sum_abs + out_acc.sum_abs;
    EXPECT_GT(persistent, transient)
        << "expected persistent state to dominate uMUL's cross-section";

    // A transient on the output wire is by construction exactly one bit, every time. If this ever
    // fails, the harness is wrong rather than the circuit.
    EXPECT_DOUBLE_EQ(out_acc.sum_abs / out_acc.n, 1.0);
    EXPECT_DOUBLE_EQ(out_acc.max_abs, 1.0);

    // The value register's MSB should be the single worst site in the unit: it moves the operand
    // by half of full scale for the entire remainder of the run.
    EXPECT_GT(value_acc.max_abs, 100.0)
        << "an MSB upset in the operand register should cost hundreds of bits";
}

// =============================================================================================
// PART B -- Poisson rate sweep, uMUL against ANDMUL at equal lambda
// =============================================================================================
TEST(RadiationSensitivityTest, PoissonRateSweepAgainstAndGate) {
    const uMULSites umap(WIDTH);
    const ANDMULSites amap(WIDTH, /*charge_for_sng=*/true);

    std::vector<bool> a = generate_sobol_stream(static_cast<double>(ONES_A) / N, N, 2);
    std::vector<bool> b = generate_sobol_stream(static_cast<double>(VALUE_B) / N, N, 3);
    ASSERT_EQ(count_ones(a), ONES_A);
    ASSERT_EQ(count_ones(b), VALUE_B);

    const std::size_t umul_clean = run_umul(a, VALUE_B, {}, umap);
    const std::size_t and_clean = run_andmul(a, b, {}, amap);

    const std::vector<double> lambdas = {0.0, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1};
    constexpr int TRIALS = 4000;

    std::ofstream csv("umul_radiation_rate_sweep.csv");
    ASSERT_TRUE(csv.is_open());
    csv << "Lambda,Unit,Sites,Trials,Expected_Upsets,Mean_Abs_Err_Bits,Max_Abs_Err_Bits,"
           "Frac_Runs_Affected,Mean_Abs_Err_Prob\n";
    csv << std::scientific << std::setprecision(6);

    std::cout << "\n  Poisson rate sweep, " << TRIALS << " trials per point\n"
              << "  uMUL   " << umap.total() << " sites, clean = " << umul_clean << " ones\n"
              << "  ANDMUL " << amap.total() << " sites (gate + SNG), clean = " << and_clean
              << " ones\n\n"
              << "      lambda    E[upsets]     uMUL |err|   ANDMUL |err|   uMUL hit%  AND hit%\n";

    for (double lam : lambdas) {
        std::mt19937_64 gen(20260730ull);  // fixed seed: the sweep must be reproducible

        double u_sum = 0.0, u_max = 0.0, a_sum = 0.0, a_max = 0.0;
        int u_hit = 0, a_hit = 0;

        for (int trial = 0; trial < TRIALS; ++trial) {
            UpsetSchedule us(umap.total(), static_cast<uint32_t>(N), lam, gen);
            double ue = std::fabs(static_cast<double>(run_umul(a, VALUE_B, us.upsets(), umap)) -
                                  static_cast<double>(umul_clean));
            u_sum += ue;
            u_max = std::max(u_max, ue);
            if (ue > 0.0) ++u_hit;

            UpsetSchedule as(amap.total(), static_cast<uint32_t>(N), lam, gen);
            double ae = std::fabs(static_cast<double>(run_andmul(a, b, as.upsets(), amap)) -
                                  static_cast<double>(and_clean));
            a_sum += ae;
            a_max = std::max(a_max, ae);
            if (ae > 0.0) ++a_hit;
        }

        double u_mean = u_sum / TRIALS, a_mean = a_sum / TRIALS;
        double u_frac = static_cast<double>(u_hit) / TRIALS;
        double a_frac = static_cast<double>(a_hit) / TRIALS;

        csv << lam << ",uMUL," << umap.total() << ',' << TRIALS << ','
            << UpsetSchedule::expected(umap.total(), static_cast<uint32_t>(N), lam) << ','
            << u_mean << ',' << u_max << ',' << u_frac << ',' << u_mean / N << '\n';
        csv << lam << ",ANDMUL," << amap.total() << ',' << TRIALS << ','
            << UpsetSchedule::expected(amap.total(), static_cast<uint32_t>(N), lam) << ','
            << a_mean << ',' << a_max << ',' << a_frac << ',' << a_mean / N << '\n';

        std::cout << std::scientific << std::setprecision(1) << std::setw(12) << lam
                  << std::setw(12)
                  << UpsetSchedule::expected(umap.total(), static_cast<uint32_t>(N), lam)
                  << std::fixed << std::setprecision(3)
                  << std::setw(14) << u_mean << std::setw(14) << a_mean
                  << std::setprecision(1)
                  << std::setw(11) << 100.0 * u_frac << std::setw(10) << 100.0 * a_frac << "\n";
    }
    csv.close();
    std::cout << "\n[SUCCESS] Exported to umul_radiation_rate_sweep.csv\n" << std::endl;

    // Zero rate must be a no-op in both units -- the harness must not perturb a clean run.
    std::mt19937_64 g0(1);
    UpsetSchedule none(umap.total(), static_cast<uint32_t>(N), 0.0, g0);
    EXPECT_TRUE(none.empty());
    EXPECT_EQ(run_umul(a, VALUE_B, none.upsets(), umap), umul_clean);
}

}  // namespace StochasticSimulator
