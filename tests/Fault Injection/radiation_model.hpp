#pragma once

/**
 * radiation_model.hpp -- a single-event-upset arrival model for whole-unit fault campaigns.
 *
 * =============================================================================================
 * PROVENANCE. Every rate constant below carries the source it came from, inline, at the point of
 * definition. Where no defensible published number exists for a configuration, this header does
 * NOT supply one -- it requires the caller to pass a rate. An invented constant with a plausible
 * exponent is worse than no constant, because it survives into a results table looking sourced.
 *
 * An earlier revision of this file DID carry invented constants (LEO_QUIET / LEO_SAA /
 * SOLAR_EVENT at 1e-6 / 5e-5 / 1e-3 upsets/bit/day). They were order-of-magnitude guesses chosen
 * to make a sweep look reasonable, traceable to no device, orbit, shielding or solar epoch. They
 * have been removed rather than relabelled.
 * =============================================================================================
 *
 * HOW A BIT FLIPS, AND WHY THE PARTICLE MATTERS
 *
 *   particle crosses silicon
 *        v   ionization along the track, ~3.6 eV per electron-hole pair in Si
 *   e-h pairs
 *        v   drift + funnelling + diffusion at a reverse-biased junction
 *   collected charge Q_coll
 *        v   compare
 *   Q_coll > Q_crit  ->  the node flips.       Q_crit ~ C_node * V_dd
 *
 * Scaling shrinks both C_node and V_dd, so Q_crit falls and sensitivity rises every node.
 *
 * The three LEO sources do this by different routes, which is why a single "radiation rate" is a
 * fiction and the model below takes a rate per environment rather than one global number:
 *
 *   TRAPPED PROTONS (South Atlantic Anomaly) -- dominant in LEO. A proton is light and usually
 *     CANNOT deposit enough charge directly. It upsets INDIRECTLY: it strikes a Si nucleus, and
 *     the heavy, short-ranged, high-LET recoil fragments deposit the charge. At deeply scaled
 *     nodes with very low Q_crit, direct proton ionization also becomes possible.
 *   GALACTIC COSMIC RAYS -- heavy ions (C, O, Fe). DIRECT ionization, characterized by linear
 *     energy transfer (LET, MeV*cm^2/mg). Low flux, high consequence. Anti-correlated with the
 *     solar cycle: GCR flux peaks at solar MINIMUM.
 *   SOLAR ENERGETIC PARTICLES -- event-driven, mostly protons plus some heavy ions.
 *
 *   Trapped electrons dominate total ionizing dose but are not normally an SEU source.
 *
 * ALPHA PARTICLES ARE NOT A SPACE-ENVIRONMENT TERM. They come from radioactive contamination in
 *   the PACKAGE -- Pb-210 -> Bi-210 -> Po-210 emitting a 5.3 MeV alpha in lead-bearing solder
 *   bumps, plus Th/U traces in mould compound, bond wire and lead frames. The rate is identical
 *   in orbit and on a bench, so it belongs in this model as a constant floor independent of
 *   orbit, never folded into the orbital rate. Conversely the terrestrial NEUTRON term (the
 *   atmospheric cascade, and the B-10(n,alpha) reaction in BPSG) DISAPPEARS in orbit, because
 *   there is no atmosphere to make the cascade. Going to LEO you drop neutrons, keep alphas
 *   unchanged, and add trapped protons plus GCR.
 *
 * WHY "DAY / NIGHT" IS THE WRONG AXIS. LEO upset rate varies GEOGRAPHICALLY, not diurnally. The
 *   South Atlantic Anomaly is a fixed region where the inner proton belt dips to low altitude,
 *   and an orbit crosses it several times a day, so the rate is bursty on the ~90 minute orbital
 *   period. See SAA_UPSET_FRACTION below.
 *
 * =============================================================================================
 * THE COMPUTATIONAL TRICK. A real per-bit-per-cycle rate is around 1e-19. Rolling a die per site
 * per cycle at that rate is not an experiment; it is a slow way to print "no error". Two things
 * make it tractable:
 *
 *   1. DRAW THE COUNT, THEN PLACE IT. Upsets over a run are Binomial(sites * cycles, lambda).
 *      Draw that once and scatter them over the (site, cycle) grid. Exact -- this IS the
 *      per-site-per-cycle process, not an approximation of it -- at cost O(upsets).
 *
 *   2. SWEEP LAMBDA, EXTRAPOLATE. To first order the expected error is linear in lambda:
 *          E[|error|]  ~=  lambda * SUM over (site, cycle) of |error given one upset there|
 *      That sum is a property of the CIRCUIT, not the orbit. Characterize it once exhaustively,
 *      and every orbit, shielding and clock-rate question becomes arithmetic. The Monte Carlo
 *      exists to find where the linearity BREAKS, which is where multi-bit effects begin.
 * =============================================================================================
 */

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

namespace StochasticSimulator {
namespace Radiation {

// ---- Unit conversion ------------------------------------------------------------------------

/** Datasheet units (upsets per bit per day) into the per-bit, per-cycle rate this model wants. */
inline double per_bit_cycle(double upsets_per_bit_day, double clock_hz) {
    return upsets_per_bit_day / (clock_hz * 86400.0);
}

/**
 * FIT into upsets/bit/day. 1 FIT = 1 failure per 1e9 device-hours, the convention used by
 * JEDEC JESD89A and by vendor reliability reports (e.g. AMD/Xilinx UG116, which quotes
 * configuration-memory and BlockRAM soft error rates in FIT/Mb corrected to New York City sea
 * level). Divide a FIT/Mb figure by 1.048576e6 first to get FIT/bit.
 */
inline double per_bit_day_from_fit(double fit_per_bit) {
    return fit_per_bit * 1.0e-9 * 24.0;
}

// ---- Published reference points ---------------------------------------------------------------
//
// These are the only rate numbers in this file, and each states its exact configuration. None of
// them is a LEO figure, because no LEO number could be verified from an accessible source while
// writing this. Do NOT silently reuse a GEO rate for a LEO study: the environments differ in
// kind, not just magnitude -- LEO is shielded from much of the GCR flux by the geomagnetic field
// but gains the trapped-proton contribution of the SAA, which GEO does not have.

/**
 * Xilinx Kintex UltraScale+ (16 nm FinFET), CREME96 prediction from measured heavy-ion data.
 * Configuration: GEOSTATIONARY orbit, SOLAR MINIMUM, 100 mils (0.1 inch) aluminium shielding.
 *
 * SOURCE: D. S. Lee et al., "Single-Event Characterization of 16 nm FinFET Xilinx UltraScale+
 *         Devices with Heavy Ion and Neutron Irradiation," IEEE NSREC/RADECS 2018.
 *         https://www.osti.gov/biblio/1570815   (values also tabulated in the open-access review
 *         at https://pmc.ncbi.nlm.nih.gov/articles/PMC11360524/ )
 *
 * The FLIP-FLOP number is the relevant one for uMUL and SuMUL: both units are made of registers,
 * not of configuration cells or block RAM. Note it is ~2000x the configuration-RAM rate -- which
 * resource you are modelling matters more than which orbit.
 */
constexpr double GEO_SOLARMIN_FLIPFLOP_UPSETS_PER_BIT_DAY = 2.08e-8;
constexpr double GEO_SOLARMIN_BLOCKRAM_UPSETS_PER_BIT_DAY = 5.29e-9;
constexpr double GEO_SOLARMIN_CONFIGRAM_UPSETS_PER_BIT_DAY = 9.18e-12;

/**
 * Fraction of LEO single-event upsets that occur inside the South Atlantic Anomaly, dominated by
 * trapped protons. The practical consequence: an orbit-AVERAGED lambda badly misrepresents the
 * physics, because the upsets are not spread evenly over the orbit -- they arrive in bursts a few
 * times a day. Modelling a mission as one uniform rate will understate correlated,
 * closely-spaced upsets and therefore understate multi-upset behaviour.
 *
 * SOURCE: Wu et al., "Analyzing SEU Rate in LEO Satellites During the Space Weather Event of
 *         May 2024," Space Weather, 2025. https://doi.org/10.1029/2025SW004608
 */
constexpr double SAA_UPSET_FRACTION = 0.90;

/**
 * Package alpha emissivity, in alphas per cm^2 per hour. ORBIT-INDEPENDENT: this is contamination
 * in the package, so the flux is the same in LEO as on a bench. Qualified ultra-low-alpha (ULA)
 * materials cut package-attributable soft error rate by 2-4 orders of magnitude relative to
 * untreated materials, which is what moved the dominant terrestrial source to cosmic neutrons in
 * the mid-1990s.
 *
 * These are EMISSIVITIES, not upset rates. Converting to upsets/bit/day needs the device's alpha
 * cross-section, which is part-specific -- see alpha_upsets_per_bit_day() below.
 *
 * SOURCE: JEDEC JESD89A, "Measurement and Reporting of Alpha Particle and Terrestrial Cosmic
 *         Ray-Induced Soft Errors in Semiconductor Devices."
 *         http://gpc.pnpi.nrcki.ru/images/files/docs/JEDEC_Standart.pdf
 *         Also: Kumar et al., "Soft Error Issue and Importance of Low Alpha Solders,"
 *         https://www.ipme.ru/e-journals/RAMS/no_23413/07_23413_kumar.pdf
 */
constexpr double ALPHA_EMISSIVITY_UNTREATED_PER_CM2_HR = 100.0;  // untreated package materials
constexpr double ALPHA_EMISSIVITY_ULA_PER_CM2_HR = 0.002;        // qualified ULA grade

/**
 * CREME96 is the standard tool for turning an orbit plus a measured device cross-section into an
 * upset rate, and its "worst week" and "worst day" scenarios are built from the measured October
 * 1989 solar particle events. There is deliberately no SOLAR_EVENT constant in this header: the
 * upset enhancement during an event is not a property of the event alone, it is the convolution
 * of the event's particle spectrum with YOUR device's cross-section curve. A part with a high LET
 * threshold can be nearly indifferent to a proton-rich event that devastates a low-Q_crit part.
 * Run CREME96 or SPENVIS for your device and orbit and pass the result in.
 *
 * SOURCE: A. J. Tylka et al., "CREME96: A Revision of the Cosmic Ray Effects on Micro-Electronics
 *         Code," IEEE Trans. Nucl. Sci. 44(6), 1997.
 */

// ---- Deriving a rate from device data --------------------------------------------------------

/**
 * First-order upset rate from a saturated cross-section and an integral particle flux. This is
 * the crude version -- the correct calculation integrates the Weibull cross-section curve
 * sigma(LET) against the differential LET spectrum, which is what CREME96 does. Use this only to
 * sanity-check a CREME96 result or to explore sensitivity, never as the published number.
 *
 * @param sigma_sat_cm2_per_bit  Saturated cross-section from accelerator testing.
 * @param integral_flux_per_cm2_s  Flux of particles above the device's LET threshold.
 */
inline double upsets_per_bit_day_from_cross_section(double sigma_sat_cm2_per_bit,
                                                    double integral_flux_per_cm2_s) {
    return sigma_sat_cm2_per_bit * integral_flux_per_cm2_s * 86400.0;
}

/**
 * The package-alpha contribution, kept separate from the orbital term because it does not vary
 * with orbit. Emissivity is per cm^2 per hour; the cross-section is the device's alpha-specific
 * one, which is NOT the same as its heavy-ion saturated cross-section (alpha LET is only about
 * 0.5-1.5 MeV*cm^2/mg, so only low-Q_crit cells respond at all).
 */
inline double alpha_upsets_per_bit_day(double emissivity_per_cm2_hr,
                                       double alpha_sigma_cm2_per_bit) {
    return emissivity_per_cm2_hr * alpha_sigma_cm2_per_bit * 24.0;
}

// ---- The arrival process ----------------------------------------------------------------------

struct Upset {
    uint32_t site;
    uint32_t cycle;
};

/**
 * One run's worth of upsets. Construction draws Binomial(sites * cycles, lambda) STRIKES and
 * places them uniformly; each strike then flips one or more sites depending on the multi-cell
 * settings. Sorted by cycle so a simulation can walk them in step with the clock. Repeats are
 * allowed and are physical: a second strike on the same flip-flop flips it back.
 *
 * MULTI-CELL UPSETS. One particle can flip several ADJACENT cells, by a track crossing multiple
 * wells or by charge sharing. The fraction rises with scaling, with LET, and sharply with grazing
 * incidence angle. MCU matters here specifically because the operand register's bits are ordered
 * by significance: an MCU flipping bits 4 and 5 together is a single-particle event worth
 * -48/1024 of operand, which an independent-single-bit model can only reach by coincidence.
 *
 * There is no default MCU fraction, because a defensible one is device- and angle-specific and
 * comes from a test campaign. mcu_probability = 0 (the default) is an explicit SBU-only model,
 * and the results should say so.
 */
class UpsetSchedule {
public:
    struct Options {
        double mcu_probability = 0.0;   // P(a strike flips more than one site)
        uint32_t mcu_max_extra = 1;     // additional adjacent sites flipped when it does
    };

    UpsetSchedule(uint32_t sites, uint32_t cycles, double lambda, std::mt19937_64& gen,
                  Options opts = Options{})
        : sites_(sites), cycles_(cycles) {
        const double site_cycles = static_cast<double>(sites) * static_cast<double>(cycles);
        if (lambda <= 0.0 || sites == 0 || cycles == 0) return;

        // Binomial is exact here. It stays well-behaved at both ends of a wide lambda sweep,
        // where a Poisson approximation would drift once lambda is no longer small.
        std::binomial_distribution<long long> how_many(
            static_cast<long long>(site_cycles), std::min(lambda, 1.0));
        long long strikes = how_many(gen);
        if (strikes <= 0) return;

        std::uniform_int_distribution<uint32_t> which_site(0, sites - 1);
        std::uniform_int_distribution<uint32_t> which_cycle(0, cycles - 1);
        std::uniform_real_distribution<double> unit(0.0, 1.0);

        list_.reserve(static_cast<std::size_t>(strikes));
        for (long long i = 0; i < strikes; ++i) {
            uint32_t s = which_site(gen);
            uint32_t c = which_cycle(gen);
            list_.push_back(Upset{s, c});
            ++strike_count_;

            // A multi-cell strike flips physically ADJACENT sites in the SAME cycle. Adjacency in
            // site index is the model's stand-in for adjacency in silicon; that is only as good
            // as the site map's ordering, which for the register files here is bit significance.
            if (opts.mcu_probability > 0.0 && unit(gen) < opts.mcu_probability) {
                for (uint32_t extra = 1; extra <= opts.mcu_max_extra; ++extra) {
                    uint32_t neighbour = s + extra;
                    if (neighbour < sites) {
                        list_.push_back(Upset{neighbour, c});
                        ++mcu_bits_;
                    }
                }
            }
        }
        std::sort(list_.begin(), list_.end(),
                  [](const Upset& a, const Upset& b) { return a.cycle < b.cycle; });
    }

    const std::vector<Upset>& upsets() const { return list_; }
    std::size_t count() const { return list_.size(); }       // total flipped sites
    std::size_t strikes() const { return strike_count_; }    // particles
    std::size_t mcu_extra_bits() const { return mcu_bits_; }
    bool empty() const { return list_.empty(); }
    uint32_t sites() const { return sites_; }
    uint32_t cycles() const { return cycles_; }

    /** Expected number of STRIKES for these dimensions at this rate -- the analytic check. */
    static double expected(uint32_t sites, uint32_t cycles, double lambda) {
        return static_cast<double>(sites) * static_cast<double>(cycles) * lambda;
    }

private:
    uint32_t sites_;
    uint32_t cycles_;
    std::vector<Upset> list_;
    std::size_t strike_count_ = 0;
    std::size_t mcu_bits_ = 0;
};

// ---- Site maps --------------------------------------------------------------------------------

/**
 * uMUL's complete radiation cross-section, for a unit of the given width.
 *
 *      RANGE            SITE                    LIFETIME     COUNT
 *      [0,   w)         value register bit      permanent    w
 *      [w,  2w)         rng_index register bit  permanent    w
 *      [2w, 3w)         RNG output bus wire     one cycle    w
 *      3w               in_0 / enable wire      one cycle    1
 *      3w + 1           output wire             one cycle    1
 *                                                            ----
 *                                                            3w + 2
 *
 * Sites within the two register ranges are ordered by BIT SIGNIFICANCE, which is what makes the
 * MCU adjacency model above meaningful -- and also what makes an MCU worse than two independent
 * upsets, since adjacent significant bits compound.
 *
 * The Sobol direction vectors are not sites: they are constants, and in hardware they are ROM or
 * tied-off wiring rather than flip-flops. Add them if your fault model covers configuration
 * memory as well as datapath -- and note that configuration cells have a very different
 * per-bit rate from flip-flops (see the Lee 2018 numbers above, which differ by ~2000x).
 */
struct uMULSites {
    unsigned width;

    explicit uMULSites(unsigned w) : width(w) {}

    uint32_t total() const { return 3u * width + 2u; }
    uint32_t value_base() const { return 0u; }
    uint32_t index_base() const { return width; }
    uint32_t bus_base() const { return 2u * width; }
    uint32_t in0_wire() const { return 3u * width; }
    uint32_t out_wire() const { return 3u * width + 1u; }

    bool is_value(uint32_t s) const { return s < width; }
    bool is_index(uint32_t s) const { return s >= width && s < 2u * width; }
    bool is_bus(uint32_t s)   const { return s >= 2u * width && s < 3u * width; }
    bool is_in0(uint32_t s)   const { return s == in0_wire(); }
    bool is_out(uint32_t s)   const { return s == out_wire(); }

    /** Persistent sites keep their damage; the rest glitch for exactly one cycle. */
    bool is_persistent(uint32_t s) const { return s < 2u * width; }

    const char* class_name(uint32_t s) const {
        if (is_value(s)) return "value";
        if (is_index(s)) return "rng_index";
        if (is_bus(s))   return "rng_bus";
        if (is_in0(s))   return "in_0_wire";
        return "out_wire";
    }
};

/**
 * The plain AND-gate multiplier's cross-section. It is STATELESS, so every site is a wire and
 * every upset costs exactly one output bit and is then gone.
 *
 *      0   in_0 wire        1   in_1 wire        2   output wire
 *
 * FAIRNESS WARNING. Three sites is the gate ALONE. uMUL's in_1 is a loaded register, but
 * ANDMUL's in_1 is a stream that something has to generate -- an SNG with its own RNG state and
 * threshold register, which is the same hardware uMUL carries inside G. Comparing 3 sites against
 * 3w + 2 credits ANDMUL with an operand that materialises from nowhere. Use include_sng = true to
 * charge it for the generator, which is the comparison that means something at the system level.
 */
struct ANDMULSites {
    unsigned sng_width;    // width of the SNG that produces in_1; 0 if you are not charging for it
    bool include_sng;

    explicit ANDMULSites(unsigned w = 0, bool charge_for_sng = false)
        : sng_width(w), include_sng(charge_for_sng) {}

    uint32_t total() const { return include_sng ? 3u + 2u * sng_width : 3u; }
    uint32_t in0_wire() const { return 0u; }
    uint32_t in1_wire() const { return 1u; }
    uint32_t out_wire() const { return 2u; }

    bool is_in0(uint32_t s) const { return s == 0u; }
    bool is_in1(uint32_t s) const { return s == 1u; }
    bool is_out(uint32_t s) const { return s == 2u; }
    bool is_sng(uint32_t s) const { return s >= 3u; }

    const char* class_name(uint32_t s) const {
        if (is_in0(s)) return "in_0_wire";
        if (is_in1(s)) return "in_1_wire";
        if (is_out(s)) return "out_wire";
        return "sng";
    }
};

}  // namespace Radiation
}  // namespace StochasticSimulator
