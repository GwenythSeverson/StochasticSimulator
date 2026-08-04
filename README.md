# StochasticSimulator
A C++ based stochastic computing simulator. 

Organization in modules is a C++ representation of the architecture of each kind of specified part. If an implimentation of a specific adder is missing, add it in the adder files. testing files live under tests and are named as the main testing point. 

GUnit/Gtests will be implimented, with a current plan to create structs that hold trial records, bit stream records, and maybe a function that records a snapshot at every step to put into a CSV to generate grpahs with in MATLAB. 

## Setup
 install CMake https://cmake.org/download/

## Tests
```powershell
# 1. Configure the build environment
cmake -B build -S .

# 2. Compile the simulator and test targets
cmake --build build --config Release

# 3. Run the verification test suite
.\build\Release\stochastic_computer.exe
```

## File explanations/ progress

Bit Stream Generator (bsg) folder-
    lfsr- Linear Feedback Shift Register
        - 2 mode galois lfsr, widths 7..16 (StreamLength::Length_128 .. Length_65536)
            - mode 0 (default seed): auto random for multi run accuracy tests and avoided correlation bais
            - mode hardware seed: deterministic mode for debugging, fault injection, and early termination tests
    sng- Stochastic Number Generator
        - the hardware comparator itself: generate_bit(p, random_val, max_scale), truncating
          threshold with a strict ">". everything else in the project matches this convention
          so LFSR and Sobol streams stay directly comparable.
    sobol- low discrepancy generator, drop in replacement for the lfsr
        - progressively stratified instead of pseudo random, so truncating a stream early still
          leaves a good estimate of p instead of a random walk error. this is what makes early
          termination work.
        - 19 tabulated dimensions, up to 32 bits wide. give every independent operand its OWN
          dimension or the streams are perfectly correlated and any multiplier/ZCE number is junk.
          asking for a dimension a narrow generator cant tell apart throws instead of silently
          handing back a duplicate stream.
        - random access (value_at / jump_to) so fault injection runs can be replayed exactly.
    fault_injector- stream level bit flipping: one specific bit, all bits, N random bits, and
        clumped MBU strikes at fixed or random locations (single clump or multi clump for
        solar flare style events). TODO in the file: a poisson/rate based version. the per unit
        register level injection lives in the modules themselves now, not here.

general_functions- shared helpers: calculate_probability, count_ones, generate_valid_stream
    (enforces the 0.05% error bound), format_bit_vector, can_terminate_early, and Julie Hsiao's
    Zero Correlation Error (calculate_zce_window).

modules folder-
    adder- MUX scaled/weighted adder. stateless, one select bit per cycle. counter based adder
        still not implemented.
    uSADD- uGEMM scaled adder (parallel counter + accumulator, Fig 3c). STATEFUL.
        - out = (p0 + ... + pn-1)/n exactly. no RNG, no estimator, no select stream. the only
          healthy error is the leftover residue, and its always a deficit.
        - the accumulator holds a RESIDUE not an operand, so an upset there is bounded and one
          time. this is the sharp contrast with uMUL and its the whole poster thesis.
        - the collapse: the answer depends only on TOTAL ones across the inputs, not where they
          sit, which is what makes it exhaustively testable.
    multiplier- plain AND gate. completely stateless, unlike the adders and dividers.
        only correct if the two streams are uncorrelated - AND(s,s) = p, not p^2.
    SuMUL- uGEMM Fig 3(a) with a SLIDING WINDOW counter estimating p1. STATEFUL.
        - kills correlation error but pays warm up, since the 2^width window boots half full.
          width is a real tradeoff (warm up vs precision); width_for_stream_length lands the
          window near 2*sqrt(N).
        - window cells are SELF SCRUBBING under fault, an upset shifts off the end on its own.
        - currently commented out of the functional test build in CMakeLists for speed.
    uMUL- uGEMM Fig 3(a) LOADED form, matching the synthesized uMUL_uni.sv. STATEFUL.
        - in_0 is a stream and the enable, in_1 is a binary value loaded into a register. no
          estimator anywhere so no warm up and no correlation error, and it beats the AND gate
          even on the AND gate's best case (measured RMSE 0.0007 vs 0.0011 at N=1024).
        - p1 = 1.0 is NOT representable, the register tops out at (2^width - 1)/2^width. thats
          the RTL's behavior, not a bug.
        - both registers are PERSISTENT, nothing self scrubs. an MSB flip in the value register
          halves or doubles the operand for the whole run.
    divider- Gaines division (GDIV) unipolar, up/down counter feedback loop / ADDIE.
        - counter depth is the knob (settling vs precision), auto sizes near sqrt(N).
        - unipolar so X must be <= Y, otherwise the counter just pins at the top. throws on bad
          input now instead of returning 0.0, since 0.0 is a legitimate quotient.

tests-
    Unit Tests- fast and structural, whole section is ~0.2s. run these constantly.
        - test_adders- truth table of the mux select signal.
            - TODO test weighted/scaled addition functionality with preset streams
        - test_multipliers- truth table of the and gate multiplier.
            - TODO add preset stream test with loop demo
        - test_dividers- convergence of the saturating counter divider. seeded, and the
          tolerances come from a measured px,py sweep at 3 seeds rather than guesses.
        - test_SuMUL / test_uMUL- split deliberately into STRUCTURAL (exact circuit assertions,
          catches a rewiring bug) and ACCURACY (statistical, explicit LFSR seeds so theyre not
          flaky). test_uMUL also carries the CharacterizeAccuracy sweep that regenerates the
          numbers quoted in uMUL.hpp.
        - test_uSADD- exact arithmetic, register sizing, the collapse property, and the fault
          injection hooks.
    Functional Tests- accuracy sweeps and characterization. slow, and several write CSVs into the
      working directory.
        - multiplier/ANDmul/test_multiplier_accuracy- **includes correlation error in results
            - creates csv for matlab graphs and later analysis
        - multiplier/SuMUL/test_SuMUL_accuracy- exhaustive, every one of the 1023 distinct
          Length_1024 LFSR streams for a target pair. -> sumul_exhaustive_trials.csv
        - multiplier/uMUL/test_uMUL_behavior- exhaustive bit level sweep, every representable
          operand pair at N=1024. ~41s, -> umul_behavior_exhaustive.csv (29 MB)
        - multiplier/uMUL/test_uMUL_convergence- early termination at three operating points
          plus an all pairs sweep. ~107s, -> umul_convergence_{low,mid,high,allpairs}.csv
        - adder/uSADD/test_uSADD_behavior- exhaustive accuracy (all 257x257 count pairs at
          N=256), uSADD vs MUX on identical operands under three correlation regimes, and early
          termination. -> usadd_behavior.csv, usadd_vs_mux.csv, usadd_early_termination.csv
        - gen_fault_tests- the old general fault test. lives here but belongs with fault
          injection by purpose. dont touch.
        - the .m files next to each test are the matlab side (early termination readers,
          divider initialization comparison, uMUL behavior reader).
    Fault Injection- single event upset campaigns, the slow CSV heavy ones.
        - radiation_model.hpp- the arrival model + site maps. every rate constant carries its
          published source inline, and where no defensible number exists it makes you pass one
          (the old invented LEO constants were removed, not relabelled). draws
          Binomial(sites*cycles, lambda) once and places the strikes, which is exact and cheap.
          has optional multi cell upset support. site maps for uMUL and the AND gate, including
          the fairness note about charging the AND gate for the SNG that feeds it.
        - multiplier/test_uMUL_radiation- whole unit campaign, uMUL vs the AND gate, faults land
          on registers and wires not just the stream. -> umul_radiation_sensitivity.csv,
          umul_radiation_rate_sweep.csv
        - divider/test_divider_fault2- exhaustive 32 bit divider sweep, ~1s.
          -> 32bit_exhaustive_divider_trials.csv
        - currently commented out of the build: adder/test_adder_fault2,
          multiplier/test_multiplier_fault2, multiplier/multiplier_injections, experiment2.
          re-enable them individually in CMakeLists as needed.

ECE_SPARK2026_poster_code- the SPARK 2026 poster. four units, one operand pair each, N=256,
    exhaustive fault enumeration. x axis is BITS FLIPPED everywhere, no orbit and no upset rate
    attached to anything on purpose.
    - Addition/MUX_add and Addition/uSADD_add, Multiplication/AND_mul and Multiplication/uMUL_mul.
      each folder has its .cpp, its CSVs, its PNGs, its plot_*.m, and a HOW_THIS_DATA_WAS_MADE.txt
      with the trial ledger and the honest limitations. read those before answering questions at
      the poster.
    - these write their CSVs next to their own source via __FILE__, not to the working directory,
      so the poster folder stays self contained.
    - POSTER_PLAN.txt is the thesis, the measured numbers, which five figures go up, and a
      "what not to claim" list.
    - regenerate with:
      `--gtest_filter=AndMulPoster.*:uMULPoster.*:MuxAddPoster.*:uSADDAddPoster.*` then run the
      four plot_*.m scripts in matlab.

keepdata- saved results worth keeping, separate from whatever the last run happened to dump.
    Adder Data-
        32 bit exhaustive adder trials + the fault analysis graph, and notes on the fault
        intensity error distribution
    Divider Data-
        CNT_starting_point_data-
            did 2 trial runs of .25/.5 with the counter initialized at different spots (0 & 16)
            and decided that the warmup wasnt incredibly effected. the
            starting_counter_doesnt_matter folder is the 64 and 256 bit comparison that settled it
        fault testing gen behavior- 32 bit exhaustive divider trials
    Multiplier Data/AND Mul-
        basic mult behavior graphs and the csvs they were pulled from, the exhaustive and ZCE
        trial data, the .mat, and bitflip_mult_graph_explanation.txt for reading the fault graphs
