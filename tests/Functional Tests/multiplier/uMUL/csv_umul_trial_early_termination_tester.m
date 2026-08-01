%% csv_umul_trial_early_termination_tester.m
%
% uMUL behaviour / convergence check. Loads ALL THREE convergence CSVs written by
% test_uMUL_convergence.cpp and lays every one of the six charts out in a single 2x3 figure, so
% the LOW / MID / HIGH operating points can be compared side by side.
%
%     umul_convergence_low.csv    0.5000 x 0.2002 = 0.1001    LOW output
%     umul_convergence_mid.csv    0.5000 x 0.8750 = 0.4375    MID output
%     umul_convergence_high.csv   0.9375 x 0.9600 = 0.9000    HIGH output
%
% ===========================================================================================
% THIS IS AN EXHAUSTIVE SAMPLE, NOT A SURVEY. The primary columns (Mean_Est, Mean_AbsError,
% Min_Est, Max_Est, Std_Est, P01_Est, P99_Est) are computed EXACTLY over every one of the
% ~10^307 ways in_0's ones can be arranged -- not sampled from them.
%
% That is possible because the output ones count after k enabled cycles is
%       out(k) = #{ j < k : sobol[j] < value }
% which depends only on k and never on WHERE in_0's ones sit. So truncating at N, the only thing
% that varies across arrangements is how many of in_0's ones fell in the first N cycles, and that
% integer is hypergeometrically distributed. Weighting out(k)/N by that pmf and summing over
% k = 0..N is the exact expectation over all arrangements, in at most 1025 terms.
%
% The Seed_* columns are the OLD 1023-LFSR-seed sample, kept so the sampling bias is visible
% rather than assumed away. Both are plotted: the exhaustive result solid, the seed sample as a
% thin dotted line. Where they separate, the seed sweep was lying.
% ===========================================================================================
%
% LAYOUT
%   Top row     BEHAVIOUR   mean estimate vs N, semilog-x, against the true product.
%   Bottom row  CONVERGENCE mean absolute error vs N, LOG-LOG, against a 1/sqrt(N) reference.
%
% WHY THE BOTTOM ROW IS LOG-LOG: a 1/sqrt(N) decay is a STRAIGHT LINE of slope -1/2 on log-log
% and a curve on semilog-x. The dashed reference makes the comparison direct -- where the measured
% curve bends BELOW it, uMUL is converging faster than a random estimator, which is the Sobol
% low-discrepancy prefix doing its job.
%
% THE SHADED BAND IS P01..P99, NOT MIN..MAX. Over all 10^307 arrangements the true support
% extremes are reachable but astronomically improbable (they need every one of in_0's ones to fall
% inside the truncation window), so a min-max band would be nearly full-scale and say nothing.
% The 1st-to-99th percentile band is what an actual run will sit inside.
%
% WHAT TO LOOK FOR, and it differs by level. Short truncations do not measure the multiplier, they
% measure how well G's Sobol PREFIX can express p1 in the few enabled cycles that have happened.
% The first eight Sobol points at width 10 are 0, 512, 768, 256, 384, 896, 640, 128, so:
%     LOW  (value 205): only j = 0 and j = 7 clear it -> 2/8 = 0.25 against a target of 0.20.
%     MID  (value 896): j = 5 is the only miss -> 7/8 = 0.875, exactly the target, from the moment
%                       eight enabled cycles exist. Locks on and never moves.
%     HIGH (value 983): all eight clear it -> 8/8 = 1.0 against 0.96. Reads p0 until the first
%                       miss arrives, then steps down.
% None of that is warm-up -- there is no estimator converging on anything -- it is quantization.
% The DimAvg_* columns average the exact result over 8 independent Sobol dimensions, which is what
% an ARRAY of these units looks like and separates "how uMUL behaves" from "what dimension 1 does".

clear; clc; close all;

%% ------------------------------------------------------------------
%% 1. CONFIGURATION -- all three files, plotted together
%% ------------------------------------------------------------------
csvNames = [ "umul_convergence_low.csv", ...
             "umul_convergence_mid.csv", ...
             "umul_convergence_high.csv" ];

levelColours = [0.20 0.45 0.85;    % LOW   blue
                0.47 0.67 0.19;    % MID   green
                0.85 0.33 0.10];   % HIGH  orange

% Walk up from this script until we find the folder containing 'keepdata'.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = scriptDir;
while ~isfolder(fullfile(repoRoot, 'keepdata'))
    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
        error("Could not locate the repository root (no 'keepdata' folder found above %s).", scriptDir);
    end
    repoRoot = parentDir;
end

%% ------------------------------------------------------------------
%% 2. LOAD ALL THREE
%% ------------------------------------------------------------------
nLevels = numel(csvNames);
D = struct([]);

for c = 1:nLevels
    csvName = csvNames(c);
    candidates = [ fullfile(repoRoot, 'keepdata', 'Multiplier Data', csvName), ...
                   fullfile(repoRoot, csvName) ];
    csvFile = "";
    for k = 1:numel(candidates)
        if isfile(candidates(k)), csvFile = candidates(k); break; end
    end
    if csvFile == ""
        error("Could not find %s in either:\n  %s\n  %s", csvName, candidates(1), candidates(2));
    end

    T = readtable(csvFile, "Delimiter", ",", "TextType", "string");

    [~, stem] = fileparts(csvName);
    D(c).level     = upper(erase(stem, "umul_convergence_"));
    D(c).N         = double(T.N);
    D(c).mean      = double(T.Mean_Est);          % EXACT over all arrangements
    D(c).absErr    = double(T.Mean_AbsError);     % EXACT
    D(c).p01       = double(T.P01_Est);
    D(c).p99       = double(T.P99_Est);
    D(c).seedMean  = double(T.Seed_Mean_Est);     % 1023-seed sample, for comparison
    D(c).seedAbs   = double(T.Seed_Mean_AbsError);
    D(c).dimAbs    = double(T.DimAvg_Mean_AbsError);
    D(c).mcMean    = double(T.MC_Mean_Est);       % random-shuffle validation of the exact result
    D(c).ideal     = double(T.Ideal_Product(1));
    D(c).log10Arr  = double(T.Log10_Arrangements(1));
    D(c).seeds     = double(T.Seed_Samples(1));
    D(c).dims      = double(T.Dims(1));
    D(c).mcTrials  = double(T.MC_Trials(1));
    D(c).operation = string(T.Operation(1));

    biasPct = 100 * max(abs(D(c).seedAbs - D(c).absErr) ./ max(D(c).absErr, eps));
    mcGap   = max(abs(D(c).mcMean - D(c).mean));
    % Single-quoted concatenation: ["a" "b"] would build a string ARRAY, and fprintf would read
    % its first element as a file identifier rather than as the format string.
    fprintf(['[INFO] %-5s exact over 10^%.0f arrangements | %d-dim avg | ' ...
             'worst seed bias %.1f%% | MC(%d) agrees to %.1e\n'], ...
            D(c).level, D(c).log10Arr, D(c).dims, biasPct, D(c).mcTrials, mcGap);
end

%% ------------------------------------------------------------------
%% 3. ONE FIGURE, SIX PANELS
%% ------------------------------------------------------------------
fig = figure("Name", "uMUL Behaviour and Convergence -- LOW / MID / HIGH", ...
             "Position", [60, 60, 1720, 940], "Color", "w");
tl = tiledlayout(fig, 2, nLevels, "TileSpacing", "compact", "Padding", "compact");

title(tl, "uMUL: output behaviour and convergence, EXHAUSTIVE over every in_0 arrangement", ...
      "FontWeight", "bold", "FontSize", 14);
subtitle(tl, sprintf(['exact over all 10^{%.0f} arrangements (hypergeometric, not sampled), ' ...
                      'cross-checked against %d random shuffles   |   ' ...
                      'dotted = old %d-seed sample   |   in_1 is a loaded binary value'], ...
                     D(1).log10Arr, D(1).mcTrials, D(1).seeds), "FontSize", 10);

%% ---- TOP ROW: behaviour at three operating points ----
for c = 1:nLevels
    ax = nexttile(tl, c);
    col = levelColours(min(c, size(levelColours,1)), :);
    N = D(c).N;

    fill(ax, [N; flipud(N)], [D(c).p01; flipud(D(c).p99)], col, ...
         "FaceAlpha", 0.15, "EdgeColor", "none");
    hold(ax, "on");
    plot(ax, N, D(c).mean, "-o", "LineWidth", 2.0, "MarkerSize", 6, ...
         "Color", col, "MarkerFaceColor", col);
    plot(ax, N, D(c).seedMean, ":", "LineWidth", 1.3, "Color", [0.35 0.35 0.35]);
    yline(ax, D(c).ideal, "--r", "LineWidth", 1.6, ...
          "Label", sprintf("%.4f", D(c).ideal), "LabelHorizontalAlignment", "left");
    hold(ax, "off");

    set(ax, "XScale", "log");
    grid(ax, "on");
    xticks(ax, N); xticklabels(ax, string(N));
    ylim(ax, [0 1]);
    xlabel(ax, "Stream Length N [log]");
    if c == 1, ylabel(ax, "Mean Output Probability"); end
    % TRIAL COUNTS, stated on every panel as a title line rather than a floating box: a box
    % large enough to be legible collides with the data in the HIGH panel, where the curve sits
    % at 0.94. Any wobble in these curves is NOT a sampling artefact -- the arrangement average
    % is exact and the Monte Carlo confirms it -- so a reader can attribute the kinks to the
    % Sobol prefix, which is the circuit itself.
    title(ax, {sprintf("%s  --  %s", D(c).level, extractOperands(D(c).operation)), ...
               sprintf("ALL 10^{%.0f} arrangements (exact)  |  MC %s  |  %d Sobol dims  |  %d seeds", ...
                       D(c).log10Arr, addThousandsSeparator(D(c).mcTrials), ...
                       D(c).dims, D(c).seeds)});
    ax.Title.FontSize = 10;
    ax.Toolbar.Visible = "off";

    if c == 1
        legend(ax, ["P01-P99 band (exact)", "mean (exact)", "1023-seed sample", "true product"], ...
               "Location", "southeast", "FontSize", 8);
    end
end

%% ---- BOTTOM: ONE large panel -- the WHOLE UNIT, every operand against every operand ----
%
% The three panels above are three points in a space of 1,049,600 operand pairs. This one sweeps
% all of them, each averaged exactly over all ~10^307 in_0 arrangements, and reports the unit's
% convergence rather than one operating point's.
%
% Note the shape: smooth and monotone. The kinks visible in the panels above are OPERAND-SPECIFIC
% Sobol prefix structure -- where one particular threshold happens to fall among the first few
% low-discrepancy points -- and they average away completely across the operand space.

apCsv = "umul_convergence_allpairs.csv";
apCandidates = [ fullfile(repoRoot, 'keepdata', 'Multiplier Data', apCsv), ...
                 fullfile(repoRoot, apCsv) ];
apFile = "";
for k = 1:numel(apCandidates)
    if isfile(apCandidates(k)), apFile = apCandidates(k); break; end
end
if apFile == ""
    error("Could not find %s. Run --gtest_filter=*ExhaustiveOverEveryOperandPair* first.", apCsv);
end

A = readtable(apFile, "Delimiter", ",", "TextType", "string");
aN    = double(A.N);
aMean = double(A.Mean_AbsError);
aMin  = double(A.Min_AbsError);
aMax  = double(A.Max_AbsError);
aP05  = double(A.P05_AbsError);
aP95  = double(A.P95_AbsError);
aMed  = double(A.Median_AbsError);
nPairs = double(A.Pairs(1));

fprintf("[INFO] ALL-PAIRS  %s operand pairs | mean err N=%d -> %.5f, N=%d -> %.5f\n", ...
        addThousandsSeparator(nPairs), aN(1), aMean(1), aN(end), aMean(end));

ax = nexttile(tl, nLevels + 1, [1 nLevels]);

% P05-P95 band: the range an arbitrary operand pair actually falls in. The true MINIMUM is
% exactly zero at every length and is deliberately not drawn -- the k1 = 0 and k0 = 0 rails emit
% nothing and match an ideal product of zero perfectly, so a "best bound" line would be a
% degenerate flat zero that a log axis cannot show anyway.
% EVERY termination point is present (N = 1..1024), so this is a dense curve. Markers and tick
% labels go on the powers of two only; 1024 of either would be an unreadable solid band.
fill(ax, [aN; flipud(aN)], [aP05; flipud(aP95)], [0.30 0.30 0.75], ...
     "FaceAlpha", 0.18, "EdgeColor", "none");
hold(ax, "on");
loglog(ax, aN, aMax,  "-", "LineWidth", 1.5, "Color", [0.75 0.15 0.15]);
loglog(ax, aN, aMean, "-", "LineWidth", 2.5, "Color", [0.15 0.30 0.75]);
loglog(ax, aN, aMed,  "-.", "LineWidth", 1.3, "Color", [0.15 0.30 0.75]);

pow2 = 2.^(0:floor(log2(max(aN))));
isP2 = ismember(aN, pow2);
loglog(ax, aN(isP2), aMean(isP2), "o", "MarkerSize", 7, ...
       "MarkerFaceColor", [0.15 0.30 0.75], "MarkerEdgeColor", "none");

% 1/sqrt(N) reference, anchored at N = 32 so the comparison is about SHAPE not offset.
anchor = find(aN == 32, 1);
if isempty(anchor), anchor = 1; end
refC = aMean(anchor) * sqrt(aN(anchor));
loglog(ax, aN, refC ./ sqrt(aN), "--", "LineWidth", 1.5, "Color", [0.45 0.45 0.45]);
hold(ax, "off");

set(ax, "XScale", "log", "YScale", "log");
grid(ax, "on");
xticks(ax, pow2); xticklabels(ax, string(pow2));
xlim(ax, [1 max(aN)]);
xlabel(ax, "Bit Stream Length N [log]");
ylabel(ax, "Mean Absolute Error [log]");
title(ax, sprintf(['WHOLE UNIT -- every operand against every operand: %s pairs, ' ...
                   'each exact over all 10^{%.0f} arrangements'], ...
                  addThousandsSeparator(nPairs), D(1).log10Arr));
legend(ax, ["P05-P95 across operand pairs", "worst pair", "MEAN across all pairs", ...
            "median pair", "powers of two", "1/\surdN reference"], ...
       "Location", "southwest", "FontSize", 9);

ax.Toolbar.Visible = "off";

% Explicit dark Color: the default text colour is light in this theme, which on a pale background
% box renders as unreadable grey-on-grey.
text(ax, 0.985, 0.94, ...
     sprintf(['no sampling anywhere on this curve\n' ...
              'true minimum is exactly 0 at every N\n' ...
              '(k_0 = 0 and k_1 = 0 rails are error-free)']), ...
     "Units", "normalized", "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
     "FontSize", 9, "Color", [0.10 0.10 0.10], ...
     "BackgroundColor", [1 1 1], "EdgeColor", [0.55 0.55 0.55], "Margin", 5);

outPng = "umul_behaviour_and_convergence_all_levels.png";
exportgraphics(fig, outPng, "Resolution", 150);
fprintf("[SUCCESS] Combined 6-panel chart saved as %s\n", outPng);

%% ------------------------------------------------------------------
%% Helpers
%% ------------------------------------------------------------------
function s = extractOperands(operation)
    tok = regexp(operation, '([\d.]+\s*x\s*[\d.]+\s*=\s*[\d.]+)', 'tokens', 'once');
    if isempty(tok), s = operation; else, s = string(tok{1}); end
end

% 1049600 -> "1,049,600"
function s = addThousandsSeparator(v)
    s = sprintf("%d", round(v));
    s = regexprep(s, '(\d)(?=(\d{3})+$)', '$1,');
end
