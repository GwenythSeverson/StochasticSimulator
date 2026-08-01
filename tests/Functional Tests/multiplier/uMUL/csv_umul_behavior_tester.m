%% csv_umul_behavior_tester.m
%
% uMUL WHOLE-UNIT BEHAVIOUR. Four panels built from the two exhaustive sweeps, which cover
% perpendicular axes and together describe the unit with nothing sampled anywhere:
%
%   umul_convergence_allpairs.csv   (test_uMUL_convergence.cpp :: ExhaustiveOverEveryOperandPair)
%       every operand pair x every stream ARRANGEMENT x every truncation length.
%       1,049,600 pairs, each averaged EXACTLY over all ~10^307 arrangements of in_0 by
%       hypergeometric weighting rather than by sampling them.
%
%   umul_behavior_exhaustive.csv    (test_uMUL_behavior.cpp)
%       every operand pair at FULL length, one row per pair, error in BITS of the 1024-cycle
%       output. 1,049,600 rows, ~29 MB. This is genuine brute force: 1.07 billion multiply()
%       calls, one real simulation per pair.
%
% WHY BOTH. The convergence sweep answers "how fast does the unit converge, on average, over
% everything" -- that is panel 1, the early-termination curve. It cannot answer "WHICH operand
% pairs are the bad ones", because it aggregates them away. The behaviour sweep keeps every pair
% separate, so panels 2-4 can show where in the operand plane the error actually lives.
%
% PANELS
%   1  Early-termination average error across ALL pairs and ALL arrangements. Log-log, with the
%      P05-P95 band across operand pairs and the worst pair. THE headline chart.
%   2  Heatmap of signed error in bits over the whole k0 x k1 operand plane at full length.
%   3  Distribution of |error| at full length, over all 1,049,600 pairs.
%   4  |error| against the ideal product, to show whether accuracy depends on output magnitude.
%
% NOTE ON THE TWO OPERAND RANGES. k0 runs 0..1024 (1025 values) because it is a TALLY of ones in
% a 1024-bit stream; k1 runs 0..1023 (1024 values) because it is the contents of a 10-BIT
% REGISTER. Storing 1024 would need an eleventh bit, which uMUL_uni.sv's `iB` does not have, so
% p1 = 1.0 is not an expressible operand. The grid is asymmetric by construction, not by slip.

clear; clc; close all;

%% ------------------------------------------------------------------
%% 1. LOCATE THE CSVs
%% ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = scriptDir;
while ~isfolder(fullfile(repoRoot, 'keepdata'))
    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
        error("Could not locate the repository root (no 'keepdata' folder found above %s).", scriptDir);
    end
    repoRoot = parentDir;
end

function p = findCsv(repoRoot, name)
    candidates = [ fullfile(repoRoot, 'keepdata', 'Multiplier Data', name), ...
                   fullfile(repoRoot, name) ];
    p = "";
    for k = 1:numel(candidates)
        if isfile(candidates(k)), p = candidates(k); return; end
    end
    error(['Could not find %s.\n' ...
           'Run:  stochastic_computer.exe --gtest_filter=*UMULConvergence*:*Behaviour*'], name);
end

convFile = findCsv(repoRoot, "umul_convergence_allpairs.csv");
behFile  = findCsv(repoRoot, "umul_behavior_exhaustive.csv");

%% ------------------------------------------------------------------
%% 2. LOAD
%% ------------------------------------------------------------------
fprintf("[INFO] reading %s\n", convFile);
A = readtable(convFile, "Delimiter", ",", "TextType", "string");
aN     = double(A.N);
aMean  = double(A.Mean_AbsError);
aMax   = double(A.Max_AbsError);
aP05   = double(A.P05_AbsError);
aP95   = double(A.P95_AbsError);
aMed   = double(A.Median_AbsError);
nPairs = double(A.Pairs(1));

% readmatrix rather than readtable: 1.05M rows of pure numerics load several times faster and
% we do not need the column names once the order is known.
%   columns: K0, K1, Out_Ones, Ideal_Ones, Err_Bits
fprintf("[INFO] reading %s (this is the ~29 MB one)\n", behFile);
tRead = tic;
M = readmatrix(behFile);
fprintf("[INFO] %d rows loaded in %.1f s\n", size(M,1), toc(tRead));

K0   = M(:,1);
K1   = M(:,2);
errB = M(:,5);
ideal = (K0 ./ 1024) .* (K1 ./ 1024);

nK0 = max(K0) + 1;   % 1025
nK1 = max(K1) + 1;   % 1024
assert(size(M,1) == nK0 * nK1, "row count does not match the k0 x k1 grid");

% The C++ sweep writes k1 in the OUTER loop and k0 in the inner, so k0 varies fastest and a
% column-major reshape lands each column on one k1. E(i,j) is then k0 = i-1, k1 = j-1.
E = reshape(errB, nK0, nK1);

fprintf("[INFO] full-length: RMS %.4f bits | max |err| %.4f bits | within 1 bit %.2f %%\n", ...
        sqrt(mean(errB.^2)), max(abs(errB)), 100*mean(abs(errB) <= 1));

%% ------------------------------------------------------------------
%% 3. FIGURE
%% ------------------------------------------------------------------
fig = figure("Name", "uMUL Whole-Unit Behaviour", "Position", [50, 50, 1700, 980], "Color", "w");
tl = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "compact");
title(tl, "uMUL whole-unit behaviour -- every operand pair, every stream arrangement", ...
      "FontWeight", "bold", "FontSize", 15);
subtitle(tl, sprintf(['%s operand pairs   |   each averaged exactly over all 10^{307} in_0 ' ...
                      'arrangements   |   nothing sampled'], addThousandsSeparator(nPairs)), ...
         "FontSize", 10);

%% ---- PANEL 1: early-termination average error, whole unit ----
% EVERY termination point is present (N = 1..1024), so this is a dense curve rather than ten
% sampled points. Markers are drawn only on the powers of two -- 1024 markers would be a solid
% band -- and the tick labels likewise, since 1024 of those would be unreadable.
ax1 = nexttile(tl, 1);
fill(ax1, [aN; flipud(aN)], [aP05; flipud(aP95)], [0.30 0.30 0.75], ...
     "FaceAlpha", 0.18, "EdgeColor", "none");
hold(ax1, "on");
loglog(ax1, aN, aMax,  "-", "LineWidth", 1.4, "Color", [0.75 0.15 0.15]);
loglog(ax1, aN, aMean, "-", "LineWidth", 2.4, "Color", [0.15 0.30 0.75]);
loglog(ax1, aN, aMed,  "-.", "LineWidth", 1.2, "Color", [0.15 0.30 0.75]);

pow2 = 2.^(0:floor(log2(max(aN))));
isP2 = ismember(aN, pow2);
loglog(ax1, aN(isP2), aMean(isP2), "o", "MarkerSize", 6, ...
       "MarkerFaceColor", [0.15 0.30 0.75], "MarkerEdgeColor", "none");

% 1/sqrt(N) reference anchored at N = 32, so the comparison is about SHAPE not offset. On log-log
% a random estimator is a straight line of slope -1/2; bending below it is the Sobol prefix.
anchor = find(aN == 32, 1); if isempty(anchor), anchor = 1; end
refC = aMean(anchor) * sqrt(aN(anchor));
loglog(ax1, aN, refC ./ sqrt(aN), "--", "LineWidth", 1.4, "Color", [0.45 0.45 0.45]);
hold(ax1, "off");

set(ax1, "XScale", "log", "YScale", "log"); grid(ax1, "on");
xticks(ax1, pow2); xticklabels(ax1, string(pow2));
xlim(ax1, [1 max(aN)]);
xlabel(ax1, "Bit Stream Length N  (early termination point) [log]");
ylabel(ax1, "Mean Absolute Error [log]");
title(ax1, sprintf("1.  Early-termination average error -- ALL pairs, ALL arrangements, all %d termination points", ...
                   numel(aN)));
legend(ax1, ["P05-P95 across operand pairs", "worst pair", "MEAN across all pairs", ...
             "median pair", "powers of two", "1/\surdN reference"], ...
       "Location", "southwest", "FontSize", 8);
ax1.Toolbar.Visible = "off";

%% ---- PANEL 2: the operand plane at full length ----
ax2 = nexttile(tl, 2);
imagesc(ax2, [0 nK1-1], [0 nK0-1], E);
set(ax2, "YDir", "normal");
lim = max(abs(errB));
clim(ax2, [-lim lim]);
colormap(ax2, divergingMap(256));
cb = colorbar(ax2);
cb.Label.String = "signed error (bits of 1024)";
xlabel(ax2, "k_1  --  loaded binary operand (0..1023)");
ylabel(ax2, "k_0  --  ones in the in_0 stream (0..1024)");
title(ax2, sprintf("2.  Full-length error over the whole operand plane  (max |err| = %.2f bits)", lim));
ax2.Toolbar.Visible = "off";

%% ---- PANEL 3: distribution of |error| at full length ----
ax3 = nexttile(tl, 3);
absErr = abs(errB);
histogram(ax3, absErr, 0:0.1:ceil(max(absErr)*10)/10, ...
          "FaceColor", [0.15 0.30 0.75], "EdgeColor", "none");
set(ax3, "YScale", "log"); grid(ax3, "on");
xline(ax3, 1.0, "--r", "LineWidth", 1.8, "Label", "1 bit", "LabelOrientation", "horizontal");
xlabel(ax3, "|error| at full length (bits)");
ylabel(ax3, "operand pairs [log]");
title(ax3, sprintf("3.  Error distribution at N = 1024  (%.2f %% within one bit)", ...
                   100*mean(absErr <= 1)));
ax3.Toolbar.Visible = "off";

%% ---- PANEL 4: does accuracy depend on the output magnitude? ----
% Binned over the ideal product. If uMUL were biased toward one end of its range this is where it
% would show; a flat mean says the unit is even-handed across its whole output span.
ax4 = nexttile(tl, 4);
nBins = 64;
edges = linspace(0, 1, nBins+1);
[~, ~, bin] = histcounts(ideal, edges);
keep = bin > 0;
centres = (edges(1:end-1) + edges(2:end)) / 2;
meanByBin = accumarray(bin(keep), absErr(keep), [nBins 1], @mean, NaN);
maxByBin  = accumarray(bin(keep), absErr(keep), [nBins 1], @max,  NaN);

plot(ax4, centres, maxByBin, "-", "LineWidth", 1.5, "Color", [0.75 0.15 0.15]);
hold(ax4, "on");
plot(ax4, centres, meanByBin, "-", "LineWidth", 2.4, "Color", [0.15 0.30 0.75]);
hold(ax4, "off");
grid(ax4, "on");
xlabel(ax4, "ideal product  p_0 \times p_1");
ylabel(ax4, "|error| at full length (bits)");
title(ax4, "4.  Error vs output magnitude  (binned over the ideal product)");
legend(ax4, ["worst pair in bin", "mean over bin"], "Location", "north", "FontSize", 8);
ax4.Toolbar.Visible = "off";

outPng = "umul_behaviour_whole_unit.png";
exportgraphics(fig, outPng, "Resolution", 150);
fprintf("[SUCCESS] saved %s\n", outPng);

%% ------------------------------------------------------------------
%% Helpers
%% ------------------------------------------------------------------

% Blue-white-red, centred on zero. MATLAB ships no diverging map, and a sequential one (parula)
% would hide the sign of the error, which is the most informative thing in panel 2.
function cmap = divergingMap(n)
    half = floor(n/2);
    lower = [linspace(0.13,1,half)', linspace(0.30,1,half)', linspace(0.65,1,half)'];
    upper = [linspace(1,0.70,n-half)', linspace(1,0.10,n-half)', linspace(1,0.10,n-half)'];
    cmap = [lower; upper];
end

% 1049600 -> "1,049,600"
function s = addThousandsSeparator(v)
    s = sprintf("%d", round(v));
    s = regexprep(s, '(\d)(?=(\d{3})+$)', '$1,');
end
