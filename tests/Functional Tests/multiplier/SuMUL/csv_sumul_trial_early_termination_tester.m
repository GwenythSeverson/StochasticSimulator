%% csv_sumul_trial_early_termination_tester.m
%
% SuMUL behaviour / sanity check. Reads sumul_exhaustive_trials.csv, written by
% test_SuMUL_accuracy.cpp, which sweeps EVERY stream arrangement for one target probability pair
% (1023 seeds per operand -> 1023^2 = 1,046,529 arrangements) and aggregates them in C++.
%
% The CSV carries one row per truncation length, plus the operation that produced it, so the
% figure titles below state exactly what was computed instead of hard-coding it here.
%
% Two figures:
%   1. Mean estimate vs N, against the true product. This is the BEHAVIOUR plot -- the counter
%      warm-up is visible as the curve climbing onto the reference line.
%   2. Mean absolute error vs N. This is the CONVERGENCE plot, matching the AND-gate script.
%
% Read fig 1 before fig 2. |error| folds the sign away, so it shows a spurious bump wherever the
% mean estimate crosses the target -- fig 1 shows what is actually happening.

clear; clc; close all;

%% ------------------------------------------------------------------
%% 1. CONFIGURATION
%% ------------------------------------------------------------------
csvName = "sumul_exhaustive_trials.csv";

% Walk up from this script until we find the folder containing 'keepdata', rather than counting
% a fixed number of '..' steps. This script lives four levels deep
% (tests/Functional Tests/multiplier/SuMUL/), and searching for a landmark means moving it again
% costs nothing.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = scriptDir;
while ~isfolder(fullfile(repoRoot, 'keepdata'))
    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
        error("Could not locate the repository root (no 'keepdata' folder found above %s).", scriptDir);
    end
    repoRoot = parentDir;
end

% The C++ test writes the CSV to whatever directory it was run from (usually the repo root);
% you then move it into keepdata. Accept either location.
candidates = [ fullfile(repoRoot, 'keepdata', 'Multiplier Data', csvName), ...
               fullfile(repoRoot, csvName) ];
csvFile = "";
for k = 1:numel(candidates)
    if isfile(candidates(k))
        csvFile = candidates(k);
        break;
    end
end
if csvFile == ""
    error("Could not find %s in either:\n  %s\n  %s", csvName, candidates(1), candidates(2));
end
fprintf("[INFO] Reading %s\n", csvFile);

%% ------------------------------------------------------------------
%% 2. LOAD CSV
%% ------------------------------------------------------------------
T = readtable(csvFile, "Delimiter", ",", "TextType", "string");

truncLens    = double(T.N);
meanEst      = double(T.Mean_Est);
meanAbsError = double(T.Mean_AbsError);
minEst       = double(T.Min_Est);
maxEst       = double(T.Max_Est);

ideal_target_prob = double(T.Ideal_Product(1));
numArrangements   = double(T.Total_Arrangements(1));
operation         = string(T.Operation(1));

% Everything the run needs to identify itself, assembled once and reused in both titles.
subtitleStr = sprintf("%s   |   all %s stream arrangements", ...
                      operation, addThousandsSeparator(numArrangements));

fprintf("[INFO] %s\n", operation);
fprintf("[INFO] %s arrangements, N = %d..%d\n", ...
        addThousandsSeparator(numArrangements), min(truncLens), max(truncLens));
fprintf("[INFO] mean est: N=%d -> %.4f,  N=%d -> %.4f  (ideal %.4f)\n", ...
        truncLens(1), meanEst(1), truncLens(end), meanEst(end), ideal_target_prob);

%% ------------------------------------------------------------------
%% 3. PLOT 1: BEHAVIOUR -- where the output actually sits vs stream length
%% ------------------------------------------------------------------
figure("Name", "SuMUL Behaviour", "Position", [150, 150, 900, 540]);

% Min/max envelope across every arrangement, drawn as a shaded band behind the mean.
fill([truncLens; flipud(truncLens)], [minEst; flipud(maxEst)], [0.20 0.45 0.85], ...
     "FaceAlpha", 0.15, "EdgeColor", "none");
hold on;
semilogx(truncLens, meanEst, "-o", "LineWidth", 2.2, "MarkerSize", 8, ...
         "Color", [0.20 0.45 0.85], "MarkerFaceColor", [0.20 0.45 0.85]);
yline(ideal_target_prob, "--r", "LineWidth", 1.8, ...
      "Label", sprintf("true product %.4f", ideal_target_prob));
hold off;

set(gca, "XScale", "log", "XDir", "normal");
grid on;
xticks(truncLens); xticklabels(string(truncLens));
ylim([0 1]);
xlabel("Bit Stream Length (N bits) [Log Scale]");
ylabel("Mean Output Probability");
title({"SuMUL Output vs Stream Length", subtitleStr});
legend(["min-max envelope", "mean estimate", "true product"], "Location", "southeast");

saveas(gcf, "sumul_behaviour_vs_length.png");
fprintf("[SUCCESS] Behaviour chart saved as sumul_behaviour_vs_length.png\n");

%% ------------------------------------------------------------------
%% 4. PLOT 2: CONVERGENCE -- mean absolute error
%% ------------------------------------------------------------------
figure("Name", "SuMUL Overall Average Error Summary", "Position", [200, 120, 900, 540]);

semilogx(truncLens, meanAbsError, "-^", "LineWidth", 2.2, "MarkerSize", 8, ...
         "Color", [0.85 0.33 0.10], "MarkerFaceColor", [0.85 0.33 0.10]);
grid on;
set(gca, "XDir", "normal");
xticks(truncLens); xticklabels(string(truncLens));

xlabel("Bit Stream Length (N bits) [Log Scale]");
ylabel("Mean Absolute Error");
title({"SuMUL Convergence: Mean Abs. Error", subtitleStr});

saveas(gcf, "sumul_overall_mean_error_summary.png");
fprintf("[SUCCESS] Global overview chart saved as sumul_overall_mean_error_summary.png\n");

%% ------------------------------------------------------------------
%% Helper: 1046529 -> "1,046,529"
%% ------------------------------------------------------------------
function s = addThousandsSeparator(v)
    s = sprintf("%d", round(v));
    s = regexprep(s, '(\d)(?=(\d{3})+$)', '$1,');
end
