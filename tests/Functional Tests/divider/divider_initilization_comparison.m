%% divider_initialization_comparison.m
%
% Compares the first 64 bits of the divider warm-up phase between two setups:
% 1) Counter initialized at 0  (From the CNTat0 folder)
% 2) Counter initialized at 16 (From the CNTat16 folder)

clear; clc; close all;

%% ------------------------------------------------------------------
%% 1. CONFIGURATION (FILE PATHS MAPPED TO YOUR DIRECTORY)
%% ------------------------------------------------------------------
% 1. Start inside: 'StochasticSimulator/tests/Functional Tests/'
scriptDir = fileparts(mfilename('fullpath'));

% 2. Step up out of 'Functional Tests' into 'tests/'
testsDir = fullfile(scriptDir, '..');

% 3. Step up out of 'tests' into the repository root ('StochasticSimulator/')
rawRepoRoot = fullfile(testsDir, '..');

% 4. FORCE Windows/OneDrive to resolve '..' links into a true absolute path
[status, info] = fileattrib(rawRepoRoot);
if ~status
    error("Could not resolve repository root directory structure.");
end
repoRoot = info.Name; % This is now a clean path with NO '..' or relative links

% 5. Step down into the target data directories
dataRootDir = fullfile(repoRoot, 'keepdata','Divider data', 'CNT_starting_point_data');

file_CNT0  = fullfile(dataRootDir, '.25_.5_CNTat0_graphs_and_files', 'divider_accuracy_trials1.csv');
file_CNT16 = fullfile(dataRootDir, '.25_.5_CNTat16_graphs_and_files', 'divider_accuracy_trials.csv');


warmupWindow = 256; 
ideal_target_prob = 0.5000;

%% ------------------------------------------------------------------
%% 2. PROCESS FILE 1: COUNTER INITIALIZED AT 0
%% ------------------------------------------------------------------
fprintf("[INFO] Processing Counter = 0 data...\n");
meanProb_CNT0 = getEarlyRunningAverage(file_CNT0, warmupWindow);

%% ------------------------------------------------------------------
%% 3. PROCESS FILE 2: COUNTER INITIALIZED AT 16
%% ------------------------------------------------------------------
fprintf("[INFO] Processing Counter = 16 data...\n");
meanProb_CNT16 = getEarlyRunningAverage(file_CNT16, warmupWindow);

%% ------------------------------------------------------------------
%% 4. GRAPHICAL COMPARISON PLOT
%% ------------------------------------------------------------------
figure("Name", "Counter Initialization Impact: First 256 Bits", "Position", [200, 200, 850, 500]);

% Plot CNT = 0 Trace (Blue)
plot(1:warmupWindow, meanProb_CNT0, ".-", "LineWidth", 2.0, "MarkerSize", 10, ...
     "Color", [0.00 0.45 0.74], "DisplayName", "Counter Init = 0");
hold on;

% Plot CNT = 16 Trace (Orange)
plot(1:warmupWindow, meanProb_CNT16, ".-", "LineWidth", 2.0, "MarkerSize", 10, ...
     "Color", [0.85 0.33 0.10], "DisplayName", "Counter Init = 16 (Midpoint)");

% Target line (Red dashed)
yline(ideal_target_prob, "--r", "Target Probability (0.5000)", ...
      "LineWidth", 1.5, "LabelHorizontalAlignment", "left", "FontSize", 10, ...
      "HandleVisibility", "off");

grid on;
box on;

xlim([1 warmupWindow]);
ylim([0 1.0]);
xticks([1, 10:10:60, warmupWindow]);

xlabel("Bit Position (Clock Cycle)");
ylabel("Mean Running Probability");
title("Warm-Up Transient: Counter Init 0 vs. 16 Comparison");
legend("Location", "southeast", "FontSize", 10);

% Save the comparison graph
saveas(gcf, "counter_init_comparison_256bits.png");
fprintf("[SUCCESS] Comparison plot saved as counter_init_comparison_256bits.png\n");


%% ------------------------------------------------------------------
%% LOCAL HELPER FUNCTION TO PROCESS A SINGLE FILE
%% ------------------------------------------------------------------
function meanRunningProb = getEarlyRunningAverage(filePath, windowSize)
    opts = detectImportOptions(filePath, "Delimiter", ",", "TextType", "string");
    opts = setvartype(opts, "string");   
    T = readtable(filePath, opts);

    varNames = string(T.Properties.VariableNames);
    cleanNames = lower(erase(varNames, ["_", " "]));

    trialCol   = varNames(cleanNames == "trial");
    targetCol  = varNames(cleanNames == "vectoroutput"); 

    if isempty(trialCol) || isempty(targetCol)
        error("Could not find required columns in %s", filePath);
    end

    trialValues = double(T.(trialCol));
    numTrials   = numel(trialValues);

    earlyProbs = zeros(numTrials, windowSize);

    for tIdx = 1:numTrials
        rowIdx = find(trialValues == tIdx, 1);
        if isempty(rowIdx)
            continue;
        end
        
        fullStream = parseBitVector(T.(targetCol)(rowIdx));
        earlyStream = fullStream(1:windowSize);
        
        % Compute cumulative running probability
        earlyProbs(tIdx, :) = cumsum(earlyStream) ./ (1:windowSize);
    end
    
    meanRunningProb = mean(earlyProbs, 1, "omitnan");
end

%% Helper function: parse a CSV-stored bit vector string
function bits = parseBitVector(strVal)
    s = char(strVal);
    s = erase(s, ["[", "]"]);     
    parts = split(string(s), ",");
    parts = strtrim(parts);
    bits = double(parts)';        
end