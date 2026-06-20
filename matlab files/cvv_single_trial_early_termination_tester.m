%% cvv_multi_trial_early_termination_tester.m
%
% Please note that this is AI-Generated code and may need adjustments to work
% correctly in your environment. I kept a given structure but chnaged quite a bit so 
% there might be some weird artifacts leftover. - Gwenyth 6/19/26
%
% Loads a CSV containing Constant Value Vector (CVV) bit streams and analyzes
% multiple trial rows to evaluate early termination behavior across subplots.
% Also aggregates data across all trials to plot the mean absolute error.

clear; clc; close all;

%% ------------------------------------------------------------------
%% 1. CONFIGURATION
%% ------------------------------------------------------------------
csvFile         = "multiplier_accuracy_trials.csv";  
trialsToAnalyze = 1:10;                    % Supports any array range (e.g., 1:1000)
truncLens       = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];  % Kept sorted ascending, you can chnage. 

% Ideal target probability used when the stream was generated
ideal_target_prob = 0.2500;

%% ------------------------------------------------------------------
%% 2. LOAD CSV AND PREPARE DATA
%% ------------------------------------------------------------------
opts = detectImportOptions(csvFile, "Delimiter", ",", "TextType", "string");
opts = setvartype(opts, "string");   
T = readtable(csvFile, opts);

varNames = string(T.Properties.VariableNames);
cleanNames = lower(erase(varNames, ["_", " "]));

trialCol   = varNames(cleanNames == "trial");
targetCol  = varNames(cleanNames == "vectoroutput"); 

if isempty(trialCol) || isempty(targetCol)
    error("Could not find expected columns (Trial, Vector_Output) in %s.", csvFile);
end

trialValues = double(T.(trialCol));

%% ------------------------------------------------------------------
%% 3. DATA PRE-ALLOCATION FOR SUMMARY
%% ------------------------------------------------------------------
numTrials = numel(trialsToAnalyze);
nLens     = numel(truncLens);

% Matrices to save calculation states across all trials
allEstProbs = zeros(numTrials, nLens);
allAbsErrors = zeros(numTrials, nLens);

%% ------------------------------------------------------------------
%% 4. LOOP THROUGH EACH REQUESTED TRIAL & CALC STATS
%% ------------------------------------------------------------------
for tIdx = 1:numTrials
    currentTrial = trialsToAnalyze(tIdx);
    
    % Find trial row
    rowIdx = find(trialValues == currentTrial, 1);
    if isempty(rowIdx)
        warning("Trial %d not found in %s. Skipping.", currentTrial, csvFile);
        continue;
    end
    
    % Extract and parse target stream
    target_str = T.(targetCol)(rowIdx);
    streamOut  = parseBitVector(target_str);
    
    for i = 1:nLens
        N = truncLens(i);
        if N > numel(streamOut)
            error("Truncation length %d exceeds stream length %d in Trial %d", N, numel(streamOut), currentTrial);
        end
        
        calculatedProb = sum(streamOut(1:N)) / N;
        allEstProbs(tIdx, i) = calculatedProb;
        allAbsErrors(tIdx, i) = abs(calculatedProb - ideal_target_prob);
    end
end

%% ------------------------------------------------------------------
%% 5. COMPUTE COMPREHENSIVE AVERAGES
%% ------------------------------------------------------------------
% Find average absolute error at every truncation step across all trials
meanAbsError = mean(allAbsErrors, 1, "omitnan");

%% ------------------------------------------------------------------
%% 6. PLOT 1: COMPREHENSIVE OVERALL SUMMARY GRAPH
%% ------------------------------------------------------------------
figure("Name", "Overall Average Error Summary", "Position", [150, 150, 800, 500]);

% Using loglog or semilogx highlights the exponential convergence properties
semilogx(truncLens, meanAbsError, "-^", "LineWidth", 2.2, "MarkerSize", 8, ...
         "Color", [0.85 0.33 0.10], "MarkerFaceColor", [0.85 0.33 0.10]);
grid on;
set(gca, "XDir", "normal");
xticks(truncLens);
xticklabels(string(truncLens));

xlabel("Bit Stream Length (N bits) [Log Scale]");
ylabel("Mean Absolute Error");
title(sprintf("Overall Convergence Layout: Mean Abs. Error across %d Trials", numTrials));

% Save summary chart
saveas(gcf, "overall_mean_error_summary.png");
fprintf("[SUCCESS] Global overview chart saved as overall_mean_error_summary.png\n");

%% ------------------------------------------------------------------
%% 7. PLOT 2: INDIVIDUAL TRIAL SUBPLOTS
%% ------------------------------------------------------------------
% Keeps individual subplots visual, clean, and uncrowded from the summary engine
nCols = ceil(sqrt(numTrials));
nRows = ceil(numTrials / nCols);

figure("Name", "Individual CVV Trial Behavior Panels", ...
       "Position", [200, 100, 350 * nCols, 300 * nRows]);

for tIdx = 1:numTrials
    currentTrial = trialsToAnalyze(tIdx);
    
    subplot(nRows, nCols, tIdx);
    semilogx(truncLens, allEstProbs(tIdx, :), "-o", "LineWidth", 1.5, "MarkerSize", 5, ...
             "MarkerFaceColor", [0.20 0.45 0.85]);
    hold on;
    yline(ideal_target_prob, "--r", "LineWidth", 1.2);
    hold off;
    
    set(gca, "XDir", "normal");    
    xticks(truncLens);
    xticklabels(string(truncLens));
    grid on;
    ylim([0 1]); 
    
    xlabel("Bit Stream Length");
    ylabel("Est. Prob.");
    title(sprintf("Trial %d", currentTrial));
end

sgtitle(sprintf("Individual CVV Stream Realizations (Target P = %.4f)", ideal_target_prob));

% Save final multi-panel figure
saveas(gcf, "multi_trial_cvv_individual_panels.png");
fprintf("[SUCCESS] Individual panels saved as multi_trial_cvv_individual_panels.png\n");

%% ------------------------------------------------------------------
%% Helper function: parse a CSV-stored bit vector string
%% ------------------------------------------------------------------
function bits = parseBitVector(strVal)
    s = char(strVal);
    s = erase(s, ["[", "]"]);     
    parts = split(string(s), ",");
    parts = strtrim(parts);
    bits = double(parts)';        
end