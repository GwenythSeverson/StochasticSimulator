%% cvv_divider_trial_early_termination_tester.m
%
% Loads a CSV containing divider Constant Value Vector (CVV) bit streams,
% analyzes early termination behavior across all 1000 trials, and aggregates
% data to plot the mean absolute error. For individual visualization, it 
% randomly selects 50 trials to plot as subplots.

clear; clc; close all;

%% ------------------------------------------------------------------
%% 1. CONFIGURATION
%% ------------------------------------------------------------------
csvFile         = "divider_accuracy_trials.csv";  
truncLens       = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];  
numSubplotPlots = 50;  % Number of random individual subplots to generate

% Ideal target probability for the divider: Pa / Pb = 0.2500 / 0.5000
ideal_target_prob = 0.5000;

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
totalTrialsInFile = numel(trialValues);
trialsToAnalyze = 1:totalTrialsInFile;

%% ------------------------------------------------------------------
%% 3. DATA PRE-ALLOCATION FOR SUMMARY
%% ------------------------------------------------------------------
numTrials = numel(trialsToAnalyze);
nLens     = numel(truncLens);

% Matrices to save calculation states across all trials
allEstProbs = zeros(numTrials, nLens);
allAbsErrors = zeros(numTrials, nLens);
% Running probability for warm-up analysis
streamLength = truncLens(end);      % 1024 bits
allRunningProb = zeros(numTrials, streamLength);

%% ------------------------------------------------------------------
%% 4. LOOP THROUGH EACH REQUESTED TRIAL & CALC STATS
%% ------------------------------------------------------------------
fprintf("[INFO] Parsing %d trials from C++ simulation output...\n", numTrials);
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
    % Running probability across the entire bitstream
runningProb = cumsum(streamOut) ./ (1:numel(streamOut));
allRunningProb(tIdx,:) = runningProb;

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
%% 6. PLOT 1: DIVIDER CONVERGENCE (LOG vs LINEAR X-AXIS)
%% ------------------------------------------------------------------

figure("Name", "Divider Convergence Comparison", ...
       "Position", [150, 150, 1200, 500]);

tiledlayout(1,2);

%% -----------------------------
%% (A) LOG-SCALED X-AXIS
%% -----------------------------
nexttile;

semilogx(truncLens, meanAbsError, "-^", ...
    "LineWidth", 2.2, ...
    "MarkerSize", 8, ...
    "Color", [0.00 0.45 0.74], ...
    "MarkerFaceColor", [0.00 0.45 0.74]);

grid on;
xticks(truncLens);
xticklabels(string(truncLens));

xlabel("Bit Stream Length (N bits) [Log Scale]");
ylabel("Mean Absolute Error");

title("Log-Scale View");

%% -----------------------------
%% (B) LINEAR X-AXIS (UNSCALED)
%% -----------------------------
nexttile;

plot(truncLens, meanAbsError, "-^", ...
    "LineWidth", 2.2, ...
    "MarkerSize", 8, ...
    "Color", [0.85 0.33 0.10], ...
    "MarkerFaceColor", [0.85 0.33 0.10]);

grid on;

xticks(truncLens);
xticklabels(string(truncLens));

xlabel("Bit Stream Length (N bits) [Linear Scale]");
ylabel("Mean Absolute Error");

title("Linear (Unscaled) View");

%% -----------------------------
%% SAVE FIGURE
%% -----------------------------
saveas(gcf, "divider_mean_error_dual_view.png");

fprintf("[SUCCESS] Dual-view convergence plot saved as divider_mean_error_dual_view.png\n");
%% ------------------------------------------------------------------
%% 7. PLOT 2: INDIVIDUAL TRIAL SUBPLOTS (RANDOM 50 SELECTION)
%% ------------------------------------------------------------------
% Select a random subset of 50 trials from the available dataset
rng("shuffle"); % Ensures different selections across fresh sessions
if numTrials >= numSubplotPlots
    subplotIndices = randperm(numTrials, numSubplotPlots);
else
    subplotIndices = 1:numTrials; % Fallback if file has fewer than 50 rows
    numSubplotPlots = numTrials;
end

nCols = ceil(sqrt(numSubplotPlots));
nRows = ceil(numSubplotPlots / nCols);

figure("Name", "Randomized Individual Divider Trial Behavior Panels", ...
       "Position", [200, 100, 250 * nCols, 200 * nRows]);

for pIdx = 1:numSubplotPlots
    tIdx = subplotIndices(pIdx);
    currentTrial = trialsToAnalyze(tIdx);
    
    subplot(nRows, nCols, pIdx);
    semilogx(truncLens, allEstProbs(tIdx, :), "-o", "LineWidth", 1.2, "MarkerSize", 4, ...
             "Color", [0.4660 0.6740 0.1880], "MarkerFaceColor", [0.4660 0.6740 0.1880]);
    hold on;
    yline(ideal_target_prob, "--r", "LineWidth", 1.0);
    hold off;
    
    set(gca, "XDir", "normal", "FontSize", 8);    
    xticks(truncLens);
    % Only show minimal x labels to clear screen real estate
    if pIdx > (numSubplotPlots - nCols)
        xticklabels(string(truncLens));
    else
        xticklabels([]);
    end
    grid on;
    ylim([0 1]); 
    
    title(sprintf("Tr. %d", currentTrial), "FontSize", 9);
end

sgtitle(sprintf("Randomized Sample of %d Divider Stream Realizations (Target P = %.4f)", ...
    numSubplotPlots, ideal_target_prob), "FontSize", 12, "FontWeight", "bold");

% Save final multi-panel figure
saveas(gcf, "multi_trial_divider_individual_panels.png");
fprintf("[SUCCESS] 50 randomized individual panels saved as multi_trial_divider_individual_panels.png\n");
%% ------------------------------------------------------------------
%% 8. PLOT 3: DIVIDER WARM-UP BEHAVIOR Empirical probability is the observed fraction of 1s in generated bitstream, computed directly from simulation data rather than theoretical assumptions.
%% ------------------------------------------------------------------

meanRunningProb = mean(allRunningProb,1,"omitnan");

figure("Name","Divider Warm-Up Behavior", ...
       "Position",[250 250 850 500]);

plot(1:streamLength, meanRunningProb, ...
    "LineWidth",2.2, ...
    "Color",[0.85 0.33 0.10]);

hold on;
yline(ideal_target_prob,"--b","Target Probability", ...
      "LineWidth",1.8);

grid on;
box on;

xlabel("Bit Position");
ylabel("Average Running Probability");

title(sprintf("Divider Warm-Up Behavior Across %d Trials",numTrials));

xlim([1 streamLength]);
ylim([0 1]);

saveas(gcf,"divider_warmup_behavior.png");

fprintf("[SUCCESS] Warm-up behavior plot saved as divider_warmup_behavior.png\n");
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