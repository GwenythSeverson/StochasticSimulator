%% Stochastic Computing Adder (MUX) Exhaustive Fault-Injection Analysis
% This script parses the exhaustive 32-bit adder fault injection CSV, generates
% a hardware-constraint focused graph suite, and saves results to .mat.

clear; clc; close all;

csvFile = '32bit_exhaustive_adder_trials.csv';
if ...
   ~isfile(csvFile)
    error('CRITICAL: %s not found. Please run the C++ simulator first.', csvFile);
end

%% 1. Ingest Data via Optimized Table Import
fprintf('Ingesting exhaustive adder CSV file... ');
tic;

% 1. Detect options while forcing the comma delimiter upfront
opts = detectImportOptions(csvFile, 'Delimiter', ',', 'VariableNamingRule', 'preserve');

% 2. Only load the columns we actually need for plotting
neededVars = {'CountA', 'CountB', 'ExpectedCleanOnes', 'BitsFlipped', ...
              'TotalError', 'PrecisionError', 'BitFlipError'};
opts.SelectedVariableNames = neededVars;

% 3. Explicitly enforce numeric types for these columns
opts = setvartype(opts, neededVars, 'double');

% 4. Treat skipped columns as string so the bracketed lists "[1, 0...]" don't break the parser
allVars = opts.VariableNames;
skippedVars = setdiff(allVars, neededVars);
opts = setvartype(opts, skippedVars, 'string'); 

% 5. Safely ingest all 3,276,800 rows
data = readtable(csvFile, opts);

% Rename 'ExpectedCleanOnes' to 'CleanOnes' to seamlessly match your plotting code
data.Properties.VariableNames{'ExpectedCleanOnes'} = 'CleanOnes';

fprintf('Done! Loaded %d rows in %.2f seconds.\n', height(data), toc);
%% 2. Save Workspace to .mat File
matFile = 'Stochastic_Adder_Exhaustive_Data.mat';
fprintf('Saving workspace variables to %s... ', matFile);
save(matFile, 'data');
fprintf('Done!\n');
%% 3. Setup Figure Layout (3x3 Grid)
figure('Name', 'Stochastic Adder Exhaustive Fault Analysis Suite', ...
       'Units', 'Normalized', 'Position', [0.05, 0.05, 0.9, 0.9]);
t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title(t, '32-Bit Exhaustive Stochastic Adder (MUX) Hardware Constraint Suite', ...
    'FontSize', 14, 'FontWeight', 'Bold');

%% Graph 1: Error Scaling & Degradation
nexttile;
[G1, BitsFlipped] = findgroups(data.BitsFlipped);
meanAbsTotalError = splitapply(@mean, abs(data.TotalError), G1);
stdAbsTotalError  = splitapply(@std, abs(data.TotalError), G1);

errorbar(BitsFlipped, meanAbsTotalError, stdAbsTotalError, 'b-o', ...
    'LineWidth', 1.2, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
grid on; xlim([0 33]);
xlabel('Number of Bits Flipped'); ylabel('Mean Absolute Total Error');
title('1. Total Error Scaling & Degradation');

%% Graph 2: Fault Sign Bias
nexttile;
[G2, cA] = findgroups(data.CountA);
meanRawTotalError = splitapply(@mean, data.TotalError, G2);

plot(cA, meanRawTotalError, 'r-s', 'LineWidth', 1.2, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
hold on; plot([0 31], [0 0], 'k--', 'LineWidth', 1);
grid on; xlim([0 31]);
xlabel('Count A (Faulted Stream)'); ylabel('Mean Raw Total Error');
title('2. Fault Sign Bias');

%% Graph 3: Error Variance and Bounds
nexttile; 
selected_flips = [2, 4, 8, 16, 24, 32];
filtered_data = data(ismember(data.BitsFlipped, selected_flips), :);

boxchart(categorical(filtered_data.BitsFlipped), filtered_data.TotalError, ...
    'BoxFaceColor', [0.2 0.2 0.2], 'MarkerStyle', 'none');
grid on;
xlabel('Number of Bits Flipped'); ylabel('Total Error Distribution Range');
title('3. Total Error Variance & Bounds');

%% Graph 4: Precision Error vs. Bit Flip Error (Line Plot)
nexttile;
meanAbsPrecision = splitapply(@mean, abs(data.PrecisionError), G1);
meanAbsBitFlip   = splitapply(@mean, abs(data.BitFlipError), G1);

plot(BitsFlipped, meanAbsPrecision, 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'MarkerSize', 5);
hold on;
plot(BitsFlipped, meanAbsBitFlip, 'm-v', 'LineWidth', 1.5, 'MarkerFaceColor', 'm', 'MarkerSize', 5);
grid on; xlim([0 33]);
xlabel('Number of Bits Flipped'); ylabel('Mean Absolute Error');
legend('Precision (Quantization) Error', 'Bit Flip (Fault) Error', 'Location', 'NorthWest');
title('4. Precision vs. Bit Flip Error Scaling');

%% Graph 5: Resulting Bitstream Fault Sign Bias
nexttile;
[G6, cleanOnes] = findgroups(data.CleanOnes);
meanRawErrorOut = splitapply(@mean, data.TotalError, G6);

plot(cleanOnes, meanRawErrorOut, 'm-o', 'LineWidth', 1.2, 'MarkerFaceColor', 'm', 'MarkerSize', 4);
hold on; plot([0 31], [0 0], 'k--', 'LineWidth', 1);
grid on; xlim([0 31]);
xlabel('Expected Clean Output Ones Count'); ylabel('Mean Raw Total Error');
title('5. Resulting Bitstream Fault Sign Bias');

%% Graph 6: Input Stream A Trial Distribution
nexttile;
[cA_Counts, cA_Values] = groupcounts(data.CountA);
bar(cA_Values, cA_Counts, 'FaceColor', [0.8 0.4 0.4]);
grid on; xlim([-1 32]);
xlabel('Count A (Input Ones Count)'); ylabel('Number of Trials');
title('6. Input Stream A Trial Distribution');

%% Graph 7: REPLACED - Hardware Precision Error Trend Suite (Fail-Safe Mapping)
nexttile([1, 3]);

% Define the specific target values of CountB we want to inspect
target_B = [0, 8, 16, 24, 31];
colors = ['r', 'g', 'b', 'm', 'k'];
labels = {'Count B = 0', 'Count B = 8', 'Count B = 16', 'Count B = 24', 'Count B = 31'};

hold on;
for i = 1:length(target_B)
    % 1. Extract ONLY the rows matching the current Count B target
    subTable = data(data.CountB == target_B(i), :);
    
    % 2. Group this clean slice by CountA to compress the 100 organizations
    [G_sub, unique_cA] = findgroups(subTable.CountA);
    mean_abs_prec = splitapply(@mean, abs(subTable.PrecisionError), G_sub);
    
    % 3. Plot the true, unscrambled trend line
    plot(unique_cA, mean_abs_prec, 'Color', colors(i), 'LineWidth', 2, ...
         'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', colors(i));
end
hold off;

% Format the trend graph
grid on; xlim([0 31]); ylim([0 1.2]);
xlabel('Count A (Input Ones Count)');
ylabel('Mean Absolute Quantization Error');
legend(labels, 'Location', 'NorthWest');
title('7. Hardware Quantization Error Trend Profiles (Parabolic Variance Limits)');
%% 4. Export Figure to PNG
fig = gcf; 
fprintf('Exporting layout... ');
exportgraphics(fig, 'Stochastic_Adder_Exhaustive_Fault_Analysis.png', 'Resolution', 300, 'BackgroundColor', 'current');
fprintf('Done! Saved as PNG.\n');