%% Stochastic Computing Multiplier Fault-Injection Analysis
clear; clc; close all;

csvFile = '32bit_ZCE_multiplier_trials.csv';
if ~isfile(csvFile)
    error('CRITICAL: %s not found. Please run the C++ simulator first.', csvFile);
end

%% 0. High-Performance Text Scanning Ingestion
fprintf('Scanning file structurally... ');
tic;

% Open the file for reading
fid = fopen(csvFile, 'r');

% Skip the very first header line
fgetl(fid);

% CREATE A BULLETPROOF FORMAT SPECIFIER:
% We match the exact 12 columns written by C++. 
% %f = number, %*q = completely skip a quoted string block (our messy vectors!)
% Column layout: Trial(f), cA(f), cB(f), CleanOnes(f), Org(f), Flips(f), VecA(*q), VecB(*q), VecA_F(*q), VecOut(*q), FaultedOnes(f), Error(f)
formatSpec = '%f %f %f %f %f %f %*q %*q %*q %*q %f %f';

% Parse the raw text file directly into numeric vectors using commas as delimiters
rawCells = textscan(fid, formatSpec, 'Delimiter', ',', 'TreatAsEmpty', {'NA','NaN'});
fclose(fid);

% Put the data into a clean, tight table
data = table(rawCells{2}, rawCells{3}, rawCells{6}, rawCells{8}, ...
    'VariableNames', {'CountA', 'CountB', 'BitsFlipped', 'Error'});
%% Data Pruning: Remove Deterministic Extremes (0 and 32)
% This removes any trial where CountA or CountB is exactly 0 or 32
isExtreme = (data.CountA == 0 | data.CountA == 32 | ...
             data.CountB == 0 | data.CountB == 32);
data = data(~isExtreme, :); 

fprintf('Pruned deterministic extremes. Remaining active trials: %d\n', height(data));
fprintf('Done! Successfully loaded %d valid rows in %.2f seconds.\n', height(data), toc);

%% --- SANITY CHECK DIAGNOSTIC ---
if height(data) == 0 || any(isnan(data.Error))
    warning('Data might still contain issues. Checking first row elements:');
    disp(head(data, 3));
end
disp(groupcounts(data, 'CountA'));
%% Setup Figure Layout (2x2 Grid)
figure('Name', 'Stochastic Multiplier Fault Analysis Suite', ...
       'Units', 'Normalized', 'Position', [0.05, 0.05, 0.9, 0.85]);
t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title(t, '32-Bit ZCE Stochastic Multiplier Fault Injection Suite', ...
    'FontSize', 14, 'FontWeight', 'Bold');

%% 1. Error Scaling vs. Fault Intensity
nexttile;
[G1, BitsFlipped] = findgroups(data.BitsFlipped);
meanAbsError = splitapply(@mean, abs(data.Error), G1);
stdAbsError  = splitapply(@std, abs(data.Error), G1);

errorbar(BitsFlipped, meanAbsError, stdAbsError, 'b-o', ...
    'LineWidth', 1.2, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
grid on; xlim([0 33]);
xlabel('Number of Bits Flipped'); ylabel('Mean Absolute Error (MAE)');
title('1. Error Scaling & Degradation');

%% 2. Asymmetric Fault Sign Bias vs. Stream Density
nexttile;
[G2, cA] = findgroups(data.CountA);
meanRawError = splitapply(@mean, data.Error, G2);

plot(cA, meanRawError, 'r-s', 'LineWidth', 1.2, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
hold on; plot([0 32], [0 0], 'k--', 'LineWidth', 1);
grid on; xlim([0 32]);
xlabel('Count A (Faulted Stream)'); ylabel('Mean Raw Error (Faulted - Clean)');
title('2. Fault Sign Bias');

%% 3. Error Variance Spread
nexttile; 
selected_flips = [2, 4, 8, 16, 24, 32];
filtered_data = data(ismember(data.BitsFlipped, selected_flips), :);

boxchart(categorical(filtered_data.BitsFlipped), filtered_data.Error, ...
    'BoxFaceColor', [0.2 0.2 0.2], 'MarkerStyle', 'x');
grid on;
xlabel('Number of Bits Flipped'); ylabel('Error Distribution Range');
title('3. Error Variance and Bounds');

%% 4. Mathematical Coverage Verification Grid
nexttile;
scatter(data.CountA, data.CountB, 20, 'Filled', ...
    'MarkerFaceColor', [0 .5 .5], 'MarkerFaceAlpha', 0.4);
grid on; xlim([-1 33]); ylim([-1 33]);
xlabel('Count A'); ylabel('Count B');
title('4. ZCE Resolution Grid');
%% 5. Export Figure to PNG
% Get the current figure handle and save it with high resolution
fig = gcf; 
exportgraphics(fig, 'Stochastic_Multiplier_Fault_Analysis.png', 'Resolution', 300);
fprintf('Figure successfully exported and saved as PNG!\n');