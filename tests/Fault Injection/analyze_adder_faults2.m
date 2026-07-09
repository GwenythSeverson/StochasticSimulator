%% Stochastic Computing Adder (MUX) Exhaustive Fault-Injection Analysis
% This script parses the exhaustive 32-bit adder fault injection CSV, generates
% the same graph suite as the multiplier analysis, and saves results to .mat.

clear; clc; close all;

csvFile = '32bit_exhaustive_adder_trials.csv';
if ~isfile(csvFile)
    error('CRITICAL: %s not found. Please run the C++ simulator first.', csvFile);
end

%% 1. Ingest Data via Text Scanning (High Performance)
fprintf('Ingesting exhaustive adder CSV file... ');
tic;

fid = fopen(csvFile, 'r');
fgetl(fid); % Skip header line

% Column format: Trial(f), cA(f), cB(f), CleanOnes(f), Org(f), Flips(f),
% VecA(*q), VecB(*q), VecA_F(*q), VecOut(*q), FaultedOnesCount(f),
% TotalError(f), PrecisionError(f), BitFlipError(f)
formatSpec = '%f %f %f %f %f %f %*q %*q %*q %*q %f %f %f %f';
rawCells = textscan(fid, formatSpec, 'Delimiter', ',', 'TreatAsEmpty', {'NA','NaN'});
fclose(fid);

% Parse non-skipped cells into table
data = table(rawCells{2}, rawCells{3}, rawCells{4}, rawCells{6}, rawCells{8}, rawCells{9}, rawCells{10}, ...
    'VariableNames', {'CountA', 'CountB', 'CleanOnes', 'BitsFlipped', 'TotalError', 'PrecisionError', 'BitFlipError'});

fprintf('Done! Loaded %d rows in %.2f seconds.\n', height(data), toc);

%% 2. Save Workspace to .mat File
matFile = 'Stochastic_Adder_Exhaustive_Data.mat';
fprintf('Saving workspace variables to %s... ', matFile);
save(matFile);
fprintf('Done!\n');

%% 3. Setup Figure Layout (3x3 Grid)
figure('Name', 'Stochastic Adder Exhaustive Fault Analysis Suite', ...
       'Units', 'Normalized', 'Position', [0.05, 0.05, 0.9, 0.9]);
t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title(t, '32-Bit Exhaustive Stochastic Adder (MUX) Fault Injection Suite', ...
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
    'BoxFaceColor', [0.2 0.2 0.2], 'MarkerStyle', 'x');
grid on;
xlabel('Number of Bits Flipped'); ylabel('Total Error Distribution Range');
title('3. Total Error Variance & Bounds');

%% Graph 4: Exhaustive ZCE Resolution Grid
nexttile;
scatter(data.CountA, data.CountB, 20, 'Filled', ...
    'MarkerFaceColor', [0 .5 .5], 'MarkerFaceAlpha', 0.1);
grid on; xlim([-1 32]); ylim([-1 32]);
xlabel('Count A'); ylabel('Count B');
title('4. Exhaustive 32x32 Resolution Grid');

%% Graph 5: Precision Error vs. Bit Flip Error (Line Plot)
nexttile;
meanAbsPrecision = splitapply(@mean, abs(data.PrecisionError), G1);
meanAbsBitFlip   = splitapply(@mean, abs(data.BitFlipError), G1);

plot(BitsFlipped, meanAbsPrecision, 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'MarkerSize', 5);
hold on;
plot(BitsFlipped, meanAbsBitFlip, 'm-v', 'LineWidth', 1.5, 'MarkerFaceColor', 'm', 'MarkerSize', 5);
grid on; xlim([0 33]);
xlabel('Number of Bits Flipped'); ylabel('Mean Absolute Error');
legend('Precision (Quantization) Error', 'Bit Flip (Fault) Error', 'Location', 'NorthWest');
title('5. Precision vs. Bit Flip Error Scaling');

%% Graph 6: Resulting Bitstream Fault Sign Bias
nexttile;
[G6, cleanOnes] = findgroups(data.CleanOnes);
meanRawErrorOut = splitapply(@mean, data.TotalError, G6);

plot(cleanOnes, meanRawErrorOut, 'm-o', 'LineWidth', 1.2, 'MarkerFaceColor', 'm', 'MarkerSize', 4);
hold on; plot([0 31], [0 0], 'k--', 'LineWidth', 1);
grid on; xlim([0 31]);
xlabel('Expected Clean Output Ones Count'); ylabel('Mean Raw Total Error');
title('6. Resulting Bitstream Fault Sign Bias');

%% Graph 7: Clean Output Trial Distribution
nexttile([1, 2]);
[cleanOnesCounts, cleanOnesValues] = groupcounts(data.CleanOnes);
bar(cleanOnesValues, cleanOnesCounts, 'FaceColor', [0.4 0.6 0.8]);
grid on; xlim([-1 32]);
xlabel('Expected Clean Output Ones Count'); ylabel('Number of Trials');
title('7. Clean Output Trial Distribution');

%% Graph 8: Input Stream A Trial Distribution
nexttile;
[cA_Counts, cA_Values] = groupcounts(data.CountA);
bar(cA_Values, cA_Counts, 'FaceColor', [0.8 0.4 0.4]);
grid on; xlim([-1 32]);
xlabel('Count A (Input Ones Count)'); ylabel('Number of Trials');
title('8. Input Stream A Trial Distribution');

%% 4. Export Figure to PNG
fig = gcf; 
exportgraphics(fig, 'Stochastic_Adder_Exhaustive_Fault_Analysis.png', 'Resolution', 300, 'BackgroundColor', 'current');
fprintf('Figure successfully exported and saved as PNG!\n');
