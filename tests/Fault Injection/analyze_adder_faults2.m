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
t = tiledlayout(3,4,'TileSpacing','Compact','Padding','Compact');
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

%% Graph 6: Input Stream A & B Trial Distributions (Grouped)
nexttile;
[cA_Counts, cA_Values] = groupcounts(data.CountA);
[cB_Counts, cB_Values] = groupcounts(data.CountB);

% Render side-by-side grouped bars for both inputs
bar(cA_Values, [cA_Counts, cB_Counts], 'grouped');
grid on; xlim([-1 32]);
xlabel('Input Ones Count'); ylabel('Number of Trials');
legend('Stream A (Count A)', 'Stream B (Count B)', 'Location', 'SouthOutside', 'Orientation', 'Horizontal');
title('6. Input Stream A & B Trial Distributions');



%% Graph 7: Trial-Level Error Scatter
nexttile;

scatter(data.BitsFlipped, ...
        abs(data.TotalError), ...
        3, ...                 % Marker size
        '.', ...
        'MarkerEdgeAlpha',0.10);

grid on;
box on;

xlabel('Number of Bits Flipped');
ylabel('|Total Error|');

title('8. Trial-Level Error Scatter');

xlim([0 33]);




%% Graph 11: Improved Global Error Distribution

nexttile;

% Use more bins to show distribution shape
histogram(data.TotalError,...
          150,...
          'Normalization','pdf',...
          'EdgeColor','k');

hold on;

xline(0,...
      'r--',...
      'LineWidth',1.5);

grid on;
box on;

xlabel('Total Error');
ylabel('Probability Density');

title('11. Global Signed Error Distribution');


%% Graph 12: Fault Severity Error Distributions

nexttile;

selected_flips = [2 8 16 32];

hold on;

for k = 1:length(selected_flips)

    temp = data(data.BitsFlipped == selected_flips(k),:);

    histogram(temp.TotalError,...
              80,...
              'Normalization','pdf',...
              'DisplayStyle','stairs',...
              'LineWidth',1.5);

end

xline(0,...
      'r--',...
      'LineWidth',1);

grid on;
box on;

xlabel('Total Error');
ylabel('Probability Density');

legend('2 Bits','8 Bits','16 Bits','32 Bits');

title('Error Distribution by Fault Intensity');

%% Graph: Error Contribution Per Flipped Bit
% Tests whether each injected fault contributes equal error

nexttile;

% Group data by number of flipped bits
[G_fault, flipCount] = findgroups(data.BitsFlipped);

% Mean absolute error at each fault level
meanAbsError = splitapply(@mean, abs(data.TotalError), G_fault);

% Normalize by number of faults
errorPerBit = meanAbsError ./ flipCount;

% Avoid zero fault division
errorPerBit(flipCount == 0) = NaN;


plot(flipCount,...
     errorPerBit,...
     'b-o',...
     'LineWidth',2,...
     'MarkerSize',5,...
     'MarkerFaceColor','b');


grid on;
box on;

xlabel('Number of Bits Flipped');
ylabel('Mean |Error| / Flipped Bit');

title('Per-Fault Error Contribution');

xlim([0 33]);
%% Export Figure

fig = gcf;

fprintf('Exporting layout... ');

exportgraphics(fig,...
    'Stochastic_Adder_Exhaustive_Fault_Analysis.png',...
    'Resolution',300,...
    'BackgroundColor','current');

fprintf('Done! Saved as PNG.\n');