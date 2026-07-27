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
%% 3. Setup Figure Layout (chnaging grid as needed  Grid)
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


%% Graph 12: Fault Severity Error Distributions (True Discrete PMF)
nexttile;

selected_flips = [2, 8, 16, 32];
colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'}; % Distinct, high-contrast colors
hold on;

for k = 1:length(selected_flips)
    % 1. Extract the data subset
    temp = data(data.BitsFlipped == selected_flips(k), :);
    
    % 2. Calculate unique error values and their exact probabilities
    % This avoids artificial binning stairs entirely!
    [unique_errors, ~, idx] = unique(temp.TotalError);
    probabilities = accumarray(idx, 1) / numel(temp.TotalError);
    
    % 3. Plot the envelope (showing the clean, underlying bell curves)
    plot(unique_errors, probabilities, '-', 'Color', colors{k}, 'LineWidth', 2);
    
    % 4. Add discrete markers to represent the true discrete state steps of SC
    plot(unique_errors, probabilities, '.', 'Color', colors{k}, 'MarkerSize', 8);
end

% Draw the zero-error axis reference line
xline(0, 'r--', 'LineWidth', 1.2);

grid on;
box on;
xlabel('Total Error (Ones Count Deviation)');
ylabel('True State Probability');
legend('2 Bits Envelope', '2 Bits States', ...
       '8 Bits Envelope', '8 Bits States', ...
       '16 Bits Envelope', '16 Bits States', ...
       '32 Bits Envelope', '32 Bits States', ...
       'Location', 'NorthEast');
title('True Discrete Error Distribution by Fault Intensity');
xlim([-16 16]); % Focus on the active physical bounds

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
%% Isolation Graph: Degradation and Jaggedness at High Fault Densities (24-32 Flips)
% Run this snippet in your workspace containing your loaded 'data' table.

figure('Name', 'High-Fault Jaggedness Isolation (24 to 32 Flips)', ...
       'Units', 'Normalized', 'Position', [0.1, 0.1, 0.8, 0.7]);

% 1. Filter data for the high-fault regime (even steps from 24 to 32)
target_flips = [24, 26, 28, 30, 32];
colors = {[0.4660, 0.6740, 0.1880], ... % Light Green (24)
          [0.3010, 0.7450, 0.9330], ... % Cyan-Blue (26)
          [0.8500, 0.3250, 0.0980], ... % Orange-Red (28)
          [0.4940, 0.1840, 0.5560], ... % Purple (30)
          [0.6350, 0.0780, 0.1840]};    % Dark Red (32)

% 2. Create a 1x2 visualization layout
t_iso = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
title(t_iso, 'Systemic Coherence & Comb-Filtering in High-Fault SC MUX Adders', ...
    'FontSize', 14, 'FontWeight', 'Bold');

%% Left Panel: Waterfall PMF Curves (Showing the transition from smooth to jagged)
ax1 = nexttile(t_iso);
hold on;

for k = 1:length(target_flips)
    % Extract the specific fault intensity subset
    temp = data(data.BitsFlipped == target_flips(k), :);
    
    % Compute exact discrete probabilities (PMF)
    [unique_errors, ~, idx] = unique(temp.TotalError);
    probabilities = accumarray(idx, 1) / numel(temp.TotalError);
    
    % Offset each line vertically slightly to create a clean pseudo-3D waterfall view
    y_offset = (k-1) * 0.015; 
    
    % Plot the continuous envelope
    plot(unique_errors, probabilities + y_offset, '-', ...
         'Color', colors{k}, 'LineWidth', 2);
     
    % Plot the discrete state markers
    plot(unique_errors, probabilities + y_offset, '.', ...
         'Color', colors{k}, 'MarkerSize', 6);
end

xline(0, 'k--', 'LineWidth', 1.2, 'Alpha', 0.5);
grid on; box on;
xlabel('Total Error (Ones Count Deviation)');
ylabel('Probability + Vertical Offset');
title('A. Evolution of State Probabilities (24 \rightarrow 32 Flips)');

% Add custom labels to the offset curves
legend_labels = cellfun(@(x) sprintf('%d Flips', x), num2cell(target_flips), 'UniformOutput', false);
legend(legend_labels, 'Location', 'NorthWest', 'Box', 'off');
xlim([-16 16]);

%% Right Panel: The "Entropy Collapse" Metric (Quantifying the Jaggedness)
% As the curves get more jagged, the variance of adjacent states increases.
% We calculate the mean absolute difference between neighboring probability states.
ax2 = nexttile(t_iso);

jaggedness_metric = zeros(length(target_flips), 1);

for k = 1:length(target_flips)
    temp = data(data.BitsFlipped == target_flips(k), :);
    [unique_errors, ~, idx] = unique(temp.TotalError);
    probabilities = accumarray(idx, 1) / numel(temp.TotalError);
    
    % Calculate the "roughness" (first derivative of the PMF)
    % Higher values mean huge spikes and steep drops (jagged), lower means smooth transitions.
    jaggedness_metric(k) = mean(abs(diff(probabilities)));
end

bar(target_flips, jaggedness_metric, 0.5, 'FaceColor', [0.15 0.15 0.15], 'EdgeColor', 'k');
hold on;
plot(target_flips, jaggedness_metric, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');

grid on; box on;
xlabel('Number of Bits Flipped');
ylabel('PMF Roughness Index (Mean |\Delta P|)');
title('B. Quantification of Structural Comb Filtering');
xticks(target_flips);

% Link axes for easy zooming
linkaxes([ax1, ax2], 'y');
%% Export Figure

fig = gcf;

fprintf('Exporting layout... ');

exportgraphics(fig,...
    'Stochastic_Adder_Exhaustive_Fault_Analysis.png',...
    'Resolution',300,...
    'BackgroundColor','current');

fprintf('Done! Saved as PNG.\n');