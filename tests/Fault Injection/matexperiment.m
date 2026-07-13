%% Exhaustive 32-Bit Stochastic Computing Adder Simulation Framework
% Consolidates generation, MUX simulation, fault injection, and data visualization
% with native, auto-scaling textual data cards (Zero overlapping text).

clear; clc; close all;

%% 1. Configuration & Global Hyperparameters
STREAM_LEN = 32;               % Length of the stochastic bitstreams
NUMS_TO_SWEEP = 0:31;          % The 32 unique numbers represented by ones counts
ORGANIZATIONS_PER_PAIR = 15;   % Number of random layouts per pair
MAX_FLIPS = 32;                % Maximum bit flips injected into Stream A

total_numerical_pairs = length(NUMS_TO_SWEEP)^2;
total_trials = total_numerical_pairs * ORGANIZATIONS_PER_PAIR * MAX_FLIPS;

fprintf('Executing simulation loops... ');

% Preallocate arrays
countA_arr = zeros(total_trials, 1);
countB_arr = zeros(total_trials, 1);
flips_arr  = zeros(total_trials, 1);
ideal_val_arr        = zeros(total_trials, 1);
clean_ones_arr       = zeros(total_trials, 1);
faulted_ones_arr     = zeros(total_trials, 1);
total_error_arr      = zeros(total_trials, 1);
precision_error_arr  = zeros(total_trials, 1);
bit_flip_error_arr   = zeros(total_trials, 1);

%% 2. Run Simulation Loop
tic;
idx = 1;

for cA = NUMS_TO_SWEEP
    for cB = NUMS_TO_SWEEP
        ideal_val = (cA + cB) / 2; 
        
        for org = 1:ORGANIZATIONS_PER_PAIR
            vecA_orig = [true(1, cA), false(1, STREAM_LEN - cA)];
            vecA_orig = vecA_orig(randperm(STREAM_LEN));
            
            vecB_orig = [true(1, cB), false(1, STREAM_LEN - cB)];
            vecB_orig = vecB_orig(randperm(STREAM_LEN));
            
            vecS = [true(1, STREAM_LEN/2), false(1, STREAM_LEN/2)];
            vecS = vecS(randperm(STREAM_LEN));
            
            % Clean Hardware MUX
            vecOut_clean = zeros(1, STREAM_LEN);
            for i = 1:STREAM_LEN
                if vecS(i), vecOut_clean(i) = vecA_orig(i);
                else,       vecOut_clean(i) = vecB_orig(i); end
            end
            clean_ones = sum(vecOut_clean);
            
            % Sweep Faults
            for flips = 1:MAX_FLIPS
                vecA_faulted = vecA_orig;
                flip_indices = randperm(STREAM_LEN, flips);
                vecA_faulted(flip_indices) = ~vecA_faulted(flip_indices);
                
                % Faulted Hardware MUX
                vecOut_faulted = zeros(1, STREAM_LEN);
                for i = 1:STREAM_LEN
                    if vecS(i), vecOut_faulted(i) = vecA_faulted(i);
                    else,       vecOut_faulted(i) = vecB_orig(i); end
                end
                faulted_ones = sum(vecOut_faulted);
                
                % Record Metrics
                countA_arr(idx) = cA;
                countB_arr(idx) = cB;
                flips_arr(idx)  = flips;
                ideal_val_arr(idx)       = ideal_val;
                clean_ones_arr(idx)      = clean_ones;
                faulted_ones_arr(idx)    = faulted_ones;
                total_error_arr(idx)     = faulted_ones - ideal_val;
                precision_error_arr(idx) = clean_ones - ideal_val;
                bit_flip_error_arr(idx)  = faulted_ones - clean_ones;
                
                idx = idx + 1;
            end
        end
    end
end

data = table(countA_arr, countB_arr, flips_arr, ideal_val_arr, ...
             clean_ones_arr, faulted_ones_arr, total_error_arr, precision_error_arr, bit_flip_error_arr, ...
             'VariableNames', {'CountA', 'CountB', 'BitsFlipped', 'IdealVal', ...
                               'CleanOnes', 'FaultedOnes', 'TotalError', 'PrecisionError', 'BitFlipError'});

% Compute global diagnostics
global_mae_total     = mean(abs(data.TotalError));
global_mae_precision = mean(abs(data.PrecisionError));
global_mae_bitflip   = mean(abs(data.BitFlipError));

fprintf('Done in %.2f seconds.\n', toc);

%% 3. Clean Diagnostics Layout (No Floating Textboxes)
figure('Name', 'Stochastic Adder Diagnostic Suite', 'Units', 'Normalized', 'Position', [0.05, 0.05, 0.9, 0.88]);
t = tiledlayout(2, 2, 'TileSpacing', 'Loose', 'Padding', 'Normal');
title(t, '32-Bit Stochastic Adder Architecture Fault-Injection Manifest', 'FontSize', 16, 'FontWeight', 'Bold');

% --- PANEL 1: General Fault Sign Bias ---
nexttile;
[G1, cA] = findgroups(data.CountA);
meanRawTotalError = splitapply(@mean, data.TotalError, G1);

plot(cA, meanRawTotalError, 'r-s', 'LineWidth', 1.2, 'MarkerFaceColor', 'r', 'MarkerSize', 4); hold on;
mean_flips = mean(1:MAX_FLIPS);
expected_fault_error = (mean_flips / (2 * STREAM_LEN)) * (STREAM_LEN - 2 * cA);
plot(cA, expected_fault_error, 'k--', 'LineWidth', 2);
grid on; xlim([0 31]);

title('1. Aggregate Fault Sign Bias vs. Theory', 'FontSize', 12, 'FontWeight', 'Bold');
subtitle({sprintf('Data Pool: %d Total Trials Across Full Sweep', total_trials), ...
          sprintf('Mathematical Model: E[Error] = (f / 2N) * (N - 2k_A)  |  Global MAE = %.3f', global_mae_total)}, ...
          'FontSize', 9, 'FontAngle', 'italic');
xlabel('Count A (Initial Stream Density of Corrupted Input)'); 
ylabel('Mean Raw Total Error (Output minus Ideal)');
legend('Empirical Hardware Sim', 'Theoretical Model', 'Location', 'SouthWest');

% --- PANEL 2: Selected Fault Levels ---
nexttile; hold on;
target_flips = [4, 8, 16, 24, 32];
colors = lines(length(target_flips));

for j = 1:length(target_flips)
    f_val = target_flips(j);
    sub_data = data(data.BitsFlipped == f_val, :);
    [G2, cA_sub] = findgroups(sub_data.CountA);
    mean_err_sub = splitapply(@mean, sub_data.BitFlipError, G2);
    plot(cA_sub, mean_err_sub, '-o', 'Color', colors(j,:), 'LineWidth', 1.5, ...
         'MarkerFaceColor', colors(j,:), 'MarkerSize', 3, 'DisplayName', sprintf('%d Flips', f_val));
end
plot([0 31], [0 0], 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
grid on; xlim([0 31]);

title('2. Error Trajectories Across Selected Fault Levels', 'FontSize', 12, 'FontWeight', 'Bold');
subtitle({'Data Pool: Isolated Subsets of Specific Fault Intensities', ...
          'Architectural Behavior: Error paths pivot symmetrically at 16 ones (50% density), where E[Err] = 0.'}, ...
          'FontSize', 9, 'FontAngle', 'italic');
xlabel('Initial Stream A Density (Ones Count from 0 to 31)'); 
ylabel('Mean Pure Bit-Flip Error (Output minus Clean)');
legend('Location', 'SouthWest');

% --- PANEL 3: Error Components Scaling ---
nexttile;
[G3, unique_flips] = findgroups(data.BitsFlipped);
mean_abs_precision = splitapply(@mean, abs(data.PrecisionError), G3);
mean_abs_bitflip   = splitapply(@mean, abs(data.BitFlipError), G3);

plot(unique_flips, mean_abs_precision, 'g-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'g'); hold on;
plot(unique_flips, mean_abs_bitflip, 'm-v', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
plot([1 32], [0.5 0.5], 'r:', 'LineWidth', 2); 
grid on; xlim([1 32]); ylim([0 max(mean_abs_bitflip)+1]);

title('3. Error Roots: Resolution Limit vs. Fault Progression', 'FontSize', 12, 'FontWeight', 'Bold');
subtitle({sprintf('Data Pool: Grouped Strictly by Injected Bit Flips (X-axis)  |  Quantization MAE = %.3f', global_mae_precision), ...
          sprintf('Hardware Insight: Quantization error is static; Environmental Bit-Flip error (MAE = %.3f) scales linearly.', global_mae_bitflip)}, ...
          'FontSize', 9, 'FontAngle', 'italic');
xlabel('Fault Intensity (Total Number of Random Bits Flipped in Stream A)'); 
ylabel('Mean Absolute Error Magnitude (Ones Count)');
legend('Quantization (Precision) Error', 'Bit Flip (Fault) Error', 'Worst-Case Quantization Limit (\pm0.5)', 'Location', 'NorthWest');

% --- PANEL 4: Output State Vulnerability Envelope ---
nexttile;
[G4, expected_out] = findgroups(data.CleanOnes);
mean_raw_err = splitapply(@mean, data.TotalError, G4);
std_raw_err  = splitapply(@std, data.TotalError, G4);

fill([expected_out; flipud(expected_out)], [mean_raw_err + std_raw_err; flipud(mean_raw_err - std_raw_err)], ...
     [1 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
plot(expected_out, mean_raw_err, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
plot([0 31], [0 0], 'k--', 'LineWidth', 1);
grid on; xlim([0 31]);

title('4. Output State Vulnerability Envelope', 'FontSize', 12, 'FontWeight', 'Bold');
subtitle({'Data Pool: Grouped by Target Clean Output State Value', ...
          'Hardware Insight: Shaded standard deviation (\sigma) envelope displays the localized effect of MUX logic masking.'}, ...
          'FontSize', 9, 'FontAngle', 'italic');
xlabel('Expected Clean Output Value (Ones Count Target)'); 
ylabel('Mean Raw Total Error Deviation');
legend('Empirical Hardware Variance Window (\sigma)', 'Mean Hardware Profile', 'Location', 'SouthWest');