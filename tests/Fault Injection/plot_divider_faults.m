%% plot_divider_error_isolation.m
% =========================================================================
% STOCHASTIC DIVIDER ERROR ISOLATION PIPELINE
% Classifies and plots fault metrics by physical hardware zones
% to isolate exactly where errors occur in the architecture.
% =========================================================================
clear; clc; close all;

fprintf('Reading dataset: divider_exhaustive_data.csv...\n');
opts = detectImportOptions('divider_exhaustive_data.csv');
data = readtable('divider_exhaustive_data.csv', opts);

stream_len = 32;
flip_axis = 1:stream_len;

%% =========================================================================
%% FIGURE 1: Absolute Fault Tracking for Core Profiles (Preserved)
%% =========================================================================
figure('Position', [100, 100, 800, 400]);
idx_mid  = (data.CountA == 16 & data.CountB == 16);
idx_low  = (data.CountA == 2  & data.CountB == 2);
idx_high = (data.CountA == 30 & data.CountB == 30);

hold on;
plot(data.Flips(idx_mid),  abs(data.Error(idx_mid)),  '-o', 'LineWidth', 2, 'Color', [0.12, 0.47, 0.71], 'DisplayName', 'Midscale (0.5 / 0.5)');
plot(data.Flips(idx_low),  abs(data.Error(idx_low)),  '-^', 'LineWidth', 2, 'Color', [0.20, 0.63, 0.17], 'DisplayName', 'Low Edge (0.05 / 0.05)');
plot(data.Flips(idx_high), abs(data.Error(idx_high)), '-s', 'LineWidth', 2, 'Color', [0.89, 0.10, 0.11], 'DisplayName', 'High Edge (0.95 / 0.95)');

title('Absolute Output Error Progression vs. Injected Bit Flips', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Number of Flipped Bits in Stream A', 'FontSize', 10);
ylabel('Absolute Error |FaultedOut - CleanOut|', 'FontSize', 10);
grid on; legend('Location', 'NorthWest');
hold off;

%% =========================================================================
%% FIGURE 2: Error Isolation by Functional Operating Zones
%% =========================================================================
figure('Position', [150, 150, 800, 400]);
hold on;

% Isolate trials into three distinct physical architectural zones:
% Zone 1: Fraction Zone (Numerator is significantly smaller than Denominator)
idx_fraction = (data.CountA < data.CountB);
% Zone 2: Balanced Point Zone (Numerator is equal or very close to Denominator)
idx_balanced = (data.CountA == data.CountB);
% Zone 3: Saturation Zone (Numerator is greater than Denominator, forcing clean output to 1.0)
idx_saturate = (data.CountA > data.CountB);

err_fraction = zeros(stream_len, 1);
err_balanced = zeros(stream_len, 1);
err_saturate = zeros(stream_len, 1);

for f = 1:stream_len
    err_fraction(f) = mean(abs(data.Error(idx_fraction & data.Flips == f)));
    err_balanced(f) = mean(abs(data.Error(idx_balanced & data.Flips == f)));
    err_saturate(f) = mean(abs(data.Error(idx_saturate & data.Flips == f)));
end

plot(flip_axis, err_fraction, '-o', 'LineWidth', 2, 'Color', [0.49, 0.18, 0.56], 'DisplayName', 'Fractional Zone (X < Y)');
plot(flip_axis, err_balanced, '-x', 'LineWidth', 2, 'Color', [0.12, 0.47, 0.71], 'DisplayName', 'Balanced Zone (X = Y)');
plot(flip_axis, err_saturate, '-s', 'LineWidth', 2, 'Color', [0.89, 0.10, 0.11], 'DisplayName', 'Saturated Zone (X > Y)');

title('Error Isolation across Physical Operational Zones', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Bit-Flip Injection Depth', 'FontSize', 10);
ylabel('Mean Absolute Error', 'FontSize', 10);
grid on; legend('Location', 'NorthWest');
hold off;

%% =========================================================================
%% FIGURE 3: Error Isolation by Denominator Sparsity (Y Density)
%% =========================================================================
figure('Position', [200, 200, 800, 400]);
hold on;

% Isolate trials based entirely on how "full" the denominator stream is
% Sparse Denominator: CountB <= 4 (Very few 1s, highly unstable division)
idx_sparse_Y = (data.CountB <= 4 & data.CountB > 0);
% Dense Denominator: CountB >= 24 (Packed with 1s, counter pulls down strongly)
idx_dense_Y  = (data.CountB >= 24);

err_sparse_Y = zeros(stream_len, 1);
err_dense_Y  = zeros(stream_len, 1);

for f = 1:stream_len
    err_sparse_Y(f) = mean(abs(data.Error(idx_sparse_Y & data.Flips == f)));
    err_dense_Y(f)  = mean(abs(data.Error(idx_dense_Y & data.Flips == f)));
end

plot(flip_axis, err_sparse_Y, '-v', 'LineWidth', 2, 'Color', [0.85, 0.32, 0.10], 'DisplayName', 'Sparse Denominator (CountB \leq 4)');
plot(flip_axis, err_dense_Y,  '-d', 'LineWidth', 2, 'Color', [0.20, 0.63, 0.17], 'DisplayName', 'Dense Denominator (CountB \geq 24)');

title('Error Isolation: Denominator Density Vulnerability', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Bit-Flip Injection Depth', 'FontSize', 10);
ylabel('Mean Absolute Error', 'FontSize', 10);
grid on; legend('Location', 'NorthWest');
hold off;

fprintf('[SUCCESS] New error isolation graphs rendered.\n');