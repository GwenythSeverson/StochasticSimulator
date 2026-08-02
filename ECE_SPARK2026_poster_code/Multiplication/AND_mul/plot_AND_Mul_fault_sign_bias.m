%% plot_AND_Mul_fault_sign_bias.m                     ECE SPARK 2026 -- AND multiplier
%
% Reads AND_Mul_sign_bias.csv (written by AND_Mul_sign_bias.cpp in this folder) and saves
% AND_Mul_poster_fault_sign_bias.png here. Nothing outside this folder is touched.
%
% ------------------------------------------------------------------------------------------
% WHAT IT SHOWS. The average number of output 1s ADDED or SUBTRACTED by faults, against the
% ones-count of the faulted input stream. Every other AND_mul figure reports error MAGNITUDE;
% this one reports DIRECTION, which is why it is worth its own panel.
%
% THE LINE CROSSES ZERO AT HALF DENSITY, and that is the result:
%     a SPARSE stream (few 1s) is mostly 0s, so a flip usually ADDS a one   -> positive
%     a DENSE stream (many 1s) is mostly 1s, so a flip usually REMOVES one -> negative
%     at 16 of 32 the two are equally likely and the bias is exactly zero
% So the fault bias of an AND-gate multiplier is not a property of the gate. It is a property
% of the OPERAND it happens to be carrying: the data sets the sign, not the hardware.
%
% ------------------------------------------------------------------------------------------
% IT IS A STRAIGHT LINE, AND THAT IS THE POINT OF REGENERATING IT
%
% The original panel (analyze_multiplier_faults2.m, "Graph 2") was drawn from
% 32bit_exhaustive_multiplier_trials.csv, which sweeps only a few sampled arrangements per
% operand pair and only part of the flip range. It wobbled by a few tenths of a bit around the
% trend, and that wobble was sampling noise, not physics.
%
% AND_Mul_sign_bias.cpp computes the same quantity EXHAUSTIVELY -- every Count A, every Count B,
% every flip count 0 to 32, every combination of flips within each count, and every arrangement
% of the two streams, each carrying its exact multiplicity. Done that way the answer is
%
%       mean 1s added/subtracted  =  8 - CountA/2
%
% exactly. The .cpp asserts it to 1e-9 and measures a worst deviation of 5e-14, i.e. floating
% point. Any curvature on this figure would be a bug.
%
% ONE THING THE EXHAUSTIVE VERSION FIXED. An intermediate version averaged over the campaign's
% canonical "zero correlation error" overlap, round(cA*cB/32), rather than over all
% arrangements. Rounding the overlap to a whole number of 1s sends every half-integer up, which
% biases every odd Count B by -1 and put a real 0.24-bit kink at Count A = 16 -- right where
% this figure makes its point. Averaging over arrangements removes the rounding entirely, and
% only then is the line exact. See the .cpp header.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
csvFile   = fullfile(scriptDir, 'AND_Mul_sign_bias.csv');
if ~isfile(csvFile)
    error(['CRITICAL: %s not found.\n' ...
           'Run:  stochastic_computer.exe --gtest_filter=AndMulSignBias.*'], csvFile);
end

%% POSTER STYLE -- identical to the other scripts in this folder set.
STY = struct( ...
    'bg',    [0.62 0.62 0.62], ...   % the PLOT PANEL only
    'figbg', [1.00 1.00 1.00], ...   % figure surround + exported margin: WHITE
    'fg',    [0.10 0.10 0.10], ...   % axes, ticks, labels, titles
    'grid',  [1.00 1.00 1.00], ...   % divider lines, WHITE on the grey panel
    'tick',  14, ...                 % tick labels
    'label', 15, ...                 % axis labels
    'title', 15);                    % titles
MARGIN_PX = 40;
BLUE = [0.10 0.28 0.85];
RED  = [0.90 0.20 0.20];

%% ---------------------------------------------------------------------------------------
%% Load and collapse to one point per Count A
%% ---------------------------------------------------------------------------------------
T = readtable(csvFile, 'Delimiter', ',');
% One row per (Count A, flip count). Averaging over flip count uniformly is what "from 0
% faulted to all 32 faulted" means: every flip level counts once.
[G, cA] = findgroups(T.CountA);
meanDelta = splitapply(@mean, T.MeanOnesDelta, G);

fprintf('Loaded %d rows, %d Count A values (%d to %d), %d flip levels each\n', ...
        height(T), numel(cA), min(cA), max(cA), height(T)/numel(cA));
fprintf('Mean 1s delta runs %+.4f at Count A = %d  to  %+.4f at Count A = %d\n', ...
        meanDelta(1), cA(1), meanDelta(end), cA(end));
fprintf('Worst deviation from the exact line 8 - CountA/2: %.3e\n', ...
        max(abs(meanDelta - (8 - cA/2))));
zero_at = cA(abs(meanDelta) < 1e-9);
if ~isempty(zero_at)
    fprintf('Bias is exactly zero at Count A = %d (half density)\n', zero_at(1));
end

%% ---------------------------------------------------------------------------------------
%% FIGURE
%% ---------------------------------------------------------------------------------------
fig = figure('Name', 'AND Multiplier -- Fault Sign Bias', ...
             'Units', 'Normalized', 'Position', [0.14, 0.12, 0.54, 0.72]);
ax = axes(fig);
hold(ax, 'on');

% Red dashed at zero, matching the other figures where the red dashed line is always "the
% answer you wanted" -- here that is zero signed error, i.e. no bias.
yline(ax, 0, '--', 'Color', RED, 'LineWidth', 2.2, ...
      'Label', 'no bias', 'LabelHorizontalAlignment', 'right', ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 12);
plot(ax, cA, meanDelta, '-s', 'Color', BLUE, 'LineWidth', 2.4, ...
     'MarkerFaceColor', BLUE, 'MarkerEdgeColor', BLUE, 'MarkerSize', 6);
hold(ax, 'off');

grid(ax, 'on');
ax.XMinorGrid = 'off';  ax.YMinorGrid = 'off';
ax.XMinorTick = 'off';  ax.YMinorTick = 'off';

% X runs 0 to 32 and so does the DATA -- the sweep now covers a fully dense stream, so there is
% no empty space at the right-hand end the way there was when the campaign stopped at 31.
xlim(ax, [0 32]);
xticks(ax, 0:4:32);

% Y in output 1s, signed and symmetric about zero, because the whole content of this figure is
% that the sign flips. Explicit + / - on every tick.
lim  = ceil(max(abs(meanDelta)));
step = max(1, round(lim / 4));
tk   = -lim:step:lim;
lab  = arrayfun(@(v) sprintf('%+d', v), tk, 'UniformOutput', false);
lab{tk == 0} = '0';
yticks(ax, tk);  yticklabels(ax, lab);
ylim(ax, [-lim lim]);

% The two halves named where they happen. The line runs top-left to bottom-right, so the
% top-right and bottom-left corners are empty and the labels sit clear of it.
text(ax, 31, lim * 0.86, '1s ADDED', ...
     'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Color', STY.fg);
text(ax, 1, -lim * 0.86, '1s SUBTRACTED', ...
     'HorizontalAlignment', 'left', 'FontSize', 14, 'FontWeight', 'bold', 'Color', STY.fg);

xlabel(ax, 'Ones in the Faulted Input Stream  (of 32)', 'FontSize', 15);
ylabel(ax, {'1s Added or Subtracted on Average', 'once all faults are applied'}, 'FontSize', 15);
title(ax, {'AND-Gate Multiplier -- Fault Sign Bias', ...
           'exhaustive: every operand pair, every flip count 0-32, every combination'}, ...
      'FontSize', 15);

style_axes(ax, STY);
set(fig, 'Color', STY.figbg);
drawnow;

png = fullfile(scriptDir, 'AND_Mul_poster_fault_sign_bias.png');
exportgraphics(fig, png, 'Resolution', 300, 'BackgroundColor', 'current');
pad_png(png, MARGIN_PX);
fprintf('\nFigure saved to %s\n', png);

% ---------------------------------------------------------------------------------------
% Helpers, identical to the other scripts in this folder set.
% ---------------------------------------------------------------------------------------
function style_axes(ax, STY)
    ax.Color      = STY.bg;
    ax.XColor     = STY.fg;
    ax.YColor     = STY.fg;
    ax.GridColor  = STY.grid;
    ax.GridAlpha  = 0.9;   % near-opaque, so the dividers read as clean white
    ax.FontSize   = STY.tick;
    ax.LineWidth  = 1.0;
    ax.Box        = 'on';
    ax.Title.Color  = STY.fg;  ax.Title.FontSize  = STY.title;
    ax.XLabel.Color = STY.fg;  ax.XLabel.FontSize = STY.label;
    ax.YLabel.Color = STY.fg;  ax.YLabel.FontSize = STY.label;
end

function pad_png(file, px)
    img = imread(file);
    if size(img, 3) == 1, img = repmat(img, 1, 1, 3); end
    bg = img(1, 1, :);
    [h, w, ~] = size(img);
    out = repmat(bg, h + 2*px, w + 2*px, 1);
    out(px+1:px+h, px+1:px+w, :) = img;
    imwrite(out, file);
end
