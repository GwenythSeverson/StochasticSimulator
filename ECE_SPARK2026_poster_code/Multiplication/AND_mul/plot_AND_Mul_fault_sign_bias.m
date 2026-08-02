%% plot_AND_Mul_fault_sign_bias.m                     ECE SPARK 2026 -- AND multiplier
%
% Redraws ONE panel -- "Fault Sign Bias" -- from the 3x3 suite in
% tests/Fault Injection/multiplier/analyze_multiplier_faults2.m, on its own and in the poster
% palette, so it can be mounted beside the other AND_mul figures.
%
%     32bit_exhaustive_multiplier_trials.csv  ->  AND_Mul_poster_fault_sign_bias.png
%
% ------------------------------------------------------------------------------------------
% WHAT IT SHOWS. Mean RAW (signed, not absolute) total error against the ones-count of the
% faulted input stream. Every other figure in this folder reports magnitude; this one reports
% DIRECTION, and that is the whole reason it is worth a panel of its own.
%
% THE LINE CROSSES ZERO IN THE MIDDLE, and that is the result:
%     a SPARSE stream (few ones) has mostly zeros to flip, so a random flip usually ADDS a one
%       and the output error is POSITIVE
%     a DENSE stream (many ones) has mostly ones to flip, so a flip usually REMOVES one and the
%       error is NEGATIVE
%     at half density the two are equally likely and the bias vanishes
% So the fault bias of an AND-gate multiplier is not a property of the gate -- it is a property
% of the OPERAND it happens to be carrying. The sign is set by the data, not by the hardware.
%
% NOTE THE STREAM LENGTH. This CSV is the 32-BIT exhaustive campaign, not the 256-bit sweep the
% rest of this folder uses. It is a different experiment at a different resolution, and the two
% should not be read as points on one curve. It is here because the sign result does not depend
% on length and it is the cleanest place to see it.
%
% ------------------------------------------------------------------------------------------
% WHICH ERROR COLUMN, AND WHY IT IS NOT THE ONE THE ORIGINAL PANEL USED
%
% The CSV carries three error columns, and test_multiplier_fault2.cpp defines them as:
%
%     ideal_val       = countA * countB / 32          the real-valued product, unquantized
%     PrecisionError  = ExpectedCleanOnes - ideal_val rounding only; nothing to do with faults
%     BitFlipError    = FaultedOnes - ExpectedCleanOnes   <- ones ADDED or REMOVED by the flips
%     TotalError      = FaultedOnes - ideal_val       = BitFlipError + PrecisionError
%
% The original 3x3 panel plotted TOTALERROR. That is a perfectly good quantity -- it is the
% error against the true mathematical answer -- but it is NOT "ones added or removed", because
% it also carries the quantisation error of rounding the ideal product to a whole number of
% ones. Since this figure's y axis is labelled in output BITS ADDED AND REMOVED, it has to plot
% the column that actually means that, which is BITFLIPERROR.
%
% HOW MUCH DIFFERENCE DOES IT MAKE? The precision term is small next to the fault term -- at
% most 0.25 bits against a range of +/-8 -- so the overall shape is identical either way. But it
% is NOT negligible where this figure makes its point: at Count A = 16 the mean TotalError is
% +0.006 while the mean BitFlipError is -0.244, so the two cross zero one step apart --
% TotalError between 16 and 17, BitFlipError between 15 and 16. Plotting TotalError would put
% the "no bias" point in the wrong place by exactly the amount the rounding contributes.
%
% Set METRIC below to 'TotalError' to reproduce the original panel instead.
%
% ORIGINAL SOURCE: analyze_multiplier_faults2.m, "Graph 2: Fault Sign Bias". Changes here: the
% leading "2." is dropped from the title (it is no longer part of a numbered suite), the styling
% matches the poster set, the y axis is expressed in signed output bits, and the metric is
% BitFlipError as argued above.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
% AND_mul/ -> Multiplication/ -> ECE_SPARK2026_poster_code/ -> repo root
repoRoot  = fullfile(scriptDir, '..', '..', '..');

% Prefer the archived copy under keepdata/, which is the one that is meant to stay put; fall
% back to the working-directory copy the C++ suite writes.
candidates = { ...
    fullfile(repoRoot, 'keepdata', 'Multiplier Data', 'AND Mul', '32bit_exhaustive_multiplier_trials.csv'), ...
    fullfile(repoRoot, '32bit_exhaustive_multiplier_trials.csv')};
csvFile = '';
for c = candidates
    if isfile(c{1}), csvFile = c{1}; break; end
end
if isempty(csvFile)
    error(['CRITICAL: 32bit_exhaustive_multiplier_trials.csv not found in keepdata/ or at the ' ...
           'repo root.\nRun the C++ fault-injection suite first.']);
end

%% POSTER STYLE -- identical to the other four scripts in this folder set.
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

% 'BitFlipError' = ones added/removed by the flips alone -- what the y axis claims to show.
% 'TotalError'   = that plus the quantisation error, which is what the original panel plotted.
METRIC = 'BitFlipError';

%% ---------------------------------------------------------------------------------------
%% Load. Only three columns are needed, so the four bitstring columns are skipped outright --
%% textscan's %*q is what keeps this fast on a file with 32-element vectors in every row.
%% ---------------------------------------------------------------------------------------
fprintf('Ingesting %s ... ', csvFile);
tic;
fid = fopen(csvFile, 'r');
fgetl(fid);                       % header
% Trial, cA, cB, CleanOnes, Org, Flips, VecA, VecB, VecA_F, VecOut, FaultedOnes, TotalErr, ...
raw = textscan(fid, '%f %f %f %f %f %f %*q %*q %*q %*q %f %f %f %f', ...
               'Delimiter', ',', 'TreatAsEmpty', {'NA', 'NaN'});
fclose(fid);
countA = raw{2};
switch METRIC
    case 'BitFlipError', metricCol = raw{10};   % FaultedOnes - ExpectedCleanOnes
    case 'TotalError',   metricCol = raw{8};    % FaultedOnes - ideal (includes rounding)
    otherwise, error('METRIC must be ''BitFlipError'' or ''TotalError''.');
end
fprintf('%d rows in %.2f s  (metric: %s)\n', numel(countA), toc, METRIC);

[G, cA] = findgroups(countA);
meanRawTotalError = splitapply(@mean, metricCol, G);

fprintf('Count A spans %d to %d; mean raw error runs %+.3f to %+.3f\n', ...
        min(cA), max(cA), meanRawTotalError(1), meanRawTotalError(end));
% Where the bias changes sign, reported rather than assumed.
cross = find(sign(meanRawTotalError(1:end-1)) ~= sign(meanRawTotalError(2:end)), 1);
if ~isempty(cross)
    fprintf('Sign change between Count A = %d and %d\n', cA(cross), cA(cross+1));
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
plot(ax, cA, meanRawTotalError, '-s', 'Color', BLUE, 'LineWidth', 2.2, ...
     'MarkerFaceColor', BLUE, 'MarkerEdgeColor', BLUE, 'MarkerSize', 6);
hold(ax, 'off');

grid(ax, 'on');
ax.XMinorGrid = 'off';  ax.YMinorGrid = 'off';
ax.XMinorTick = 'off';  ax.YMinorTick = 'off';

% X SPANS THE WHOLE 32-BIT STREAM, 0 to 32, so the reader sees the operand at its real scale
% rather than stopping at whatever the sweep happened to reach.
%
% TWO THINGS ABOUT THE ENDS, both deliberate:
%   0 IS KEPT, not started at 1. Count A = 0 is an all-zeros stream and it is a genuine measured
%     point -- the largest positive bias on the whole figure, in fact -- so starting the axis at
%     1 would drop the single most informative sample.
%   32 IS EMPTY. The C++ sweep runs `for countA = 0; countA < STREAM_LEN`, so a fully dense
%     stream is never built. The axis still runs out to 32 to show the full operand range; the
%     line simply stops one short of it, which is the truth about this campaign.
xlim(ax, [0 32]);
xticks(ax, 0:4:32);

% Y IN OUTPUT BITS, SIGNED AND SYMMETRIC. TotalError is already a count of output ones
% (faulted minus clean), not a fraction, so no conversion is needed -- but the axis has to SAY
% so, and it has to say which direction is which. Ticks are forced symmetric about zero and
% given explicit + / - signs, because the entire content of this figure is that the sign flips.
lim  = ceil(max(abs(meanRawTotalError)));
step = max(1, round(lim / 4));
tk   = -lim:step:lim;
lab  = arrayfun(@(v) sprintf('%+d', v), tk, 'UniformOutput', false);
lab{tk == 0} = '0';
yticks(ax, tk);  yticklabels(ax, lab);
ylim(ax, [-lim lim]);

% The two halves named where they happen. The data runs top-left to bottom-right, so the
% top-right and bottom-left corners are empty and the labels sit clear of the line.
text(ax, 31, lim * 0.86, '1s ADDED', ...
     'HorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold', 'Color', STY.fg);
text(ax, 1, -lim * 0.86, '1s SUBTRACTED', ...
     'HorizontalAlignment', 'left', 'FontSize', 14, 'FontWeight', 'bold', 'Color', STY.fg);

xlabel(ax, 'Ones in the Faulted Input Stream  (of 32)', 'FontSize', 15);
ylabel(ax, {'1s Added or Subtracted on Average', 'once all faults are applied'}, 'FontSize', 15);
% Subtitle kept to two short clauses. The "sparse biases high, dense biases low" reading is
% already written into the plot itself by the two in-axes labels, so repeating it here only
% made the line overflow both edges of the panel.
title(ax, {'AND-Gate Multiplier -- Fault Sign Bias', ...
           sprintf('32-bit exhaustive campaign   |   %s trials', ...
                   commafy(numel(countA)))}, 'FontSize', 15);

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
function s = commafy(n)
    s = fliplr(regexprep(fliplr(sprintf('%d', round(n))), '(\d{3})(?=\d)', '$1,'));
end

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
