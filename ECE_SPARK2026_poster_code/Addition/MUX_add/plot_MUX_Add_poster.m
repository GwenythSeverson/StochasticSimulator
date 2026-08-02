%% plot_MUX_Add_poster.m                              ECE SPARK 2026 -- MUX scaled adder
%
% Companion to AND_mul/plot_AND_Mul_poster.m. Same axes, same marks, same conventions, so the
% multiplication and addition figures read as one set.
%
%     MUX_Add_poster_warmup.csv  -> plot 1, zero-fault early termination
%     MUX_Add_poster_faults.csv  -> plot 2, fault response 0..36 flipped bits
%
% MUX_Add_poster_sensitivity.csv is also written by the .cpp but is NOT plotted here. It is the
% exact single-flip damage for all 24 site classes, and it is reported as a table in
% HOW_THIS_DATA_WAS_MADE.txt section 4 instead -- the numbers are more use written out than as
% a bar chart, since there are only two of them (0 and 1/256).
%
% ------------------------------------------------------------------------------------------
% THE UNIT. out = sel ? a : b. With a 0.5 select stream this is the standard stochastic scaled
% adder, (a + b)/2. THREE 256-bit operand streams, so 768 flippable bits -- the select stream is
% as real and as exposed as the data, and this sweep includes it.
%
% A single flip on this surface has only two possible effects: nothing at all, or exactly
% +/- 1/256. Which one depends on the cycle:
%     a flip   matters only when sel = 1 (otherwise b is being routed)
%     b flip   matters only when sel = 0
%     sel flip matters only when a and b DISAGREE
% Exactly half of each stream is therefore inert, so 384 of the 768 sites are live and every one
% of them is worth the same 1/256. Nothing on this surface is more dangerous than anything else,
% which is the flat-risk profile a fully unary unit has and uMUL does not.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
warmCsv = fullfile(scriptDir, 'MUX_Add_poster_warmup.csv');
faultCsv = fullfile(scriptDir, 'MUX_Add_poster_faults.csv');
for fchk = {warmCsv, faultCsv}
    if ~isfile(fchk{1})
        error(['CRITICAL: %s not found.\n' ...
               'Run:  stochastic_computer.exe --gtest_filter=MuxAddPoster.*'], fchk{1});
    end
end

TARGET     = 0.5;
SHOW_FLIPS = [0 1 2 4 8 16 32];
TAIL       = 0.05;

%% POSTER STYLE -- identical in all four scripts in this folder set.
%
% MEDIUM GREY, NOT THE MATLAB DEFAULT DARK. Mounted on a light poster, a near-black panel reads
% as a hole punched in the page and pulls the eye away from the text around it. A medium grey
% sits quietly next to white card while still separating the plot from the background.
%
% Everything else follows: on mid-grey, white text loses contrast, so the axes, labels and titles
% go NEAR-BLACK and the grid goes light, reading as a soft guide rather than as data. The data
% colours are unchanged apart from the blue, darkened from [0.20 0.40 0.95] so it does not go
% muddy against grey in print.
%
% Font sizes are set HERE and applied in style_axes, which overrides whatever each xlabel/title
% call asked for. One place to tune, and it keeps the four scripts identical.
STY = struct( ...
    'bg',    [0.62 0.62 0.62], ...   % the PLOT PANEL only
    'figbg', [1.00 1.00 1.00], ...   % figure surround + exported margin: WHITE
    'fg',    [0.10 0.10 0.10], ...   % axes, ticks, labels, titles
    'grid',  [1.00 1.00 1.00], ...   % divider lines, WHITE on the grey panel
    'legbg', [0.72 0.72 0.72], ...   % legend panel, a touch lighter so it lifts off the axes
    'tick',  14, ...                 % tick labels
    'label', 15, ...                 % axis labels
    'title', 15);                    % titles
MARGIN_PX = 40;                      % border around every exported PNG; 300 DPI ~ 3.4 mm

W = readtable(warmCsv,  'Delimiter', ',');
F = readtable(faultCsv, 'Delimiter', ',');
BLUE = [0.10 0.28 0.85];
RED  = [0.90 0.20 0.20];

% Equal-tailed middle 90%: 5th and 95th percentiles, nothing forced symmetric. Written the same
% way in every script in this poster set so one estimator is used throughout.
function [lo, hi] = tail_band(vals, probs, tail)
    [vals, ord] = sort(vals);
    probs = probs(ord);  probs = probs / sum(probs);
    lo = vals(find(cumsum(probs) >= tail, 1));
    hi = vals(find(cumsum(probs, 'reverse') >= tail, 1, 'last'));
end

flipLevels = unique(F.BitsFlipped);
nf = numel(flipLevels);
mn = zeros(nf,1); mx = zeros(nf,1); mu = zeros(nf,1);
bLo = zeros(nf,1); bHi = zeros(nf,1);
for k = 1:nf
    sel = F.BitsFlipped == flipLevels(k);
    v = F.OutputFraction(sel);  p = F.Probability(sel);  p = p / sum(p);
    mn(k) = min(v);  mx(k) = max(v);  mu(k) = sum(v .* p);
    [bLo(k), bHi(k)] = tail_band(v, p, TAIL);
end

fprintf('MUX adder: %d outcome rows across %d flip levels, 768-bit operand surface\n', ...
        height(F), nf);
fprintf('\n  flips     mean      min       max      5th pct      95th pct\n');
for k = 1:nf
    if ismember(flipLevels(k), SHOW_FLIPS)
        fprintf('  %5d  %8.6f  %8.6f  %8.6f  %11.6f  %12.6f\n', ...
                flipLevels(k), mu(k), mn(k), mx(k), bLo(k), bHi(k));
    end
end

%% ---------------------------------------------------------------------------------------
%% FIGURE 1 -- zero-fault early termination
%% ---------------------------------------------------------------------------------------
fig1 = figure('Name', 'MUX Adder -- Zero-Fault Warm-Up', ...
              'Units', 'Normalized', 'Position', [0.10, 0.14, 0.52, 0.70]);
ax1 = axes(fig1);
plot(ax1, W.N, W.MeanAbsError, '-o', 'LineWidth', 2.4, 'MarkerSize', 8, ...
     'Color', [0.10 0.28 0.85], 'MarkerFaceColor', [0.10 0.28 0.85]);
set(ax1, 'XScale', 'log');  grid(ax1, 'on');
ax1.XMinorGrid = 'off';  ax1.YMinorGrid = 'off';
ax1.XMinorTick = 'off';  ax1.YMinorTick = 'off';
xticks(ax1, W.N); xticklabels(ax1, string(W.N));
xlim(ax1, [min(W.N)*0.85, max(W.N)*1.15]);
xlabel(ax1, 'Bit Stream Length  (early termination point) [log]', 'FontSize', 13);
ylabel(ax1, 'Mean Absolute Error', 'FontSize', 13);
title(ax1, {'MUX Scaled Adder -- Zero-Fault Warm-Up', ...
            sprintf(['(0.5 + 0.5)/2 = 0.5   |   exhaustive over all 10^{%.0f} decorrelated ' ...
                     'arrangements'], W.Log10Arrangements(1))}, 'FontSize', 13);
ax1.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% FIGURE 2 -- fault response
%% ---------------------------------------------------------------------------------------
fig2 = figure('Name', 'MUX Adder -- Fault Response', ...
              'Units', 'Normalized', 'Position', [0.14, 0.10, 0.56, 0.72]);
ax2 = axes(fig2);
hold(ax2, 'on');

% Zero has no place on a log axis, so the 0-flip stick is parked at x = 0.5 and labelled "0".
% Its position is cosmetic; its height is the real measurement.
xpos = SHOW_FLIPS;  xpos(xpos == 0) = 0.5;
capBlue = 1.055;  capRed = 1.075;

for s = 1:numel(SHOW_FLIPS)
    k = find(flipLevels == SHOW_FLIPS(s), 1);
    if isempty(k), continue; end
    x = xpos(s);
    plot(ax2, [x x], [mn(k) mx(k)], '-', 'Color', BLUE, 'LineWidth', 1.6);
    plot(ax2, [x/capBlue x*capBlue], [mn(k) mn(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);
    plot(ax2, [x/capBlue x*capBlue], [mx(k) mx(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);
    if bLo(k) > mn(k) || bHi(k) < mx(k)
        plot(ax2, [x/capRed x*capRed], [bLo(k) bLo(k)], '-', 'Color', RED, 'LineWidth', 2.4);
        plot(ax2, [x/capRed x*capRed], [bHi(k) bHi(k)], '-', 'Color', RED, 'LineWidth', 2.4);
    end
end
selShow = ismember(flipLevels, SHOW_FLIPS);
hMean = plot(ax2, xpos, mu(selShow), 'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', BLUE, 'MarkerEdgeColor', BLUE, 'LineWidth', 1.0);
yline(ax2, TARGET, 'r--', 'LineWidth', 2.2, 'Label', 'intended result = 0.5', ...
      'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'FontSize', 11);
hold(ax2, 'off');

set(ax2, 'XScale', 'log');  grid(ax2, 'on');
ax2.XMinorGrid = 'off';  ax2.YMinorGrid = 'off';
ax2.XMinorTick = 'off';  ax2.YMinorTick = 'off';
xlim(ax2, [0.38 max(SHOW_FLIPS) * 2.6]);
ylim(ax2, [0.25 0.75]);
xticks(ax2, xpos);  xticklabels(ax2, string(SHOW_FLIPS));
xlabel(ax2, 'Bits Flipped  (over the 768-bit operand surface, a + b + select) [log]', 'FontSize', 12);
ylabel(ax2, 'Fraction of Ones in the Final Answer', 'FontSize', 13);
title(ax2, {'MUX Scaled Adder -- Fault Response', ...
            'three operand streams, 768 flippable bits, select INCLUDED'}, 'FontSize', 13);
hold(ax2, 'on');
pB = plot(ax2, NaN, NaN, '-', 'Color', BLUE, 'LineWidth', 1.8);
pR = plot(ax2, NaN, NaN, '-', 'Color', RED,  'LineWidth', 2.4);
hold(ax2, 'off');
lg = legend(ax2, [pB pR hMean], {'full reach, caps = measured extremes', ...
        sprintf('middle %d%% of trials, 5th-95th pct (dashes)', round((1-2*TAIL)*100)), ...
        'trial-weighted mean'}, 'Location', 'northwest', 'FontSize', 11);
lg.Box = 'on';  ax2.Toolbar.Visible = 'off';

for ax = [ax1 ax2]
    style_axes(ax, STY);
end
set([fig1 fig2], 'Color', STY.figbg);
style_legend(lg, STY);
ax1.YAxis.TickLabelFormat = '%.2f';
ax2.YAxis.TickLabelFormat = '%.2f';
drawnow;

png1 = fullfile(scriptDir, 'MUX_Add_poster_warmup.png');
png2 = fullfile(scriptDir, 'MUX_Add_poster_faults.png');
exportgraphics(fig1, png1, 'Resolution', 300, 'BackgroundColor', 'current');
exportgraphics(fig2, png2, 'Resolution', 300, 'BackgroundColor', 'current');
for f = {png1, png2}, pad_png(f{1}, MARGIN_PX); end
fprintf('\nFigure 1 saved to %s\n', png1);
fprintf('Figure 2 saved to %s\n', png2);

% Apply the poster palette and font sizes to one axes. Title and axis labels do NOT follow
% XColor/YColor in MATLAB, so each is set by hand; miss them and they stay light-on-grey and
% vanish. Font sizes are set here too, overriding the per-call values.
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

% Legends do not inherit the axes colours; left alone they stay dark-on-dark and vanish on grey.
function style_legend(lg, STY)
    set(lg, 'Color', STY.legbg, 'TextColor', STY.fg, 'EdgeColor', STY.fg);
end

% exportgraphics crops flush to the content, so titles and tick labels sit hard against the PNG
% edge and look cramped once the image is dropped onto a poster. This adds a quiet border in
% whatever colour the image already has at its corner, so it follows the palette automatically.
% Done afterwards rather than via exportgraphics' Padding option, which only offers 'tight' or
% 'figure' -- neither is a controllable number of pixels.
function pad_png(file, px)
    img = imread(file);
    if size(img, 3) == 1, img = repmat(img, 1, 1, 3); end
    bg = img(1, 1, :);                       % corner pixel = the exported background colour
    [h, w, ~] = size(img);
    out = repmat(bg, h + 2*px, w + 2*px, 1);
    out(px+1:px+h, px+1:px+w, :) = img;
    imwrite(out, file);
end
