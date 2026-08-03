%% plot_AND_Mul_poster.m                              ECE SPARK 2026 -- AND multiplier
%
% Reads the two CSVs written by AND_Mul_poster.cpp in this folder and saves the figure here as
% a PNG. Nothing outside this folder is touched.
%
%     AND_Mul_poster_warmup.csv  -> plot 1, zero-fault early termination
%     AND_Mul_poster_faults.csv  -> plot 2, fault response 0..36 flipped bits
%
% ------------------------------------------------------------------------------------------
% PLOT 1 -- ZERO-FAULT EARLY TERMINATION
%   Mean absolute error against stream length, over the zero-error trials: pairs of 256-bit 0.5
%   streams that the AND gate puts at exactly 0.25. Truncating the stream early is the only
%   thing costing accuracy here, since there is no fault. Same shape as the ANDmul convergence
%   plot in keepdata, drawn from this run's own trials.
%
% PLOT 2 -- FAULT RESPONSE
%   X  bits flipped, at 0, 1, 2, 4, 8, 16 and 32, on a LOG axis. Zero is placed at x = 0.5 and
%      labelled "0", since log(0) does not exist -- its position is cosmetic, its height real.
%   Y  the fraction of ones the finished answer carries. 0.25 sits in the middle with a red
%      dashed line across it -- that is the intended result.
%   BLUE STICK   the full reach of the error in both directions. Every cap is a real measured
%                run, and the reach is complete rather than sampled: at 36 flips the stick runs
%                0.109375 to 0.390625, which is every flip landing where it removes a one, and
%                every flip landing where it adds one. Those two tips carry probabilities near
%                1e-21, so 20,000 random trials would never land on them -- they are found by
%                constructing the layout and running it, which is why the caps are exact.
%   RED DASHES   the EQUAL-TAILED middle 90%: the 5th and 95th percentiles, with 5% of trials
%                cut from each end. Nothing forces it symmetric -- it comes out symmetric here
%                because this distribution IS symmetric, which is a measurement rather than an
%                assumption. Drawn only where it is a proper subset of the blue: at 0, 1 and 2
%                flips there are only 1, 3 and 5 reachable outcomes, so 90% cannot be separated
%                from 100% and the dashes would sit on the caps.
%   BLUE DOT     the trial-weighted mean.
%
% ------------------------------------------------------------------------------------------
% WHERE THE NUMBERS COME FROM. Every value plotted here was measured by running real 256-bit
% streams through the gate -- roughly 750,000 simulations, reported in the figure 2 subtitle and
% carried per flip level in the LevelRealTrials column. Nothing is evaluated from a formula.
%
% The combinatorics supply only the WEIGHT on each measured outcome. The literal sweep -- every
% 256-bit form of 0.5, paired with every other, times every layout of up to 36 flips across 512
% positions -- is about 1e300 trials and cannot be run. But an AND gate cannot see WHERE a bit is,
% so all of those collapse onto a handful of outcomes per flip level whose trial counts are exact.
%
% THAT WEIGHTING IS ITSELF CHECKED AGAINST RANDOM TRIALS. The CSV carries MC_Probability beside
% Probability: 20,000 genuinely random trials per flip level, drawn with no steering toward any
% outcome. This script prints the worst gap between the two, which runs under 0.007 -- about what
% sampling error alone gives at 20,000 draws. See the .cpp header for the full argument.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
warmCsv  = fullfile(scriptDir, 'AND_Mul_poster_warmup.csv');
faultCsv = fullfile(scriptDir, 'AND_Mul_poster_faults.csv');
for fchk = {warmCsv, faultCsv}
    if ~isfile(fchk{1})
        error(['CRITICAL: %s not found.\n' ...
               'Run:  stochastic_computer.exe --gtest_filter=AndMulPoster.*'], fchk{1});
    end
end

TARGET     = 0.25;
SHOW_FLIPS = [0 1 2 4 8 16 32];  % the x positions requested, on a log axis
TAIL       = 0.05;               % cut this much off EACH end -> red band is the middle 90%
MARGIN_PX  = 40;                 % border added around every exported PNG, at 300 DPI ~ 3.4 mm

%% POSTER PALETTE -- shared by all four scripts in this folder set.
%
% MEDIUM GREY, NOT THE MATLAB DEFAULT DARK. Mounted on a light poster, a near-black panel reads
% as a hole punched in the page and pulls the eye away from the text around it. A medium grey
% sits quietly next to white card while still separating the plot from the background.
%
% Everything else follows from that choice: on mid-grey, white text loses contrast, so the axes,
% labels and titles go NEAR-BLACK, and the grid goes light so it reads as a soft guide rather
% than as data. The saturated data colours (blue / red / orange) are unchanged -- they carry
% enough contrast against grey and keeping them means the figures still match each other.
% Font sizes are set HERE and applied in style_axes, which overrides whatever each xlabel/title
% call asked for. One place to tune, and it keeps the four scripts identical.
STY = struct( ...
    'bg',    [0.80 0.80 0.80], ...   % the PLOT PANEL only -- light grey
    'figbg', [1.00 1.00 1.00], ...   % figure surround + exported margin: WHITE
    'fg',    [0.10 0.10 0.10], ...   % axes, ticks, labels, titles
    'grid',  [1.00 1.00 1.00], ...   % divider lines, WHITE on the grey panel
    'legbg', [0.93 0.93 0.93], ...   % legend panel, near-white so it lifts off the light panel
    'tick',  14, ...                 % tick labels
    'label', 15, ...                 % axis labels
    'title', 15);                    % titles

%% ---------------------------------------------------------------------------------------
%% Load
%% ---------------------------------------------------------------------------------------
W = readtable(warmCsv,  'Delimiter', ',');
F = readtable(faultCsv, 'Delimiter', ',');

% Real gate runs behind the fault figure. LevelRealTrials is constant within a flip level and
% repeated on each of that level's rows, so take one value per level and add them up.
[lvl, iFirst] = unique(F.BitsFlipped);
realRuns = sum(F.LevelRealTrials(iFirst));
mcPerLevel = F.MC_Trials(1);
% Worst disagreement between the enumerated share of trials and what the random trials measured.
mcGap = 0;
for k = 1:numel(lvl)
    sel = F.BitsFlipped == lvl(k);
    pe = F.Probability(sel); pe = pe / sum(pe);
    mcGap = max(mcGap, max(abs(pe - F.MC_Probability(sel))));
end

fprintf('Warm-up : %d truncation lengths, exact over 10^%.0f arrangements\n', ...
        height(W), W.Log10Arrangements(1));
fprintf('Faults  : %d outcome rows across %d flip levels\n', ...
        height(F), numel(lvl));
fprintf('Trials  : %s real gate simulations (%s random per flip level)\n', ...
        commafy(realRuns), commafy(mcPerLevel));
fprintf('Check   : worst enumerated-vs-random probability gap = %.4f\n', mcGap);

% 748214 -> "748,214". Trial counts on a poster are unreadable without the separators.
function s = commafy(n)
    s = fliplr(regexprep(fliplr(sprintf('%d', round(n))), '(\d{3})(?=\d)', '$1,'));
end

%% EQUAL-TAILED band: cut TAIL of the trials off each end, keep the middle, per flip level.
%
% The lower edge is the 5th percentile, the upper edge the 95th, so 45% of the trials sit between
% each edge and the median. The band is NOT forced to be symmetric.
%
% THIS REPLACED TWO EARLIER RULES, and the history is worth keeping because both were wrong in
% instructive ways.
%
%   (1) NARROWEST INTERVAL HOLDING 90%. Subtly wrong for this data. At 2 flips the centre three
%       outcomes hold 87.6% -- just under target -- and extending one step either way gives
%       93.79% (low) or 93.83% (high) for the SAME span. An effectively exact tie, resolved by
%       whichever direction the scan reached first, which planted a fake leftward skew.
%
%   (2) SYMMETRIC ABOUT THE MEAN, i.e. the narrowest [mean - d, mean + d] holding 90%. That has
%       no tie to break, and for the AND gate it is harmless because this distribution really is
%       symmetric. But it BUILDS IN the symmetry rather than measuring it, so it cannot be
%       trusted to report asymmetry when asymmetry is there -- and in the uMUL twin it is: a
%       value-register hit puts an outcome at 0.000 while the high side only reaches 0.375, and
%       the symmetric rule dragged the high edge out to match the low one, describing a symmetry
%       the hardware does not have.
%
% Equal-tailed assumes nothing about shape. Written identically in both scripts so the two
% figures use one estimator rather than two. For the AND gate the numbers barely move -- which
% is itself the useful check that this distribution's symmetry is real and not an artifact of
% how the band was drawn.
function [lo, hi] = tail_band(vals, probs, tail)
    [vals, ord] = sort(vals);
    probs = probs(ord);
    probs = probs / sum(probs);
    lo = vals(find(cumsum(probs) >= tail, 1));                    % 5th percentile
    hi = vals(find(cumsum(probs, 'reverse') >= tail, 1, 'last')); % 95th percentile
end

flipLevels = unique(F.BitsFlipped);
nf = numel(flipLevels);
mn = zeros(nf,1); mx = zeros(nf,1); mu = zeros(nf,1);
bLo = zeros(nf,1); bHi = zeros(nf,1); cov = zeros(nf,1);
for k = 1:nf
    sel = F.BitsFlipped == flipLevels(k);
    v = F.OutputFraction(sel);
    p = F.Probability(sel);
    p = p / sum(p);
    mn(k) = min(v);  mx(k) = max(v);  mu(k) = sum(v .* p);
    [bLo(k), bHi(k)] = tail_band(v, p, TAIL);
    % True coverage of the drawn band, so the "90%" claim on the figure is checked rather than
    % asserted. Discrete outcomes mean the edges land ON a value, so coverage runs a little over.
    cov(k) = sum(p(v >= bLo(k) & v <= bHi(k)));
end

fprintf('\n  flips     mean      min       max      5th pct       95th pct    coverage\n');
for k = 1:nf
    if ismember(flipLevels(k), SHOW_FLIPS)
        fprintf('  %5d  %8.6f  %8.6f  %8.6f  %11.6f  %13.6f  %9.4f\n', ...
                flipLevels(k), mu(k), mn(k), mx(k), bLo(k), bHi(k), cov(k));
    end
end

%% ---------------------------------------------------------------------------------------
%% FIGURE 1 -- zero-fault early termination            -> AND_Mul_poster_warmup.png
%%
%% The two plots are separate figures and separate PNGs. They answer different questions on
%% different x axes -- stream length versus fault intensity -- and side by side one of them
%% always ends up squeezed. Apart is also what a poster wants: each can be placed and scaled
%% on its own.
%% ---------------------------------------------------------------------------------------
fig1 = figure('Name', 'AND Multiplier -- Zero-Fault Warm-Up', ...
              'Units', 'Normalized', 'Position', [0.10, 0.14, 0.52, 0.70]);

ax1 = axes(fig1);
plot(ax1, W.N, W.MeanAbsError, '-o', 'LineWidth', 2.4, 'MarkerSize', 8, ...
     'Color', [0.10 0.28 0.85], 'MarkerFaceColor', [0.10 0.28 0.85]);
set(ax1, 'XScale', 'log');
grid(ax1, 'on');
% A log axis turns MINOR grid lines on by default -- dotted verticals at every 2x, 3x, 4x
% subdivision -- which reads as background noise behind the data. Majors only.
ax1.XMinorGrid = 'off';  ax1.YMinorGrid = 'off';
ax1.XMinorTick = 'off';  ax1.YMinorTick = 'off';
xticks(ax1, W.N); xticklabels(ax1, string(W.N));
xlim(ax1, [min(W.N)*0.85, max(W.N)*1.15]);
xlabel(ax1, 'Bit Stream Length  (early termination point) [log]', 'FontSize', 13);
ylabel(ax1, 'Mean Absolute Error', 'FontSize', 13);
% TITLE IS THE OPERATION AND THE TRIAL COUNT, NOTHING ELSE. The poster already labels each panel
% above it ("AND Mul Fault Behavior" and so on), so repeating the unit name inside the image is
% duplication. What the image has to carry is the arithmetic it ran and how much of it it ran.
%
% NOTE this panel quotes ARRANGEMENTS, not realRuns. realRuns comes from the fault CSV and counts
% the fault campaign only; the warm-up is a separate sweep whose trial basis is the exhaustive
% arrangement count in its own file. Quoting the fault campaign's number here would be wrong.
title(ax1, sprintf('0.5 \\times 0.5 = 0.25   |   exhaustive over all 10^{%.0f} arrangements', ...
                   W.Log10Arrangements(1)), 'FontSize', 15);
ax1.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% FIGURE 2 -- fault response                          -> AND_Mul_poster_faults.png
%% ---------------------------------------------------------------------------------------
fig2 = figure('Name', 'AND Multiplier -- Fault Response', ...
              'Units', 'Normalized', 'Position', [0.14, 0.10, 0.56, 0.72]);
ax2 = axes(fig2);
hold(ax2, 'on');

% ONLY the requested levels are drawn. An earlier version filled in every level 0..36 as faded
% sticks to close the gap between 16 and 32; on a log axis that gap disappears on its own, and
% the infill just read as clutter.
%
% ZERO ON A LOG AXIS. log(0) is undefined, so the 0-flip stick is placed at x = 0.5 -- one octave
% left of 1, where the next tick would fall -- and labelled "0". Its position is cosmetic; its
% height is the real measurement, and it is the only tick whose spacing is not literal.
xpos = SHOW_FLIPS;
xpos(xpos == 0) = 0.5;

% Whisker style rather than filled bars: a thin blue spine from the worst under-shoot to the
% worst over-shoot, a short cap on each end, and a red dash at each edge of the 90% interval.
% The filled version put a lot of ink on what is really four numbers per level, and at the low
% flip counts the red block swallowed the blue entirely.
%
% Cap widths are MULTIPLICATIVE so they render the same visual width at every tick on a log axis;
% an additive half-width would look tiny at 1 and enormous at 32. Kept small -- at the previous
% 1.14/1.22 each cap spanned ~5% of the axis and read as a bar rather than a tick.
capBlue = 1.055;
capRed  = 1.075;   % 90% dashes a touch wider so they read as separate marks
BLUE = [0.10 0.28 0.85];
RED  = [0.90 0.20 0.20];

for s = 1:numel(SHOW_FLIPS)
    k = find(flipLevels == SHOW_FLIPS(s), 1);
    if isempty(k), continue; end
    x = xpos(s);

    % spine
    plot(ax2, [x x], [mn(k) mx(k)], '-', 'Color', BLUE, 'LineWidth', 1.6);
    % caps at the true extremes
    plot(ax2, [x/capBlue x*capBlue], [mn(k) mn(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);
    plot(ax2, [x/capBlue x*capBlue], [mx(k) mx(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);

    % 90% dashes ONLY where they say something the blue caps do not.
    %
    % At 0, 1 and 2 flips there are just 1, 3 and 5 reachable outcomes, and 90% of the trials
    % cannot be carved out of a set that small -- the centred band swallows the whole range and
    % the red lands exactly on the blue. Drawing it there implies a distinction that does not
    % exist and makes the low end look broken. Suppressed instead, and noted in the caption.
    if bLo(k) > mn(k) || bHi(k) < mx(k)
        plot(ax2, [x/capRed x*capRed], [bLo(k) bLo(k)], '-', 'Color', RED, 'LineWidth', 2.4);
        plot(ax2, [x/capRed x*capRed], [bHi(k) bHi(k)], '-', 'Color', RED, 'LineWidth', 2.4);
    end
end

selShow = ismember(flipLevels, SHOW_FLIPS);
hMean = plot(ax2, xpos, mu(selShow), 'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', BLUE, 'MarkerEdgeColor', BLUE, 'LineWidth', 1.0);
% Label on the RIGHT: on the left it sat on top of the 0- and 1-flip sticks.
yline(ax2, TARGET, 'r--', 'LineWidth', 2.2, ...
      'Label', 'intended result = 0.25', 'LabelHorizontalAlignment', 'right', ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 11);
hold(ax2, 'off');

set(ax2, 'XScale', 'log');
grid(ax2, 'on');
ax2.XMinorGrid = 'off';  ax2.YMinorGrid = 'off';
ax2.XMinorTick = 'off';  ax2.YMinorTick = 'off';
% Right margin widened from 1.45x to 2.6x purely to park the "intended result = 0.25" label
% clear of the 32-flip whisker; at 1.45x the text sat on top of that stick's red dashes.
xlim(ax2, [0.38 max(SHOW_FLIPS) * 2.6]);
% Y bounded to 0.10 - 0.40, i.e. 0.25 +/- 0.15. This is SHARED with uMUL's operand-intact figure
% so the two can be mounted side by side and compared by eye without a scale trap. It also frames
% the data tightly: the widest plotted level (32 flips) spans 0.125 to 0.375, and even the full
% f = 36 extremes at 0.109375 / 0.390625 sit inside it.
ylim(ax2, [0.10 0.40]);
xticks(ax2, xpos);
xticklabels(ax2, string(SHOW_FLIPS));
xlabel(ax2, 'Bits Flipped  (across both 256-bit streams) [log]', 'FontSize', 13);
ylabel(ax2, 'Final Answer', 'FontSize', 15);
% Same single-line title as ax1: operation, then how many real gate runs stand behind it.
% commafy(realRuns), not a hardcoded number -- the count moves if the campaign is re-run.
title(ax2, sprintf('0.5 \\times 0.5 = 0.25   |   %s real gate simulations', ...
                   commafy(realRuns)), 'FontSize', 15);

% Legend proxies. hold(ax2,'off') ran above, and a bare plot() after hold off CLEARS the axes --
% with patch() proxies that went unnoticed, but as line objects it wiped the whole panel and left
% legend() with dead handles. Re-enable hold before drawing them.
hold(ax2, 'on');
pBlue = plot(ax2, NaN, NaN, '-', 'Color', BLUE, 'LineWidth', 1.8);
pRed  = plot(ax2, NaN, NaN, '-', 'Color', RED,  'LineWidth', 2.4);
hold(ax2, 'off');
lg = legend(ax2, [pBlue pRed hMean], ...
       {'full error reach, caps = measured extremes', ...
        sprintf('middle %d%% of trials, 5th-95th pct (dashes)', round((1-2*TAIL)*100)), ...
        'trial-weighted mean'}, 'Location', 'northwest', 'FontSize', 13);
lg.Box = 'on';
ax2.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% Light touch-ups, then export each figure on its own
%% ---------------------------------------------------------------------------------------
for ax = [ax1 ax2]
    style_axes(ax, STY);
end
set([fig1 fig2], 'Color', STY.figbg);
% The legend has to be restyled explicitly: it does not inherit the axes colours, so left alone
% it stays dark-on-dark from the default theme and becomes unreadable on grey.
style_legend(lg, STY);
ax1.YAxis.TickLabelFormat = '%.2f';
ax2.YAxis.TickLabelFormat = '%.2f';
drawnow;

png1 = fullfile(scriptDir, 'AND_Mul_poster_warmup.png');
png2 = fullfile(scriptDir, 'AND_Mul_poster_faults.png');
exportgraphics(fig1, png1, 'Resolution', 300, 'BackgroundColor', 'current');
exportgraphics(fig2, png2, 'Resolution', 300, 'BackgroundColor', 'current');
pad_png(png1, MARGIN_PX);
pad_png(png2, MARGIN_PX);
fprintf('\nFigure 1 saved to %s\n', png1);
fprintf('Figure 2 saved to %s\n', png2);

% exportgraphics crops flush to the content, so titles and tick labels sit hard against the PNG
% edge and look cramped the moment the image is dropped onto a poster next to anything else.
% This re-opens each export and adds a quiet border in whatever colour the image already has at
% its corner, so it works with the dark default theme and with a white one without being told
% which. Done afterwards rather than via exportgraphics' own Padding option because that one
% only offers 'tight' or 'figure' -- neither is a controllable number of pixels.
% Apply the poster palette and font sizes to one axes. Title and axis labels do NOT follow
% XColor/YColor in MATLAB, so each is set by hand; miss them and they stay light-on-grey and
% vanish. Font sizes are set here too, overriding the per-call values, so there is exactly one
% place to tune them.
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

function style_legend(lg, STY)
    set(lg, 'Color', STY.legbg, 'TextColor', STY.fg, 'EdgeColor', STY.fg);
end

function pad_png(file, px)
    img = imread(file);
    if size(img, 3) == 1, img = repmat(img, 1, 1, 3); end
    bg = img(1, 1, :);                       % corner pixel = the exported background colour
    [h, w, ~] = size(img);
    out = repmat(bg, h + 2*px, w + 2*px, 1);
    out(px+1:px+h, px+1:px+w, :) = img;
    imwrite(out, file);
end
