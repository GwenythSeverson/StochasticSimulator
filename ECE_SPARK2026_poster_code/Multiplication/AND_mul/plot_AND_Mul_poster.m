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
%   RED DASHES   the CENTRED interval holding at least 90% of trials -- a symmetric window grown
%                outward from the mean. Drawn only where it is a proper subset of the blue: at
%                0, 1 and 2 flips there are only 1, 3 and 5 reachable outcomes, so 90% cannot be
%                separated from 100% and the dashes would sit on the caps.
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
BAND_MASS  = 0.90;               % centred red band covers at least this much of the trials

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

%% CENTRED interval containing at least BAND_MASS of the trials, per flip level.
%
% Grows a symmetric window outward from the outcome nearest the mean until it holds 90% of the
% trials. Symmetric by construction, which is the point: the band is meant to say "this is where
% the middle 90% of answers land", and a lopsided one invites the reader to see a bias that is
% not there.
%
% THIS REPLACED A "NARROWEST INTERVAL" RULE, which was subtly wrong for this data. At 2 flips the
% centre three outcomes hold 87.6% -- just under target -- and extending one step either way gives
% 93.79% (low side) or 93.83% (high side) for the SAME span. An exact tie, resolved by whichever
% the scan happened to reach first, which planted a fake leftward skew on the figure. Growing
% symmetrically has no tie to break.
% 748214 -> "748,214". Trial counts on a poster are unreadable without the separators.
function s = commafy(n)
    s = fliplr(regexprep(fliplr(sprintf('%d', round(n))), '(\d{3})(?=\d)', '$1,'));
end

function [lo, hi] = centred_band(vals, probs, mass)
    [vals, ord] = sort(vals);
    probs = probs(ord);
    probs = probs / sum(probs);
    n = numel(vals);
    mu = sum(vals .* probs);
    [~, c] = min(abs(vals - mu));      % start at the outcome closest to the mean
    r = 0;
    while true
        i = max(1, c - r);
        j = min(n, c + r);
        if sum(probs(i:j)) >= mass || (i == 1 && j == n), break; end
        r = r + 1;
    end
    lo = vals(max(1, c - r));
    hi = vals(min(n, c + r));
end

flipLevels = unique(F.BitsFlipped);
nf = numel(flipLevels);
mn = zeros(nf,1); mx = zeros(nf,1); mu = zeros(nf,1);
bLo = zeros(nf,1); bHi = zeros(nf,1);
for k = 1:nf
    sel = F.BitsFlipped == flipLevels(k);
    v = F.OutputFraction(sel);
    p = F.Probability(sel);
    p = p / sum(p);
    mn(k) = min(v);  mx(k) = max(v);  mu(k) = sum(v .* p);
    [bLo(k), bHi(k)] = centred_band(v, p, BAND_MASS);
end

fprintf('\n  flips     mean      min       max     90%% band low   90%% band high\n');
for k = 1:nf
    if ismember(flipLevels(k), SHOW_FLIPS)
        fprintf('  %5d  %8.6f  %8.6f  %8.6f  %13.6f  %14.6f\n', ...
                flipLevels(k), mu(k), mn(k), mx(k), bLo(k), bHi(k));
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
     'Color', [0.15 0.30 0.75], 'MarkerFaceColor', [0.15 0.30 0.75]);
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
title(ax1, {'AND-Gate Multiplier -- Zero-Fault Warm-Up', ...
            sprintf('0.5 \\times 0.5 = 0.25   |   exhaustive over all 10^{%.0f} arrangements', ...
                    W.Log10Arrangements(1))}, 'FontSize', 14);
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
BLUE = [0.20 0.40 0.95];
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
ylim(ax2, [0 0.5]);
xticks(ax2, xpos);
xticklabels(ax2, string(SHOW_FLIPS));
xlabel(ax2, 'Bits Flipped  (across both 256-bit streams) [log]', 'FontSize', 13);
ylabel(ax2, 'Fraction of Ones in the Final Answer', 'FontSize', 13);
% Subtitle kept to two clauses. A third ("bounds are measured extremes") pushed the line past
% both edges of the axes at this font size; that claim lives in the legend entry instead.
title(ax2, {'AND-Gate Multiplier -- Fault Response', ...
            sprintf('0.5 \\times 0.5 = 0.25   |   %s real gate simulations', ...
                    commafy(realRuns))}, 'FontSize', 14);

% Legend proxies. hold(ax2,'off') ran above, and a bare plot() after hold off CLEARS the axes --
% with patch() proxies that went unnoticed, but as line objects it wiped the whole panel and left
% legend() with dead handles. Re-enable hold before drawing them.
hold(ax2, 'on');
pBlue = plot(ax2, NaN, NaN, '-', 'Color', BLUE, 'LineWidth', 1.8);
pRed  = plot(ax2, NaN, NaN, '-', 'Color', RED,  'LineWidth', 2.4);
hold(ax2, 'off');
lg = legend(ax2, [pBlue pRed hMean], ...
       {'full error reach, caps = measured extremes', ...
        sprintf('centred %d%% of trials (dashes)', round(BAND_MASS*100)), ...
        'trial-weighted mean'}, 'Location', 'northwest', 'FontSize', 11);
lg.Box = 'on';
ax2.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% Light touch-ups, then export each figure on its own
%% ---------------------------------------------------------------------------------------
% Default theme left alone -- forcing white made the panels heavier than they needed to be.
% Only the font size and tick format are nudged, so the axes keep whatever look MATLAB gives them.
for ax = [ax1 ax2]
    ax.FontSize = 12;
    ax.LineWidth = 1.0;
    ax.Box = 'on';
end
ax1.YAxis.TickLabelFormat = '%.2f';
ax2.YAxis.TickLabelFormat = '%.2f';
drawnow;

png1 = fullfile(scriptDir, 'AND_Mul_poster_warmup.png');
png2 = fullfile(scriptDir, 'AND_Mul_poster_faults.png');
exportgraphics(fig1, png1, 'Resolution', 300, 'BackgroundColor', 'current');
exportgraphics(fig2, png2, 'Resolution', 300, 'BackgroundColor', 'current');
fprintf('\nFigure 1 saved to %s\n', png1);
fprintf('Figure 2 saved to %s\n', png2);
