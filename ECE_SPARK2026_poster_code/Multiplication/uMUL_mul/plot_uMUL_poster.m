%% plot_uMUL_poster.m                                    ECE SPARK 2026 -- uMUL multiplier
%
% Twin of AND_mul/plot_AND_Mul_poster.m. Same axes, same marks, same conventions, so the two
% pairs of figures can be placed side by side on the poster and read against each other.
%
%     uMUL_poster_warmup.csv  -> plot 1, zero-fault early termination
%     uMUL_poster_faults.csv  -> plot 2, fault response 0..36 flipped bits
%
% ------------------------------------------------------------------------------------------
% PLOT 1 -- ZERO-FAULT EARLY TERMINATION
%   Mean absolute error against stream length, exhaustive over all C(256,128) = 5.8e75
%   arrangements of the in_0 stream. Truncating early is the only thing costing accuracy, since
%   there is no fault. Compare directly against the AND gate's plot 1: uMUL sits about 40%
%   lower from L = 8 onward.
%
%   THE CURVE IS NOT MONOTONE, AND THAT IS REAL. L = 4 (0.031) beats L = 8 (0.091). uMUL's
%   generator walks a Sobol sequence in Gray-code order, whose partial sums are exact at some
%   lengths and lopsided at others; at L = 4 every reachable prefix count lands on exactly
%   0.25, which is a genuine alignment rather than a fluke of sampling. From L = 8 onward the
%   curve is clean and monotone. See HOW_THIS_DATA_WAS_MADE.txt section 5.
%
% PLOTS 2 AND 3 -- FAULT RESPONSE, SPLIT BY WHETHER THE OPERAND REGISTER SURVIVED
%   Common axes. X is bits flipped over the 264-bit operand surface (256-bit in_0 stream +
%   8-bit value register) at 0, 1, 2, 4, 8, 16 and 32, on a LOG axis; zero is drawn at x = 0.5
%   and labelled "0", since log(0) does not exist -- its position is cosmetic, its height real.
%   Y is the fraction of ones the finished answer carries, against an intended 0.25.
%
%   PLOT 2, OPERAND INTACT   no register bit struck; only the in_0 stream took damage.
%       BLUE STICK  full reach, caps at the measured extremes
%       RED DASHES  the equal-tailed middle 90% (5th / 95th percentile), not forced symmetric
%       BLUE DOT    trial-weighted mean, pinned at 0.25 at every flip level
%       Vertical scale is 0.10-0.40, SHARED with the AND gate's fault figure so the two can be
%       mounted side by side. On plot 3's full 0-0.70 axis this panel would be a flat line.
%
%   PLOT 3, OPERAND DAMAGED  at least one of the 8 register bits struck.
%       ORANGE STICK  full reach, caps at the measured extremes
%       ORANGE DOT    trial-weighted mean -- also 0.25, which is the surprise; see the long note
%                     at that figure. Its MEDIAN is 0.27, so the mean is not "typical".
%       Vertical scale is the full 0-0.70.
%
%   WHY THEY ARE SEPARATE PLOTS. Merged into one distribution these two populations make a
%   MIXTURE, and a percentile over a mixture jumps discontinuously the moment the mixing weight
%   crosses the tail fraction. It did, at f = 14, by 0.164 -- an artifact of the statistic, not
%   of the hardware. Conditioned, each population's band edge moves by at most 1/256 per flip
%   level, which is the output quantisation floor. See HOW_THIS_DATA_WAS_MADE.txt section 6e.
%
%   READ THE PAIR AGAINST THE AND GATE. At ONE flip the AND gate spans 0.246 to 0.254 and has
%   no second population at all, because all 512 of its operand bits are worth the same 1/256.
%   uMUL's intact population matches it; its damaged population spans 0.000 to 0.375 from the
%   very first flip. uMUL fails rarely and catastrophically; the AND gate fails constantly and
%   negligibly. That contrast is the poster.
%
% PLOT 4 -- THE 8-BIT OPERAND REGISTER, EXPANDED
%   Where the damage in plot 3 comes from: each register bit's single-flip effect, against the
%   stream bit's +/- 1/256 for scale. Seven bits push the answer up, one takes it to zero.
%
% ------------------------------------------------------------------------------------------
% WHERE THE NUMBERS COME FROM. Every value plotted here was measured by clocking real 256-bit
% streams through a real UnaryMultiplier -- roughly 800,000 simulations, reported in the
% figure 2 subtitle. Nothing is evaluated from a formula.
%
% The combinatorics supply only the WEIGHT on each measured outcome. uMUL's output depends only
% on (value, ones-count of in_0) and never on where those ones sit, so the whole 5.8e75-strong
% arrangement space collapses onto a pair of integers. That weighting is itself checked against
% MC_Probability -- 20,000 genuinely random trials per flip level -- and this script prints the
% worst gap between the two. See the .cpp header and the .txt for the full argument.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
warmCsv  = fullfile(scriptDir, 'uMUL_poster_warmup.csv');
faultCsv = fullfile(scriptDir, 'uMUL_poster_faults.csv');
for fchk = {warmCsv, faultCsv}
    if ~isfile(fchk{1})
        error(['CRITICAL: %s not found.\n' ...
               'Run:  stochastic_computer.exe --gtest_filter=uMULPoster.*'], fchk{1});
    end
end

TARGET     = 0.25;
SHOW_FLIPS = [0 1 2 4 8 16 32];  % the x positions, on a log axis
TAIL       = 0.05;               % cut this much off EACH end -> red band is the middle 90%

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

fprintf('Warm-up : %d truncation lengths, exact over 10^%.0f in_0 arrangements\n', ...
        height(W), W.Log10Arrangements(1));
fprintf('Faults  : %d outcome rows across %d flip levels, 264-bit operand surface\n', ...
        height(F), numel(lvl));
fprintf('Trials  : %s real gate simulations (%s random per flip level)\n', ...
        commafy(realRuns), commafy(mcPerLevel));
fprintf('Check   : worst enumerated-vs-random probability gap = %.4f\n', mcGap);

% 748214 -> "748,214". Trial counts on a poster are unreadable without the separators.
function s = commafy(n)
    s = fliplr(regexprep(fliplr(sprintf('%d', round(n))), '(\d{3})(?=\d)', '$1,'));
end

%% CENTRED interval containing at least BAND_MASS of the trials, per flip level.
%
% EQUAL-TAILED: cut TAIL of the trials off each end and keep what is left. The lower edge is the
% 5th percentile, the upper edge the 95th, so 45% of the trials sit between each edge and the
% median. The band is NOT forced to be symmetric, and for uMUL it must not be.
%
% WHY THIS REPLACED A SYMMETRIC-ABOUT-THE-MEAN BAND. The previous rule grew one window
% [mean - d, mean + d] until it held 90%. That is only meaningful if the two sides of the
% distribution are equally reachable, and uMUL's are not: flipping value bit 7 zeroes the operand
% and puts an outcome at 0.000, while the far side only reaches 0.375 at one flip. Forcing the
% window symmetric made the low edge drag the high edge out with it, so the band reported reach
% on the high side that no comparable share of trials actually occupies -- it described a
% symmetry the hardware does not have. At 32 flips the symmetric rule gave 0.055 to 0.445;
% equal-tailed gives the real, lopsided answer.
%
% The zero outcome is not hidden by this, it is correctly placed: at one flip it is 1 site in
% 264 = 0.38% of trials, far inside the 5% tail, so it belongs outside a 90% band and inside the
% blue spine -- which is exactly where the figure now puts it.
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
bLo = zeros(nf,1); bHi = zeros(nf,1); pReg = zeros(nf,1); cov = zeros(nf,1);
% Conditional statistics, computed once here and used by BOTH figure 2 and figure 4.
%   c* = operand INTACT  (no register bit struck)
%   d* = operand DAMAGED (at least one struck).  NaN at f = 0, where that population is empty.
cLo = zeros(nf,1); cHi = zeros(nf,1); cMn = zeros(nf,1); cMx = zeros(nf,1); cMu = zeros(nf,1);
dLo = nan(nf,1);   dHi = nan(nf,1);   dMn = nan(nf,1);   dMx = nan(nf,1);   dMu = nan(nf,1);
pDam = zeros(nf,1);
for k = 1:nf
    sel = F.BitsFlipped == flipLevels(k);
    v = F.OutputFraction(sel);
    p = F.Probability(sel);
    p = p / sum(p);
    mn(k) = min(v);  mx(k) = max(v);  mu(k) = sum(v .* p);
    pReg(k) = sum(F.ProbRegisterHit(sel));   % share of trials whose value register was hit
    [bLo(k), bHi(k)] = tail_band(v, p, TAIL);
    % True coverage of the drawn band, so the "90%" claim on the figure is checked rather than
    % asserted. Discrete outcomes mean the edges land ON a value, so coverage runs a little over.
    cov(k) = sum(p(v >= bLo(k) & v <= bHi(k)));

    % ---- split by whether the 8-bit operand register survived --------------------------------
    % ProbRegisterClean/Hit are JOINT masses (share of ALL trials at this f), so each is
    % renormalised by its own total to become a conditional distribution.
    pc = F.ProbRegisterClean(sel);
    pd = F.ProbRegisterHit(sel);
    pDam(k) = sum(pd) / (sum(pc) + sum(pd));
    if sum(pc) > 0
        i = pc > 0;
        [cLo(k), cHi(k)] = tail_band(v(i), pc(i), TAIL);
        cMn(k) = min(v(i));  cMx(k) = max(v(i));  cMu(k) = sum(v(i) .* pc(i)) / sum(pc(i));
    end
    if sum(pd) > 0
        i = pd > 0;
        [dLo(k), dHi(k)] = tail_band(v(i), pd(i), TAIL);
        dMn(k) = min(v(i));  dMx(k) = max(v(i));  dMu(k) = sum(v(i) .* pd(i)) / sum(pd(i));
    end
end

fprintf('\n  flips     mean      min       max      5th pct       95th pct    coverage   P(reg hit)\n');
for k = 1:nf
    if ismember(flipLevels(k), SHOW_FLIPS)
        fprintf('  %5d  %8.6f  %8.6f  %8.6f  %11.6f  %13.6f  %9.4f  %10.4f\n', ...
                flipLevels(k), mu(k), mn(k), mx(k), bLo(k), bHi(k), cov(k), pReg(k));
    end
end

%% ---------------------------------------------------------------------------------------
%% FIGURE 1 -- zero-fault early termination            -> uMUL_poster_warmup.png
%% ---------------------------------------------------------------------------------------
fig1 = figure('Name', 'uMUL -- Zero-Fault Warm-Up', ...
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
title(ax1, {'uMUL -- Zero-Fault Warm-Up', ...
            sprintf('0.5 \\times 0.5 = 0.25   |   exhaustive over all 10^{%.0f} in_0 arrangements', ...
                    W.Log10Arrangements(1))}, 'FontSize', 14);
ax1.Toolbar.Visible = 'off';


%% ---------------------------------------------------------------------------------------
%% Shared whisker geometry, used by figures 2 and 3
%%
%% ZERO ON A LOG AXIS. log(0) is undefined, so the 0-flip stick is placed at x = 0.5 -- one octave
%% left of 1, where the next tick would fall -- and labelled "0". Its position is cosmetic; its
%% height is the real measurement, and it is the only tick whose spacing is not literal.
%%
%% Cap widths are MULTIPLICATIVE so they render the same visual width at every tick on a log axis;
%% an additive half-width would look tiny at 1 and enormous at 32.
%% ---------------------------------------------------------------------------------------
xpos = SHOW_FLIPS;
xpos(xpos == 0) = 0.5;
capBlue = 1.055;
capRed  = 1.075;   % 90% dashes a touch wider so they read as separate marks
BLUE   = [0.10 0.28 0.85];
RED    = [0.90 0.20 0.20];
ORANGE = [0.95 0.55 0.10];

%% ---------------------------------------------------------------------------------------
%% FIGURE 2 -- FAULT RESPONSE, OPERAND INTACT        -> uMUL_poster_faults_intact.png
%%
%% Trials where NONE of the 8 value-register bits was struck -- only the 256-bit in_0 stream took
%% damage. This is the ordinary case: at one flip it is 97% of all trials.
%%
%% NOTE THE EXPANDED VERTICAL SCALE. This panel runs 0.15 to 0.35, not 0 to 0.70 like its damaged
%% twin. Plotted on the damaged figure's axis this whole curve would be a flat line, which is
%% itself the point -- but it would show nothing. The two figures are deliberately NOT on a shared
%% y axis, and both say so in their subtitles.
%% ---------------------------------------------------------------------------------------
fig2 = figure('Name', 'uMUL -- Fault Response, Operand Intact', ...
              'Units', 'Normalized', 'Position', [0.10, 0.12, 0.54, 0.72]);
ax2 = axes(fig2);
hold(ax2, 'on');

for s = 1:numel(SHOW_FLIPS)
    k = find(flipLevels == SHOW_FLIPS(s), 1);
    if isempty(k), continue; end
    x = xpos(s);
    % full reach, caps at the measured extremes
    plot(ax2, [x x], [cMn(k) cMx(k)], '-', 'Color', BLUE, 'LineWidth', 1.6);
    plot(ax2, [x/capBlue x*capBlue], [cMn(k) cMn(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);
    plot(ax2, [x/capBlue x*capBlue], [cMx(k) cMx(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);
    % middle 90%, suppressed where it is the whole reach (too few outcomes to carve 90% from)
    if cLo(k) > cMn(k) || cHi(k) < cMx(k)
        plot(ax2, [x/capRed x*capRed], [cLo(k) cLo(k)], '-', 'Color', RED, 'LineWidth', 2.4);
        plot(ax2, [x/capRed x*capRed], [cHi(k) cHi(k)], '-', 'Color', RED, 'LineWidth', 2.4);
    end
end
selShow = ismember(flipLevels, SHOW_FLIPS);
hMean = plot(ax2, xpos, cMu(selShow), 'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', BLUE, 'MarkerEdgeColor', BLUE, 'LineWidth', 1.0);
yline(ax2, TARGET, 'r--', 'LineWidth', 2.2, 'Label', 'intended result = 0.25', ...
      'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'FontSize', 11);
hold(ax2, 'off');

set(ax2, 'XScale', 'log');
grid(ax2, 'on');
ax2.XMinorGrid = 'off';  ax2.YMinorGrid = 'off';
ax2.XMinorTick = 'off';  ax2.YMinorTick = 'off';
xlim(ax2, [0.38 max(SHOW_FLIPS) * 2.6]);
% Y bounded to 0.10 - 0.40, i.e. 0.25 +/- 0.15, SHARED with the AND gate's fault-response figure.
% That is the whole point of this panel: mounted beside AND_Mul_poster_faults.png on the same
% scale, it shows that uMUL's ordinary case -- 97% of single upsets -- behaves like the AND gate.
% The comparison only works if the axes match, so do not retune one without the other.
ylim(ax2, [0.10 0.40]);
xticks(ax2, xpos);  xticklabels(ax2, string(SHOW_FLIPS));
xlabel(ax2, 'Bits Flipped  (over the 264-bit operand surface) [log]', 'FontSize', 13);
ylabel(ax2, 'Fraction of Ones in the Final Answer', 'FontSize', 13);
title(ax2, {'uMUL -- Fault Response, OPERAND INTACT', ...
            'no register bit struck   |   scale 0.10 to 0.40, shared with the AND gate'}, ...
      'FontSize', 13);

hold(ax2, 'on');
pBlue = plot(ax2, NaN, NaN, '-', 'Color', BLUE, 'LineWidth', 1.8);
pRed  = plot(ax2, NaN, NaN, '-', 'Color', RED,  'LineWidth', 2.4);
hold(ax2, 'off');
lg2 = legend(ax2, [pBlue pRed hMean], ...
       {'full reach, caps = measured extremes', ...
        sprintf('middle %d%% of trials, 5th-95th pct (dashes)', round((1-2*TAIL)*100)), ...
        'trial-weighted mean'}, 'Location', 'northwest', 'FontSize', 11);
lg2.Box = 'on';
ax2.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% FIGURE 3 -- FAULT RESPONSE, OPERAND DAMAGED       -> uMUL_poster_faults_damaged.png
%%
%% Trials where AT LEAST ONE of the 8 value-register bits was struck. Broad from ONE flip onward,
%% and the floor is 0.000 the whole way across -- when the operand dies it can die completely, and
%% that is as true at f = 1 as at f = 36.
%%
%% THE MEAN SITS ON 0.25, WHICH IS THE SURPRISE. Measured 0.250000 at f = 1, drifting only to
%% 0.249974 by f = 36: even the corrupted-operand trials are unbiased on average. It is exact at
%% f = 1 and checkable by hand -- one struck register bit gives eight equally likely operands
%% 0, 192, 160, 144, 136, 132, 130, 129, whose answers are 0, 96, 80, 72, 68, 66, 65, 65 out of
%% 256. Those sum to 512, so the mean is 64/256 = 0.25 exactly. The single catastrophic zero is
%% worth precisely as much as the seven upward pushes combined.
%%
%% DO NOT READ THE MEAN AS TYPICAL. This population's MEDIAN is 0.266 to 0.273, well above its
%% mean, because most strikes push the answer up and a few pull it all the way to zero. Mean 0.25,
%% median 0.27, floor 0.000 is the fat-tail signature quantified in section 6c of the .txt, and it
%% is exactly why an error bar is the wrong summary for uMUL.
%% ---------------------------------------------------------------------------------------
fig3 = figure('Name', 'uMUL -- Fault Response, Operand Damaged', ...
              'Units', 'Normalized', 'Position', [0.14, 0.10, 0.54, 0.72]);
ax3 = axes(fig3);
hold(ax3, 'on');

for s = 1:numel(SHOW_FLIPS)
    k = find(flipLevels == SHOW_FLIPS(s), 1);
    if isempty(k) || isnan(dMn(k)), continue; end   % f = 0 has no damaged population
    x = xpos(s);
    plot(ax3, [x x], [dMn(k) dMx(k)], '-', 'Color', ORANGE, 'LineWidth', 1.6);
    plot(ax3, [x/capBlue x*capBlue], [dMn(k) dMn(k)], '-', 'Color', ORANGE, 'LineWidth', 1.8);
    plot(ax3, [x/capBlue x*capBlue], [dMx(k) dMx(k)], '-', 'Color', ORANGE, 'LineWidth', 1.8);
end
hMeanD = plot(ax3, xpos, dMu(selShow), 'o', 'MarkerSize', 7, ...
              'MarkerFaceColor', ORANGE, 'MarkerEdgeColor', ORANGE, 'LineWidth', 1.0);
yline(ax3, TARGET, 'r--', 'LineWidth', 2.2, 'Label', 'intended result = 0.25', ...
      'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'FontSize', 11);
hold(ax3, 'off');

set(ax3, 'XScale', 'log');
grid(ax3, 'on');
ax3.XMinorGrid = 'off';  ax3.YMinorGrid = 'off';
ax3.XMinorTick = 'off';  ax3.YMinorTick = 'off';
xlim(ax3, [0.38 max(SHOW_FLIPS) * 2.6]);
ylim(ax3, [0 0.70]);
xticks(ax3, xpos);  xticklabels(ax3, string(SHOW_FLIPS));
xlabel(ax3, 'Bits Flipped  (over the 264-bit operand surface) [log]', 'FontSize', 13);
ylabel(ax3, 'Fraction of Ones in the Final Answer', 'FontSize', 13);
title(ax3, {'uMUL -- Fault Response, OPERAND DAMAGED', ...
            'at least one of the 8 register bits struck   |   full vertical scale, 0 to 0.70'}, ...
      'FontSize', 13);

hold(ax3, 'on');
pOr = plot(ax3, NaN, NaN, '-', 'Color', ORANGE, 'LineWidth', 1.8);
hold(ax3, 'off');
lg3 = legend(ax3, [pOr hMeanD], ...
       {'full reach, caps = measured extremes', ...
        'trial-weighted mean (the median is higher -- see note)'}, ...
       'Location', 'northwest', 'FontSize', 11);
lg3.Box = 'on';
ax3.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% FIGURE 4 -- THE 8-BIT OPERAND REGISTER, EXPANDED  -> uMUL_poster_register.png
%%
%% Where the damage comes from. Each of the 8 value-register bits as the change it makes to the
%% answer when it alone is flipped, against the stream bit's +/- 1/256 for scale.
%%
%% Bits 0 through 6 rise geometrically -- +0.004, +0.004, +0.008, +0.016, +0.031, +0.062, +0.125
%% -- because setting bit b adds 2^b to an operand of 128. Bit 7 goes the other way and goes all
%% the way: value = 128 has exactly ONE bit set, so flipping bit 7 clears the operand outright,
%% the strict ">" comparator can never fire again, and the answer is a hard zero. One upset there
%% is worth 64 stream upsets.
%%
%% SEVEN OF THE EIGHT BITS PUSH UP AND ONLY ONE PUSHES DOWN, but that one outweighs the other
%% seven combined -- which is why figure 3's damaged population has a mean of exactly 0.25 and a
%% median well above it. Every bar is a measured gate run, from uMUL_poster_sensitivity.csv.
%% ---------------------------------------------------------------------------------------
S = readtable(fullfile(scriptDir, 'uMUL_poster_sensitivity.csv'), 'Delimiter', ',');

fig4 = figure('Name', 'uMUL -- Operand Register Sensitivity', ...
              'Units', 'Normalized', 'Position', [0.18, 0.16, 0.52, 0.64]);
ax4 = axes(fig4);

regDelta = zeros(8,1);
for b = 0:7
    r = find(string(S.Site) == sprintf("value_bit_%d", b), 1);
    regDelta(b+1) = S.Delta(r);
end
streamDelta = max(abs(S.Delta(startsWith(string(S.Site), "stream"))));

hold(ax4, 'on');
bar(ax4, 0:7, regDelta, 0.72, 'FaceColor', RED, 'EdgeColor', STY.fg, 'LineWidth', 0.6);
% The stream bit's damage, for scale. On this axis it is almost invisible, which is the message.
yline(ax4, +streamDelta, '--', 'Color', BLUE, 'LineWidth', 1.4);
yline(ax4, -streamDelta, '--', 'Color', BLUE, 'LineWidth', 1.4, ...
      'Label', 'any one of the 256 stream bits = \pm1/256', 'LabelHorizontalAlignment', 'left', ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 10, 'Color', BLUE);
yline(ax4, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0);
for b = 0:7
    d = regDelta(b+1);
    if b == 7
        % BELOW the bar, not above: above puts red text inside a red bar. Right-aligned to the
        % axis edge, since centring on the last bar runs the text off the panel.
        text(ax4, 7.62, d - 0.020, 'operand \rightarrow 0', 'HorizontalAlignment', 'right', ...
             'FontSize', 11, 'Color', RED, 'FontWeight', 'bold');
    else
        text(ax4, b, d + 0.010, sprintf('%+.3f', d), 'HorizontalAlignment', 'center', ...
             'FontSize', 12, 'Color', STY.fg);
    end
end
hold(ax4, 'off');
xlim(ax4, [-0.7 7.7]);
ylim(ax4, [-0.30 0.20]);
xticks(ax4, 0:7);
grid(ax4, 'on');  ax4.XMinorGrid = 'off';  ax4.YMinorGrid = 'off';
xlabel(ax4, 'Value-register bit   (value = 128, so only bit 7 is set)', 'FontSize', 13);
ylabel(ax4, 'Change in the Answer', 'FontSize', 13);
title(ax4, {'uMUL -- The 8-Bit Operand Register, Expanded', ...
            'ONE upset in bit 7 is worth 64 stream upsets'}, 'FontSize', 13);
ax4.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% Light touch-ups, then export each figure on its own
%% ---------------------------------------------------------------------------------------
% Same palette and font sizes as the AND-gate script, so all four PNGs look like one set.
for ax = [ax1 ax2 ax3 ax4]
    style_axes(ax, STY);
end
set([fig1 fig2 fig3 fig4], 'Color', STY.figbg);
style_legend(lg2, STY);
style_legend(lg3, STY);
ax1.YAxis.TickLabelFormat = '%.2f';
ax2.YAxis.TickLabelFormat = '%.2f';
ax3.YAxis.TickLabelFormat = '%.2f';
drawnow;

png1 = fullfile(scriptDir, 'uMUL_poster_warmup.png');
png2 = fullfile(scriptDir, 'uMUL_poster_faults_intact.png');
png3 = fullfile(scriptDir, 'uMUL_poster_faults_damaged.png');
png4 = fullfile(scriptDir, 'uMUL_poster_register.png');
exportgraphics(fig1, png1, 'Resolution', 300, 'BackgroundColor', 'current');
exportgraphics(fig2, png2, 'Resolution', 300, 'BackgroundColor', 'current');
exportgraphics(fig3, png3, 'Resolution', 300, 'BackgroundColor', 'current');
exportgraphics(fig4, png4, 'Resolution', 300, 'BackgroundColor', 'current');
for f = {png1, png2, png3, png4}, pad_png(f{1}, MARGIN_PX); end
fprintf('\nFigure 1 saved to %s\n', png1);
fprintf('Figure 2 saved to %s\n', png2);
fprintf('Figure 3 saved to %s\n', png3);
fprintf('Figure 4 saved to %s\n', png4);

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
% edge and look cramped the moment the image is dropped onto a poster next to anything else.
% This adds a quiet border in whatever colour the image already has at its corner, so it follows
% the palette automatically. Done afterwards rather than via exportgraphics' own Padding option,
% which only offers 'tight' or 'figure' -- neither is a controllable number of pixels.
function pad_png(file, px)
    img = imread(file);
    if size(img, 3) == 1, img = repmat(img, 1, 1, 3); end
    bg = img(1, 1, :);                       % corner pixel = the exported background colour
    [h, w, ~] = size(img);
    out = repmat(bg, h + 2*px, w + 2*px, 1);
    out(px+1:px+h, px+1:px+w, :) = img;
    imwrite(out, file);
end
