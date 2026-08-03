%% plot_uMUL_poster.m                                    ECE SPARK 2026 -- uMUL multiplier
%
% Twin of AND_mul/plot_AND_Mul_poster.m. Same axes, same marks, same conventions, so the two
% pairs of figures can be placed side by side on the poster and read against each other.
%
%     uMUL_poster_warmup.csv           -> plot 1, zero-fault early termination
%     uMUL_poster_faults_stream.csv    -> plot 2, REGISTER HARDENED, f = 0..36 over 256 sites
%     uMUL_poster_faults_all.csv       -> plot 3, NOTHING HARDENED,  f = 0..36 over 264 sites
%     uMUL_poster_sensitivity.csv      -> plot 4, per-bit single-flip delta
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
% PLOTS 2 AND 3 -- A BEFORE-AND-AFTER ON ONE DESIGN DECISION: HARDEN THE REGISTER, OR NOT
%
%   THE TWO PANELS ARE THE SAME EXPERIMENT WITH ONE CHANGE -- whether the 8 value-register bits
%   are in the target set. Same operand, same flip counts, same log x axis, same marks, so
%   everything that differs between them is caused by those 8 bits.
%
%   PLOT 2, REGISTER HARDENED   surface = the 256 in_0 stream bits only.
%       X is bits flipped at 0, 1, 2, 4, 8, 16, 32 on a LOG axis; zero is drawn at x = 0.5 and
%       labelled "0", since log(0) does not exist -- its position is cosmetic, its height real.
%       BLUE STICK  full reach, caps at the measured extremes
%       RED DASHES  the equal-tailed middle 90% (5th / 95th percentile), not forced symmetric
%       BLUE DOT    trial-weighted mean, pinned at 0.25 at every flip level
%       Vertical scale 0.10-0.40, SHARED with the AND gate's fault figure -- this panel damages
%       the same KIND of storage the AND gate has, so the two are directly comparable.
%       HARD-BOUNDED at 0.25 +/- f/256, asserted in the .cpp. 32 flips move it 0.1875 to 0.3125.
%
%   PLOT 3, NOTHING HARDENED    surface = all 264 bits, flips spread UNIFORMLY over stream AND
%       register with no site favoured. This is the panel that describes a part in a radiation
%       field, where an upset does not respect module boundaries.
%       ORANGE STICK  full reach, caps at the measured extremes
%       ORANGE DOT    trial-weighted mean -- still 0.25 at every flip count, which is the surprise
%       Vertical scale is the full 0-0.70. NO 90% band; see below.
%       UNBOUNDED FROM ONE UPSET. The floor is a hard 0.000 at EVERY flip count from 1 to 36,
%       asserted in the .cpp: whatever f is, some layout puts one flip on value bit 7 and the
%       rest anywhere, and clearing the only set bit of value = 128 makes the operand 0.
%
%   THE COMPARISON IS THE POSTER:  32 flips with the register safe  ->  0.1875 to 0.3125
%                                   1 flip  with nothing safe       ->  0.0000 to 0.3750
%   256 bits that are cheap to lose and 8 that are not. The AND gate has no plot 3 at all,
%   because all 512 of its operand bits are worth the same 1/256 -- it fails constantly and
%   negligibly, uMUL fails rarely and catastrophically.
%
%   HOW OFTEN, NOT JUST HOW BAD. The zero floor on plot 3 carries P(bit 7 struck) = f/264
%   exactly -- ONE site out of 264 -- which the .cpp asserts and this script prints: 0.4% at one
%   flip, 12% at 32. Quote that number alongside the floor or the panel overstates the risk.
%
%   NO 90% BAND ON PLOT 3, DELIBERATELY. That distribution is a mixture of a tight cluster at
%   0.25 and a cluster near zero, and a quantile over a mixture steps discontinuously the instant
%   the mixing weight crosses the tail fraction -- f/264 passes 5% at f = 13.2. An earlier
%   version showed a 0.164 cliff there that is in the statistic, not the hardware. See
%   HOW_THIS_DATA_WAS_MADE.txt section 6e.
%
% PLOT 4 -- THE 8-BIT OPERAND REGISTER, EXPANDED
%   Where plot 3's floor comes from: each register bit's single-flip damage MAGNITUDE, against
%   one stream bit's 1/256 for scale. A geometric ladder, 1x to 64x. Absolute value, so the
%   ladder holds for any operand -- only the SIGN depends on what the register holds.
%
% ------------------------------------------------------------------------------------------
% WHERE THE NUMBERS COME FROM. Every value plotted here was measured by clocking real 256-bit
% streams through a real UnaryMultiplier -- about 740,000 simulations for the stream campaign and
% a further 181,000 for the register campaign, each reported in ITS OWN figure title and neither
% borrowed from the other. Nothing is evaluated from a formula.
%
% The combinatorics supply only the WEIGHT on each measured outcome. uMUL's output depends only
% on (value, ones-count of in_0) and never on where those ones sit, so the whole 5.8e75-strong
% arrangement space collapses onto a pair of integers. That weighting is itself checked against
% MC_Probability -- 20,000 genuinely random trials per flip level -- and this script prints the
% worst gap between the two. See the .cpp header and the .txt for the full argument.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
warmCsv   = fullfile(scriptDir, 'uMUL_poster_warmup.csv');
streamCsv = fullfile(scriptDir, 'uMUL_poster_faults_stream.csv');
allCsv    = fullfile(scriptDir, 'uMUL_poster_faults_all.csv');
for fchk = {warmCsv, streamCsv, allCsv}
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
    'bg',    [0.80 0.80 0.80], ...   % the PLOT PANEL only -- light grey
    'figbg', [1.00 1.00 1.00], ...   % figure surround + exported margin: WHITE
    'fg',    [0.10 0.10 0.10], ...   % axes, ticks, labels, titles
    'grid',  [1.00 1.00 1.00], ...   % divider lines, WHITE on the grey panel
    'legbg', [0.93 0.93 0.93], ...   % legend panel, near-white so it lifts off the light panel
    'tick',  14, ...                 % tick labels
    'label', 15, ...                 % axis labels
    'title', 15);                    % titles
MARGIN_PX = 40;                      % border around every exported PNG; 300 DPI ~ 3.4 mm

%% ---------------------------------------------------------------------------------------
%% Load
%% ---------------------------------------------------------------------------------------
W  = readtable(warmCsv,   'Delimiter', ',');
FS = readtable(streamCsv, 'Delimiter', ',');   % register HARDENED  -- 256 sites
FR = readtable(allCsv,    'Delimiter', ',');   % nothing hardened   -- 264 sites

% Real gate runs behind each fault figure, counted SEPARATELY -- the two campaigns share nothing,
% so each panel quotes only its own. LevelRealTrials is constant within a flip level and repeated
% on each of that level's rows, so take one value per level and add them up; summing the column
% directly would multiply by the row count.
[lvlS, iS] = unique(FS.BitsFlipped);
[lvlR, iR] = unique(FR.BitsFlipped);
runsStream = sum(FS.LevelRealTrials(iS));
runsAll    = sum(FR.LevelRealTrials(iR));
mcPerLevel = FS.MC_Trials(1);

% Worst disagreement between the enumerated share of trials and what the random trials measured,
% checked in both campaigns independently.
gapOf = @(T, lv) max(arrayfun(@(L) ...
    max(abs(T.Probability(T.BitsFlipped == L) / sum(T.Probability(T.BitsFlipped == L)) ...
          - T.MC_Probability(T.BitsFlipped == L))), lv));
mcGapS = gapOf(FS, lvlS);
mcGapR = gapOf(FR, lvlR);

fprintf('Warm-up   : %d truncation lengths, exact over 10^%.0f in_0 arrangements\n', ...
        height(W), W.Log10Arrangements(1));
fprintf('Hardened  : %d outcome rows, %d flip levels, 256 sites, %s real gate simulations\n', ...
        height(FS), numel(lvlS), commafy(runsStream));
fprintf('Unhardened: %d outcome rows, %d flip levels, 264 sites, %s real gate simulations\n', ...
        height(FR), numel(lvlR), commafy(runsAll));
fprintf('Random    : %s trials per flip level in each campaign\n', commafy(mcPerLevel));
fprintf('Check     : worst enumerated-vs-random gap  hardened %.4f   unhardened %.4f\n', ...
        mcGapS, mcGapR);

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

% ---------------------------------------------------------------------------------------
% Per-flip-level summary of one campaign. Written once and called twice, because the two
% campaigns are now structurally identical -- a flip count, a set of outcomes, a weight on each.
% The old version needed two different code paths here only because it was conditioning one
% joint population; with separate sweeps that asymmetry is gone.
% ---------------------------------------------------------------------------------------
function S = summarise(T, TAIL)
    S.lvl = unique(T.BitsFlipped);
    n = numel(S.lvl);
    [S.mn, S.mx, S.mu, S.lo, S.hi, S.cov] = deal(zeros(n,1));
    for k = 1:n
        sel = T.BitsFlipped == S.lvl(k);
        v = T.OutputFraction(sel);
        p = T.Probability(sel);  p = p / sum(p);
        S.mn(k) = min(v);  S.mx(k) = max(v);  S.mu(k) = sum(v .* p);
        [S.lo(k), S.hi(k)] = tail_band(v, p, TAIL);
        % True coverage of the drawn band, so the "90%" claim on the figure is checked rather
        % than asserted. Discrete outcomes mean the edges land ON a value, so coverage runs over.
        S.cov(k) = sum(p(v >= S.lo(k) & v <= S.hi(k)));
    end
end

SS = summarise(FS, TAIL);   % register hardened -- 256 sites
SR = summarise(FR, TAIL);   % nothing hardened  -- 264 sites

fprintf('\nREGISTER HARDENED (256 sites, register never struck)\n');
fprintf('  flips     mean      min       max      5th pct       95th pct    coverage\n');
for k = 1:numel(SS.lvl)
    if ismember(SS.lvl(k), SHOW_FLIPS)
        fprintf('  %5d  %8.6f  %8.6f  %8.6f  %11.6f  %13.6f  %9.4f\n', ...
                SS.lvl(k), SS.mu(k), SS.mn(k), SS.mx(k), SS.lo(k), SS.hi(k), SS.cov(k));
    end
end

% P(bit 7 struck) = f/264 -- the mass sitting on the zero floor. Printed rather than plotted,
% because a 90% band over this population steps discontinuously (see the figure 3 header).
fprintf('\nNOTHING HARDENED (264 sites, flips spread evenly over stream + register)\n');
fprintf('  flips     mean      min       max    P(reg hit)  P(bit 7)   f/264\n');
for k = 1:numel(SR.lvl)
    if ismember(SR.lvl(k), SHOW_FLIPS)
        sel = FR.BitsFlipped == SR.lvl(k);
        fprintf('  %5d  %8.6f  %8.6f  %8.6f  %10.4f  %8.4f  %8.4f\n', ...
                SR.lvl(k), SR.mu(k), SR.mn(k), SR.mx(k), ...
                sum(FR.ProbRegisterHit(sel)), sum(FR.ProbMSBLost(sel)), SR.lvl(k) / 264);
    end
end
% The headline comparison, printed so it is checked every run rather than remembered.
fprintf('\nHardened,   worst reach at 32 flips : %.6f to %.6f\n', ...
        SS.mn(SS.lvl == 32), SS.mx(SS.lvl == 32));
fprintf('Unhardened, worst reach at  1 flip  : %.6f to %.6f\n', ...
        SR.mn(SR.lvl == 1), SR.mx(SR.lvl == 1));
% THE ASSERTION THE PANEL RESTS ON: the floor is a hard zero at EVERY flip count, not just high
% ones. If this ever prints a nonzero the "you can zero out from anywhere" claim is dead.
fprintf('Unhardened, floor across f = 1..36  : max = %.6f (should be 0)\n', ...
        max(SR.mn(SR.lvl >= 1)));

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
ylabel(ax1, 'Mean Absolute Error', 'FontSize', 15);
% TITLE IS THE OPERATION AND THE TRIAL BASIS, NOTHING ELSE -- the poster labels the panel above
% it. Arrangements, not realRuns: realRuns counts the fault campaign, this is the warm-up sweep.
title(ax1, sprintf('0.5 \\times 0.5 = 0.25   |   exhaustive over all 10^{%.0f} in_0 arrangements', ...
                   W.Log10Arrangements(1)), 'FontSize', 15);
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
%% FIGURE 2 -- FAULT RESPONSE, STREAM ONLY           -> uMUL_poster_faults_intact.png
%%
%% Its own campaign: f flips land on the 256-bit in_0 stream and the value register is never
%% touched. Not a filtered slice of a joint sweep -- see the header.
%%
%% THE PANEL IS HARD-BOUNDED, AND THE .cpp ASSERTS IT. With the register protected, one stream
%% flip moves n1 by exactly one, so the answer cannot leave 0.25 +/- f/256. At 32 flips that is
%% 0.125 to 0.375 in the absolute worst case, and the measured extremes sit inside it.
%%
%% NOTE THE EXPANDED VERTICAL SCALE. This panel runs 0.10 to 0.40, not 0 to 0.70 like figure 3.
%% Plotted on figure 3's axis this whole curve would be a flat line, which is itself the point --
%% but it would show nothing. The two are deliberately NOT on a shared y axis.
%% ---------------------------------------------------------------------------------------
fig2 = figure('Name', 'uMUL -- Fault Response, Stream Only', ...
              'Units', 'Normalized', 'Position', [0.10, 0.12, 0.54, 0.72]);
ax2 = axes(fig2);
hold(ax2, 'on');

for s = 1:numel(SHOW_FLIPS)
    k = find(SS.lvl == SHOW_FLIPS(s), 1);
    if isempty(k), continue; end
    x = xpos(s);
    % full reach, caps at the measured extremes
    plot(ax2, [x x], [SS.mn(k) SS.mx(k)], '-', 'Color', BLUE, 'LineWidth', 1.6);
    plot(ax2, [x/capBlue x*capBlue], [SS.mn(k) SS.mn(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);
    plot(ax2, [x/capBlue x*capBlue], [SS.mx(k) SS.mx(k)], '-', 'Color', BLUE, 'LineWidth', 1.8);
    % middle 90%, suppressed where it is the whole reach (too few outcomes to carve 90% from)
    if SS.lo(k) > SS.mn(k) || SS.hi(k) < SS.mx(k)
        plot(ax2, [x/capRed x*capRed], [SS.lo(k) SS.lo(k)], '-', 'Color', RED, 'LineWidth', 2.4);
        plot(ax2, [x/capRed x*capRed], [SS.hi(k) SS.hi(k)], '-', 'Color', RED, 'LineWidth', 2.4);
    end
end
selShow = ismember(SS.lvl, SHOW_FLIPS);
hMean = plot(ax2, xpos, SS.mu(selShow), 'o', 'MarkerSize', 6, ...
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
% scale, it shows that uMUL's stream behaves like the AND gate's streams, because it IS the same
% kind of storage. The comparison only works if the axes match, so do not retune one alone.
ylim(ax2, [0.10 0.40]);
xticks(ax2, xpos);  xticklabels(ax2, string(SHOW_FLIPS));
xlabel(ax2, 'Bits Flipped  (over the 256-bit in\_0 stream) [log]', 'FontSize', 14);
ylabel(ax2, 'Final Answer', 'FontSize', 15);
% Its OWN trial count now -- runsStream, not a figure shared with the register panel. The two
% panels are also told apart by their x-axis labels and ranges, which is what replaced the
% descriptive title line the poster no longer wants.
title(ax2, sprintf('0.5 \\times 0.5 = 0.25   |   %s real gate simulations', ...
                   commafy(runsStream)), 'FontSize', 15);

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
%% FIGURE 3 -- FAULT RESPONSE, NOTHING HARDENED      -> uMUL_poster_faults_damaged.png
%%
%% Its own campaign: f flips spread UNIFORMLY over all 264 operand bits -- 256 stream and 8
%% register, no site favoured. This is the panel that describes a part sitting in a radiation
%% field, because an upset does not respect module boundaries.
%%
%% IT IS THE SAME EXPERIMENT AS FIGURE 2 WITH ONE CHANGE: the 8 register bits are in the target
%% set. Same operand, same flip counts, same x axis, same marks. Everything that differs between
%% the two panels is caused by those 8 bits, which is what makes the pair an argument for
%% selective hardening rather than two unrelated measurements.
%%
%%      figure 2, register hardened     f = 32  ->  0.1875 to 0.3125, bounded by f/256
%%      figure 3, nothing hardened      f =  1  ->  0.0000 to 0.3750, unbounded from ONE upset
%%
%% THE FLOOR IS A HARD ZERO AT EVERY FLIP COUNT FROM ONE, and that is the point of the panel.
%% Whatever f is, there is a layout that puts one flip on value bit 7 and the rest anywhere; bit 7
%% is the only set bit of value = 128, so clearing it makes the operand 0, the strict ">"
%% comparator never fires, and the answer is a bit-perfect zero regardless of what happened to the
%% stream. The .cpp asserts min == 0 at every f from 1 to 36.
%%
%%     (An earlier revision made this panel a register-ONLY sweep, f = 0..8. That answered "what
%%      if you aim f upsets at the register", which is a fault-injection question rather than a
%%      radiation one, and it had the odd artifact that damage FELL as flips rose -- value' = 0
%%      needs the subset {bit 7} ALONE, so a second register flip lifts the operand off zero.
%%      Spreading uniformly removes the artifact: the extra flips go on the stream, where they
%%      cannot undo a dead operand.)
%%
%% HOW OFTEN, NOT JUST HOW BAD. The zero floor carries P(bit 7 struck) = f/264 exactly -- ONE site
%% out of 264 -- which the .cpp asserts. That is 0.4% at one flip and 13.6% at 36. Rare and
%% catastrophic, and the .m prints the column so the "rare" half is quantified rather than implied.
%%
%% NO 90% BAND IS DRAWN HERE, deliberately. This distribution is a MIXTURE of a tight cluster at
%% 0.25 and a cluster near zero, and a quantile over a mixture jumps discontinuously the instant
%% the mixing weight crosses the tail fraction -- f/264 passes 5% at f = 13.2, and an earlier
%% version of this figure showed a 0.164 cliff there that is in the statistic, not the hardware.
%% Reach plus mean plus the printed P(bit 7) says everything a band would, without the artifact.
%% See HOW_THIS_DATA_WAS_MADE.txt section 6e.
%%
%% THE MEAN SITS ON 0.25 ANYWAY. Do not read it as typical: it is the balance point between the
%% two modes. At one register flip the eight possible operands 0, 192, 160, 144, 136, 132, 130,
%% 129 give answers 0, 96, 80, 72, 68, 66, 65, 65 out of 256, which sum to 512 -- so the single
%% catastrophic zero is worth exactly as much as the seven upward pushes combined, and 64 is not
%% one of the values on that list.
%% ---------------------------------------------------------------------------------------
fig3 = figure('Name', 'uMUL -- Fault Response, Nothing Hardened', ...
              'Units', 'Normalized', 'Position', [0.14, 0.10, 0.54, 0.72]);
ax3 = axes(fig3);
hold(ax3, 'on');

% SAME xpos, capBlue AND log axis AS FIGURE 2. The two panels are only comparable if the x
% geometry is identical, and the whole design of the pair is that comparison.
for s = 1:numel(SHOW_FLIPS)
    k = find(SR.lvl == SHOW_FLIPS(s), 1);
    if isempty(k), continue; end
    x = xpos(s);
    plot(ax3, [x x], [SR.mn(k) SR.mx(k)], '-', 'Color', ORANGE, 'LineWidth', 1.6);
    plot(ax3, [x/capBlue x*capBlue], [SR.mn(k) SR.mn(k)], '-', 'Color', ORANGE, 'LineWidth', 1.8);
    plot(ax3, [x/capBlue x*capBlue], [SR.mx(k) SR.mx(k)], '-', 'Color', ORANGE, 'LineWidth', 1.8);
end
selShowR = ismember(SR.lvl, SHOW_FLIPS);
hMeanD = plot(ax3, xpos, SR.mu(selShowR), 'o', 'MarkerSize', 7, ...
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
xlabel(ax3, 'Bits Flipped  (over all 264 operand bits, stream + register) [log]', 'FontSize', 14);
ylabel(ax3, 'Final Answer', 'FontSize', 15);
title(ax3, sprintf('0.5 \\times 0.5 = 0.25   |   %s real gate simulations', ...
                   commafy(runsAll)), 'FontSize', 15);

hold(ax3, 'on');
pOr = plot(ax3, NaN, NaN, '-', 'Color', ORANGE, 'LineWidth', 1.8);
hold(ax3, 'off');
lg3 = legend(ax3, [pOr hMeanD], ...
       {'full reach, caps = measured extremes', ...
        'trial-weighted mean'}, ...
       'Location', 'northwest', 'FontSize', 11);
lg3.Box = 'on';
ax3.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% FIGURE 4 -- THE 8-BIT OPERAND REGISTER, EXPANDED  -> uMUL_poster_register.png
%%
%% Where the damage comes from. Each of the 8 value-register bits as the MAGNITUDE of the change
%% it makes to the answer when it alone is flipped, against the stream bit's 1/256 for scale.
%%
%% --- WHY THIS PANEL IS ABSOLUTE VALUE, AND WHY THAT MAKES IT UNIVERSAL --------------------
%%
%% MAGNITUDE IS A PROPERTY OF THE BIT. SIGN IS A PROPERTY OF THE DATA. Flipping register bit b
%% changes the operand by exactly 2^b, so the answer changes by
%%
%%       |delta|  =  (2^b / 2^width) * p0        <- p0 = density of the in_0 stream
%%
%% and NOTHING in that expression depends on the operand's value. The sign does, and only the
%% sign: if bit b was CLEAR the flip sets it and the answer rises, if it was SET the flip clears
%% it and the answer falls. That is the whole of the difference.
%%
%% PLOTTING THE SIGNED VERSION MADE THE FIGURE ABOUT 128, NOT ABOUT uMUL. With value = 128 the
%% register is one-hot, so seven bars went up and exactly one went down -- and that 7-up-1-down
%% shape is an accident of the operand chosen. At value = 64 it would be bit 6 pointing down. At
%% value = 255 every bar would point down. A reader could not tell which parts of the picture
%% were the architecture and which were the example.
%%
%% ABSOLUTE VALUE REMOVES THE EXAMPLE AND LEAVES THE ARCHITECTURE. What is left is a clean
%% geometric ladder, 1/512 to 128/512, doubling at every bit, TRUE FOR EVERY OPERAND:
%%
%%       bit b   |delta| = 2^b/512      vs one stream bit = 1/256
%%           0      0.002               0.5x
%%           1      0.004               1x     <- one register bit is worth one stream bit
%%           ...
%%           7      0.250               64x    <- one register bit is worth 64 stream bits
%%
%% THE 64:1 RATIO IS THE POSTER'S NUMBER AND IT IS OPERAND-INDEPENDENT. The top register bit
%% always carries 2^7/2^8 = half the operand's full scale, and a stream bit always carries 1/256
%% of the answer, so the ratio is 2^(width-1)/... = 64 for any 8-bit value and any stream density.
%% Say "one upset in the top register bit does what 64 stream upsets do" and it is true whatever
%% number the register happens to hold.
%%
%% WHAT IS LOST BY TAKING THE ABSOLUTE VALUE, stated so the poster does not overclaim: the reader
%% can no longer see that for THIS operand bit 7 is the one that goes to zero. That fact is not
%% deleted, it is relocated -- it is figure 3's floor, and it is section 6b of the .txt. This
%% panel now answers "how much can one register bit do", and figure 3 answers "in which
%% direction, for this operand". Two questions, two panels.
%%
%% Every bar is a measured gate run, from uMUL_poster_sensitivity.csv. The magnitudes are taken
%% with abs() at plot time; the CSV keeps the signed Delta, so nothing is thrown away on disk.
%% ---------------------------------------------------------------------------------------
S = readtable(fullfile(scriptDir, 'uMUL_poster_sensitivity.csv'), 'Delimiter', ',');

fig4 = figure('Name', 'uMUL -- Operand Register Sensitivity', ...
              'Units', 'Normalized', 'Position', [0.18, 0.16, 0.52, 0.64]);
ax4 = axes(fig4);

regDelta = zeros(8,1);
for b = 0:7
    r = find(string(S.Site) == sprintf("value_bit_%d", b), 1);
    regDelta(b+1) = abs(S.Delta(r));       % <- MAGNITUDE ONLY. See the note above.
end
streamDelta = max(abs(S.Delta(startsWith(string(S.Site), "stream"))));

% THE IDEAL LADDER, (2^b / 2^width) * p0. NOT DRAWN -- the reference line, its legend and the
% explanatory annotation were all taken off the panel to leave just the bars and their labels.
% It is still computed, because the check below is what guarantees the bars mean what the .txt
% says they mean; it simply fails loudly at generation time instead of arguing on the poster.
predicted = (2.^(0:7))' / 512;

% WHERE THE LADDER IS EXACT, AND WHY IT IS NOT EXACT EVERYWHERE. Damage in OUTPUT BITS is the
% number of consumed Sobol points the operand sweeps past, and the first n1 = p0*256 points form
% a (0,k,1)-net: every dyadic interval of width 2^(8-k) holds exactly one point. So a bit-b flip,
% which moves the operand across an aligned interval of width 2^b, crosses exactly 2^b*p0 points
% -- PROVIDED that number is at least 1. Below
%       b* = log2(1/p0)
% the interval is finer than the net can resolve and the bit collapses onto the 1-output-bit
% floor. At p0 = 0.5, b* = 1, so bit 0 alone is floored. At p0 = 0.25 it would be bits 0 and 1.
% Verified exhaustively over all 256 operands and both densities in
% tests/Unit Tests/test_uMUL.cpp :: RegisterBitDamageCountsSobolPointsCrossed.
p0     = 0.5;
bstar  = log2(1 / p0);
exactB = 0:7 >= bstar;
if any(abs(regDelta(exactB) - predicted(exactB)) > 1e-12)
    warning('uMUL:registerLadder', ...
            ['Register magnitudes depart from (2^b/2^w)*p0 ABOVE the resolution limit b* = %g. ' ...
             'That ladder is what makes figure 4 operand-independent -- check the RNG width.'], ...
            bstar);
end

hold(ax4, 'on');
bar(ax4, 0:7, regDelta, 0.72, 'FaceColor', RED, 'EdgeColor', STY.fg, 'LineWidth', 0.6);
% The stream bit's damage, for scale. One line now rather than a +/- pair, because on an
% absolute-value axis there is only one magnitude to draw. Label sits RIGHT, over the tall bars,
% where there is empty panel; on the left it ran straight through bits 0-2's value text.
yline(ax4, streamDelta, '--', 'Color', BLUE, 'LineWidth', 1.6, ...
      'Label', 'any one of the 256 stream bits = 1/256', 'LabelHorizontalAlignment', 'right', ...
      'LabelVerticalAlignment', 'top', 'FontSize', 11, 'Color', BLUE);
for b = 0:7
    d = regDelta(b+1);
    % Multiples of a stream bit -- the comparison the poster actually makes. Printed on every bar
    % because the geometric doubling is the point and the low bars are too short to read by eye.
    text(ax4, b, d + 0.009, sprintf('%.3f\n%gx', d, round(d / streamDelta * 2) / 2), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', STY.fg);
end

% NO REFERENCE LINE AND NO ANNOTATION ON THIS PANEL. The bars and their labels are the whole
% figure; the reasoning behind them (damage in output bits = Sobol points crossed, and why bits 0
% and 1 share a rung below b* = log2(1/p0)) lives in HOW_THIS_DATA_WAS_MADE.txt section 6b and is
% asserted in tests/Unit Tests/test_uMUL.cpp, which is where an argument belongs rather than
% competing with the data for panel space.
hold(ax4, 'off');
xlim(ax4, [-0.7 7.7]);
ylim(ax4, [0 0.30]);
xticks(ax4, 0:7);
grid(ax4, 'on');  ax4.XMinorGrid = 'off';  ax4.YMinorGrid = 'off';
xlabel(ax4, 'Value-register bit', 'FontSize', 15);
ylabel(ax4, 'Size of the Change in the Answer', 'FontSize', 15);
% This panel is not a fault-campaign figure, so it quotes its own basis: all 8 single-bit upsets.
title(ax4, '0.5 \times 0.5 = 0.25   |   every single-bit upset in the 8-bit register', ...
      'FontSize', 15);
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
