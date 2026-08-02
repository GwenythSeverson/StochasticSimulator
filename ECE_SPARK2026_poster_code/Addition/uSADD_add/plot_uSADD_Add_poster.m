%% plot_uSADD_Add_poster.m                          ECE SPARK 2026 -- uGEMM scaled adder
%
% Companion to MUX_add/plot_MUX_Add_poster.m and the two multiplier scripts. Same axes, same
% marks, same conventions.
%
%     uSADD_Add_poster_warmup.csv  -> plot 1, zero-fault early termination
%     uSADD_Add_poster_faults.csv  -> plot 2, fault response over the 512-bit stream surface
%
% uSADD_Add_poster_state.csv is also written by the .cpp but is NOT plotted here. It holds the
% measured PC and accumulator strike data -- every bit, every strike cycle -- and the headline
% from it is a single number rather than a shape: the worst damage ONE strike on uSADD's binary
% state can do is 1 output bit, 0.0039, against uMUL's 0.25 for one strike on its operand
% register. That is a sentence, not a bar chart, so it is reported in the console summary and in
% HOW_THIS_DATA_WAS_MADE.txt rather than given a panel. The CSV is still generated and still the
% place to go for the per-bit, per-cycle detail.
%
% ------------------------------------------------------------------------------------------
% THE UNIT. A parallel counter sums the input bits; an accumulator holds the running credit and
% emits one bit whenever it reaches n. Nothing is discarded, so out = (a + b)/2 exactly, with no
% correlation requirement, no select stream, and no RNG.
%
% TWO OPERAND STREAMS, NOT THREE. The MUX next door needs a 0.5 select stream and pays for it
% with 256 extra fault sites; uSADD does not have one. Its surface is 512 stream bits plus THREE
% bits of state: a ONE-BIT accumulator and a 2-bit PC bus.
%
% WHY THE ACCUMULATOR IS ONE BIT. uGEMM takes the output from the CARRY of the accumulator, so
% the register only ever stores the residue -- for two inputs that is {0, 1}. The wide pre-drain
% sum exists inside the adder, not in a flip-flop. An earlier revision of this project sized the
% register for the sum instead, which gave it two bits; the healthy arithmetic was identical but
% the accumulator's fault cross-section was double what the real design has.
%
% WHY THAT STATE IS CHEAP TO LOSE. Both uSADD and uMUL carry binary state, and the contrast is
% the whole poster:
%     uMUL's 8-bit register holds an OPERAND. One MSB strike sets it to zero and every remaining
%       cycle of the run is wrong: 0.25 of damage from one upset, never repaired.
%     uSADD's 1-bit accumulator holds a RESIDUE -- credit not yet spent. One strike perturbs the
%       running total once, the unit keeps integrating correctly, and the damage stops there:
%       at most 1 output bit, 0.0039, whatever the strike and whenever it lands -- the same
%       as a single stream bit is worth. Its state is no more dangerous than its data.
% Same idea (put state in a binary register), opposite consequence, because of WHAT the register
% holds.
%
% ------------------------------------------------------------------------------------------
% WHY THE 1-FLIP WHISKER ONLY GOES DOWN -- the first thing anyone asks about figure 2
%
% At one flip the ONLY reachable outcomes are 127 and 128 output ones. There is no 129, so the
% whisker has nothing above 0.5 to draw. That is the floor, not a plotting bug:
%
%     uSADD emits floor(total input ones / 2), and the clean total here is 256 -- an EVEN number
%     sitting exactly on a floor boundary.
%         flip a 1 -> 0   total 255   floor(255/2) = 127   loses a whole output bit
%         flip a 0 -> 1   total 257   floor(257/2) = 128   GAINS NOTHING
%     The extra credit from the second case is real, but it is stranded in the accumulator as
%     residue rather than becoming an output one. So a single flip can take a bit away or do
%     nothing; it can never add one.
%
% IT IS A PARITY EFFECT AND IT REPEATS. The net credit change always has the same parity as the
% flip count, so an odd f always leaves the total odd and the floor always eats the half bit:
%
%     f      1     2     3     4     5     6
%     reach -1/+0 -1/+1 -2/+1 -2/+2 -3/+2 -3/+3
%     mean  .4980 .5000 .4980 .5000 .4980 .5000
%
% Every ODD flip level reaches one step further down than up and sits half an output bit low;
% every EVEN one is symmetric and exactly on target. This is the same residue that is the unit's
% only healthy error -- the accumulator can WITHHOLD output but never MANUFACTURE it, so uSADD's
% error is always a deficit. It is the one place in this poster set where a mean is not pinned on
% its target, and it is worth saying out loud rather than letting a reviewer find it.

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
warmCsv  = fullfile(scriptDir, 'uSADD_Add_poster_warmup.csv');
faultCsv = fullfile(scriptDir, 'uSADD_Add_poster_faults.csv');
stateCsv = fullfile(scriptDir, 'uSADD_Add_poster_state.csv');
for fchk = {warmCsv, faultCsv, stateCsv}
    if ~isfile(fchk{1})
        error(['CRITICAL: %s not found.\n' ...
               'Run:  stochastic_computer.exe --gtest_filter=uSADDAddPoster.*'], fchk{1});
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
S = readtable(stateCsv, 'Delimiter', ',');

BLUE   = [0.10 0.28 0.85];
RED    = [0.90 0.20 0.20];
ORANGE = [0.95 0.55 0.10];

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

fprintf('uSADD: %d outcome rows across %d flip levels, 512-bit stream surface\n', height(F), nf);
fprintf('\n  flips     mean      min       max      5th pct      95th pct\n');
for k = 1:nf
    if ismember(flipLevels(k), SHOW_FLIPS)
        fprintf('  %5d  %8.6f  %8.6f  %8.6f  %11.6f  %12.6f\n', ...
                flipLevels(k), mu(k), mn(k), mx(k), bLo(k), bHi(k));
    end
end
fprintf('\nWorst single state strike: %d output bits = %.6f\n', ...
        max(abs(S.DeltaOnes)), max(abs(S.DeltaFraction)));

%% ---------------------------------------------------------------------------------------
%% FIGURE 1 -- zero-fault early termination
%% ---------------------------------------------------------------------------------------
fig1 = figure('Name', 'uSADD -- Zero-Fault Warm-Up', ...
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
title(ax1, {'uSADD -- Zero-Fault Warm-Up', ...
            sprintf(['(0.5 + 0.5)/2 = 0.5   |   exhaustive over all 10^{%.0f} arrangements, ' ...
                     'NO decorrelation filter'], W.Log10Arrangements(1))}, 'FontSize', 13);
ax1.Toolbar.Visible = 'off';

%% ---------------------------------------------------------------------------------------
%% FIGURE 2 -- fault response over the stream surface
%% ---------------------------------------------------------------------------------------
fig2 = figure('Name', 'uSADD -- Fault Response', ...
              'Units', 'Normalized', 'Position', [0.14, 0.10, 0.56, 0.72]);
ax2 = axes(fig2);
hold(ax2, 'on');

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
% Y bounded to 0.35 - 0.65, i.e. 0.5 +/- 0.15, SHARED with the MUX fault figure. That shared
% axis is the whole point of the addition pair: on it, uSADD's whisker is visibly HALF the
% MUX's at every flip level (0.438-0.563 against 0.375-0.625 at 32 flips), which is the 1/2
% scaling attenuating input faults before they reach the output. Do not retune one without the
% other or the comparison stops being fair.
ylim(ax2, [0.35 0.65]);
xticks(ax2, xpos);  xticklabels(ax2, string(SHOW_FLIPS));
xlabel(ax2, 'Bits Flipped  (over the 512-bit stream surface, a + b) [log]', 'FontSize', 12);
ylabel(ax2, 'Fraction of Ones in the Final Answer', 'FontSize', 13);
% Subtitle leads with the operand expression, matching the AND gate's "0.5 x 0.5 = 0.25" and
% this folder's own warm-up figure, so the four fault-response panels state their arithmetic the
% same way. The "no select stream" point is still made -- it is in the surface size, 512 against
% the MUX's 768 -- and spelling it out here as well was one clause too many for the line.
title(ax2, {'uSADD -- Fault Response', ...
            '(0.5 + 0.5)/2 = 0.5   |   two operand streams, 512 flippable bits'}, ...
      'FontSize', 13);
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

png1 = fullfile(scriptDir, 'uSADD_Add_poster_warmup.png');
png2 = fullfile(scriptDir, 'uSADD_Add_poster_faults.png');
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
function pad_png(file, px)
    img = imread(file);
    if size(img, 3) == 1, img = repmat(img, 1, 1, 3); end
    bg = img(1, 1, :);                       % corner pixel = the exported background colour
    [h, w, ~] = size(img);
    out = repmat(bg, h + 2*px, w + 2*px, 1);
    out(px+1:px+h, px+1:px+w, :) = img;
    imwrite(out, file);
end
