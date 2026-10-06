% FIT_SCALING  How fast do the times of the scaling grid grow with N? For
% every time (build, simulation), algorithm and qudit dimension d, fits the
% runs with N >= FIT_MIN_N by least squares (polyfit) with a polynomial of
% degree 1 and one of degree 2 in N, and tests whether the quadratic term is
% needed. Runs no experiment: it only reads results/scaling/scaling_raw.csv,
% written by scaling_accuracy.m.
%
% Expected: the build time of the Bullock construction grows as N^2 in our
% implementation (each of its O(N) reflectors reads and updates the whole
% state vector, O(N) per reflector); those of the two QdKP-Trees grow as N.
% The simulation time should grow as N^2 for all three (O(N) gates, O(N)
% each once the simulation is sparse), so it checks the method: a fit that
% does not find the quadratic term there cannot be trusted to reject it for
% the build times.
%
% The fits:
%   - every trial is a point (N_TRIALS per size), so the residuals include
%     the spread of the runs, and the test has degrees of freedom left even
%     with few sizes;
%   - N is scaled to x = N / N_max, so the coefficients are in seconds: c2
%     and c1 are the quadratic and linear contributions at the largest N, c0
%     the fixed cost;
%   - the quadratic term is needed when the F-test of the degree-2 against
%     the degree-1 fit gives p < ALPHA (the same as a t-test on c2) and c2
%     is at least SHARE_MIN of the fitted time at N_max: a significant but
%     tiny c2 is curvature (caches, memory), not quadratic growth;
%   - the exponent p of a power law t ~ N^p, a degree-1 fit of log t against
%     log N over the same runs, is reported next to them: a significant c2
%     alone does not tell N^2 from, say, N^1.5. The fixed cost of the
%     smallest sizes pulls p down, so it underestimates the growth at the
%     largest N (the local slope between the last two sizes is in tables.txt).
% FIT_MIN_N is set per time. The build times are fitted over every size: the
% fits are in seconds, so the small sizes (a fixed cost of milliseconds)
% barely move the coefficients, but they add runs to the F-test and make p
% smaller, and they raise the largest relative residual (neither polynomial
% describes a fixed cost). The simulation times start at N = 512, since QCLAB
% simulates with dense matrices below it and one fit would mix two regimes.
% A (time, algorithm, d) with fewer than 3 sizes from FIT_MIN_N on gets no
% quadratic fit ("too few sizes").
%
% Output (written to OUT_DIR):
%   fit_scaling.txt   the report: one line per (time, algorithm, d), then the
%                     fits that disagree with EXPECTED (also on the console)
%   fit_scaling.csv   every coefficient, one row per (time, algorithm, d)
%
% How to read the report. Each line is one curve of Figure 5.2 (build) or
% 5.3 (simulation): one algorithm in the panel of one d.
%   sizes, N range  the distinct N fitted (N_TRIALS runs each) and their span
%   degree 1        R^2 and res of t = c1 x + c0
%   degree 2        c2, c1, c0 of t = c2 x^2 + c1 x + c0 (s, at N_max), the
%                   standard error of c2 in % of c2 (+-c2), and res
%   p               chance that the parabola beats the line only by luck
%   share           c2 / fitted time at N_max: how much of the largest time
%                   comes from the N^2 term
%   exponent        slope of log t against log N over the whole range
%   verdict         "quadratic term needed": p < ALPHA and share >= SHARE_MIN;
%                   "linear (small curvature)": p < ALPHA, |share| < SHARE_MIN;
%                   "linear": p >= ALPHA; "concave": share <= -SHARE_MIN;
%                   "too few sizes": fewer than 3 sizes, no parabola
%                   "(NOT AS EXPECTED)" when the verdict contradicts EXPECTED.
% What carries the conclusion:
%   - share and +-c2. A large share with a small +-c2 is quadratic growth
%     (Bullock build, d = 2..4: about 3/4 of the time at N_max, c2 known to a
%     few %). A share of a few % is linear growth with some curvature. A +-c2
%     near or above 100% means c2 cannot be told from 0, whatever its share:
%     the data cannot decide (Bullock build, d = 5: 5 sizes up to 15625 and a
%     large spread between runs, so "linear" there is "undecided", not
%     evidence of linear growth).
%   - the simulation rows are the control: all should be quadratic. If they
%     are not, the fits cannot be trusted on the build times either.
% What not to read too much into:
%   - p. The runs of a size differ by a few %, so p is tiny for any visible
%     curvature, and the small sizes of the build fits add runs with almost no
%     residual (in seconds), which lowers p further. Read it as "the curve is
%     not exactly a line", not as "the curve is quadratic".
%   - degree-1 R^2: above 0.94 even for curves that are plainly quadratic,
%     since the largest N dominate a fit in seconds.
%   - res (both fits): relative misses are largest at the smallest N, where
%     the time is a fixed cost that no polynomial in N describes and where
%     the fit in seconds puts almost no weight; for the build fits, which
%     start at N = d^2, they reach hundreds of %. Compare res between the two
%     fits of a line, not with a threshold.
%   - exponent: the fixed cost of the small sizes pulls it down (below 1 for
%     the build times over every size); it is no estimate of the growth at
%     the largest N, which the local slopes of tables.txt give.
%   - a share just above SHARE_MIN. The threshold is a convention: the
%     Householder build at d = 2 (share about 11%) is a linear curve with a
%     mild upturn over its last two sizes (t/N grows by about 9%), not
%     quadratic growth. Its cause is not determined by these runs: the
%     Givens tree, with as many gates and qudits at d = 2, has about the same
%     c2 in seconds, and a cost per gate growing with n (O(N log N) in total)
%     fits these times worse than an N^2 term.
%
% Run from anywhere:  run('experiments/fit_scaling.m')

clear; clc;

%% ------------------------------------------------------------------------
%  Configuration (edit here)
%  ------------------------------------------------------------------------

ENSEMBLE  = 'real';                      % the timed ensemble of scaling_accuracy.m
TIMES     = {'buildTime', 'simTime'};    % columns of scaling_raw.csv to fit
% fit only the runs with N >= FIT_MIN_N, per time (0: every size)
FIT_MIN_N = struct('buildTime', 0, 'simTime', 512);
ALPHA     = 0.01;                        % significance level of the quadratic term
SHARE_MIN = 0.10;                        % smallest share of c2 in the fitted time at N_max

ALGS  = {'bullock', 'qdkptree', 'qdkpgivens'};
NAMES = {'Bullock', 'Householder', 'Givens'};
% expected degree of growth of each time, per algorithm in the order of ALGS
EXPECTED = struct('buildTime', [2, 1, 1], 'simTime', [2, 2, 2]);

HERE    = fileparts(mfilename('fullpath'));
RES     = fullfile(HERE, 'results');
OUT_REL = fullfile('results', 'fits');      % printed in the log (no local paths)
OUT_DIR = fullfile(HERE, OUT_REL);
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

raw = readtable(fullfile(RES, 'scaling', 'scaling_raw.csv'), 'TextType', 'string');
raw = raw(raw.ensemble == ENSEMBLE, :);

fid = fopen(fullfile(OUT_DIR, 'fit_scaling.txt'), 'w');
out = [1, fid];                          % the report goes to the console and the file
say(out, 'Polynomial fits of the times of the scaling grid, written by fit_scaling.m on %s\n', ...
  char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
logFile = fullfile(RES, 'run_all.log');
if exist(logFile, 'file')
  L = readlines(logFile);
  say(out, 'results from: %s\n', L(1));
end
say(out, ['%s ensemble, every trial a point; x = N / N_max, so ' ...
  'c2, c1, c0 are in s (contributions at the largest N)\n'], ENSEMBLE);
say(out, ['quadratic term needed: F-test of degree 2 against degree 1 with p < %g, ' ...
  'and c2 >= %g of the fitted time at N_max\n'], ALPHA, SHARE_MIN);
say(out, ['res = largest relative residual |t - fit| / t; exponent = slope of ' ...
  'log t against log N\n']);

%% ------------------------------------------------------------------------
%  Fits
%  ------------------------------------------------------------------------

fits = {};
for tm = TIMES
  minN = FIT_MIN_N.(tm{1});
  section(out, char(sprintf('%s, runs with N >= %d (expected degree: %s)', tm{1}, minN, ...
    strjoin(compose('%s %d', string(NAMES(:)), EXPECTED.(tm{1})(:)), ', '))));
  say(out, '%-11s %2s %5s %13s | %-17s | %-41s | %7s %6s %8s  %s\n', '', '', '', '', ...
    'degree 1', 'degree 2', '', '', '', '');
  say(out, '%-11s %2s %5s %13s | %8s %8s | %8s %8s %8s %6s %8s | %7s %6s %8s  %s\n', ...
    'algorithm', 'd', 'sizes', 'N range', 'R^2', 'res', 'c2 [s]', 'c1 [s]', 'c0 [s]', ...
    '+-c2', 'res', 'p', 'share', 'exponent', 'verdict');
  for a = 1:numel(ALGS)
    for d = unique(raw.d(raw.algorithm == ALGS{a}))'
      R = raw(raw.algorithm == ALGS{a} & raw.d == d & raw.N >= minN, :);
      F = fitTimes(R.N, R.(tm{1}), ALPHA, SHARE_MIN);
      expected = EXPECTED.(tm{1})(a);
      agrees = (expected == 2 && F.verdict == "quadratic term needed") || ...
               (expected == 1 && startsWith(F.verdict, "linear"));
      say(out, '%-11s %2d %5d %13s | %8.5f %7.1f%% | %8.3g %8.3g %8.3g %5.0f%% %7.1f%% | %7.1e %5.0f%% %8.2f  %s%s\n', ...
        NAMES{a}, d, F.nSizes, sprintf('%d-%d', F.Nmin, F.Nmax), F.lin_R2, 100 * F.lin_res, ...
        F.quad_c2, F.quad_c1, F.quad_c0, 100 * F.quad_c2_se / abs(F.quad_c2), ...
        100 * F.quad_res, F.p, 100 * F.share, F.exponent, F.verdict, ...
        repmat(' (NOT AS EXPECTED)', 1, double(~agrees)));
      fits{end+1} = [table(string(tm{1}), string(ALGS{a}), d, expected, agrees, ...
        'VariableNames', {'time', 'algorithm', 'd', 'expected', 'agrees'}), ...
        struct2table(F)]; %#ok<SAGROW>
    end
  end
end

fits = vertcat(fits{:});
writetable(fits, fullfile(OUT_DIR, 'fit_scaling.csv'));

section(out, 'Against the expected degrees');
say(out, '%d of %d fits as expected\n', nnz(fits.agrees), height(fits));
for k = find(~fits.agrees)'
  say(out, '  %s %s d=%d: expected degree %d, %s (exponent %.2f)\n', fits.time(k), ...
    fits.algorithm(k), fits.d(k), fits.expected(k), fits.verdict(k), fits.exponent(k));
end

fclose(fid);
fprintf('\nWrote %s\n', fullfile(OUT_REL, 'fit_scaling.txt'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'fit_scaling.csv'));

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

function F = fitTimes(N, t, alpha, shareMin)
% Degree-1 and degree-2 least-squares fits of the times t against
% x = N / max(N), the F-test of the quadratic term, and the exponent of a
% power law t ~ N^p. Fields the data cannot support are NaN.
N = N(:);
t = t(:);
m = numel(t);
F.nSizes = numel(unique(N));
F.nRuns  = m;
F.Nmin   = min(N);
F.Nmax   = max(N);
[F.lin_c1, F.lin_c0, F.lin_R2, F.lin_res] = deal(NaN);
[F.quad_c2, F.quad_c2_se, F.quad_c1, F.quad_c0, F.quad_R2, F.quad_res] = deal(NaN);
[F.Fstat, F.p, F.share, F.exponent] = deal(NaN);
F.verdict = "too few sizes";

x   = N / F.Nmax;
tss = sum((t - mean(t)).^2);
if F.nSizes < 2 || m < 3, return; end
[p1, S1] = polyfit(x, t, 1);
F.lin_c1  = p1(1);
F.lin_c0  = p1(2);
F.lin_R2  = 1 - S1.normr^2 / tss;
F.lin_res = max(abs(t - polyval(p1, x)) ./ t);
q = polyfit(log(N), log(t), 1);
F.exponent = q(1);

if F.nSizes < 3 || m < 4, return; end
[p2, S2] = polyfit(x, t, 2);
F.quad_c2  = p2(1);
F.quad_c1  = p2(2);
F.quad_c0  = p2(3);
F.quad_R2  = 1 - S2.normr^2 / tss;
F.quad_res = max(abs(t - polyval(p2, x)) ./ t);
% standard error of c2 from the covariance of the coefficients (polyfit doc)
Rinv = S2.R \ eye(3);
C = (Rinv * Rinv') * S2.normr^2 / S2.df;
F.quad_c2_se = sqrt(C(1, 1));

% F-test of the nested fits, F(1, m - 3): upper tail through betainc, which
% needs no toolbox (fcdf is in the Statistics Toolbox)
df = m - 3;
F.Fstat = max(S1.normr^2 - S2.normr^2, 0) / (S2.normr^2 / df);
F.p = betainc(df / (df + F.Fstat), df / 2, 1 / 2);
F.share = F.quad_c2 / polyval(p2, 1);

if F.p >= alpha
  F.verdict = "linear";
elseif F.share >= shareMin
  F.verdict = "quadratic term needed";
elseif F.share <= -shareMin
  F.verdict = "concave";
else
  F.verdict = "linear (small curvature)";
end
end

function say(fids, varargin)
% fprintf the same text to every file identifier in fids
for f = fids
  fprintf(f, varargin{:});
end
end

function section(fids, title)
say(fids, '\n%s\n%s\n', title, repmat('=', 1, numel(title)));
end
