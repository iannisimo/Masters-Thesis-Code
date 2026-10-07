% FIT_SCALING  How fast do the times of the scaling grid grow with N? For
% every time (build, simulation), algorithm and qudit dimension d, fits the
% runs with N >= FIT_MIN_N by a polynomial of degree 1 and one of degree 2 in
% N, and tests whether the quadratic term is needed. Runs no experiment: it
% only reads results/scaling/scaling_raw.csv, written by scaling_accuracy.m.
%
% Expected: the build time of the Bullock construction grows as N^2 in our
% implementation (each of its O(N) reflectors reads and updates the whole
% state vector, O(N) per reflector); those of the two QdKP-Trees grow as N.
% The simulation time should grow as N^2 for all three (O(N) gates, O(N)
% each once the simulation is sparse), so it checks the method: a fit that
% does not find the quadratic term there cannot be trusted to reject it for
% the build times.
%
% The fits are done twice, with two measures of the misfit of a run:
%   - absolute, t - fit, in seconds: ordinary least squares, what polyfit
%     computes;
%   - relative, (t - fit) / t: weighted least squares with weights 1 / t^2.
% The relative fit is the one that carries the conclusion. The runs of a
% size differ by a percentage of their time, at every size, so the relative
% misfit has about the same spread everywhere, as the F-test assumes. The
% absolute misfit does not: the largest sizes decide the fit alone, and the
% runs of the small sizes enter the F-test with almost no residual, so its p
% is far from the true chance. The absolute fit is kept to show what that
% does to the verdicts (section "Absolute against relative misfit").
%
% Common to both:
%   - every trial is a point (N_TRIALS per size), so the residuals include
%     the spread of the runs;
%   - N is scaled to x = N / N_max, so the coefficients are in seconds: c2
%     and c1 are the quadratic and linear contributions at the largest N, c0
%     the fixed cost;
%   - the quadratic term is needed when the F-test of the degree-2 against
%     the degree-1 fit gives p < ALPHA (the same as a t-test on c2) and c2 is
%     positive; a significant negative c2 is a concave curve.
% FIT_MIN_N is set per time. The build times are fitted over every size: c0
% takes the fixed cost of the small sizes. The simulation times start at
% N = 512, since QCLAB simulates with dense matrices below it and one fit
% would mix two regimes. A (time, algorithm, d) with fewer than 3 sizes from
% FIT_MIN_N on gets no quadratic fit ("too few sizes").
%
% Output (written to OUT_DIR):
%   fit_scaling.txt   the report: one table per (misfit, time), the fits
%                     that disagree with EXPECTED, and the two misfits side
%                     by side (also on the console)
%   fit_scaling.csv   every coefficient, one row per (misfit, time,
%                     algorithm, d)
%   pgf/fit_<alg>_d<d>_runs.csv, pgf/fit_<alg>_d<d>_fits.csv
%                     the build-time curves of SHOW, for the figures of the
%                     appendix: the runs with their misfits, and both
%                     degree-2 fits (see exportFits)
%
% How to read the report. Each line is one curve of Figure 5.2 (build) or
% 5.3 (simulation): one algorithm in the panel of one d.
%   sizes, N range  the distinct N fitted (N_TRIALS runs each) and their span
%   degree 1        res of t = c1 x + c0
%   degree 2        c2, c1, c0 of t = c2 x^2 + c1 x + c0 (s, at N_max), the
%                   standard error of c2 in % of c2 (+-c2), and res
%   res             largest relative misfit |t - fit| / t over the runs
%   p               chance that the parabola beats the line only by luck
%   share           c2 / fitted time at N_max: how much of the largest time
%                   comes from the N^2 term
%   verdict         "quadratic term needed": p < ALPHA and c2 > 0;
%                   "concave": p < ALPHA and c2 < 0; "linear": p >= ALPHA;
%                   "too few sizes": fewer than 3 sizes, no parabola;
%                   "(NOT AS EXPECTED)" when the verdict contradicts EXPECTED.
% The verdict says whether there is an N^2 term; share says how much of the
% time it is, and +-c2 how well it is known.
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
MISFITS   = {'relative', 'absolute'};    % the first one carries the conclusion

ALGS  = {'bullock', 'qdkptree', 'qdkpgivens'};
NAMES = {'Bullock', 'Householder', 'Givens'};
% build-time curves exported for the figures of the appendix: every
% algorithm and every d of the grid
SHOW      = [repelem(ALGS(:), 4, 1), num2cell(repmat((2:5)', numel(ALGS), 1))];
% misfit range drawn in those figures, in % (their ymin/ymax in
% tex/part/exp-fits.tex): runs beyond it are exported as edge markers
MISFIT_MAX = 100;
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
say(out, ['%s ensemble, every trial a point; x = N / N_max, so c2, c1, c0 are in s ' ...
  '(contributions at the largest N)\n'], ENSEMBLE);
say(out, ['misfit: relative (t - fit) / t, which carries the conclusion, and absolute ' ...
  't - fit (polyfit), for comparison\n']);
say(out, ['quadratic term needed: F-test of degree 2 against degree 1 with p < %g, ' ...
  'and c2 > 0\n'], ALPHA);
say(out, 'res = largest relative misfit |t - fit| / t\n');

%% ------------------------------------------------------------------------
%  Fits
%  ------------------------------------------------------------------------

fits = {};
for mf = MISFITS
  relative = strcmp(mf{1}, 'relative');
  for tm = TIMES
    minN = FIT_MIN_N.(tm{1});
    section(out, char(sprintf('%s misfit, %s, runs with N >= %d (expected degree: %s)', ...
      mf{1}, tm{1}, minN, ...
      strjoin(compose('%s %d', string(NAMES(:)), EXPECTED.(tm{1})(:)), ', '))));
    say(out, '%-11s %2s %5s %13s | %-8s | %-41s | %7s %6s  %s\n', '', '', '', '', ...
      'degree 1', 'degree 2', '', '', '');
    say(out, '%-11s %2s %5s %13s | %8s | %8s %8s %8s %6s %8s | %7s %6s  %s\n', ...
      'algorithm', 'd', 'sizes', 'N range', 'res', 'c2 [s]', 'c1 [s]', 'c0 [s]', ...
      '+-c2', 'res', 'p', 'share', 'verdict');
    for a = 1:numel(ALGS)
      for d = unique(raw.d(raw.algorithm == ALGS{a}))'
        R = raw(raw.algorithm == ALGS{a} & raw.d == d & raw.N >= minN, :);
        F = fitTimes(R.N, R.(tm{1}), ALPHA, relative);
        expected = EXPECTED.(tm{1})(a);
        agrees = (expected == 2 && F.verdict == "quadratic term needed") || ...
                 (expected == 1 && F.verdict == "linear");
        say(out, '%-11s %2d %5d %13s | %7.1f%% | %8.3g %8.3g %8.3g %5.0f%% %7.1f%% | %7.1e %5.0f%%  %s%s\n', ...
          NAMES{a}, d, F.nSizes, sprintf('%d-%d', F.Nmin, F.Nmax), 100 * F.lin_res, ...
          F.quad_c2, F.quad_c1, F.quad_c0, 100 * F.quad_c2_se / abs(F.quad_c2), ...
          100 * F.quad_res, F.p, 100 * F.share, F.verdict, ...
          repmat(' (NOT AS EXPECTED)', 1, double(~agrees)));
        fits{end+1} = [table(string(mf{1}), string(tm{1}), string(ALGS{a}), d, expected, agrees, ...
          'VariableNames', {'misfit', 'time', 'algorithm', 'd', 'expected', 'agrees'}), ...
          struct2table(F)]; %#ok<SAGROW>
      end
    end
  end
end

fits = vertcat(fits{:});
writetable(fits, fullfile(OUT_DIR, 'fit_scaling.csv'));

section(out, 'Against the expected degrees');
for mf = MISFITS
  M = fits(fits.misfit == mf{1}, :);
  say(out, '%s misfit: %d of %d fits as expected\n', mf{1}, nnz(M.agrees), height(M));
  for k = find(~M.agrees)'
    say(out, '  %s %s d=%d: expected degree %d, %s\n', M.time(k), ...
      M.algorithm(k), M.d(k), M.expected(k), M.verdict(k));
  end
end

% the two misfits on the same curve: where the verdicts differ, the F-test
% of the absolute fit is the one to distrust
section(out, 'Absolute against relative misfit');
say(out, '%-9s %-11s %2s | %-31s | %-31s\n', '', '', '', 'relative', 'absolute');
say(out, '%-9s %-11s %2s | %7s %6s  %-16s | %7s %6s  %-16s\n', 'time', 'algorithm', 'd', ...
  'p', 'share', 'verdict', 'p', 'share', 'verdict');
rel = fits(fits.misfit == "relative", :);
abso = fits(fits.misfit == "absolute", :);   % same rows in the same order
for k = 1:height(rel)
  differ = rel.verdict(k) ~= abso.verdict(k);
  say(out, '%-9s %-11s %2d | %7.1e %5.0f%%  %-16s | %7.1e %5.0f%%  %-16s%s\n', ...
    rel.time(k), NAMES{strcmp(ALGS, rel.algorithm(k))}, rel.d(k), ...
    rel.p(k), 100 * rel.share(k), shortVerdict(rel.verdict(k)), ...
    abso.p(k), 100 * abso.share(k), shortVerdict(abso.verdict(k)), ...
    repmat('  <- differ', 1, double(differ)));
end

fclose(fid);
fprintf('\nWrote %s\n', fullfile(OUT_REL, 'fit_scaling.txt'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'fit_scaling.csv'));

% the curves of the appendix figure
PGF_DIR = fullfile(OUT_DIR, 'pgf');
if ~exist(PGF_DIR, 'dir'), mkdir(PGF_DIR); end
for k = 1:size(SHOW, 1)
  [alg, d] = SHOW{k, :};
  R = raw(raw.algorithm == alg & raw.d == d & raw.N >= FIT_MIN_N.buildTime, :);
  curve = fits.time == "buildTime" & fits.algorithm == alg & fits.d == d;
  name = sprintf('fit_%s_d%d', alg, d);
  exportFits(PGF_DIR, name, R.N, R.buildTime, ...
    fits(curve & fits.misfit == "relative", :), fits(curve & fits.misfit == "absolute", :), MISFIT_MAX);
  fprintf('Wrote %s_{runs,fits}.csv\n', fullfile(OUT_REL, 'pgf', name));
end

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

function F = fitTimes(N, t, alpha, relative)
% Degree-1 and degree-2 fits of the times t against x = N / max(N), on the
% relative misfit if relative is true and on the absolute one otherwise,
% and the F-test of the quadratic term. Fields the data cannot support are
% NaN.
N = N(:);
t = t(:);
m = numel(t);
F.nSizes = numel(unique(N));
F.nRuns  = m;
F.Nmin   = min(N);
F.Nmax   = max(N);
[F.lin_c1, F.lin_c0, F.lin_res] = deal(NaN);
[F.quad_c2, F.quad_c2_se, F.quad_c1, F.quad_c0, F.quad_res] = deal(NaN);
[F.Fstat, F.p, F.share] = deal(NaN);
F.verdict = "too few sizes";

x = N / F.Nmax;
if F.nSizes < 2 || m < 3, return; end
L = fitPoly(x, t, 1, relative);
F.lin_c1  = L.coef(2);
F.lin_c0  = L.coef(1);
F.lin_res = L.maxMisfit;

if F.nSizes < 3 || m < 4, return; end
Q = fitPoly(x, t, 2, relative);
F.quad_c2    = Q.coef(3);
F.quad_c1    = Q.coef(2);
F.quad_c0    = Q.coef(1);
F.quad_c2_se = Q.se(3);
F.quad_res   = Q.maxMisfit;

% F-test of the nested fits, F(1, m - 3): upper tail through betainc, which
% needs no toolbox (fcdf is in the Statistics Toolbox)
df = m - 3;
F.Fstat = max(L.rss - Q.rss, 0) / (Q.rss / df);
F.p = betainc(df / (df + F.Fstat), df / 2, 1 / 2);
F.share = F.quad_c2 / sum(Q.coef);       % the fitted time at x = 1 is c0 + c1 + c2

if F.p >= alpha
  F.verdict = "linear";
elseif F.quad_c2 > 0
  F.verdict = "quadratic term needed";
else
  F.verdict = "concave";
end
end

function S = fitPoly(x, t, degree, relative)
% Least-squares fit of t by a polynomial of the given degree in x. With
% relative true it minimizes the relative misfit sum(((t - fit) ./ t).^2):
% dividing every row of the problem by its t turns it into an ordinary
% least-squares problem. With relative false it minimizes the absolute
% misfit sum((t - fit).^2), as polyfit does.
%   coef       coefficients [c0; c1; ...], in increasing degree
%   se         their standard errors
%   rss        the minimized sum of squared misfits (relative or absolute)
%   maxMisfit  the largest relative misfit |t - fit| / t, for both
X = x .^ (0:degree);                     % columns 1, x, x^2, ...
if relative
  w = 1 ./ t;                            % row k divided by t(k)
else
  w = ones(size(t));
end
A = X .* w;
b = t .* w;
[Qa, Ra] = qr(A, 0);
S.coef = Ra \ (Qa' * b);
r = b - A * S.coef;                      % the misfits being minimized
S.rss = sum(r .^ 2);
S.maxMisfit = max(abs(t - X * S.coef) ./ t);
% covariance of the coefficients: inv(A' A) times the variance of the misfit
Rinv = Ra \ eye(degree + 1);
S.se = sqrt(diag(Rinv * Rinv') * S.rss / (numel(t) - degree - 1));
end

function exportFits(dir, name, N, t, rel, abso, misfitMax)
% Two CSVs for the appendix figure, from the runs N, t of one curve and the
% rows rel, abso of its relative and absolute fits in the table of fits:
%   <name>_runs.csv  every run: N, its time per amplitude t / N in us, and
%                    its relative misfit (t - fit) / t in % under the
%                    degree-1 fit on the relative misfit (misfit_lin) and
%                    the degree-2 fit of each misfit (misfit_rel, misfit_abs);
%                    for each of these, <col>_in is the misfit where it lies
%                    within +-misfitMax (NaN elsewhere), and <col>_hi and
%                    <col>_lo are +misfitMax and -misfitMax where it lies
%                    above or below that range (NaN elsewhere), so that the
%                    figure can mark those runs on the edge of the panel
%   <name>_fits.csv  divided by N, in us, on 200 values of N evenly spaced
%                    in log N over the sizes: the degree-1 (lin) and degree-2
%                    (rel) fits on the relative misfit, and the degree-2 fit
%                    on the absolute misfit (abs); NaN where a fit is not
%                    positive, which a log axis cannot draw
Nmax = max(N);
fitAt = @(F, n) F.quad_c0 + F.quad_c1 * (n / Nmax) + F.quad_c2 * (n / Nmax) .^ 2;
linAt = @(F, n) F.lin_c0 + F.lin_c1 * (n / Nmax);
runs = table(N, 1e6 * t ./ N, 100 * (t - linAt(rel, N)) ./ t, 100 * (t - fitAt(rel, N)) ./ t, ...
  100 * (t - fitAt(abso, N)) ./ t, 'VariableNames', {'N', 'perN', 'misfit_lin', 'misfit_rel', 'misfit_abs'});
for col = {'misfit_lin', 'misfit_rel', 'misfit_abs'}
  v = runs.(col{1});
  inside = v;              inside(abs(v) > misfitMax) = NaN;
  hi = NaN(size(v));       hi(v > misfitMax) = misfitMax;
  lo = NaN(size(v));       lo(v < -misfitMax) = -misfitMax;
  runs.([col{1} '_in']) = inside;
  runs.([col{1} '_hi']) = hi;
  runs.([col{1} '_lo']) = lo;
end
writetable(runs, fullfile(dir, [name '_runs.csv']));

n = logspace(log10(min(N)), log10(Nmax), 200)';
linN = 1e6 * linAt(rel, n) ./ n;
relN = 1e6 * fitAt(rel, n) ./ n;
absN = 1e6 * fitAt(abso, n) ./ n;
linN(linN <= 0) = NaN;
relN(relN <= 0) = NaN;
absN(absN <= 0) = NaN;
writetable(table(n, linN, relN, absN, 'VariableNames', {'N', 'lin', 'rel', 'abs'}), ...
  fullfile(dir, [name '_fits.csv']));
end

function v = shortVerdict(v)
% the verdict in a word, for the side-by-side table
if v == "quadratic term needed", v = "quadratic"; end
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
