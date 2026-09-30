% THESIS_TABLES  The tables of Chapter 5, and the checks behind the claims of
% its text, computed from the result CSVs. Runs no experiment: it only reads
% the files written by scaling_accuracy.m and d2_consistency.m.
%
% The figures of the chapter read their CSVs directly; the tables are typed
% in the .tex, so this script writes, for every table, each cell at full
% precision next to the value as it must be printed and the run it comes
% from, to cross-reference the thesis against the data.
%
% Output (written to OUT_DIR):
%   tables.txt            the report: one section per table, then the checks
%   tab_random_error.tex  body rows of tab:exp-random-error
%   tab_d2.tex            body rows of tab:exp-d2
%   tab_times.tex         body rows of tab:exp-times
% The .tex files hold the rows between \midrule and \bottomrule, as printed
% in the report, so they can be \input by the tables.
%
% Printed values: errors and distances with two significant digits,
% $m \cdot 10^{e}$ (sprintf '%.1e' rounding), exact zeros as $0$; times in
% seconds with two decimals.
%
% Run from anywhere:  run('experiments/thesis_tables.m')

clear; clc;

%% ------------------------------------------------------------------------
%  Configuration (edit here)
%  ------------------------------------------------------------------------

TOL = 1e-8;                              % pass threshold of scaling_accuracy.m

HERE    = fileparts(mfilename('fullpath'));
RES     = fullfile(HERE, 'results');
PGF     = fullfile(RES, 'scaling', 'pgf');
OUT_REL = fullfile('results', 'tables');   % printed in the log (no local paths)
OUT_DIR = fullfile(HERE, OUT_REL);
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

raw   = readtable(fullfile(RES, 'scaling', 'scaling_raw.csv'), 'TextType', 'string');
sweep = readtable(fullfile(RES, 'scaling', 'sweep_raw.csv'),   'TextType', 'string');
D2    = readtable(fullfile(RES, 'd2', 'd2_raw.csv'));

fid = fopen(fullfile(OUT_DIR, 'tables.txt'), 'w');
fprintf(fid, 'Chapter 5 tables and checks, written by thesis_tables.m on %s\n', ...
  char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
logFile = fullfile(RES, 'run_all.log');
if exist(logFile, 'file')
  L = readlines(logFile);
  fprintf(fid, 'results from: %s\n', L(1));
end
fprintf(fid, ['Each value: full precision (%%.17g), then as printed in the ' ...
  'thesis, then the run it comes from.\n']);

%% ------------------------------------------------------------------------
%  tab:exp-random-error: largest error and unitarity defect, random states
%  of the scaling grid (the data of fig:exp-error)
%  ------------------------------------------------------------------------

rows = {'bullock',    'complex', 'Bullock construction'; ...
        'bullock',    'real',    ''; ...
        'qdkptree',   'complex', 'Householder QdKP-Tree'; ...
        'qdkptree',   'real',    ''; ...
        'qdkpgivens', 'real',    'Givens QdKP-Tree'};
section(fid, 'tab:exp-random-error (source: results/scaling/pgf/error_<alg>_<ensemble>.csv)');
tex = {};
for r = 1:size(rows, 1)
  f = sprintf('error_%s_%s.csv', rows{r, 1}, rows{r, 2});
  T = readtable(fullfile(PGF, f));
  [e, ie] = max(T.error);
  [u, iu] = max(T.unitarity);
  fprintf(fid, '%s, %s: %d runs\n', rows{r, 1}, rows{r, 2}, height(T));
  fprintf(fid, '  max error      %-24s -> %-22s at %s\n', exact(e), texSci(e), gridRun(T, ie));
  fprintf(fid, '  max unitarity  %-24s -> %-22s at %s\n', exact(u), texSci(u), gridRun(T, iu));
  tex{end+1} = sprintf('    %-25s & %-7s & %-20s & %s \\\\', ...
    rows{r, 3}, rows{r, 2}, texSci(e), texSci(u)); %#ok<SAGROW>
end
writeRows(fullfile(OUT_DIR, 'tab_random_error.tex'), tex);
fprintf(fid, '\nrows (tab_random_error.tex):\n');
fprintf(fid, '%s\n', tex{:});

%% ------------------------------------------------------------------------
%  tab:exp-d2: the qudit loaders against the qubit KP-Tree at d = 2
%  ------------------------------------------------------------------------

section(fid, 'tab:exp-d2 (source: results/d2/d2_raw.csv)');
fprintf(fid, '%d runs: n = %d..%d, %d trials per n\n', height(D2), min(D2.n), max(D2.n), ...
  height(D2) / numel(unique(D2.n)));
gKP = gateCell(D2.gates_kp,     2.^D2.n - 1, '$2^n - 1$');
gG  = gateCell(D2.gates_givens, 2.^D2.n - 1, '$2^n - 1$');
gH  = gateCell(D2.gates_house,  2.^D2.n,     '$2^n$');
fprintf(fid, 'gates: KP-Tree %s, Givens %s, Householder %s (checked on every run)\n', gKP, gG, gH);

cols = {'err_kp', 'err_givens', 'err_house', 'givens_col0', 'givens_U', 'house_col0', 'house_U'};
v = struct();
for c = 1:numel(cols)
  [v.(cols{c}), i] = max(D2.(cols{c}));
  fprintf(fid, '  max %-12s %-24s -> %-22s at n=%d trial=%d\n', cols{c}, ...
    exact(v.(cols{c})), texSci(v.(cols{c})), D2.n(i), D2.trial(i));
end
tex = {sprintf('    %-21s & %-9s & %-20s & %-22s & %s \\\\', 'KP-Tree (qubit)', gKP, ...
         texSci(v.err_kp), '---', '---'), ...
       sprintf('    %-21s & %-9s & %-20s & %-22s & %s \\\\', 'Givens QdKP-Tree', gG, ...
         texSci(v.err_givens), texSci(v.givens_col0), texSci(v.givens_U)), ...
       sprintf('    %-21s & %-9s & %-20s & %-22s & %s \\\\', 'Householder QdKP-Tree', gH, ...
         texSci(v.err_house), texSci(v.house_col0), texSci(v.house_U))};
writeRows(fullfile(OUT_DIR, 'tab_d2.tex'), tex);
fprintf(fid, '\nrows (tab_d2.tex):\n');
fprintf(fid, '%s\n', tex{:});

%% ------------------------------------------------------------------------
%  tab:exp-times: build and simulation times at the largest N of each d
%  (the last points of fig:exp-build and fig:exp-sim)
%  ------------------------------------------------------------------------

ALGS  = {'bullock', 'qdkptree', 'qdkpgivens'};   % column order of the table
NAMES = {'Bullock', 'Householder', 'Givens'};
f  = dir(fullfile(PGF, 'scaling_bullock_d*.csv'));
DS = sort(cellfun(@(s) sscanf(s, 'scaling_bullock_d%d.csv'), {f.name}));

section(fid, 'tab:exp-times (source: results/scaling/pgf/scaling_<alg>_d<d>.csv, last row)');
tex = {};
top = struct('alg', {}, 'd', {}, 'N', {}, 'nGates', {}, 'b', {}, 'bmin', {}, 'bmax', {}, ...
  's', {}, 'smin', {}, 'smax', {}, 'bslope', {}, 'sslope', {});
for d = DS
  cells = cell(2, numel(ALGS));
  Ns = zeros(1, numel(ALGS));
  for a = 1:numel(ALGS)
    T = readtable(fullfile(PGF, sprintf('scaling_%s_d%d.csv', ALGS{a}, d)));
    k = height(T);
    top(end+1) = struct('alg', ALGS{a}, 'd', d, 'N', T.N(k), 'nGates', T.nGates(k), ...
      'b', T.build_med(k), 'bmin', T.build_min(k), 'bmax', T.build_max(k), ...
      's', T.sim_med(k), 'smin', T.sim_min(k), 'smax', T.sim_max(k), ...
      'bslope', T.build_slope(k), 'sslope', T.sim_slope(k)); %#ok<SAGROW>
    Ns(a) = T.N(k);
    cells{1, a} = sprintf('%.2f', T.build_med(k));
    cells{2, a} = sprintf('%.2f', T.sim_med(k));
    fprintf(fid, 'd=%d N=%-6d %-11s build median %-24s -> %-5s  sim median %-24s -> %s\n', ...
      d, T.N(k), NAMES{a}, exact(T.build_med(k)), cells{1, a}, exact(T.sim_med(k)), cells{2, a});
  end
  if any(Ns ~= Ns(1))
    fprintf(fid, '  WARNING: the largest N differs between the algorithms at d=%d: %s\n', ...
      d, mat2str(Ns));
  end
  tex{end+1} = sprintf('    %-4d & %-6d & %-7s & %-11s & %-9s & %-7s & %-11s & %s \\\\', ...
    d, Ns(1), cells{1, :}, cells{2, :}); %#ok<SAGROW>
end
writeRows(fullfile(OUT_DIR, 'tab_times.tex'), tex);
fprintf(fid, '\nrows (tab_times.tex):\n');
fprintf(fid, '%s\n', tex{:});

%% ------------------------------------------------------------------------
%  Checks behind the text
%  ------------------------------------------------------------------------

section(fid, 'Checks behind the text');
keep = {'algorithm', 'ensemble', 'd', 'n', 'N', 'trial', 'nGates', 'error', 'unitarity'};
G = raw(:, keep);    G.set = repmat("grid", height(G), 1);
S = sweep(:, keep);  S.set = repmat("equal-N", height(S), 1);
A = [G; S];

% every random run passes (error and defect below TOL; NaN fails)
pass = A.error < TOL & A.unitarity < TOL;
fprintf(fid, 'random runs: %d in the scaling grid, %d at equal N; %d fail (TOL = %g)\n', ...
  height(G), height(S), nnz(~pass), TOL);
for k = find(~pass)'
  fprintf(fid, '  FAIL %s error %s unitarity %s\n', setRun(A, k), exact(A.error(k)), exact(A.unitarity(k)));
end

% largest error and defect, over the grid alone and over every random run
for sub = {"grid", "all"}
  if sub{1} == "grid", B = A(A.set == "grid", :); else, B = A; end
  [e, ie] = max(B.error);
  [u, iu] = max(B.unitarity);
  m = max(e, u);
  fprintf(fid, 'largest over %-4s: error %s at %s\n', sub{1}, exact(e), setRun(B, ie));
  fprintf(fid, '                  unitarity %s at %s\n', exact(u), setRun(B, iu));
  fprintf(fid, '                  TOL / largest = 10^%.2f (orders of magnitude below the threshold)\n', ...
    log10(TOL / m));
end

% gate counts against Chapter 4: (N-1)/(d-1) + 1 reflectors and phase, N - 1 rotations
ref = (A.N - 1) ./ (A.d - 1) + 1;
isG = A.algorithm == "qdkpgivens";
ref(isG) = A.N(isG) - 1;
bad = find(A.nGates ~= ref);
fprintf(fid, 'gate counts: %d of %d random runs differ from the formulas of Chapter 4\n', ...
  numel(bad), height(A));
for k = bad'
  fprintf(fid, '  %s: %d gates, expected %d\n', setRun(A, k), A.nGates(k), ref(k));
end

% fig:exp-error caption: exact zeros are not drawn on the log axis
fprintf(fid, 'fig:exp-error, runs not drawn (value exactly 0):\n');
for r = 1:size(rows, 1)
  T = readtable(fullfile(PGF, sprintf('error_%s_%s.csv', rows{r, 1}, rows{r, 2})));
  fprintf(fid, '  %-10s %-7s error = 0: %d, unitarity = 0: %d\n', rows{r, 1}, rows{r, 2}, ...
    nnz(T.error == 0), nnz(T.unitarity == 0));
end

% sec:exp-correctness-d2
g = ~isnan(D2.gatewise_same_qubits);
fprintf(fid, ['d = 2: Givens vs KP-Tree, max distance on |0..0> %s, on X %s; ' ...
  'gate by gate (n <= %d): same qubits on %d of %d runs, max ||M_G - M_KP|| %s\n'], ...
  exact(max(D2.givens_col0)), exact(max(D2.givens_U)), max(D2.n(g)), ...
  nnz(D2.gatewise_same_qubits(g) == 1), nnz(g), exact(max(D2.gatewise_max_diff(g))));
fprintf(fid, ['d = 2: Householder vs KP-Tree, largest |house_col0 - err_house| %s; ' ...
  'distance on X from %s to %s\n'], exact(max(abs(D2.house_col0 - D2.err_house))), ...
  exact(min(D2.house_U)), exact(max(D2.house_U)));

% todo notes: spread of the three runs, build time per gate and last slopes at
% the largest N of each d
fprintf(fid, 'largest N of each d: spread (max/min - 1), build time per gate, slope from the previous size\n');
for t = top
  fprintf(fid, '  %-10s d=%d N=%-6d build spread %5.1f%%  sim spread %5.1f%%  build/gate %.4f ms  slopes build %.2f sim %.2f\n', ...
    t.alg, t.d, t.N, 100 * (t.bmax / t.bmin - 1), 100 * (t.smax / t.smin - 1), ...
    1e3 * t.b / t.nGates, t.bslope, t.sslope);
end

% build slope between the last two sizes on complex states (the times of the
% table and figures are on real states)
fprintf(fid, 'complex states, build slope between the last two sizes (median of the trials):\n');
C = raw(raw.ensemble == "complex", :);
for alg = unique(C.algorithm)'
  for d = unique(C.d(C.algorithm == alg))'
    Cd = C(C.algorithm == alg & C.d == d, :);
    N = unique(Cd.N);
    t = arrayfun(@(x) median(Cd.buildTime(Cd.N == x)), N);
    fprintf(fid, '  %-10s d=%d N=%-6d %.2f\n', alg, d, N(end), ...
      log(t(end) / t(end-1)) / log(N(end) / N(end-1)));
  end
end

fclose(fid);
fprintf('Wrote %s\n', fullfile(OUT_REL, 'tables.txt'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'tab_random_error.tex'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'tab_d2.tex'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'tab_times.tex'));

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

function s = texSci(x)
% x as printed in the tables: two significant digits, $m \cdot 10^{e}$, with
% the rounding of sprintf('%.1e'); $m$ alone for 1 <= |x| < 10, $0$ for 0.
if x == 0, s = '$0$'; return; end
if isnan(x), s = 'NaN'; return; end
tok = regexp(sprintf('%.1e', x), '^(-?\d\.\d)e([+-]\d+)$', 'tokens', 'once');
e = str2double(tok{2});
if e == 0
  s = sprintf('$%s$', tok{1});
else
  s = sprintf('$%s \\cdot 10^{%d}$', tok{1}, e);
end
end

function s = exact(x)
% x with enough digits to recover the double of the CSV
s = sprintf('%.17g', x);
end

function s = gridRun(T, i)
% the run of row i of a pgf error file
s = sprintf('N=%d d=%d n=%d trial=%d', T.N(i), T.d(i), T.n(i), T.trial(i));
end

function s = setRun(A, i)
% the run of row i of the combined scaling-grid and equal-N table
s = sprintf('%s %s %s N=%d d=%d n=%d trial=%d', A.set(i), A.algorithm(i), ...
  A.ensemble(i), A.N(i), A.d(i), A.n(i), A.trial(i));
end

function s = gateCell(counts, expected, formula)
% the formula if every run has that many gates, a flag otherwise
if all(counts == expected)
  s = formula;
else
  s = sprintf('MISMATCH(%d)', nnz(counts ~= expected));
end
end

function section(fid, title)
fprintf(fid, '\n%s\n%s\n', title, repmat('=', 1, numel(title)));
end

function writeRows(file, rows)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', rows{:});
fclose(fid);
end
