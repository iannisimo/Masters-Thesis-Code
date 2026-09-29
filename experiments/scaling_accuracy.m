% SCALING_ACCURACY  Runtime scaling and correctness of the state-preparation
% functions in code/ (QR_bullock.m, QdKPTree_Householder.m,
% QdKPTree_Givens.m), for Chapter 5. The qubit KP-Tree (KPTree.m) is not part
% of these tests: it is compared with the qudit loaders at d = 2 only, in
% d2_consistency.m.
%
%   1a  build time vs N = d^n at fixed d, with the local slope between
%       consecutive sizes
%   1b  simulation time vs N at fixed d, same slope
%   2a  error ||C|0> - psi||_2 vs N, on random real and complex states
%   2b  pass/fail on edge-case states (zero entries, basis states, sparse, ...),
%       and the number of gates of each circuit that are the identity
%   2c  unitarity of every circuit C (random and edge-case runs): with X a
%       random N x UNIT_K matrix with orthonormal columns and Y = C X,
%       defect = ||Y'Y - I||_2 (C itself is never formed, N x N may not fit).
%       A small ||C|0> - psi|| does not imply a unitary C (see the Givens
%       child-index bug), so both are checked.
%   3   fixed-N sweep: for every N in SWEEP_NS and every d with d^n = N,
%       n >= SWEEP_MIN_n, the same runs at equal state size
%
% Every run also records the size of its circuit (circuitStats): nGates and
% depth (greedy layering).
%
% Part 1 and 2a share the same runs: for every algorithm and every d in DS,
% n grows from MIN_n until d^n > N_MAX, or until the median build +
% simulation time of one trial exceeds TIME_BUDGET (that n is still kept).
% Every algorithm is timed on TIMING_ENSEMBLE, so that all of them get the
% same states; the other ensemble only enters the error and unitarity data.
%
% Output (written to OUT_DIR; see code/experiments/DATA.md):
%   scaling_raw.csv   one row per random trial
%   edge_raw.csv      one row per (edge case, d, n, algorithm)
%   sweep_raw.csv     one row per fixed-N trial (3)
%   sweep_summary.csv mean over the trials of each (N, d, algorithm)
%   pgf/              the same data laid out for pgfplots (see writePgfData),
%                     rewritten from scratch on every run
%   pgf/sweep/        part 3 (see writeSweepPgf)
%
% Run from anywhere:  run('code/experiments/scaling_accuracy.m')

clear; clc;

%% ------------------------------------------------------------------------
%  Configuration (edit here)
%  ------------------------------------------------------------------------

DS          = [2, 3, 4, 5];              % qudit dimensions for the scaling runs
ALGORITHMS  = {'bullock', 'qdkptree', 'qdkpgivens'};
REAL_ONLY   = {'qdkpgivens'};            % these only get real states
QUBIT_ONLY  = {};                        % these only run with d = 2
ENSEMBLES   = {'complex', 'real'};       % random-state ensembles (2a)
TIMING_ENSEMBLE = 'real';                % the times (1a, 1b, 3) of every algorithm
MIN_n       = 2;                         % smallest number of qudits
N_MAX       = 2^14;                      % never go above this state size
TIME_BUDGET = 60;                        % s; stop growing n once a trial takes longer
N_TRIALS    = 1;                         % random states per (algorithm, d, n, ensemble)
BASE_SEED   = 0;
TOL         = 1e-8;                      % error / unitarity defect above this is a FAIL
UNIT_K      = 4;                         % columns of the unitarity probe (2c)
ID_TOL      = 1e-12;                     % a gate is the identity if ||M - I||_1 < ID_TOL

% edge cases (2b): real states, each run on every (d, n) row of EDGE_DN
EDGE_DN     = [2 4; 3 3; 4 2; 5 2];
EDGE_CASES  = {'random', 'psi0_zero', 'zero_leaf_node', 'zero_subtree', ...
               'aligned_nodes', 'basis_first', 'basis_last', 'uniform', ...
               'sparse_d', 'sparse_sqrtN'};

% circuit-size columns recorded for every run (see circuitStats)
STATS_NAMES = {'nGates', 'depth'};

% fixed-N sweep (3): 729 = 3^6 = 9^3 = 27^2, 4096 = 2^12 = 4^6 = 8^4 = 16^3 = 64^2,
% 6561 = 3^8 = 9^4 = 81^2, 15625 = 5^6 = 25^3 = 125^2
SWEEP_NS    = [729, 4096, 6561, 15625];  % state sizes; every d with d^n = N is run
SWEEP_MIN_n = 2;                         % skip d with n < SWEEP_MIN_n

HERE    = fileparts(mfilename('fullpath'));
OUT_REL = fullfile('results', 'scaling');   % printed in the log (no local paths)
OUT_DIR = fullfile(HERE, OUT_REL);
addpath(fullfile(HERE, '..'));           % the functions add qclab and util themselves
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

%% ------------------------------------------------------------------------
%  1a, 1b, 2a, 2c: random states, growing n at fixed d
%  ------------------------------------------------------------------------

% warm-up (untimed): the first call of each algorithm pays for JIT/class loading
for a = 1:numel(ALGORITHMS)
  f = algorithmFunction(ALGORITHMS{a});
  [~] = f(2, 2, 0, 0);
end

rows = {};
for a = 1:numel(ALGORITHMS)
  alg = ALGORITHMS{a};
  f   = algorithmFunction(alg);
  ens = supportedEnsembles(alg, ENSEMBLES, REAL_ONLY);

  for d = DS
    if any(strcmp(alg, QUBIT_ONLY)) && d ~= 2, continue; end

    n = MIN_n;
    while d^n <= N_MAX
      N = d^n;
      trialTime = zeros(1, N_TRIALS);          % build + sim, TIMING_ENSEMBLE

      for e = 1:numel(ens)
        isComplex = strcmp(ens{e}, 'complex');
        for trial = 1:N_TRIALS
          % the same seed gives the same psi to every algorithm
          seed = BASE_SEED + trial + 1000 * find(strcmp(ens{e}, ENSEMBLES));
          fprintf('%-10s d=%-2d n=%-2d N=%-6d %-7s trial %d: ', alg, d, n, N, ens{e}, trial);
          [cir, err, bt, st] = f(d, n, double(isComplex), seed);
          unit = unitarityDefect(cir, d, n, UNIT_K);
          cs   = circuitStats(cir, n);
          if strcmp(ens{e}, TIMING_ENSEMBLE), trialTime(trial) = bt + st; end

          rows{end+1} = [{alg, d, n, N, ens{e}, trial, seed}, statsCells(cs), ...
            {bt, st, err, unit}]; %#ok<SAGROW>

          status = 'ok';
          if ~(err <= TOL && unit <= TOL), status = 'FAIL'; end
          fprintf('build %8.3fs, sim %8.3fs, err %.1e, unit %.1e  %s\n', ...
            bt, st, err, unit, status);
        end
      end

      if median(trialTime) > TIME_BUDGET
        fprintf('%-10s d=%-2d: stopping at n=%d (%.1fs per trial > %ds)\n', ...
          alg, d, n, median(trialTime), TIME_BUDGET);
        break
      end
      n = n + 1;
    end
  end
end

raw = cell2table(vertcat(rows{:}), 'VariableNames', [{'algorithm', 'd', 'n', 'N', ...
  'ensemble', 'trial', 'seed'}, STATS_NAMES, {'buildTime', 'simTime', 'error', 'unitarity'}]);
writetable(raw, fullfile(OUT_DIR, 'scaling_raw.csv'));

%% ------------------------------------------------------------------------
%  2b, 2c: edge-case states
%  ------------------------------------------------------------------------

erows = {};
for c = 1:numel(EDGE_CASES)
  for k = 1:size(EDGE_DN, 1)
    d = EDGE_DN(k, 1);
    n = EDGE_DN(k, 2);
    rng(BASE_SEED + 7919 * c + k);
    psi = edgeState(EDGE_CASES{c}, d, n);

    for a = 1:numel(ALGORITHMS)
      alg = ALGORITHMS{a};
      if any(strcmp(alg, QUBIT_ONLY)) && d ~= 2, continue; end
      fprintf('edge %-15s d=%d n=%d %-10s ', EDGE_CASES{c}, d, n, alg);
      [status, err, unit, msg, cs, nId] = runEdgeCase(algorithmFunction(alg), d, n, psi, UNIT_K, TOL, ID_TOL);
      erows{end+1} = [{EDGE_CASES{c}, d, n, alg, status}, statsCells(cs), {nId, err, unit, msg}]; %#ok<SAGROW>
      fprintf('%-5s err %.1e, unit %.1e, %d/%d identity  %s\n', status, err, unit, nId, cs.nGates, msg);
    end
  end
end

edges = cell2table(vertcat(erows{:}), 'VariableNames', [{'edgeCase', 'd', 'n', ...
  'algorithm', 'status'}, STATS_NAMES, {'nIdentity', 'error', 'unitarity', 'message'}]);
writetable(edges, fullfile(OUT_DIR, 'edge_raw.csv'));

%% ------------------------------------------------------------------------
%  3: fixed-N sweep, equal state size across d
%  ------------------------------------------------------------------------

srows = {};
for N = SWEEP_NS
  for d = dimensionsFor(N, SWEEP_MIN_n)
    n = round(log(N) / log(d));
    for a = 1:numel(ALGORITHMS)
      alg = ALGORITHMS{a};
      if any(strcmp(alg, QUBIT_ONLY)) && d ~= 2, continue; end
      isComplex = strcmp(TIMING_ENSEMBLE, 'complex');
      for trial = 1:N_TRIALS
        seed = BASE_SEED + trial;        % same psi for every algorithm
        fprintf('sweep N=%-5d d=%-3d n=%-2d %-10s %-7s trial %d: ', N, d, n, alg, TIMING_ENSEMBLE, trial);
        f = algorithmFunction(alg);
        [cir, err, bt, st] = f(d, n, double(isComplex), seed);
        unit = unitarityDefect(cir, d, n, UNIT_K);
        cs   = circuitStats(cir, n);
        srows{end+1} = [{N, d, n, alg, TIMING_ENSEMBLE, trial, seed}, statsCells(cs), ...
          {bt, st, err, unit}]; %#ok<SAGROW>
        status = 'ok';
        if ~(err <= TOL && unit <= TOL), status = 'FAIL'; end
        fprintf('%5d gates, depth %5d, build %.3fs, sim %.3fs, err %.1e, unit %.1e  %s\n', ...
          cs.nGates, cs.depth, bt, st, err, unit, status);
      end
    end
  end
end

sweep = cell2table(vertcat(srows{:}), 'VariableNames', [{'N', 'd', 'n', 'algorithm', ...
  'ensemble', 'trial', 'seed'}, STATS_NAMES, {'buildTime', 'simTime', 'error', 'unitarity'}]);
writetable(sweep, fullfile(OUT_DIR, 'sweep_raw.csv'));
sweepSummary = groupsummary(sweep, {'N', 'd', 'n', 'algorithm'}, 'mean', ...
  [STATS_NAMES, {'buildTime', 'simTime', 'error', 'unitarity'}]);
writetable(sweepSummary, fullfile(OUT_DIR, 'sweep_summary.csv'));

%% ------------------------------------------------------------------------
%  Save for pgfplots
%  ------------------------------------------------------------------------

% pgf/ holds only what this run writes: no stale files from older layouts
if exist(fullfile(OUT_DIR, 'pgf'), 'dir'), rmdir(fullfile(OUT_DIR, 'pgf'), 's'); end
writePgfData(raw, edges, ALGORITHMS, ENSEMBLES, TIMING_ENSEMBLE, EDGE_CASES, ...
  fullfile(OUT_DIR, 'pgf'));
writeSweepPgf(sweepSummary, ALGORITHMS, fullfile(OUT_DIR, 'pgf', 'sweep'));

fprintf('\nWrote %s\n', fullfile(OUT_REL, 'scaling_raw.csv'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'edge_raw.csv'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'sweep_raw.csv'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'sweep_summary.csv'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'pgf'));


%% ========================================================================
%  Helpers
%  ========================================================================

function f = algorithmFunction(alg)
switch alg
  case 'bullock',    f = @QR_bullock;
  case 'qdkptree',   f = @QdKPTree_Householder;
  case 'qdkpgivens', f = @QdKPTree_Givens;
  case 'kptree',     f = @KPTree;
  otherwise,         error('Unknown algorithm "%s"', alg);
end
end

function ens = supportedEnsembles(alg, ensembles, realOnly)
% The ensembles alg can prepare, in the configured order.
ens = ensembles;
if any(strcmp(alg, realOnly))
  ens = ensembles(strcmp(ensembles, 'real'));
end
end

function defect = unitarityDefect(cir, d, n, k)
% ||Y'Y - I|| for Y = C X, X random N x k with orthonormal columns.
% Zero for a unitary C; C itself is never formed.
N = d^n;
k = min(k, N);
[X, ~] = qr(randn(N, k) + 1i * randn(N, k), 0);
Y = cir.apply('R', 'N', n, X, 0, d);
defect = norm(Y' * Y - eye(k));
end

function psi = edgeState(name, d, n)
% Real test states that stress the special cases of the constructions.
% Random entries are Gaussian, so they carry random signs.
N = d^n;
switch name
  case 'random'                 % control: a generic real state
    psi = randn(N, 1);
  case 'psi0_zero'              % first amplitude zero: node phases b_{k0} = 0
    psi = randn(N, 1);
    psi(1) = 0;
  case 'zero_leaf_node'         % the d children of leaf node 1 (0-based) all zero
    psi = randn(N, 1);
    psi(d+1:2*d) = 0;
  case 'zero_subtree'           % the whole subtree of the root's child 1 zero
    psi = randn(N, 1);
    psi(d^(n-1)+1:2*d^(n-1)) = 0;
  case 'aligned_nodes'          % every leaf node already along e1
    psi = zeros(N, 1);
    psi(1:d:end) = randn(N/d, 1);
  case 'basis_first'            % |0...0>: nothing to do
    psi = zeros(N, 1);
    psi(1) = 1;
  case 'basis_last'             % |d-1 ... d-1>
    psi = zeros(N, 1);
    psi(N) = 1;
  case 'uniform'                % uniform superposition
    psi = ones(N, 1);
  case 'sparse_d'               % d nonzero entries at random positions
    psi = sparseState(N, d);
  case 'sparse_sqrtN'           % ~sqrt(N) nonzero entries at random positions
    psi = sparseState(N, round(sqrt(N)));
  otherwise
    error('Unknown edge case "%s"', name);
end
psi = psi / norm(psi);
end

function psi = sparseState(N, k)
psi = zeros(N, 1);
psi(randperm(N, k)) = randn(k, 1);
end

function [status, err, unit, msg, cs, nId] = runEdgeCase(f, d, n, psi, k, tol, idTol)
% status: 'pass', 'fail' (wrong state or not unitary), 'nan' (NaN in the
% result) or 'error' (the construction threw; msg holds the message).
% cs: circuitStats of the circuit (NaN fields if it was not built).
% nId: number of gates of the circuit that are the identity (NaN if not built).
msg = '';
cs  = circuitStats([], n);
nId = NaN;
try
  [cir, err] = f(d, n, 0, [], psi);
  cs   = circuitStats(cir, n);
  nId  = countIdentityGates(cir, d, idTol);
  unit = unitarityDefect(cir, d, n, k);
  if isnan(err) || isnan(unit)
    status = 'nan';
  elseif err <= tol && unit <= tol
    status = 'pass';
  else
    status = 'fail';
  end
catch ME
  status = 'error';
  err    = NaN;
  unit   = NaN;
  msg    = ME.message;
end
end

function nId = countIdentityGates(cir, d, idTol)
% Gates whose matrix (controls included) is the identity: the reflectors of
% makeHouseholder on a node with nothing to reflect, and a phase gate with
% phase 1. Forms each gate's matrix: only for the small edge-case circuits.
nId = 0;
for k = 1:cir.nbObjects
  M = cir.objects(k).matrix(d);
  if norm(M - eye(size(M)), 1) < idTol, nId = nId + 1; end
end
end

function ds = dimensionsFor(N, minN)
% All d >= 2 such that N = d^n for an integer n >= minN.
ds = [];
for d = 2:N
  n = round(log(N) / log(d));
  if n >= minN && d^n == N
    ds(end+1) = d; %#ok<AGROW>
  end
end
end

function stats = circuitStats(cir, n)
% Gate count and depth of a flat circuit of single-target gates.
% Depth: each gate is placed one layer after the latest gate that touches
% any of its qudits (controls + target). cir = [] gives NaN fields.
if isempty(cir)
  stats = struct('nGates', NaN, 'depth', NaN);
  return
end
layerOfQudit = zeros(1, n);         % last layer used on each qudit
for k = 1:cir.nbObjects
  q = cir.objects(k).qubits + 1;     % 0-based -> 1-based
  layerOfQudit(q) = max(layerOfQudit(q)) + 1;
end
stats.nGates = cir.nbObjects;
stats.depth  = max(layerOfQudit);
end

function c = statsCells(stats)
% circuitStats fields as a cell row, in the order of STATS_NAMES.
c = {stats.nGates, stats.depth};
end

function writeSweepPgf(summary, algs, pgfDir)
% Part 3, for pgfplots:
%   sweep_N<N>.csv   one row per d: d, n, ref_householder = (N-1)/(d-1) + 1,
%                    ref_givens = N - 1, and <alg>_<metric> (mean over the
%                    trials; NaN for algorithms not run at that d)
metrics = {'nGates', 'depth', 'buildTime', 'simTime', 'error', 'unitarity'};
if ~exist(pgfDir, 'dir'), mkdir(pgfDir); end

for N = unique(summary.N)'
  S  = summary(summary.N == N, :);
  d  = unique(S.d);
  T  = table(d, round(log(N) ./ log(d)), (N - 1) ./ (d - 1) + 1, repmat(N - 1, size(d)), ...
    'VariableNames', {'d', 'n', 'ref_householder', 'ref_givens'});
  for a = 1:numel(algs)
    for m = 1:numel(metrics)
      col = NaN(size(d));
      for i = 1:numel(d)
        r = S.d == d(i) & strcmp(S.algorithm, algs{a});
        if any(r), col(i) = S.(['mean_', metrics{m}])(r); end
      end
      T.([algs{a}, '_', metrics{m}]) = col;
    end
  end
  writetable(T, fullfile(pgfDir, sprintf('sweep_N%d.csv', N)));
end
end

function writePgfData(raw, edges, algs, ensembles, timingEns, edgeCases, pgfDir)
% Write the results in a layout pgfplots can read directly
% (\addplot table[x=..., y=...]{file}, col sep=comma):
%
%   scaling_<alg>_d<d>.csv   1a, 1b: one row per N, on timingEns:
%                            N, n, nGates, depth, build_med/min/max,
%                            sim_med/min/max, and build_slope, sim_slope:
%                            log(t_i/t_{i-1}) / log(N_i/N_{i-1}) of the medians
%                            from the previous size (NaN on the first row)
%   error_<alg>_<ens>.csv    2a, 2c: N, d, n, trial, error, unitarity for every
%                            random trial of that algorithm and ensemble
%   edge_table.csv           2b: one row per edge case, one column per algorithm:
%                            'pass', or '<k>/<m> <status>' when k of the m
%                            (d, n) pairs did not pass (worst status shown)
%   edge_max.csv             2b, 2c: same layout, max error and max unitarity
%                            defect over the (d, n) pairs (<alg>_err, <alg>_unit)
%   edge_identity.csv        2b: one row per (edge case, d, n): <alg>_nGates and
%                            <alg>_nIdentity, the gates that are the identity
if ~exist(pgfDir, 'dir'), mkdir(pgfDir); end

% --- 1a, 1b ------------------------------------------------------------
for a = 1:numel(algs)
  alg = algs{a};
  R   = raw(strcmp(raw.algorithm, alg) & strcmp(raw.ensemble, timingEns), :);
  for d = unique(R.d)'
    Rd = R(R.d == d, :);
    N  = unique(Rd.N);
    T  = table(N, round(log(N) / log(d)), 'VariableNames', {'N', 'n'});
    z  = zeros(size(N));
    [T.nGates, T.depth, T.build_med, T.build_min, T.build_max, ...
      T.sim_med, T.sim_min, T.sim_max] = deal(z);
    for i = 1:numel(N)
      Ri = Rd(Rd.N == N(i), :);
      T.nGates(i) = Ri.nGates(1);    T.depth(i) = Ri.depth(1);
      T.build_med(i) = median(Ri.buildTime);
      T.build_min(i) = min(Ri.buildTime); T.build_max(i) = max(Ri.buildTime);
      T.sim_med(i)   = median(Ri.simTime);
      T.sim_min(i)   = min(Ri.simTime);   T.sim_max(i)   = max(Ri.simTime);
    end
    T.build_slope = localSlope(T.N, T.build_med);
    T.sim_slope   = localSlope(T.N, T.sim_med);
    writetable(T, fullfile(pgfDir, sprintf('scaling_%s_d%d.csv', alg, d)));
  end
end

% --- 2a, 2c ------------------------------------------------------------
for a = 1:numel(algs)
  for e = 1:numel(ensembles)
    R = raw(strcmp(raw.algorithm, algs{a}) & strcmp(raw.ensemble, ensembles{e}), ...
      {'N', 'd', 'n', 'trial', 'error', 'unitarity'});
    if isempty(R), continue; end
    writetable(R, fullfile(pgfDir, sprintf('error_%s_%s.csv', algs{a}, ensembles{e})));
  end
end

% --- 2b, 2c ------------------------------------------------------------
worst = {'error', 'nan', 'fail'};              % reported in this order of severity
tbl = table(edgeCases(:), 'VariableNames', {'edgeCase'});
mx  = table(edgeCases(:), 'VariableNames', {'edgeCase'});
for a = 1:numel(algs)
  txt = strings(numel(edgeCases), 1);
  me  = NaN(numel(edgeCases), 1);
  mu  = NaN(numel(edgeCases), 1);
  for c = 1:numel(edgeCases)
    E = edges(strcmp(edges.algorithm, algs{a}) & strcmp(edges.edgeCase, edgeCases{c}), :);
    if isempty(E), txt(c) = "--"; continue; end
    bad = ~strcmp(E.status, 'pass');
    if ~any(bad)
      txt(c) = "pass";
    else
      w = worst{find(cellfun(@(s) any(strcmp(E.status, s)), worst), 1)};
      txt(c) = sprintf('%d/%d %s', nnz(bad), height(E), w);
    end
    me(c) = max(E.error);
    mu(c) = max(E.unitarity);
  end
  tbl.(algs{a}) = txt;
  mx.([algs{a}, '_err'])  = me;
  mx.([algs{a}, '_unit']) = mu;
end
writetable(tbl, fullfile(pgfDir, 'edge_table.csv'));
writetable(mx,  fullfile(pgfDir, 'edge_max.csv'));

keys = unique(edges(:, {'edgeCase', 'd', 'n'}), 'stable');
for a = 1:numel(algs)
  [g, id] = deal(NaN(height(keys), 1));
  for i = 1:height(keys)
    r = strcmp(edges.algorithm, algs{a}) & strcmp(edges.edgeCase, keys.edgeCase{i}) ...
      & edges.d == keys.d(i) & edges.n == keys.n(i);
    if any(r), g(i) = edges.nGates(r); id(i) = edges.nIdentity(r); end
  end
  keys.([algs{a}, '_nGates'])    = g;
  keys.([algs{a}, '_nIdentity']) = id;
end
writetable(keys, fullfile(pgfDir, 'edge_identity.csv'));
end

function s = localSlope(x, y)
% Slope of log y against log x between each row and the previous one.
s = [NaN; diff(log(y)) ./ diff(log(x))];
end
