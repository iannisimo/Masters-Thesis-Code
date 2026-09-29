% D2_CONSISTENCY  The qubit KP-Tree (KPTree.m) against the qudit tree loaders
% (QdKPTree_Givens.m, QdKPTree_Householder.m) at d = 2, on the same real
% states. At d = 2 the Givens QdKP-Tree should build the circuit of the
% KP-Tree, and the Householder QdKP-Tree the same state with 2 x 2
% reflectors in place of the rotations.
%
% For every n in NS and every trial (the same seeds as the real ensemble of
% scaling_accuracy.m):
%   - the preparation error of each algorithm;
%   - Givens vs KP-Tree: ||(C_G - C_KP) X||_2 with X random N x UNIT_K with
%     orthonormal columns (0 if the two unitaries coincide), the same on the
%     first column only, and, for n <= GATEWISE_MAX_n, a gate-by-gate check:
%     same qudits and largest ||M_G - M_KP||_2 over the gates;
%   - Householder vs KP-Tree: the same two distances (the first column is
%     the prepared state, so it should agree; the unitaries need not);
%   - gate counts, build and simulation times.
%
% Output (written to OUT_DIR):
%   d2_raw.csv       one row per (n, trial)
%   d2_summary.csv   one row per n: max of the distances and errors, median
%                    of the times over the trials (pgfplots-ready)
%
% Run from anywhere:  run('code/experiments/d2_consistency.m')

clear; clc;

%% ------------------------------------------------------------------------
%  Configuration (edit here)
%  ------------------------------------------------------------------------

NS             = 2:14;                   % number of qubits
N_TRIALS       = 3;                      % random real states per n
SEED0          = 2000;                   % seeds SEED0 + trial (real ensemble of scaling_accuracy.m)
UNIT_K         = 4;                      % columns of the random probe X
GATEWISE_MAX_n = 10;                     % gate-by-gate check up to this n

HERE    = fileparts(mfilename('fullpath'));
OUT_REL = fullfile('results', 'd2');        % printed in the log (no local paths)
OUT_DIR = fullfile(HERE, OUT_REL);
addpath(fullfile(HERE, '..'));           % the functions add qclab and util themselves
addpath(HERE);                           % saveShard, loadShards
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

% Jobs for run_parallel.sh, chosen by the environment variable JOB (unset:
% everything, here): list (print them, delete old shards), n:<n> (one n ->
% shards/n_<n>.mat), merge (the output files from the shards).
JOB       = strtrim(getenv('JOB'));
SHARD_DIR = fullfile(OUT_DIR, 'shards');
jobs      = arrayfun(@(n) sprintf('n:%d', n), NS, 'UniformOutput', false);
if strcmp(JOB, 'list')
  if exist(SHARD_DIR, 'dir'), rmdir(SHARD_DIR, 's'); end
  fprintf('@@JOB %s\n', jobs{:});
  return
end
if ~isempty(JOB) && ~any(strcmp(JOB, [jobs, {'merge'}]))
  error('Unknown JOB "%s" (JOB=list prints the valid ones)', JOB);
end

%% ------------------------------------------------------------------------
%  Runs
%  ------------------------------------------------------------------------

% warm-up (untimed): the first call of each algorithm pays for class loading
if ~strcmp(JOB, 'merge')
  [~] = KPTree(2, 2, 0, 0);
  [~] = QdKPTree_Givens(2, 2, 0, 0);
  [~] = QdKPTree_Householder(2, 2, 0, 0);
end

rows = {};
for n = NS
  if ~(isempty(JOB) || strcmp(JOB, sprintf('n:%d', n))), continue; end
  N = 2^n;
  for trial = 1:N_TRIALS
    seed = SEED0 + trial;                % same seed, same state for all three
    fprintf('n=%-2d N=%-6d trial %d (build/sim): ', n, N, trial);

    [cK, eK, bK, sK] = KPTree(2, n, 0, seed);
    fprintf('KP %.3f/%.3fs err %.1e | ', bK, sK, eK);
    [cG, eG, bG, sG] = QdKPTree_Givens(2, n, 0, seed);
    fprintf('Givens %.3f/%.3fs err %.1e | ', bG, sG, eG);
    [cH, eH, bH, sH] = QdKPTree_Householder(2, n, 0, seed);
    fprintf('House %.3f/%.3fs err %.1e | ', bH, sH, eH);

    rng(seed + 7919);
    [X, ~] = qr(randn(N, UNIT_K) + 1i * randn(N, UNIT_K), 0);
    e0 = zeros(N, 1); e0(1) = 1;
    YK = cK.apply('R', 'N', n, [e0, X], 0, 2);
    YG = cG.apply('R', 'N', n, [e0, X], 0, 2);
    YH = cH.apply('R', 'N', n, [e0, X], 0, 2);
    dG0 = norm(YG(:, 1) - YK(:, 1));  dG = norm(YG(:, 2:end) - YK(:, 2:end));
    dH0 = norm(YH(:, 1) - YK(:, 1));  dH = norm(YH(:, 2:end) - YK(:, 2:end));

    sameQubits = NaN; gateDiff = NaN;
    if n <= GATEWISE_MAX_n
      [sameQubits, gateDiff] = compareGates(cG, cK);
    end

    rows{end+1} = {n, N, trial, seed, cK.nbObjects, cG.nbObjects, cH.nbObjects, ...
      eK, eG, eH, dG0, dG, dH0, dH, sameQubits, gateDiff, ...
      bK, sK, bG, sG, bH, sH}; %#ok<SAGROW>

    fprintf('gates %d/%d/%d | Givens-KP: col0 %.1e, U %.1e, gatewise %s | House-KP: col0 %.1e, U %.1e\n', ...
      cK.nbObjects, cG.nbObjects, cH.nbObjects, dG0, dG, gatewiseText(sameQubits, gateDiff), dH0, dH);
  end
end

if startsWith(JOB, 'n:'), saveShard(SHARD_DIR, JOB, rows); return; end
if strcmp(JOB, 'merge'), rows = loadShards(SHARD_DIR, jobs); end

raw = cell2table(vertcat(rows{:}), 'VariableNames', {'n', 'N', 'trial', 'seed', ...
  'gates_kp', 'gates_givens', 'gates_house', 'err_kp', 'err_givens', 'err_house', ...
  'givens_col0', 'givens_U', 'house_col0', 'house_U', 'gatewise_same_qubits', ...
  'gatewise_max_diff', 'build_kp', 'sim_kp', 'build_givens', 'sim_givens', ...
  'build_house', 'sim_house'});
writetable(raw, fullfile(OUT_DIR, 'd2_raw.csv'));

%% ------------------------------------------------------------------------
%  Summary per n (max of distances and errors, median of times)
%  ------------------------------------------------------------------------

maxCols = {'err_kp', 'err_givens', 'err_house', 'givens_col0', 'givens_U', ...
  'house_col0', 'house_U', 'gatewise_max_diff'};
medCols = {'build_kp', 'sim_kp', 'build_givens', 'sim_givens', 'build_house', 'sim_house'};
S = table(NS(:), 2.^NS(:), 'VariableNames', {'n', 'N'});
for i = 1:numel(NS)
  R = raw(raw.n == NS(i), :);
  S.gates_kp(i) = R.gates_kp(1); S.gates_givens(i) = R.gates_givens(1); S.gates_house(i) = R.gates_house(1);
  S.gatewise_same_qubits(i) = min(R.gatewise_same_qubits);
  for c = maxCols, S.(c{1})(i) = max(R.(c{1})); end
  for c = medCols, S.(c{1})(i) = median(R.(c{1})); end
end
writetable(S, fullfile(OUT_DIR, 'd2_summary.csv'));

fprintf('\nWrote %s\n', fullfile(OUT_REL, 'd2_raw.csv'));
fprintf('Wrote %s\n', fullfile(OUT_REL, 'd2_summary.csv'));
if strcmp(JOB, 'merge'), rmdir(SHARD_DIR, 's'); end


%% ========================================================================
%  Helpers
%  ========================================================================

function [same, maxDiff] = compareGates(cA, cB)
% Gate by gate: do the two circuits act on the same qudits, and how far
% apart are the matrices of corresponding gates (on the qudits they span)?
same = cA.nbObjects == cB.nbObjects;
maxDiff = NaN;
if ~same, return; end
maxDiff = 0;
for k = 1:cA.nbObjects
  a = cA.objects(k);
  b = cB.objects(k);
  if ~isequal(a.qubits, b.qubits)
    same = false;
    maxDiff = NaN;
    return
  end
  maxDiff = max(maxDiff, norm(full(a.matrix(2)) - full(b.matrix(2))));
end
end

function t = gatewiseText(same, diff)
if isnan(same), t = '(skipped)';
elseif same, t = sprintf('same qubits, max diff %.1e', diff);
else, t = 'DIFFERENT gates';
end
end
