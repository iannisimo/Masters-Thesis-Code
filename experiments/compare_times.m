function compare_times(dirA, dirB)
% COMPARE_TIMES  Compare two runs of scaling_accuracy.m on the same machine,
% e.g. a serial one (run_all.sh) and a parallel one (run_parallel.sh):
%   cp -r results results_serial && ./run_parallel.sh
%   matlab -batch "compare_times('results_serial', 'results')"
% Prints the median ratio B/A of the build and the simulation times, per
% algorithm and d, over the timed runs with N >= 1000 (smaller runs are too
% short to compare); then the same for the fixed-N sweep. Also checks that
% both runs built the same circuits (gates, depth) with the same errors.
if nargin < 2, dirB = 'results'; end
compareFile(dirA, dirB, 'scaling_raw.csv', {'algorithm', 'd', 'n', 'ensemble', 'trial'}, {'algorithm', 'd'});
compareFile(dirA, dirB, 'sweep_raw.csv', {'N', 'd', 'algorithm', 'trial'}, {'algorithm', 'd'});
end

function compareFile(dirA, dirB, file, keys, groups)
A = readtable(fullfile(dirA, 'scaling', file));
B = readtable(fullfile(dirB, 'scaling', file));
if ismember('ensemble', keys)             % only the timed ensemble
  A = A(strcmp(A.ensemble, 'real'), :);
  B = B(strcmp(B.ensemble, 'real'), :);
end
vals = setdiff(A.Properties.VariableNames, keys, 'stable');
A = renamevars(A, vals, strcat(vals, '_A'));
B = renamevars(B, vals, strcat(vals, '_B'));
J = innerjoin(A, B, 'Keys', keys);

fprintf('\n%s: %d runs in both; same gates and depth: %s; largest difference: error %.1e, unitarity %.1e\n', ...
  file, height(J), yesNo(isequal(J.nGates_A, J.nGates_B) && isequal(J.depth_A, J.depth_B)), ...
  max(abs(J.error_A - J.error_B)), max(abs(J.unitarity_A - J.unitarity_B)));

if ismember('N', keys), N = J.N; else, N = J.N_A; end
J = J(N >= 1000, :);
if isempty(J), fprintf('no runs with N >= 1000 to compare the times\n'); return; end
J.build = J.buildTime_B ./ J.buildTime_A;
J.sim   = J.simTime_B   ./ J.simTime_A;
G = groupsummary(J, groups, 'median', {'build', 'sim'});
fprintf('time %s / time %s, median over N >= 1000:\n', dirB, dirA);
for i = 1:height(G)
  fprintf('  %-10s d=%-3d build x%.2f  sim x%.2f  (%d runs)\n', G.algorithm{i}, G.d(i), ...
    G.median_build(i), G.median_sim(i), G.GroupCount(i));
end
fprintf('  all: build x%.2f  sim x%.2f\n', median(J.build), median(J.sim));
end

function s = yesNo(b)
if b, s = 'yes'; else, s = 'NO'; end
end
