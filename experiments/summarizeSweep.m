function summary = summarizeSweep(sweep, algs, outDir)
% SUMMARIZESWEEP  Part 3 of scaling_accuracy.m (equal state size): the median
% over the trials of every (N, d, algorithm), written as
%   <outDir>/sweep_summary.csv   one row per (N, d, n, algorithm): GroupCount
%                                and median_<metric> for every metric
%   <outDir>/pgf/sweep/          the same values laid out for pgfplots (see
%                                writeSweepPgf)
% Called by scaling_accuracy.m on the runs it just made. With no arguments,
% it reads results/scaling/sweep_raw.csv and rewrites both from it, without
% running any experiment.
%
% Run from anywhere:  run('experiments/summarizeSweep.m')  or  summarizeSweep()

if nargin == 0
  here   = fileparts(mfilename('fullpath'));
  outDir = fullfile(here, 'results', 'scaling');
  sweep  = readtable(fullfile(outDir, 'sweep_raw.csv'), 'TextType', 'char');
  algs   = {'bullock', 'qdkptree', 'qdkpgivens'};
end

metrics = {'nGates', 'depth', 'buildTime', 'simTime', 'error', 'unitarity'};
summary = groupsummary(sweep, {'N', 'd', 'n', 'algorithm'}, 'median', metrics);
writetable(summary, fullfile(outDir, 'sweep_summary.csv'));
writeSweepPgf(summary, algs, metrics, fullfile(outDir, 'pgf', 'sweep'));
end

function writeSweepPgf(summary, algs, metrics, pgfDir)
% Part 3, for pgfplots:
%   sweep_N<N>.csv   one row per d: d, n, ref_householder = (N-1)/(d-1) + 1,
%                    ref_givens = N - 1, and <alg>_<metric> (median over the
%                    trials; NaN for algorithms not run at that d)
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
        if any(r), col(i) = S.(['median_', metrics{m}])(r); end
      end
      T.([algs{a}, '_', metrics{m}]) = col;
    end
  end
  writetable(T, fullfile(pgfDir, sprintf('sweep_N%d.csv', N)));
end
end
