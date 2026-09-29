function rows = loadShards(dir, jobs)
% LOADSHARDS  The result rows saved by saveShard for the given jobs,
% concatenated in the order of jobs (the order of a serial run). Fails if a
% shard is missing, so that a merge never mixes in an incomplete run.
rows = {};
missing = {};
for j = 1:numel(jobs)
  file = fullfile(dir, [strrep(jobs{j}, ':', '_'), '.mat']);
  if ~exist(file, 'file'), missing{end+1} = jobs{j}; continue; end %#ok<AGROW>
  s = load(file, 'rows');
  rows = [rows, s.rows]; %#ok<AGROW>
end
if ~isempty(missing)
  error('loadShards:missing', 'No shard for the jobs (not run, or failed): %s', strjoin(missing, ', '));
end
end
