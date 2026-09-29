function saveShard(dir, job, rows)
% SAVESHARD  Save the result rows of one parallel job (see run_parallel.sh)
% to <dir>/<job>.mat, with ':' in the job name replaced by '_'.
if ~exist(dir, 'dir'), mkdir(dir); end
file = [strrep(job, ':', '_'), '.mat'];
save(fullfile(dir, file), 'rows');
fprintf('\nWrote shard %s (%d rows)\n', file, numel(rows));
end
