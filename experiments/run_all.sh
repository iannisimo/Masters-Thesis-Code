#!/usr/bin/env bash
# run_all.sh: run every experiment headless and log to results/run_all.log.
#
# usage: ./run_all.sh [script.m ...]      (default: all of them, in order)
#   MATLAB=/path/to/matlab   MATLAB binary (default: matlab on PATH)
#   TIMEOUT=7200             seconds allowed per script (default 7200)
#
# A full run takes a few minutes. To keep it running after logging out:
#   nohup ./run_all.sh > /dev/null 2>&1 &     then   tail -f results/run_all.log
# Run nothing else heavy on the machine meanwhile: the build and simulation
# times are wall-clock times.
cd "$(dirname "$(realpath "$0")")" || exit 1
MATLAB=${MATLAB:-matlab}
TIMEOUT=${TIMEOUT:-7200}
SCRIPTS=("$@")
[ ${#SCRIPTS[@]} -eq 0 ] && SCRIPTS=(scaling_accuracy.m d2_consistency.m)
command -v "$MATLAB" > /dev/null || { echo "matlab not found; set MATLAB=/path/to/matlab" >&2; exit 1; }

# MATLAB statement that runs the MATLAB code $1 and then kills MATLAB: on some
# machines MATLAB randomly hangs or segfaults on exit, so it is never left to
# exit by itself. Prints @@OK if the code ran without error, else the error
# report and @@ERR; the exit code of MATLAB is meaningless (killed).
mwrap() {
  echo "try, $1; fprintf('\n@@OK\n'); catch mErr, fprintf(2, '%s\n', getReport(mErr, 'extended', 'hyperlinks', 'off')); fprintf('\n@@ERR\n'); end; system(sprintf('kill -9 %d', feature('getpid')));"
}
# exit status of a MATLAB run from its output file $1 and the exit code $2 of
# timeout: 0 on @@OK, 124 on timeout, else 1
mstatus() {
  if grep -q '^@@OK' "$1"; then echo 0; elif [ "$2" -eq 124 ]; then echo 124; else echo 1; fi
}

trap 'trap - INT TERM; pkill -P $$; pkill -f "run.'\''$PWD/"; exit 130' INT TERM   # also stop the MATLAB job
mkdir -p results
LOG=results/run_all.log
: > "$LOG"
{
  echo "host $(hostname), $(nproc) cores, $(date '+%F %T')"
  "$MATLAB" -batch "$(mwrap "fprintf('MATLAB %s\\n', version)")" < /dev/null 2>&1 | grep '^MATLAB '
} | tee -a "$LOG"

fail=0
for s in "${SCRIPTS[@]}"; do
  echo "===== $s ($(date +%T))" | tee -a "$LOG"
  # in the background + wait, so that a kill reaches the trap at once
  { timeout -k 10 "$TIMEOUT" "$MATLAB" -batch "$(mwrap "run('$PWD/$s')")" < /dev/null 2>&1; echo $? > results/.rc; } \
    | tee results/.out | grep --line-buffered -v -E '^@@|matlab: line [0-9]+: +[0-9]+ Killed' | tee -a "$LOG" &
  wait $!
  rc=$(mstatus results/.out "$(cat results/.rc)"); rm -f results/.rc results/.out
  echo "exit=$rc (${s%.m})" | tee -a "$LOG"
  [ "$rc" -ne 0 ] && fail=1
done
echo "===== all done ($(date +%T))" | tee -a "$LOG"
exit $fail
