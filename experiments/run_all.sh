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

mkdir -p results
LOG=results/run_all.log
: > "$LOG"
{
  echo "host $(hostname), $(nproc) cores, $(date '+%F %T')"
  "$MATLAB" -batch "fprintf('MATLAB %s\n', version)" < /dev/null 2>&1 | tail -1
} | tee -a "$LOG"

fail=0
for s in "${SCRIPTS[@]}"; do
  echo "===== $s ($(date +%T))" | tee -a "$LOG"
  timeout -k 10 "$TIMEOUT" "$MATLAB" -batch "run('$PWD/$s')" < /dev/null 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  echo "exit=$rc (${s%.m})" | tee -a "$LOG"
  [ "$rc" -ne 0 ] && fail=1
done
echo "===== all done ($(date +%T))" | tee -a "$LOG"
exit $fail
