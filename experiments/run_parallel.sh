#!/usr/bin/env bash
# run_parallel.sh: run the experiments as independent jobs, one MATLAB per job,
# each pinned to its own physical core with a single computational thread
# (the hyperthread siblings are left idle), then merge them into results/.
# The output files are those of run_all.sh.
#
# usage: ./run_parallel.sh [script.m ...]   (default: all of them)
#   MATLAB=/path/to/matlab   MATLAB binary (default: matlab on PATH)
#   CORES="0 1 2 ..."        logical CPUs to use (default: the first thread of
#                            every physical core, alternating between sockets)
#   TIMEOUT=7200             seconds allowed per job
#
# Log: results/run_parallel.log (one line per job), results/logs/<job>.log.
# To keep it running after logging out:  nohup ./run_parallel.sh > /dev/null 2>&1 &
# The jobs of a script: JOB=list matlab -batch "run('script.m')" (see the script).
cd "$(dirname "$(realpath "$0")")" || exit 1
MATLAB=${MATLAB:-matlab}
TIMEOUT=${TIMEOUT:-7200}
SCRIPTS=("$@")
[ ${#SCRIPTS[@]} -eq 0 ] && SCRIPTS=(scaling_accuracy.m d2_consistency.m)
command -v "$MATLAB" > /dev/null || { echo "matlab not found; set MATLAB=/path/to/matlab" >&2; exit 1; }
command -v taskset > /dev/null || { echo "taskset not found (package util-linux)" >&2; exit 1; }
trap 'kill 0' INT TERM

# One logical CPU per physical core (the first of its thread siblings),
# ordered socket 0, 1, 2, ..., 0, 1, 2, ... so that the jobs spread evenly.
physicalCores() {
  local c
  for c in /sys/devices/system/cpu/cpu[0-9]*; do
    [ -r "$c/topology/thread_siblings_list" ] || continue
    [ "$(sed 's/[,-].*//' "$c/topology/thread_siblings_list")" = "${c##*cpu}" ] || continue
    echo "$(cat "$c/topology/physical_package_id") ${c##*cpu}"
  done | sort -n -k1,1 -k2,2 | awk '{print r[$1]++, $1, $2}' | sort -n -k1,1 -k2,2 | awk '{print $3}'
}
read -r -a CPUS <<< "${CORES:-$(physicalCores | tr '\n' ' ')}"
[ ${#CPUS[@]} -gt 0 ] || { echo "no cores found; set CORES" >&2; exit 1; }

mkdir -p results/logs
rm -f results/logs/*.log
SUMMARY=results/run_parallel.log
: > "$SUMMARY"
log() { echo "$*" | tee -a "$SUMMARY"; }

# run one job of one script, pinned to one CPU; the output goes to its log
runJob() {  # cpu script job
  local cpu=$1 s=$2 job=$3 t0=$SECONDS rc
  local name="${s%.m}_${job//:/_}"
  JOB=$job timeout -k 10 "$TIMEOUT" taskset -c "$cpu" "$MATLAB" -singleCompThread \
    -batch "run('$PWD/$s')" < /dev/null > "results/logs/$name.log" 2>&1
  rc=$?
  log "$(date +%T)  cpu $cpu  $name  exit=$rc  $((SECONDS - t0))s"
  return $rc
}

# 1. the jobs of every script (JOB=list also deletes the old shards)
JOBS=()
for s in "${SCRIPTS[@]}"; do
  out=$(JOB=list taskset -c "${CPUS[0]}" "$MATLAB" -singleCompThread \
    -batch "fprintf('@@VER %s\n', version); run('$PWD/$s')" < /dev/null 2>&1)
  list=$(sed -n 's/^@@JOB //p' <<< "$out")
  [ -n "$list" ] || { echo "$out"; echo "no jobs from $s" >&2; exit 1; }
  for j in $list; do JOBS+=("$s $j"); done
  VER=$(sed -n 's/^@@VER //p' <<< "$out")
done
nJobs=${#JOBS[@]}
nCpu=${#CPUS[@]}
[ "$nCpu" -gt "$nJobs" ] && nCpu=$nJobs
log "host $(hostname), MATLAB $VER, $(date '+%F %T')"
log "$nJobs jobs on $nCpu cores (one thread each): ${CPUS[*]:0:nCpu}"

# 2. core i runs jobs i, i + nCpu, i + 2 nCpu, ... one after the other
for ((i = 0; i < nCpu; i++)); do
  (
    fail=0
    for ((j = i; j < nJobs; j += nCpu)); do
      runJob "${CPUS[i]}" ${JOBS[j]} || fail=1
    done
    exit $fail
  ) &
done
fail=0
for p in $(jobs -p); do wait "$p" || fail=1; done
if [ $fail -ne 0 ]; then
  log "===== some jobs failed, not merging (see results/logs/)"
  exit 1
fi

# 3. merge the shards of every script into the output files
for s in "${SCRIPTS[@]}"; do runJob "${CPUS[0]}" "$s" merge || fail=1; done
log "===== all done ($(date +%T))"
exit $fail
