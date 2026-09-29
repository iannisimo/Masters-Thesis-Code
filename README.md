# Qudit state preparation with QCLAB

MATLAB code for the experiments of my master's thesis: circuits that prepare a
state ψ ∈ ℂ^(dⁿ) on n qudits of dimension d, built and simulated with a
qudit extension of [QCLAB](https://github.com/iannisimo/qclab/tree/qudit).

| File | Algorithm |
|---|---|
| `QR_bullock.m` | Bullock et al., QR decomposition with d × d reflectors (real and complex states) |
| `QdKPTree_Householder.m` | Householder QdKP-Tree (real and complex states) |
| `QdKPTree_Givens.m` | Givens QdKP-Tree (real states only) |
| `KPTree.m` | qubit KP-Tree, the d = 2 baseline (real states only) |
| `util/` | helpers (tree, reflectors, circuit assembly) |
| `experiments/` | the benchmark scripts and their results |

## Requirements

- MATLAB (tested with R2025b), no toolboxes.
- The QCLAB fork, branch `qudit`, included as the submodule `qclab/`.

## Setup

```bash
git clone --recurse-submodules https://github.com/iannisimo/Masters-Thesis-Code.git
# already cloned without it:  git submodule update --init
```

## Usage

Every loader has the same signature, and adds its own paths:

```matlab
addpath('/path/to/Masters-Thesis-Code')
% d, n, IMAG (0 = real state), seed, psi_in (optional: a given state)
[circuit, err, buildTime, simTime] = QdKPTree_Householder(3, 4, 1, 42);
```

`err` is ‖C|0…0⟩ − ψ‖₂; `circuit` is a `qclab.QCircuit`.

## Experiments

```bash
cd experiments
./run_all.sh                          # both scripts, a few minutes in total
nohup ./run_all.sh > /dev/null 2>&1 & # same, keeps running after logout
```

`run_all.sh` runs `scaling_accuracy.m` and `d2_consistency.m` with
`matlab -batch` and logs to `results/run_all.log`. Set `MATLAB=/path/to/matlab`
if `matlab` is not on the PATH. The scripts can also be run on their own from
MATLAB, e.g. `run('experiments/scaling_accuracy.m')`. Their parameters (the
dimensions, largest state size, number of trials) are at the top of each script.

The times are wall-clock times, so run the experiments on an otherwise idle
machine. `experiments/DATA.md` describes every output file and column.

### In parallel

```bash
./run_parallel.sh                     # same output as run_all.sh
```

`run_parallel.sh` splits the scripts into independent jobs (one per algorithm
and d, per N and d of the equal-size sweep, per n of the d = 2 check) and runs
one MATLAB per job with `-singleCompThread`, pinned with `taskset` to its own
physical core; hyperthread siblings are left idle. The jobs are spread over the
sockets, and a final merge writes the same files as a serial run. `CORES="0 1 2"`
restricts it to the given logical CPUs. Logs: `results/run_parallel.log` (one
line per job) and `results/logs/`.

Jobs running side by side share caches and memory bandwidth. To check that
this does not change the times, compare a serial and a parallel run on the
same machine:

```bash
./run_all.sh && cp -r results results_serial && ./run_parallel.sh
matlab -batch "compare_times('results_serial', 'results')"
```

## License

MIT (see `LICENSE`), except `qclab/`, which keeps its own license.
