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

## License

MIT (see `LICENSE`), except `qclab/`, which keeps its own license.
