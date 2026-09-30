# Experiment data: what is recorded and what it is for

Two scripts write the data under `results/`, and a third derives the tables of Chapter 5 from it; `results/` holds nothing else. Rerun them in this order to regenerate it:

| Script | Writes | Chapter 5 sections |
|---|---|---|
| `scaling_accuracy.m` | `results/scaling/` | correctness (random inputs; the structured ones are run but not reported), build time, simulation time, equal state size |
| `d2_consistency.m` | `results/d2/` | consistency at d = 2, cost of the qudit generalization |
| `thesis_tables.m` | `results/tables/` | the three tables, and the checks behind the text; reads the CSVs above, runs no experiment |

`results/run_all.log` is the console output of the last serial run: one line per run, ending in `ok` or `FAIL`. A parallel run (`run_parallel.sh`) writes `results/run_parallel.log` instead (one line per job, with its core, exit code and duration) and the console output of every job to `results/logs/`; its output files are the same. The two runs differ only in the times and, at the level of 1e-16, in `unitarity` (the probe uses multithreaded linear algebra in a serial run, one thread in a parallel one).

The algorithms, as named in the files:

| File name | Algorithm | Accepts |
|---|---|---|
| `bullock` | Bullock et al. (`QR_bullock.m`) | real and complex states |
| `qdkptree` | Householder QdKP-Tree (`QdKPTree_Householder.m`) | real and complex states |
| `qdkpgivens` | Givens QdKP-Tree (`QdKPTree_Givens.m`) | real states only |
| `kp` | qubit KP-Tree (`KPTree.m`), only in `d2_consistency.m` | real states, d = 2 only |

## Input ensembles

Every state is normalized, and every seed is fixed, so a rerun prepares the same states.

| Ensemble | State | Seeds | Run by | Recorded and used for |
|---|---|---|---|---|
| **complex** | `randn + 1i*randn` | 1000 + trial | Bullock, Householder | **Correctness only.** Error and unitarity defect (Figure 5.1, Table 5.1). Its times are in `scaling_raw.csv` but are not exported or used. |
| **real** | `randn` | 2000 + trial | all three | **Correctness and every time.** All algorithms are timed on this ensemble, so the time plots compare them on the same states. |
| **real, equal N** | `randn` | trial | all three | The equal-N sweep. At a given N every d and every algorithm prepares **the same state**. |
| **structured** | 10 real inputs that reach the degenerate cases (zero entries, zero nodes, basis states, uniform, sparse) | 7919·case + k | all three | Pass/fail, error, defect, identity gates. **Run but not reported in the thesis** (see below). |
| **real, d = 2** | `randn` | 2000 + trial | KP-Tree, both QdKP-Trees | The same states as the real ensemble of the scaling grid at d = 2. |

Complex states cost more time. At N ≥ 1000, Bullock builds about 1.4× slower on them (median; 1.05–2.3×) and simulates about 1.3× slower. The Householder tree builds at the same speed and simulates about 1.1× slower. This is why every time comes from the real ensemble.

## Structured inputs: run, not reported

The structured (degenerate) inputs are still run by `scaling_accuracy.m`, and all 120 runs pass, but Chapter 5 does not mention them: they add nothing to the narrative, which rests on the random states. Their data (`edge_*.csv`, the degenerate-inputs part of `scaling_accuracy.m`) is kept as a check on the implementations. The future-work remark on identity gates in sparse states is stated in the chapter without their counts.

## The grids

| Grid | Points | Runs |
|---|---|---|
| Scaling in N | d ∈ {2, 3, 4, 5}, n = 2, 3, … while N = dⁿ ≤ 2¹⁴: 13, 7, 6, 5 sizes | 3 × 31 per (algorithm, ensemble), 465 in total |
| Degenerate inputs | 10 inputs × (d, n) ∈ {(2,4), (3,3), (4,2), (5,2)} × 3 algorithms | 120 (fixed states, run once) |
| Equal N | N = 729 (d = 3, 9, 27), 4096 (d = 2, 4, 8, 16, 64), 6561 (d = 3, 9, 81), 15625 (d = 5, 25, 125) × 3 algorithms | 3 × 42 = 126 |
| d = 2 consistency | n = 2, …, 14; gate-by-gate check up to n = 10 | 3 × 13 per algorithm |

Every random configuration runs 3 times (`N_TRIALS = 3`, trials 1, 2, 3 with the seeds above; trial 1 is the state of the earlier single-trial runs). The times in `pgf/scaling_*` are the median over the trials, with `*_min`/`*_max` for the spread, those in `d2_summary.csv` the median; `sweep_summary.csv` and `pgf/sweep/` hold the mean.

## The values

### Recorded for every run

| Value | Meaning | Use |
|---|---|---|
| `error` | ‖C\|0…0⟩ − ψ‖₂ | Does the circuit prepare the state? |
| `unitarity` | ‖Y†Y − I₄‖₂ with Y = CX, X random N×4 with orthonormal columns | Is C unitary? A small error alone does not prove it: a tree loader with wrong node values still maps \|0…0⟩ to ψ. |
| `buildTime` [s] | from the input state to the finished circuit | Classical preprocessing: linear for the trees, towards quadratic for Bullock. |
| `simTime` [s] | simulating the circuit on \|0…0⟩ | The cost of the simulator, O(N) per gate once sparse. |
| `nGates` | gates in the circuit, phase gate included | Checks the counts of Chapter 4: Bullock and Householder (N−1)/(d−1) + 1, Givens N − 1 (at d = 2, Householder 2ⁿ and Givens 2ⁿ − 1). |
| `depth` | greedy layering: each gate one layer after the latest gate on any of its qudits (controls and target) | **Recorded, not yet in the chapter.** It would show how far each circuit parallelizes as built, without QRAM. |

### Degenerate inputs only

| Value | Meaning | Use |
|---|---|---|
| `status` | `pass`, `fail`, `nan` or `error` | Every input is handled (all 120 pass). |
| `nIdentity` | gates whose matrix is the identity, phase gate included | How many gates do nothing. On sparse inputs both reflector-based loaders have many (14 and 13 of 16 at (2,4) with d nonzero entries), while the Givens tree skips empty nodes; this backs the sparse-state future-work paragraph. On the uniform state only Bullock has them (its reflectors collapse nodes ahead of their turn). Not quoted in the thesis. |

## Files

### `results/scaling/`

| File | Content |
|---|---|
| `scaling_raw.csv` | One row per random run: algorithm, d, n, N, ensemble, trial, seed, nGates, depth, buildTime, simTime, error, unitarity. Holds everything, including the unused complex-state times. |
| `edge_raw.csv` | One row per (input, d, n, algorithm): status, nGates, depth, nIdentity, error, unitarity, and the error message if the construction failed. |
| `sweep_raw.csv`, `sweep_summary.csv` | The equal-N runs; the summary is the mean over trials. |

### `results/scaling/pgf/` (read by the thesis; rewritten from scratch on every run)

| File | Content | Used in |
|---|---|---|
| `scaling_<alg>_d<d>.csv` | Real ensemble, one row per N: n, nGates, depth, build/sim median, min and max, and `build_slope`/`sim_slope`, the slope of log t against log N from the previous size | Figures 5.2 and 5.3, Table 5.4, every slope quoted in the text |
| `error_<alg>_<ens>.csv` | Every random run: N, d, n, trial, error, unitarity | Figure 5.1, Table 5.1 |
| `edge_max.csv` | Largest error and defect per input and algorithm | not used |
| `edge_table.csv` | Pass/fail per input and algorithm | not used |
| `edge_identity.csv` | nGates and nIdentity per (input, d, n) and algorithm | not used |
| `sweep/sweep_N<N>.csv` | One row per d: n, the predicted counts `ref_householder` and `ref_givens`, and for each algorithm nGates, depth, buildTime, simTime, error, unitarity | Figures 5.4 and 5.5 (section "Equal state size") |

### `results/d2/`

| File | Content | Used in |
|---|---|---|
| `d2_raw.csv`, `d2_summary.csv` | Per n: gates, error of each loader; `*_col0` = ‖(C − C_KP)\|0…0⟩‖₂ and `*_U` = ‖(C − C_KP)X‖₂ against the KP-Tree; the gate-by-gate check (`gatewise_*`, n ≤ 10); build and simulation times of all three | Table 5.3; the qubit-vs-qudit cost comparison in the simulation section |

### `results/tables/` (`thesis_tables.m`)

The figures read their CSVs directly; the tables are typed in the thesis, so they are cross-referenced against this output.

| File | Content |
|---|---|
| `tables.txt` | Every cell of Tables 5.1, 5.3 and 5.4 at full precision (`%.17g`), the value as printed in the thesis, and the run it comes from; then the checks behind the text: every random run passes, the largest error and defect and how far below 1e-8 they are, gate counts against Chapter 4, runs with an error of exactly 0 (not drawn in Figure 5.1), the d = 2 distances and gate-by-gate check, the spread, time per gate and last slope at the largest N of each d, and the build slopes on complex states. |
| `tab_random_error.tex`, `tab_d2.tex`, `tab_times.tex` | The body rows of the three tables (between `\midrule` and `\bottomrule`), exactly as printed in `tables.txt`. |

## Not recorded

- **Complex states for the Givens tree:** it only accepts real amplitudes.
- **Many trials:** 3 per configuration give a median and a range, not a confidence interval; close times can still swap order.
- **Hardware cost:** QRAM query cost, and the decomposition of multi-controlled gates into two-qudit gates, are not measured.
