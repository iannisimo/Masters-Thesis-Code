# Experiment data: what is recorded and what it is for

Two scripts write everything under `results/`; `results/` holds nothing else. Rerun both to regenerate it:

| Script | Writes | Chapter 5 sections |
|---|---|---|
| `scaling_accuracy.m` | `results/scaling/` | correctness, degenerate inputs, build time, simulation time, equal state size |
| `d2_consistency.m` | `results/d2/` | consistency at d = 2, cost of the qudit generalization |

`results/run_all.log` is the console output of the last run: one line per run, ending in `ok` or `FAIL`.

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
| **real, equal N** | `randn` | 1 + trial | all three | The equal-N sweep. At a given N every d and every algorithm prepares **the same state**. |
| **structured** | 10 real inputs that reach the degenerate cases (zero entries, zero nodes, basis states, uniform, sparse) | 7919·case + k | all three | Pass/fail, error, defect, identity gates. The chapter keeps one sentence on them (all pass) and uses the identity counts in the future-work paragraph. |
| **real, d = 2** | `randn` | 2000 + trial | KP-Tree, both QdKP-Trees | The same states as the real ensemble of the scaling grid at d = 2. |

Complex states cost more time. At N ≥ 1000, Bullock builds about 1.4× slower on them (median; 1.05–2.3×) and simulates about 1.3× slower. The Householder tree builds at the same speed and simulates about 1.1× slower. This is why every time comes from the real ensemble.

## The grids

| Grid | Points | Runs |
|---|---|---|
| Scaling in N | d ∈ {2, 3, 4, 5}, n = 2, 3, … while N = dⁿ ≤ 2¹⁴: 13, 7, 6, 5 sizes | 31 per (algorithm, ensemble), 155 in total |
| Degenerate inputs | 10 inputs × (d, n) ∈ {(2,4), (3,3), (4,2), (5,2)} × 3 algorithms | 120 |
| Equal N | N = 729 (d = 3, 9, 27), 4096 (d = 2, 4, 8, 16, 64), 6561 (d = 3, 9, 81), 15625 (d = 5, 25, 125) × 3 algorithms | 42 |
| d = 2 consistency | n = 2, …, 14; gate-by-gate check up to n = 10 | 13 per algorithm |

Every configuration runs once (`N_TRIALS = 1`). The `*_min`/`*_max` columns and the summary files become meaningful when `N_TRIALS` > 1.

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
| `status` | `pass`, `fail`, `nan` or `error` | Every input is handled (all 120 pass); stated in one sentence in the correctness section. |
| `nIdentity` | gates whose matrix is the identity, phase gate included | How many gates do nothing. On sparse inputs both reflector-based loaders have many (14 and 13 of 16 at (2,4) with d nonzero entries), while the Givens tree skips empty nodes; this backs the sparse-state future-work paragraph. On the uniform state only Bullock has them (its reflectors collapse nodes ahead of their turn). |

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
| `edge_max.csv` | Largest error and defect per input and algorithm | the bounds quoted for the degenerate inputs |
| `edge_table.csv` | Pass/fail per input and algorithm | the check behind "all pass" |
| `edge_identity.csv` | nGates and nIdentity per (input, d, n) and algorithm | the identity-reflector counts |
| `sweep/sweep_N<N>.csv` | One row per d: n, the predicted counts `ref_householder` and `ref_givens`, and for each algorithm nGates, depth, buildTime, simTime, error, unitarity | Figures 5.4 and 5.5 (section "Equal state size") |

### `results/d2/`

| File | Content | Used in |
|---|---|---|
| `d2_raw.csv`, `d2_summary.csv` | Per n: gates, error of each loader; `*_col0` = ‖(C − C_KP)\|0…0⟩‖₂ and `*_U` = ‖(C − C_KP)X‖₂ against the KP-Tree; the gate-by-gate check (`gatewise_*`, n ≤ 10); build and simulation times of all three | Table 5.3; the qubit-vs-qudit cost comparison in the simulation section |

## Not recorded

- **Complex states for the Givens tree:** it only accepts real amplitudes.
- **Repeated trials:** `N_TRIALS = 1`, so times are single runs; small differences can flip between runs.
- **Hardware cost:** QRAM query cost, and the decomposition of multi-controlled gates into two-qudit gates, are not measured.
