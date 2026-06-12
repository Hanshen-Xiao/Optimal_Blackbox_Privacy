# Optimal Black-Box Privatization for Entropic Secrets

This repository contains the code used to regenerate the experimental figures
for the paper on black-box privatization for entropic secrets.

The experiments study two scalar side-channel leakage models:

1. AES-style power leakage, modeled as aggregate Hamming weight.
2. RSA timing leakage, modeled by a black-box square-and-multiply timing
   simulator and discretized into 1 ns timing bins.

The scripts reproduce the main 2x4 figures and the history-dependent versus
mechanism-agnostic composition comparison.  Generated data and figures are
written under `code/outputs/`.

## Repository Layout

```text
.
├── README.md
├── requirements.txt
├── code/
│   ├── gaussian_comparison_figures.py
│   ├── figures_5_to_8.py
│   ├── combined_all_figures.py
│   ├── rsa_timing_empirical_figures.py
│   ├── history_vs_agnostic_composition.py
│   ├── discrete_noise_convex.py
│   ├── discrete_leakage_experiments.py
│   ├── empirical_composition.py
│   ├── rsa_timing_one_round.py
│   └── ...
├── scripts/
│   ├── run_main_figures.sh
│   └── copy_paper_figures.sh
└── paper_figures/
```

## Environment

The code was developed with Python 3.12.  A recent Python 3.10+ environment
should also work.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

The experiments use:

- `numpy`
- `scipy`
- `cvxpy`
- `pillow`

Some discrete optimization routines use the SciPy HiGHS linear-programming
solver through `scipy.optimize.linprog`; finite-alpha PRW convex programs use
CVXPY with the bundled CLARABEL solver.

## Quick Reproduction

From the repository root, run:

```bash
bash scripts/run_main_figures.sh
```

This generates:

- `code/outputs/combined/power_experiments_2x4.png`
- `code/outputs/rsa_timing_empirical/rsa_timing_experiments_2x4.png`
- `code/outputs/history_vs_agnostic/history_vs_agnostic_composition.png`

The RSA timing script uses 100,000 empirical samples per case and can take
substantially longer than the closed-form power experiments.  The power
experiments are usually much faster because their leakage distribution is
available analytically.

To copy the regenerated paper figures into `paper_figures/`, run:

```bash
bash scripts/copy_paper_figures.sh
```

This creates:

- `paper_figures/power_main.png`
- `paper_figures/timing_main.png`
- `paper_figures/history_vs_agnostic_composition.png`

These are the figure filenames used in the LaTeX paper.

## Figure-by-Figure Commands

Run all commands from `code/`.

```bash
cd code
```

### Power Panels (a)--(d)

```bash
python gaussian_comparison_figures.py
```

Outputs:

- `outputs/gaussian_comparison/gaussian_comparison_results.csv`
- `outputs/gaussian_comparison/gaussian_comparison_zoom_40_80_results.csv`
- `outputs/gaussian_comparison/fig1_4_combined_1x4_sigma025_30_zoom4080_v7.png`

These panels compare Fano mutual information, PAC-`alpha`, PRW with
`alpha in {16,64,256}`, and the MAP ground truth under Gaussian noise.

### Power Panels (e)--(h)

```bash
python figures_5_to_8.py
```

Outputs:

- `outputs/figures_5_to_8/figures_5_to_8_results.csv`
- `outputs/figures_5_to_8/fig5_8_combined_1x4_v5.png`

These panels study entropy, constrained randomization, and history-dependent
composition.  The optimized discrete mechanisms use integer perturbation
support radius 2048 and minimize expected absolute perturbation.

### Combined Power 2x4 Figure

```bash
python combined_all_figures.py
```

Outputs:

- `outputs/combined/power_experiments_2x4.png`

This script reads the CSV files produced by `gaussian_comparison_figures.py`
and `figures_5_to_8.py`.

### RSA Timing 2x4 Figure

```bash
python rsa_timing_empirical_figures.py
```

Outputs:

- `outputs/rsa_timing_empirical/rsa_timing_empirical_results.csv`
- `outputs/rsa_timing_empirical/rsa_timing_experiments_2x4.png`

The timing leakage is treated as a black-box simulator.  The script samples
100,000 exponent strings for each case, evaluates the timing proxy, rounds the
output to 1 ns bins, and computes all accounting quantities from the empirical
finite distribution.  No confidence correction is applied in these plots.

### History-Dependent vs. Mechanism-Agnostic Composition

```bash
python history_vs_agnostic_composition.py
```

Outputs:

- `outputs/history_vs_agnostic/history_vs_agnostic_composition.csv`
- `outputs/history_vs_agnostic/history_vs_agnostic_composition.png`

This comparison uses Gaussian perturbations, PRW order `alpha=256`, and
composition lengths up to 100.  The mechanism-agnostic curve uses the balanced
Holder allocation from the agnostic composition theorem.

## Experiment Parameters

The main parameters are hard-coded in the scripts for reproducibility:

- PRW orders in the main upper-row comparisons: `16`, `64`, `256`.
- PAC-`alpha` order grid: see `alpha_grid()` in
  `gaussian_comparison_figures.py` and `rsa_timing_empirical_figures.py`.
- Gaussian-noise x-axis for power panels (a)--(d): `sigma = 0.25, 0.50, ..., 30`.
- Zoom range for power panels (a)--(d): `sigma = 40, 41, ..., 80`.
- RSA timing samples per case: `100_000`.
- RSA timing seed: `20260606`.
- Discrete perturbation radius: `2048`.
- History-dependent Gaussian composition order: `alpha=256`.
- Composition targets: full identification at `2^-240` and `2^-220`.
- History-vs-agnostic comparison horizon: up to `100` rounds.

## LP Used for Optimized Discrete Perturbations

For a finite leakage domain and finite perturbation support
`U = {u_1, ..., u_K}`, the optimized perturbation is represented by a
probability vector `q`, where `q_j = Pr[N = u_j]`.

For full identification, group secrets by leakage atom `b` and let
`m_b = max_x Pr[x]` among secrets with leakage `b`.  The MAP posterior success
after additive perturbation is

```text
sum_y max_b m_b q_{y-b}.
```

Introducing epigraph variables `t_y` gives the LP:

```text
minimize      sum_j |u_j| q_j
subject to    q_j >= 0
              sum_j q_j = 1
              sum_j u_j q_j = 0        # only for zero-mean noise
              t_y >= m_b q_{y-b}       # for all feasible y,b
              sum_y t_y <= target
```

For non-negative noise, the support is `[0, 2048]` and the zero-mean equality is
omitted.  For a second-moment utility objective, replace the cost `|u_j|` by
`u_j^2`.

The arbitrary-task version replaces `m_b q_{y-b}` by the task-weighted linear
score for each adversarial report.  See Appendix B of the paper for the formal
notation.

## Notes on Runtime and Determinism

The power experiments are deterministic up to numerical integration and solver
tolerances.  RSA timing uses fixed random seeds, so repeated runs on the same
software stack should produce the same empirical samples and nearly identical
figures.  Small numerical differences can occur across CVXPY, CLARABEL, SciPy,
or BLAS versions.

## Legacy and Diagnostic Scripts

The repository also includes scripts used during development:

- `discrete_leakage_experiments.py`
- `convex_replot_experiments.py`
- `power_short_experiments.py`
- `power_focus_experiments.py`
- `power_nonuniform_experiments.py`
- `rsa_timing_one_round.py`
- `rsa_timing_composition_sweep.py`
- `rsa_confidence_sample_requirements.py`
- `confidence_composition.py`
- `generate_learning_workflow.py`

These are not required for the final main figures, but they regenerate earlier
diagnostics, confidence/sample-complexity checks, and intermediate plots used
while developing the experiments.

