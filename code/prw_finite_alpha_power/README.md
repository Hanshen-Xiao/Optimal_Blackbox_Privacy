# Finite-alpha PRW power experiments

This folder contains a self-contained reimplementation of the Gaussian-noise
power/Hamming-weight experiments using finite-alpha PRW accounting.

The key distinction from `../figures_5_to_8.py` is the composition rule.  Here,
for a total posterior-success target `rho` and prior success `q`, the code uses
the finite-alpha Bernoulli divergence budget

```text
D_alpha(Ber(rho) || Ber(q))
```

and splits its logarithm evenly across `T` rounds.  For each round, the Gaussian
variance `sigma^2` is found by solving the PRW integral

```text
alpha * log int (sum_h p_H(h) phi_{sigma^2}(y-h)^alpha)^{1/alpha} dy
```

against the per-round budget.  The lower-row main calibration uses
`alpha=256` so the near-prior targets remain numerically finite, while the
alpha-sensitivity panels compare `alpha in {16,64,256}`.  This follows the old
notebook profile described
in the anonymous repository README: finite alpha, optimal reference distribution
`W_t`, and calibration in `sigma^2`.

Run:

```bash
/Users/hsxiao/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 run_power_prw_finite_alpha.py
```

Outputs are written to `outputs/`.
