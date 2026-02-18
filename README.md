# Optimal Universal Black-Box Privacy Preservation and Tight Adversarially-Adaptive Composition for Entropic Secret

This repository provides simulation code of the paper:

> **Optimal Universal Black-Box Privacy Preservation and Tight Adversarially-Adaptive Composition for Entropic Secret**

It contains two simulation tasks:

- RSA timing leakage
- AES S-box Hamming-weight leakage

---

# 1. Installation

## Requirements

- Python == 3.12
- NumPy
- SciPy
- Matplotlib
- PyCryptodome

Install dependencies:

```bash
pip install numpy scipy matplotlib pycryptodome
```

# 2. RSA timing leakage simulation

## Run the simulation

```bash
python rsa_simulation.py
```

The script will:

1. Sample a secret\
2. Run T adaptive rounds\
3. Calibrate σ² each round\
4. Output the transcript

# 3. AES S-box Hamming-weight leakage simulation

## Run the simulation

```
Run **all cells** in quickway.ipynb in order.
```

The notebook will:

- Initialize the keyspace
- Define AES S-box and Hamming-weight leakage
- Sample a secret from the prior
- Perform adaptive multi-round interaction
- Calibrate σ² per round via binary search
- Compute optimal reference distribution \( W_t \)
- Output the transcript and calibrated noise levels

# 4. Key Parameters

Parameter Meaning

---

alpha - Rényi order
T - Number of rounds
m - Sample size
r_t - Per-round α budget
gamma_t - Failure probability
R_t - Worst-case fallback
T_max - Leakage truncation bound
