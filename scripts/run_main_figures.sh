#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT/code"

echo "[1/5] Power Gaussian-comparison panels (a)-(d)"
python gaussian_comparison_figures.py

echo "[2/5] Power entropy/mechanism/composition panels (e)-(h)"
python figures_5_to_8.py

echo "[3/5] Combined power 2x4 figure"
python combined_all_figures.py

echo "[4/5] RSA timing 2x4 figure"
python rsa_timing_empirical_figures.py

echo "[5/5] History-dependent vs. mechanism-agnostic composition"
python history_vs_agnostic_composition.py

echo "Done. Main generated figures:"
echo "  $ROOT/code/outputs/combined/power_experiments_2x4.png"
echo "  $ROOT/code/outputs/rsa_timing_empirical/rsa_timing_experiments_2x4.png"
echo "  $ROOT/code/outputs/history_vs_agnostic/history_vs_agnostic_composition.png"

