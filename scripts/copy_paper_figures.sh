#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ROOT/paper_figures"

cp "$ROOT/code/outputs/combined/power_experiments_2x4.png" \
  "$ROOT/paper_figures/power_main.png"
cp "$ROOT/code/outputs/rsa_timing_empirical/rsa_timing_experiments_2x4.png" \
  "$ROOT/paper_figures/timing_main.png"
cp "$ROOT/code/outputs/history_vs_agnostic/history_vs_agnostic_composition.png" \
  "$ROOT/paper_figures/history_vs_agnostic_composition.png"

echo "Copied paper figures to $ROOT/paper_figures"

