# Tukey Adaptive Elastic Net

> Robust sparse regression via Tukey's biweight loss + adaptive elastic net penalty, fitted by proximal AdaGrad and tuned by robust BIC.

![R](https://img.shields.io/badge/R-%3E%3D4.0-276DC3?logo=r&logoColor=white)
![Status](https://img.shields.io/badge/status-research-orange)
![Methods](https://img.shields.io/badge/competitors-7%20methods-blueviolet)
![Configs](https://img.shields.io/badge/configurations-81-informational)

---

## Overview

This repository provides the full simulation study accompanying the manuscript on the **Tukey-AdEnet** estimator — a robust, penalized regression method designed for sparse linear models under outlier contamination.

Classical penalized estimators (Lasso, Elastic Net, adaptive variants) break down when responses or design points are corrupted. **Tukey-AdEnet** replaces the squared loss with Tukey's redescending biweight, applies adaptive elastic net penalties, and selects tuning parameters via a **robust BIC (RBIC)** criterion — remaining consistent under contamination while retaining the variable-selection properties of the elastic net.

**Key properties:**
- Redescending influence function — outliers are down-weighted, not just shrunk
- Adaptive weights from an initial robust fit produce oracle-consistent selection
- Proximal AdaGrad optimizer handles the non-convex, non-smooth objective
- RBIC over a 2-D `(λ₁, λ₂)` grid avoids cross-validation under contamination
- Scales to high-dimensional regimes (`p > n`)

---

## Table of Contents

- [Installation](#installation)
- [Quick Check](#quick-check)
- [Simulation Design](#simulation-design)
- [Usage](#usage)
  - [Full Manuscript Run](#full-manuscript-run)
  - [Pilot Run](#pilot-run)
  - [Method Subset](#method-subset)
  - [Boxplots](#boxplots)
- [Project Structure](#project-structure)
- [Output Columns](#output-columns)
- [Method Details](#method-details)
- [Competing Methods](#competing-methods)

---

## Installation

Install all R package dependencies in one step:

```sh
Rscript scripts/install_packages.R
```

**Required packages:** `glmnet`, `rqPen`, `robustHD`, `robustbase`

> If a package is missing, `run_simulation.R` skips that method by default. For manuscript-grade runs pass `--missing_action=stop` to fail fast on any missing dependency.

---

## Quick Check

Verify the estimator and RBIC tuning work before committing to a full run:

```sh
Rscript scripts/smoke_test.R
```

---

## Simulation Design

The study follows the linear model **y = Xβ + ε** across a full factorial grid of 81 configurations:

| Factor | Levels |
|---|---|
| AR(1) correlation `ρ` | 0.30, 0.60, 0.80 |
| Dimensional regime | `ζ₁₂` (p<n), `ζ₂₃` (p≈n), `ζ₅₆` (p>n) |
| Contamination scenario | Clean, response only, response + design |
| Active set | `s = 3 × ⌊p/9⌋` nonzero coefficients |
| Replications | 200 per configuration |

---

## Usage

### Full Manuscript Run

```sh
Rscript scripts/run_simulation.R \
  --reps=200 \
  --missing_action=stop \
  --output_dir=results
```

> **Note:** The full grid (200 reps × 81 configs × 7 methods × 2-D RBIC) is computationally intensive. Plan accordingly or run on a cluster.

### Pilot Run

A smaller run for testing and exploration:

```sh
Rscript scripts/run_simulation.R \
  --reps=5 \
  --scenarios=response_design \
  --regimes=zeta_2_3 \
  --n_values=500 \
  --rhos=0.30 \
  --n_lambda1=8 \
  --lambda2_factors=0,0.1,0.5 \
  --max_iter=500 \
  --output_dir=results
```

### Method Subset

Run only a specific subset of methods:

```sh
Rscript scripts/run_simulation.R \
  --methods=AdL,AdEnet,Tukey-AdL,Tukey-AdEnet \
  --reps=10
```

### Boxplots

Generate MSPE boxplots from a completed results CSV:

```sh
Rscript scripts/make_boxplots.R \
  --raw=results/comparison_raw_YYYYMMDD_HHMMSS.csv \
  --output_dir=figures
```

---

## Project Structure

```
.
├── R/
│   ├── tukey_adenet.R      # Tukey loss, gradient, proximal AdaGrad, RBIC
│   ├── simulate_data.R     # Data-generating mechanisms and simulation grid
│   ├── metrics.R           # MSPE and variable-selection metrics (C, IC)
│   └── competitors.R       # Wrappers for all seven comparison methods
└── scripts/
    ├── run_simulation.R    # Command-line simulation runner
    ├── install_packages.R  # One-step dependency installer
    ├── smoke_test.R        # Quick sanity check
    └── make_boxplots.R     # MSPE boxplot generator
```

---

## Output Columns

The runner writes two CSV files to `--output_dir`:

| File | Description |
|---|---|
| `comparison_raw_*.csv` | One row per method per replication |
| `comparison_summary_*.csv` | Manuscript-style averages by method and configuration |

**Column reference:**

| Column | Description |
|---|---|
| `C` | True zero coefficients correctly estimated as zero |
| `IC` | True nonzero coefficients incorrectly estimated as zero (false negatives) |
| `MSPE` | `(β̂ − β)ᵀ Σ (β̂ − β)` with AR(1) `Σ` |
| `lambda1` | RBIC-selected L1 tuning parameter |
| `lambda2` | RBIC-selected L2 tuning parameter |
| `criterion` | BIC or RBIC value at the selected model |
| `converged` | Whether proximal AdaGrad met the stopping tolerance |

---

## Method Details

**Loss function** — Tukey biweight score:

```
ψ(u) = u · (1 − (u/d)²)² · 𝟙(|u| ≤ d)
```

**Coordinate-wise proximal AdaGrad update:**

```
u_j    = β_j − η_j · ∇_j
β_j    = sign(u_j) · max(|u_j| − η_j λ₁ w_j, 0) / (1 + η_j λ₂)
```

where `w_j` are adaptive weights derived from an initial robust fit, and `η_j` is the per-coordinate AdaGrad learning rate.

**Initialization:** Uses `lmrob` (from `robustbase`) when `p < n`; falls back to ridge initialization for high-dimensional settings so all regimes remain runnable without strict package requirements.

---

## Competing Methods

| Method | Description |
|---|---|
| `AdL` | Adaptive Lasso |
| `AdEnet` | Adaptive Elastic Net |
| `LAD-Lasso` | Least Absolute Deviations Lasso |
| `Tukey-AdL` | Tukey loss + Adaptive Lasso |
| `S-LTS` | S-estimator / Least Trimmed Squares |
| `R-LARS` | Robust LARS |
| `Tukey-AdEnet` | **This work** |

---

*Correspondence: open an issue for questions about the simulation code.*
