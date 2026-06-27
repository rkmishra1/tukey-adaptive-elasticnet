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
- Redescending influence function — outliers are down-weighted *to zero* beyond breakdown point `d`
- Adaptive weights from an initial robust fit produce oracle-consistent selection
- Proximal AdaGrad optimizer handles the non-convex, non-smooth objective
- RBIC over a 2-D `(λ₁, λ₂)` grid avoids cross-validation under contamination
- Scales to high-dimensional regimes (`p > n`)

---

## Fitting Pipeline

<p align="center">
  <img src="docs/figures/pipeline.png" width="840" alt="Tukey-AdEnet fitting pipeline"/>
  <br><em>Figure 1 — Four-stage pipeline: robust initialisation (lmrob / ridge) → adaptive weights → RBIC grid search → proximal AdaGrad iterations</em>
</p>

---

## Why Tukey's Biweight?

The biweight loss **hard-zeros the influence of any residual beyond the tuning constant `d`** — unlike OLS (unbounded) or Huber (bounded but non-redescending). This gives the estimator a positive breakdown point even under leverage contamination.

<p align="center">
  <img src="docs/figures/loss_influence.png" width="840" alt="Loss and influence function comparison"/>
  <br><em>Figure 2 — Left: Tukey biweight loss stays bounded and flat for large residuals. Right: influence function redescends exactly to zero at |r| = d (dashed verticals), providing hard resistance to extreme outliers.</em>
</p>

The Tukey score function:

```
ψ(u) = u · (1 − (u/d)²)² · 𝟙(|u| ≤ d)
```

---

## Penalisation & Tuning

Adaptive weights `ŵ_j = 1/|β̃_j|` concentrate the L1 penalty on noise variables, shrinking them to exact zeros while leaving signal variables lightly penalised. The 2-D `(λ₁, λ₂)` pair is chosen by minimising RBIC over a warm-started grid.

<p align="center">
  <img src="docs/figures/penalty_rbic.png" width="840" alt="Regularisation path and RBIC surface"/>
  <br><em>Figure 3 — Left: coefficient paths (noise variables zero out early under heavy adaptive weights). Right: RBIC surface with selected (λ₁*, λ₂*) marked.</em>
</p>

**Coordinate-wise proximal AdaGrad update:**

```
u_j  = β_j − η_j · ∇_j
β_j  = sign(u_j) · max(|u_j| − η_j λ₁ ŵ_j, 0) / (1 + η_j λ₂)
```

---

## Simulation Results

### Performance across 7 methods (ζ₂₃ regime, ρ = 0.60, response + design contamination, δ = 10%)

<p align="center">
  <img src="docs/figures/simulation_metrics.png" width="860" alt="Simulation performance metrics"/>
  <br><em>Figure 4 — Tukey-AdEnet (blue) leads on all three metrics: highest correct zeros (C ↑), fewest false negatives (IC ↓), lowest median MSPE (↓).</em>
</p>

### Simulation Design

The study follows **y = Xβ + ε** across a full factorial grid of 81 configurations:

| Factor | Levels |
|---|---|
| AR(1) correlation `ρ` | 0.30, 0.60, 0.80 |
| Dimensional regime | `ζ₁₂` (p<n), `ζ₂₃` (p≈n), `ζ₅₆` (p>n) |
| Contamination scenario | Clean, response only, response + design |
| Active set | `s = 3 × ⌊p/9⌋` nonzero coefficients |
| Replications | 200 per configuration |

### Output Metrics

| Column | Description |
|---|---|
| `C` | True zero coefficients correctly estimated as zero ↑ |
| `IC` | True nonzero coefficients incorrectly zeroed (false negatives) ↓ |
| `MSPE` | `(β̂ − β)ᵀ Σ (β̂ − β)` with AR(1) `Σ` ↓ |

---

## Installation

```sh
Rscript scripts/install_packages.R
```

**Required packages:** `glmnet`, `rqPen`, `robustHD`, `robustbase`

> Missing packages cause that method to be skipped by default. Pass `--missing_action=stop` for manuscript-grade runs.

---

## Quick Check

```sh
Rscript scripts/smoke_test.R
```

---

## Usage

### Full Manuscript Run

```sh
Rscript scripts/run_simulation.R \
  --reps=200 \
  --missing_action=stop \
  --output_dir=results
```

> The full grid (200 reps × 81 configs × 7 methods × 2-D RBIC) is computationally intensive — plan for cluster use.

### Pilot Run

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

```sh
Rscript scripts/run_simulation.R \
  --methods=AdL,AdEnet,Tukey-AdL,Tukey-AdEnet \
  --reps=10
```

### Boxplots

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
│   ├── simulate_data.R     # DGP and simulation grid
│   ├── metrics.R           # MSPE, C, IC
│   └── competitors.R       # Wrappers for all seven comparison methods
├── scripts/
│   ├── run_simulation.R    # CLI simulation runner
│   ├── install_packages.R  # One-step dependency installer
│   ├── smoke_test.R        # Sanity check
│   └── make_boxplots.R     # MSPE boxplot generator
└── docs/figures/           # Figures embedded in this README
```

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
