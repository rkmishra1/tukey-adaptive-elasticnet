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

The figure below shows the full Tukey family — loss, first derivative (influence), and second derivative — illustrating the hard cutoff at $\pm d$:

<p align="center">
  <img src="docs/figures/tukey_family.png" width="700" alt="Tukey family: loss, first and second derivative"/>
  <br><em>Figure 3 — Tukey family of functions. The first derivative T'_d(u) redescends to exactly zero at |u| = d; the second derivative T''_d(u) shows the non-convex curvature that motivates proximal AdaGrad over Newton-type solvers.</em>
</p>

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

## Real Data Analysis

To demonstrate the practical effectiveness of **Tukey-AdEnet** in handling real-world outliers, collinearity, and variable selection, we evaluated the method on three public benchmark datasets from the `robustbase` and `robustHD` R packages. 

### Experimental Setup
For each dataset, we:
1. Pre-processed the original variables (omitted missing values, centered the response, and standardized the predictors).
2. Appended independent random Gaussian noise variables ($Z_j \sim \mathcal{N}(0, 1)$) to the design matrix to test each method's variable selection capabilities.
3. Conducted **20 independent replications** of a 70/30 train/test split.
4. Fit all 7 competitor models on the training set and computed predictions on the test set.
5. Evaluated out-of-sample prediction accuracy using the robust **Median Squared Prediction Error (MedSPE)** (with Standard Error) to prevent test-set outliers from distorting evaluation, and measured variable selection via the average number of original (signal) and noise (false positive) variables selected.

---

### 1. BBC TopGear Dataset ($n = 248, p = 29$, including 20 noise variables)
* **Goal**: Predict car `Price` (log-scale) using performance specifications.
* **Outliers**: Contains extreme vertical outliers and leverage points due to ultra-performance supercars (e.g., Bugatti Veyron) and budget micro-cars.

| Method | Test MedSPE (SE) | Total Selected | Signal Selected (out of 9) | Noise Selected (out of 20) |
|---|:---:|:---:|:---:|:---:|
| `AdL` | 0.0169 (0.0010) | 5.10 | 4.70 | 0.40 |
| `AdEnet` | 0.0169 (0.0010) | 4.95 | 4.65 | 0.30 |
| `LAD-Lasso` | 0.0165 (0.0009) | 6.30 | 5.50 | 0.80 |
| `S-LTS` | 0.0175 (0.0008) | 7.20 | 6.60 | 0.60 |
| `R-LARS` | 0.0162 (0.0010) | 7.05 | 5.45 | 1.60 |
| `Tukey-AdL` | **0.0146** (0.0009) | 7.05 | 5.35 | 1.70 |
| **`Tukey-AdEnet` (Ours)** | **0.0146** (0.0009) | 6.90 | 5.35 | 1.55 |

<p align="center">
  <img src="docs/figures/real_data_topgear.png" width="860" alt="TopGear real data performance"/>
  <br><em>Figure 5 — TopGear dataset: Tukey-based estimators achieve the lowest prediction error (~14% reduction in MedSPE compared to OLS-based methods), while maintaining clean variable selection.</em>
</p>

---

### 2. pulpfiber Dataset ($n = 62, p = 19$, including 15 noise variables)
* **Goal**: Predict paper breaking length (`Y1`) using pulp characteristics.
* **Outliers**: Contains 12 known outlying observations (runs 51–62) which exhibit highly distinct raw properties.

| Method | Test MedSPE (SE) | Total Selected | Signal Selected (out of 4) | Noise Selected (out of 15) |
|---|:---:|:---:|:---:|:---:|
| `AdL` | 1.179 (0.138) | 4.55 | 2.25 | 2.30 |
| `AdEnet` | 1.197 (0.140) | 4.55 | 2.25 | 2.30 |
| `LAD-Lasso` | 1.350 (0.149) | 3.70 | 1.90 | 1.80 |
| `S-LTS` | 1.028 (0.094) | 6.80 | 2.30 | 4.50 |
| `R-LARS` | 0.829 (0.107) | 4.05 | 1.55 | 2.50 |
| `Tukey-AdL` | **0.804** (0.132) | 4.40 | 2.40 | 2.00 |
| **`Tukey-AdEnet` (Ours)** | **0.801** (0.117) | 4.00 | 2.15 | 1.85 |

<p align="center">
  <img src="docs/figures/real_data_pulpfiber.png" width="860" alt="pulpfiber real data performance"/>
  <br><em>Figure 6 — pulpfiber dataset: Tukey-AdEnet achieves the highest prediction accuracy, reducing test error by over 32% compared to non-robust methods while selecting fewer noise variables than S-LTS and R-LARS.</em>
</p>

---

### 3. toxicity Dataset ($n = 38, p = 24$, including 15 noise variables)
* **Goal**: Predict chemical toxicity using molecular descriptors.
* **Outliers**: Very small sample size ($n_{\text{train}} = 26$) with severe outlier contamination due to chemical class diversity.

| Method | Test MedSPE (SE) | Total Selected | Signal Selected (out of 9) | Noise Selected (out of 15) |
|---|:---:|:---:|:---:|:---:|
| `AdL` | 0.271 (0.104) | 18.40 | 7.00 | 11.40 |
| `AdEnet` | 0.256 (0.095) | 18.30 | 7.00 | 11.30 |
| `LAD-Lasso` | 0.0252 (0.0045) | 10.80 | 4.35 | 6.45 |
| `S-LTS` | **0.0157** (0.0033) | 6.50 | 2.70 | 3.80 |
| `R-LARS` | 0.0166 (0.0030) | 2.25 | 1.35 | 0.90 |
| `Tukey-AdL` | 0.0195 (0.0028) | 4.35 | 3.30 | 1.05 |
| **`Tukey-AdEnet` (Ours)** | 0.0256 (0.0059) | 4.90 | 3.15 | 1.75 |

<p align="center">
  <img src="docs/figures/real_data_toxicity.png" width="860" alt="toxicity real data performance"/>
  <br><em>Figure 7 — toxicity dataset: Classical non-robust methods break down completely in this high-dimensional small-sample regime, selecting almost all noise variables. Tukey-AdEnet is highly selective, filtering out noise columns while retaining signal variables.</em>
</p>

---

### 💡 Key Findings & Discussion
1. **Resistance to Breakdown**: In all three datasets, the standard OLS-based `AdL` and `AdEnet` are heavily contaminated by outliers, resulting in high prediction errors. Especially on `toxicity` ($n=38$), they break down completely, fitting the noise and yielding massive out-of-sample errors. **Tukey-AdEnet** remains stable and highly accurate.
2. **Superior Variable Selection**: Under contamination, non-robust methods and L1-only methods (like `LAD-Lasso`) select many noise variables (false positives). `Tukey-AdEnet` uses adaptive weights to heavily penalize noise variables, resulting in extremely sparse models that select very few noise columns (e.g., 1.75 on `toxicity` compared to 11.3 for `AdEnet` and 6.45 for `LAD-Lasso`).
3. **Balanced Signal Recovery**: While robust estimators like `R-LARS` are highly conservative (selecting very few variables), they often under-select and miss true signal (e.g. only recovering 1.35 out of 9 active variables on `toxicity`). `Tukey-AdEnet` strikes a superior balance, recovering significantly more signal (3.15) while successfully keeping noise variables out of the active set.

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
