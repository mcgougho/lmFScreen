# `lmFScreen`

This repository contains the `lmFScreen` R package and the code to reproduce the figures in the paper *"Inference after F-screening in Linear Regression."*

---

## `lmFScreen` Package

The `lmFScreen` package enables valid post-selection inference following an overall F-test in linear regression, using tools from selective inference.

Given a design matrix `X` and a response vector `y`, the package:

1. Performs an **overall F-test** of the global null hypothesis  
   H₀: β₁ = β₂ = ... = βₚ = 0
2. **If (and only if)** this test is rejected, conducts follow-up tests on specified coefficients (`test_cols`)
3. Returns:
   - Selective **p-values**
   - Selective **confidence intervals**
   - Selective **point estimates**
   - Naive (unadjusted) counterparts for comparison

This procedure ensures **valid inference** conditional on rejection of the overall null hypothesis.

---

## Figures

The `figures/` directory contains R scripts to reproduce all plots from the paper.

### 1. Power and Type I Error (Global and Local Null)
- **Files:** `power_t1error_global_null.R`, `power_t1error_local_null.R`
- **Description:**  
  Demonstrates Type I error control of the selective p-value under the null and compares the power of the selective procedure with a sample-splitting approach.

---

### 2. Confidence Interval Coverage and Width
- **Files:** `CI_globalnull.R`, `CI_localnull.R`
- **Description:**  
  Compares the selective confidence intervals to naive intervals in terms of:
  - Coverage (selective vs. non-selective)
  - Width (one plot with the value of the true coefficient on the horizontal axis and one plot with the number of observations as the horizontal axis)

---

### 3. Multiple Testing Adjustments
- **File:** `multiple_corrections.R`
- **Description:**  
  Shows that traditional multiple testing corrections (Bonferroni, Scheffé) do **not** control Type I error after F-screening.

---

### 4. Selective vs. Naive Point Estimates
- **File:** `point_estimates.R`
- **Description:**  
  Plots the selective vs. non-selective point estimates for the first coefficient, under the global null.

---
