
## What is F-screening?

F-screening refers to a two-step procedure in linear regression where:

1. An **overall F-test** is used to evaluate the global null hypothesis:
H₀: β₁ = β₂ = ... = βₚ = 0
2. **if (and only if)** this test is rejected, follow-up tests are performed on coefficients of interest.

While intuitive and widely used, this approach introduces **selection bias** — standard inference procedures no longer guarantee valid Type I error control, confidence interval nominal coverage, or unbiased point estimates. The `lmFScreen` package addresses this by providing tools for **conditional selective inference**, which accounts for the conditional nature of the analysis.

---

## Installation

To install the `lmFScreen` package from GitHub, run the following in your R console: `devtools::install_github("mcgougho/lmFScreen")`.

---

## Prospective inference with `lmFScreen` function

The `lmFScreen` function enables valid inference conditional on rejection of an overall F-test in linear regression, using tools from conditional selective inference.

Given a design matrix `X` and a response vector `y`, the package:

1. Performs an overall F-test of the global null hypothesis  
2. If this test is rejected, conducts follow-up tests on specified coefficients (`test_cols`) 
3. Returns:
   - Selective **p-values**
   - Selective **confidence intervals**
   - Selective **point estimates**
   - Naive (unadjusted) counterparts for comparison

This procedure ensures **valid inference** conditional on rejection of the overall null hypothesis. 

See the tutorial: [Using lmFScreen](articles/lmFScreen.html)

---

## Retrospective inference with `psel_retro` function

The `psel_retro()` function allows for valid selective p-values in retrospective settings, where only summary statistics (not raw data) are available. It requires:
1. Sample size n
2. Number of predictors p
3. R-squared
4. Residual standard error (RSE)
5. F-statistic for the coefficient of interest

This is particularly useful when analyzing published results (e.g., from an applied paper’s Table 1) where individual-level data is unavailable.

See the tutorial for how to use it: [Using psel_retro](articles/psel_retro.html)

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
