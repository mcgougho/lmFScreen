

** If you are viewing this on GitHub, see to the following website for a more readable version **
[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mcgougho.github.io/lmFScreen/)




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

