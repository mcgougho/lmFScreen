

** If you are viewing this on GitHub, see to the following website for a more readable version **
[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mcgougho.github.io/lmFScreen/)




## What is F-screening?

Suppose that we have and $n$-vector $Y$ containing a quantitative response for $n$ obserations and an $n\times p$ design matrix $X$ containing $p$ covariates/features for each of $n$ observations. Out interest lies in the linear model $Y = X\beta +\epsilon$. 

A common data analysis pipeline, which we refer to as "F-screening" is as follows:

1. Test the "overall" null hypothesis, 
$H_0: \beta_1 = \beta_2 = ... = \beta_p = 0$
using an F-test.
2. **if (and only if)** this test is rejected, conduct inference on $\beta_j$ for some coefficient $j$.

Typically, Steps 1 and 2 are carried out using, e.g., the command `lm(Y~X)`. (However, as will be explained in the next paragraph, this analysis pipeline is problematic.)

While F-screening is intuitive and widely used, carrying out Step 2 using a “standard” approach -- for instance, testing $H_{0}^j: \beta_j=0$ using the t-test output by `lm(Y~X)` -- is problematic, as it does not account for the fact that we conduct inference on $\beta_j$ only if we rejected the “overall” null hypothesis in Step 1. Consequently, the standard t-test for $H_{0}^j: \beta_j=0$ will not control the Type 1 error, standard confidence intervals for $\beta_j$ will not attain the nominal coverage, and even the point estimate for $\beta_j$ output by `lm(Y~X)` will be biased. 

The `lmFScreen` package provides a valid inferential toolbox for conducting Step 2 in the F-screening procedure, by accounting for the fact that we conduct inference on $beta_j$ in Step 1 only if we reject the overall null hypothesis in Step 1. 

This inferential toolbox falls under the larger umbrella of “conditional selective inference”, where the “selection event” that we are conditioning on is rejection of the overall null hypothesis in Step 1.


---

## Installation

To install the `lmFScreen` package from GitHub, run the following in your R console: `devtools::install_github("mcgougho/lmFScreen")`.

---

## Prospective inference with `lmFScreen` function

The `lmFScreen` function enables valid inference conditional on rejection of an overall F-test in linear regression, using tools from conditional selective inference.

Given a design matrix `X` and a response vector `y`, the package:

1. Tests overall F-test of the overall null hypothesis $H_0: \beta_1 = \beta_2 = ... = \beta_p = 0$ using `lm(y~x)`
2. If this test is rejected, conduct inference on specified coefficients (`test_cols`) 
3. Returns:
   - Selective **p-values** for $H_0^j:\beta_j=0$
   - Selective **confidence intervals** for $\beta_j$
   - Selective **point estimates** for $\beta_j$
   - Standard (unadjusted) counterparts for comparison arising from `lm(y~x)`

Unlike the standard output arising from `lm(y~x), the selective p-values, confidence intervals, and points estimates are valid *conditional on rejection of the overall null hypothesis.*

A tutorial for how to use this function in the vignettes section of the FScreen website [![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mcgougho.github.io/lmFScreen/): [Using lmFScreen](articles/lmFScreen.html)

---

## Retrospective inference with `psel_retro` function

The `psel_retro()` function allows for valid selective p-values in retrospective settings, where only the output of `summary(lm(y~x))` (not the raw data `x` and `y`) are available. It requires:
1. Sample size n
2. Number of predictors p
3. R-squared (i.e `summary(lm(y~x))$r.squared`)
4. Residual standard error (RSE) (i.e. `summary(lm(y~x))$sigma`)
5. F-statistic for the coefficient of interest (i.e. `summary(lm(y~x))$coefficients[j, "t value"]`, where `j-1` is the index of the coefficient of interest)

This is particularly useful when analyzing published results (e.g., from an applied paper’s Table 1) where individual-level data is unavailable.

A tutorial for how to use thsi function is in vignettes section of the FScreen website [![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mcgougho.github.io/lmFScreen/): [Using psel_retro](articles/psel_retro.html)

---

