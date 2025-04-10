% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectivePvalue.R
\name{psel_retro}
\alias{psel_retro}
\title{Retrospective Selective P-Value Based on Summary Statistics}
\usage{
psel_retro(
  n,
  p,
  R_squared,
  RSE,
  Fstat,
  sigma_sq = NULL,
  alpha_ov = 0.05,
  B = 1e+06,
  min_select = 1000,
  max_attempts = 100
)
}
\arguments{
\item{n}{Sample size (number of observations).}

\item{p}{Number of predictors used in F-screening (excluding the intercept).}

\item{R_squared}{R-squared from the fitted linear model.}

\item{RSE}{Residual standard error from the fitted model.}

\item{Fstat}{Observed F-statistic for the follow-up hypothesis test of beta_i.}

\item{sigma_sq}{Optional estimate of the noise variance. If NULL, uses debiased estimate that accounts for selection.}

\item{alpha_ov}{Significance level for the overall F-test. Default is 0.05.}

\item{B}{Number of Monte Carlo samples per iteration. Default is 1,000,000.}

\item{min_select}{Minimum number of samples satisfying the selection condition. Default is 1,000.}

\item{max_attempts}{Maximum number of iterations for passing selection criterion before giving up. Default is 100.}
}
\value{
A numeric value representing the estimated selective p-value.
If no selected samples are obtained after max_attempts, the function returns NA and issues a warning.
}
\description{
Computes a selective p-value for a single regression coefficient (beta_i) using summary statistics from a linear model:
R-squared, residual standard error (RSE), and an F-statistic for the test of H_0: beta_i.
This function is useful in retrospective settings where the raw data is unavailable.
}
\details{
The selective p-value is estimated as the conditional probability that a test statistic W / Z exceeds the observed value,
where W follows a chi-squared distribution with 1 degree of freedom, and Z follows a chi-squared distribution with n - p - 1 degrees of freedom.
The selection condition is based on the rejection threshold for the overall F-test.

This method can be applied in settings where only summary statistics are available, such as published studies.
}
\examples{
data(mtcars)
mod <- lm(mpg ~ wt + hp, data = mtcars)
rse <- summary(mod)$sigma
r2 <- summary(mod)$r.squared
t_wt <- summary(mod)$coefficients["wt", "t value"]
psel_retro(n=nrow(mtcars), p=2, R_squared=r2, RSE=rse, Fstat=t_wt^2)

}
