#' Estimate Variance for Selective Inference
#'
#' This function computes an estimate of the variance parameter `sigma_sq`
#' for selective inference using a rescaled residual sum of squares (RSS).
#'
#' @param X A matrix of predictor variables.
#' @param Y A vector of response variables.
#' @param alpha_ov The significance level for the overall F-test.
#'
#' @return A numeric value representing the estimated variance `sigma_sq`.
#'
#' @details
#' The function follows these steps:
#' 1. Computes a scaling factor using Monte Carlo simulation from chi-squared distributions.
#' 2. Computes the residual sum of squares (RSS) from a linear model fit.
#' 3. Uses the scaling factor to adjust the RSS and obtain the variance estimate.
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10)
#' Y <- rnorm(10)
#' sigma_sq_est <- get_variance_estimate(X, Y, alpha_ov = 0.05)
#' print(sigma_sq_est)
#'
#' @export
get_variance_estimate <- function(X,Y,alpha_ov){
  n <- dim(X)[1]
  p <- dim(X)[2]
  cc <- (p/(n-p))*qf(1-alpha_ov, df1=p, df2=(n-p))
  Zs <- rchisq(100000, n-p)
  As <- rchisq(100000, p)
  scaling <- mean(Zs[As >= cc*Zs])
  lm_mod <- lm(Y~X)
  rss <- sum(residuals(lm_mod)^2)
  sigma_sq <- as.numeric(rss/scaling)
}
