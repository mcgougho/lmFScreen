#' Estimate Variance for Selective Inference
#'
#' Computes an estimate of the error variance to be used in selective inference,
#' based on a rescaled residual sum of squares from a linear regression model.
#'
#' @param X A numeric matrix of predictor variables with n rows and p columns.
#' @param Y A numeric response vector of length n.
#' @param alpha_ov The significance level used in the overall F-test.
#'
#' @return A numeric value representing the estimated variance used for selective inference.
#'
#' @details
#' This function estimates the variance of the error term under the assumption that
#' inference is performed conditional on rejection of the overall F-test.
#'
#'
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
