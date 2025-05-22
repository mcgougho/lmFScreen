#' Construct the Selective P-Value Function for a Single Coefficient
#'
#' This function returns a function that computes the selective p-value for a given value of b (the hypothesized value of beta1),
#' conditional on rejection of the overall F-test. The p-value is estimated via Monte Carlo integration, accounting for the selection event.
#'
#' @param X A numeric matrix of predictors with n rows and p columns. The first column corresponds to the coefficient of interest.
#' @param y A numeric response vector of length n.
#' @param sigma_sq The variance of the error term. Must be specified in advance.
#' @param yPy Optional value of the quadratic form y' P_X y. If NULL, it is computed internally.
#' @param rss Optional residual sum of squares. If NULL, it is computed internally using y and yPy.
#' @param alpha_ov Significance level for the overall F-test. Default is 0.05.
#' @param B Number of Monte Carlo samples drawn per iteration. Default is 10,000.
#' @param min_select Minimum number of selected samples required to estimate the p-value. Default is equal to B.
#' @param B_max Maximum number of Monte Carlo samples drawn in any iteration. Default is 10 million.
#' @param verbose Logical flag indicating whether to print progress messages during sampling.
#'
#' @return A function that takes a numeric argument b and returns the selective p-value corresponding to that value of beta1. If the selection condition
#' is not met by the observed data, the returned function will always return NA.
#'
#' @details
#' The returned function estimates the conditional probability that a test statistic T(b) exceeds the observed value, given that the data passes the F-test threshold.
#' Sampling proceeds until at least min_select Monte Carlo draws satisfy the selection condition. If this cannot be achieved within the given sample limits, the function returns NA.
#'
get_pselb <- function(X, y, sigma_sq, yPy = NULL, rss = NULL, alpha_ov = 0.05, B = 10000, min_select = B, B_max = 10e6, verbose = FALSE) {

  n <- dim(X)[1]
  p <- dim(X)[2]
  cc <- (p / (n - p)) * qf(1 - alpha_ov, df1 = p, df2 = (n - p))
  if(is.null(rss)) {
    U <- svd(X)$u
    yPy <- as.numeric(sum((t(U) %*% y)^2))
    rss <- as.numeric(sum(y^2) - yPy)
  }


  # Construct the orthonormal basis using SVD of X excluding the first column.
  X.1 <- X[, -1]
  U1 <- svd(X.1)$u
  u.1 <- X[, 1] - U1 %*% t(U1) %*% X[, 1]
  U.1 <- u.1 / sqrt(sum(u.1^2))

  # Compute projections of y.
  d <- sum((t(U1) %*% y)^2)
  v <- t(U.1) %*% y

  # Compute the F-distribution quantile for selection.
  F.quantile <- qf(1 - alpha_ov, p, (n - p))
  cc <- p / (n - p) * F.quantile

  # Compute terms for later calculations.
  term1 <- t(v) %*% v
  mat1 <- t(X[, 1]) %*% U.1
  term2 <- mat1 %*% v
  term3 <- t(X[, 1]) %*% (U.1 %*% t(U.1)) %*% X[, 1]
  term4 <- t(U.1) %*% X[, 1]

  # Define the function that computes the selective p-value for a given b.
  pselb <- function(b) {

    # Compute the test statistic numerator.
    a <- as.numeric(term1 - 2 * b * term2 + (b^2) * term3)

    # Compute the mean of the projected beta1 effect.
    meann <- term4 * b

    # Generate Monte Carlo samples.
    V_samples_selected <- c()
    Z_samples_selected <- c()
    select <- 0
    num_iter <- 0

    #if(verbose){cat("We are inside psel-B \n")}
    while (length(Z_samples_selected) < min_select){
      num_iter <- num_iter + 1
      if(verbose & num_iter > 1){cat("Sampling iteration number: ", num_iter, "\n")}
      V_samples <- rep(meann,B) + rnorm(B, 0, sqrt(sigma_sq))
      V2_samples <- V_samples^2
      Z_samples <-  rchisq(B, (n - p))

      # Selection criterion for hypothesis testing.
      select <- (V2_samples - cc * sigma_sq * Z_samples >= -d)
      V_samples_selected <- c(V_samples_selected, V_samples[select])
      Z_samples_selected <- c(Z_samples_selected, Z_samples[select])
      # estimate number of additional samples needed
      p_select = mean(select)
      B <- (1.1*min_select - length(Z_samples_selected)) / p_select
      B <- min(B, B_max) # not too large
      if(verbose & num_iter > 1){cat("Number selected so far: ", length(V_samples_selected), "\n")}
    }
    V_samples <- V_samples_selected
    Z_samples <- Z_samples_selected
    n_select <- length(V_samples)


    # Debugging messages for Monte Carlo estimates.
    #cat("Value of b is:", b, "\n")
    #cat("Value of B is:", B, "\n")
    #cat("Mean of V2_samples:", mean(V2_samples), "\n")
    #cat("Mean of select indicator:", mean(select), "\n")

    # Compute the numerator and denominator for p-value estimation.

    V2_samples <- V_samples^2
    numerator <- V2_samples - 2 * b * mat1 %*% V_samples + rep(b^2 * term3, n_select)
    numerator <- as.vector(numerator)
    denominator <- sigma_sq * Z_samples

    # Compute the p-value as the proportion of selected samples that satisfy the test statistic.
    pval <- mean((numerator / denominator) >= (a) / rss)
    return(pval)
  }

  return(pselb)
}



#' Retrospective Selective P-Value Based on Summary Statistics
#'
#' This function is useful for conducting valid retrospective F-screening.
#' Suppose that we have access to the outputs of a least squares linear regression model, such as from the output of summary(lm(y~X)),
#' and we want to conduct a test of the significance of a single regression coefficient (beta_j) that accounts for the rejection of an overall F-test.
#' Then this function can provide a selective p-value for beta_j based on of only a few summary statistics.
#' The arguments of this function include R-squared, residual standard error (RSE), and an F-statistic for the test of H_0: beta_j.
#' This function is especially useful in settings where the raw data is unavailable, such as published studies.
#'
#' @param n Sample size (number of observations).
#' @param p Number of predictors used in F-screening (excluding the intercept).
#' @param R_squared R-squared from the fitted linear model.
#' @param RSE Residual standard error from the fitted model.
#' @param tstat Observed t-statistic for the follow-up hypothesis test of beta_j.
#' @param sigma_sq Optional estimate of the noise variance. If NULL, uses debiased estimate that accounts for selection.
#' @param alpha_ov Significance level for the overall F-test. Default is 0.05.
#' @param B Number of Monte Carlo samples per iteration. Default is 1,000,000.
#' @param min_select Minimum number of samples satisfying the selection condition. Default is 1,000.
#' @param max_attempts Maximum number of iterations for passing selection criterion before giving up. Default is 100.
#'
#' @return A numeric value representing the estimated selective p-value.
#' If no selected samples are obtained after max_attempts, the function returns NA and issues a warning.
#'
#' @examples
#' data(mtcars)
#' mod <- lm(mpg ~ wt + hp, data = mtcars)
#' rse <- summary(mod)$sigma
#' r2 <- summary(mod)$r.squared
#' t_hp <- summary(mod)$coefficients["hp", "t value"]
#' psel_retro(n=nrow(mtcars), p=2, R_squared=r2, RSE=rse, tstat=t_hp)
#' result <- lmFScreen(mpg ~ wt + hp, data = mtcars)
#' result[["selective pvalues"]][2]
#' # the retrospective and prospective p-values coincide (up to Monte Carlo error)
#'
#' @export
psel_retro <- function(n, p, R_squared, RSE, tstat, sigma_sq = NULL, alpha_ov=0.05, B=1000000, min_select = 1000, max_attempts = 100){
  # compute needed quantities
  Fstat <- tstat^2
  a <- Fstat/(n-p-1)
  F.quantile <- qf(1-alpha_ov, p, (n-p-1))
  d <- RSE^2 * Fstat - RSE^2 * (n-p-1) * R_squared / (1-R_squared)
  cc <- F.quantile * p / (n-p-1)
  # get variance estimate
  if (is.null(sigma_sq)){
    Zs <- rchisq(100000, n-p-1)
    As <- rchisq(100000, p)
    scaling <- mean(Zs[As >= cc*Zs])
    RSS <- RSE^2 * (n-p-1)
    sigma_sq <- RSS/scaling
  }
  # compute psel
  total_W <- total_Z <- NULL
  total_select <- 0
  attempt <- 0
  # make sure total_select is large enough
  while (total_select < min_select && attempt < max_attempts) {
    attempt <- attempt + 1
    Z <- rchisq(B, n - p - 1)
    W <- rchisq(B, 1)
    select <- (W - cc * Z >= d / sigma_sq)
    total_W <- c(total_W, W[select])
    total_Z <- c(total_Z, Z[select])
    total_select <- length(total_W)
  }

  if (total_select == 0) {
    warning("No selected samples even after max_attempts. Returning NA.")
    return(NA)
  }

  pval <- mean((total_W / total_Z) >= a)
  return(pval)
}


