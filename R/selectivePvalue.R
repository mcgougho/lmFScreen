#' Compute the Selective P-Value Function for Beta1
#'
#' This function returns a function `pselb` that computes the selective p-value
#' for a given value of `b` (corresponding to `beta1`). The function uses
#' a Monte Carlo approach to evaluate the probability that the estimated coefficient
#' satisfies a selection criterion.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of response variables.
#' @param sigma_sq The variance parameter.
#' @param rss The residual sum of squares from a model fit (default: NULL, computed internally if missing).
#' @param alpha_ov The significance level for the overall F-test (default: 0.05).
#' @param B The number of Monte Carlo samples used for estimation.
#' @param min_select The minimum number of selected samples required (default: `B`).
#' @param seed A seed value for reproducibility (default: 123456).
#'
#' @return A function `pselb(b)` that computes the selective p-value for a given `b`.
#' If the overall F-test fails, the function returns `NA`.
#'
#' @details
#' The function performs the following steps:
#' 1. Computes the overall F-statistic for the regression model and checks whether it
#'    passes the selection criterion based on the quantile of an F-distribution.
#' 2. Constructs an orthonormal basis using the singular value decomposition (SVD)
#'    of `X` excluding the first column.
#' 3. Defines a Monte Carlo estimation approach using random samples from normal
#'    and chi-squared distributions to approximate the conditional probability.
#' 4. Returns a function `pselb(b)` that evaluates the p-value for any given `b`.
#'
#'
#' @export
get_pselb <- function(X, y, sigma_sq, yPy = NULL, rss = NULL, alpha_ov = 0.05, B = 10000, min_select = B, B_max = 10e6, seed = NULL, verbose = FALSE) {

  if(is.null(seed)) {
    seed <- rnorm(1, 0, 10000)
    set.seed(seed)
  }

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



#' Retrospective Selective P-Value for a Single Coefficient
#'
#' Computes the selective p-value for a linear regression coefficient from R squared, RSE
#' and the post-hoc F-statistic of interest, conditional on having passed an overall F-test.
#'
#' @param n Sample size (number of observations).
#' @param p Number of predictors of interest (non-intercept) in the model.
#' @param R_squared R-squared value from the overall fitted model.
#' @param RSE Residual standard error from the overall fitted model.
#' @param Fstat Observed F-statistic for the post-hoc.
#' @param sigma_sq Optional variance parameter. If NULL, it will be estimated using debiased estimate.
#' @param alpha_ov Significance level for the overall F-test (default: 0.05).
#' @param B Number of Monte Carlo samples per attempt (default: 10000).
#' @param min_select Minimum number of selected samples required to compute the p-value (default: 100).
#' @param max_attempts Maximum number of sampling attempts if the selection condition is rare (default: 10).
#'
#' @return A scalar value: the estimated selective p-value. If no selected samples are obtained
#' after `max_attempts`, returns NA with a warning.
#'
#'
#' @export
psel_retro <- function(n, p, R_squared, RSE, Fstat, sigma_sq = NULL, alpha_ov=0.05, B=1000000, min_select = 1000, max_attempts = 100){
  # compute needed quantities
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


