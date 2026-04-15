#' Construct the Selective P-Value Function for a Single Coefficient
#'
#' This function returns a function that computes the selective p-value for a given value of b (the hypothesized value of beta1),
#' conditional on rejection of the overall F-test. The p-value is estimated via Monte Carlo integration, accounting for the selection event.
#'
#' @param X A numeric matrix of predictors with n rows and p columns. The first column corresponds to the coefficient of interest.
#' @param y A numeric response vector of length n.
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
get_pselb <- function(X, y, alpha_ov = 0.05, B = 1000, min_select = B, max_draws = 10e6, verbose = FALSE) {
  # some checks
  if (length(y) != nrow(X)) {
    stop("Length of y must match the number of rows in X.")
  }
  if (any(is.na(y)) || any(is.na(X))) {
    stop("Input data contains NA values. Please handle missing data before calling this function.")
  }
  if (alpha_ov <= 0 || alpha_ov >= 1) {
    stop("Significance level alpha_ov must be between 0 and 1.")
  }

  n <- dim(X)[1]
  p <- dim(X)[2]

  # Construct the orthonormal basis using SVD of X excluding the first column.
  X.1 <- X[, -1]
  qrX.1 <- qr(X.1)
  r <- qrX.1$rank
  k <- n-r
  Qr <- qr.Q(qrX.1, complete = FALSE)       # Qr %*% t(Qr) = P_{X_{-1}}
  V <- qr.Q(qrX.1, complete = TRUE)[, (r+1):n, drop = FALSE]  # V %*% t(V) = I - P_{X_{-1}}
  XV <- as.numeric(crossprod(V, X[,1]))
  lambda <- sum(XV^2)

  # Compute the F-distribution quantile for selection.
  F.quantile <- qf(1 - alpha_ov, p, (n - p))
  cc <- p / (n - p) * F.quantile

  # Compute terms for later calculations.
  XtX <- crossprod(X)
  XtX_inv <- solve(XtX)
  Xty0 <- crossprod(X, y)
  x1 <- X[,1]
  x1_norm2 <- sum(x1^2)
  XtX_inv_11 <- XtX_inv[1,1]
  a0 <- sum(y^2)
  y0x1 <- sum(y * x1)
  d0 <- Qr %*% as.numeric(crossprod(Qr, y))
  dx <- Qr %*% as.numeric(crossprod(Qr, x1))
  dd0 <- sum(d0^2)
  d0dx <- sum(d0 * dx)
  ddx <- sum(dx^2)
  term20 <- sum(d0 * x1)
  term2x <- sum(dx * x1)
  term1 <- sum(x1^2)

  # Define the function that computes the selective p-value for a given b.
  pselb <- function(b) {

    yb <- y - x1 * b
    Xtyb <- Xty0 - crossprod(X, x1) * b
    beta_hat <- XtX_inv %*% Xtyb
    rss <- sum(yb^2) - as.numeric(crossprod(beta_hat, Xtyb))
    sigma2_hat <- rss / (n - p)
    Fobs <- beta_hat[1]^2 / (sigma2_hat * XtX_inv_11)

    # Corresponding fixed quantities
    a <- a0 - 2 * b * y0x1 + b^2 * term1
    d <- d0 - b * dx
    dd <- dd0 - 2 * b * d0dx + b^2 * ddx
    sqrt_a_dd <- sqrt(a - dd)
    term2 <- term20 - b * term2x

    n_sel <- 0
    n_joint <- 0
    draws <- 0
    iter <- 0

    #if(verbose){cat("We are inside psel-B \n")}
    while (n_sel < min_select && draws < max_draws){
      iter <- iter + 1
      if(verbose & iter > 1){cat("Sampling chunk: ", iter, "\n")}
      m <- as.integer(min(B, max_draws - draws))
      u1 <- rnorm(m)
      v <- rchisq(m, df = k - 1)
      XjVG <- sqrt(lambda) * sqrt_a_dd * u1 / sqrt(u1^2 + v)

      # check selection
      yPy <- dd + 1/lambda * XjVG^2 + 2*b*(term2 + XjVG) + b^2 * term1
      YtY <- a + 2*b*(term2 + XjVG) + b^2 * term1
      Fsel <- (yPy) / (YtY-yPy) # back in Y coordinates, not Y tilde (selection is on Y)
      check_selection <- (Fsel >= cc)

      # check test statistic
      Fval <- ((n-p) * XjVG^2 / lambda) / (a - dd - (1/lambda) * XjVG^2)
      check_Fval <- (Fval >= Fobs)

      # track results
      n_sel <- n_sel + sum(check_selection)
      n_joint <- n_joint + sum(check_selection & check_Fval)
      draws <- draws + m

      if (verbose && iter > 1) cat("Number selected so far: ", n_sel, "\n")
    }

    # Debugging messages for Monte Carlo estimates.
    #cat("Value of b is:", b, "\n")
    #cat("Value of B is:", B, "\n")

    # Compute the numerator and denominator for p-value estimation.
    if (n_sel == 0) return(NA_real_)
    p_sel <- n_joint / n_sel
    return(p_sel)
  }

  return(pselb)
}




#' Retrospective Selective P-Value Based on Summary Statistics
#'
#' This function is useful for conducting valid retrospective F-screening as defined in the 2025 paper "Valid F-screening in linear regression" by McGough, Witten, and Kessler (arxiv preprint: [https://arxiv.org/abs/2505.23113](https://arxiv.org/abs/2505.23113)).
#' Suppose that we have access to the outputs of an "overall" least squares linear regression model, such as from the output of summary(lm(y~X)),
#' and we want to conduct a test of the significance of a single regression coefficient (beta_j) that accounts for the rejection of the "overall" F-test.
#' Then this function can provide a selective p-value for beta_j based on of only a few summary statistics.
#' The arguments of this function include R-squared and residual standard error (RSE) from the overall model (e.g. from summary(lm(y~X))), and a t-statistic for the test of H_0: beta_j=0.
#' This function is especially useful in settings where the raw data is unavailable, such as published studies.
#'
#' @param n Sample size (number of observations).
#' @param p Number of predictors used in the "overall" least squares linear model (excluding the intercept).
#' @param R_squared R-squared from the "overall" fitted least squares linear model (e.g. from summary(lm(y~X))).
#' @param RSE Residual standard error from the "overall" fitted least squares linear model (e.g. from summary(lm(y~X))).
#' @param tstat Observed t-statistic for the post hoc hypothesis test of beta_j.
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
  # some checks
  if (n <= p) {
    stop("Sample size n must be greater than number of predictors p.")
  }
  if (R_squared < 0 || R_squared > 1) {
    stop("R-squared must be between 0 and 1.")
  }
  if (p < 1) {
    stop("Number of predictors p must be at least 1.")
  }
  if (RSE <= 0) {
    stop("Residual standard error RSE must be positive.")
  }
  if (alpha_ov <= 0 || alpha_ov >= 1) {
    stop("Significance level alpha_ov must be between 0 and 1.")
  }
  # compute needed quantities
  Fstat <- tstat^2
  r <- Fstat/(n-p-1)
  F.quantile <- qf(1-alpha_ov, p, (n-p-1))
  d <- RSE^2 * Fstat - RSE^2 * (n-p-1) * R_squared / (1-R_squared)
  a <- RSE^2 * (n-p-1)/(1-R_squared)
  cc <- F.quantile * p / (n-p-1)
  # get variance estimate
  if (is.null(sigma_sq)){
    Zs <- rchisq(100000, n-p-1)
    As <- rchisq(100000, p)
    RSS <- RSE^2 * (n-p-1)
  }
  # compute psel
  n_sel <- 0
  n_joint <- 0
  attempt <- 0
  # make sure total_select is large enough
  while (n_sel < min_select && attempt < max_attempts) {
    attempt <- attempt + 1
    Z <- rchisq(B, n - p - 1)
    W <- rchisq(B, 1)
    ratio <- W / (W+Z)
    check_sel <- ((-d + (a+d)*ratio) / (a + d - (a+d)*ratio) >= cc)
    check_Fval <- (((a+d)*ratio) / (a + d - (a+d)*ratio) >= r)
    n_sel <- n_sel + sum(check_sel)
    n_joint <- n_joint + sum(check_sel & check_Fval)
  }
  pval <- (n_joint / n_sel)
  return(pval)
}





