#' Compute the Selective Maximum Likelihood Estimate for beta1
#'
#' This function numerically approximates the conditional maximum likelihood estimate (MLE) of a single regression
#' coefficient using a Monte Carlo approximation to the selective likelihood,
#' conditional on passing the overall F-test.
#'
#' @param X A numeric matrix of predictor variables (n x p), with the first column corresponding to `beta1`.
#' @param y A numeric response vector of length n.
#' @param alpha_ov The significance level for the overall F-test (used in defining the selection region).
#' @param interval A numeric vector of length 2 giving the search interval for `beta1` (default: `c(-10, 10)`).
#' @param B The number of Monte Carlo samples used to approximate the likelihood (default: `1e6`).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta}{The MLE of `beta1` (a single coefficient) under the selective likelihood.}
#'   \item{max_likelihood}{The maximum log-likelihood value achieved at the optimal `beta1`.}
#' }
#'
compute_MLE <- function(X, y, beta_init, sigma_sq_init, alpha_ov, interval_beta = c(-10,10), interval_sigma_sq = c(-10,10), B = 1000000) {
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

  # Create the likelihood function based on the current data and parameters.
  # Note: compute_likelihood_function returns a function that computes the likelihood for a given beta1.
  lik_fun <- compute_likelihood_function(X, y, alpha_ov = alpha_ov, B = B)

  # Define a negative likelihood function to be minimized.
  neg_likelihood <- function(theta) {
    beta1 <- theta[1]
    log_sigma_sq <- theta[2]
    if (beta1 < interval_beta[1] || beta1 > interval_beta[2]) return(Inf)
    if (log_sigma_sq < log(interval_sigma_sq[1]) || log_sigma_sq > log(interval_sigma_sq[2])) return(Inf)
    -lik_fun(theta)
  }

  init <- c(beta1 = beta_init, log_sigma_sq = log(sigma_sq_init))

  # Use the optimize function to minimize the negative likelihood (thus maximizing the original likelihood).
  opt_result <- optim(init, neg_likelihood, method ="Nelder-Mead")

  list(
    beta1 = opt_result$par[1],
    sigma_sq = exp(opt_result$par[2]),
    max_likelihood = lik_fun(opt_result$par),
    optim = opt_result
  )
}



#' Construct the Selective Log-Likelihood Function for beta1
#'
#' This function constructs a log-likelihood function for a single regression
#' coefficient (`beta1`), conditional on selection via the overall F-test.
#' The log-likelihood is approximated using Monte Carlo integration.
#'
#' @param X A numeric matrix of predictors (n x p), with `beta1` corresponding to the first column.
#' @param y A numeric response vector of length n.
#' @param alpha_ov The significance level for the overall F-test (default: `0.05`).
#' @param B The number of Monte Carlo samples (default: `1e6`).
#'
#' @return A function of one argument `beta1` that returns the approximate
#'         log-likelihood, conditional on selection.
#'         If the observed data fails the selection condition, the function returns `-Inf`.
#'
#'
compute_likelihood_function <- function(X, y, alpha_ov = 0.05, B = 1000000) {
  n <- dim(X)[1]
  p <- dim(X)[2]

  if (length(y) != nrow(X)) {
    stop("Length of y must match the number of rows in X.")
  }
  if (any(is.na(y)) || any(is.na(X))) {
    stop("Input data contains NA values. Please handle missing data before calling this function.")
  }
  if (alpha_ov <= 0 || alpha_ov >= 1) {
    stop("Significance level alpha_ov must be between 0 and 1.")
  }

  cc <- (p / (n - p)) * qf(1 - alpha_ov, df1 = p, df2 = (n - p))

  # Computing P_diff = P_X - P_{X_{-1}}
  # Exclude the first column of X for projection calculations.
  X.1 <- X[ , -1]
  # Perform SVD on the remaining columns to get an orthogonal matrix.
  U1 <- svd(X.1)$u
  # Compute the component of the first column orthogonal to the subspace spanned by X.1.
  u.1 <- X[ , 1] - U1 %*% t(U1) %*% X[ , 1]
  # Normalize the orthogonal vector.
  U.1 <- u.1 / sqrt(sum(u.1^2))
  # Create the projection matrix onto the orthogonal direction.
  P_diff <- U.1 %*% t(U.1)
  # Project y
  w_obs <- as.numeric(t(U.1) %*% y)

  # Computing I-P_X
  # QR decomposition
  qrX <- qr(X)
  # Get the basis for the column space
  Q_full <- qr.Q(qrX, complete = TRUE)
  # Column space basis for X (spanned by first p columns of Q)
  U2 <- Q_full[, 1:p, drop = FALSE]
  # Column space basis for I-P_X
  U_perp <- Q_full[, (p+1):ncol(Q_full), drop = FALSE]
  # project y
  z_obs <- as.vector(t(U_perp) %*% y)

  # Check a condition that involves a global variable 'cc'.
  # If the condition fails, print a message and return -Inf to indicate a bad parameter.
  d <- sum((t(U1) %*% y)^2)
  condition <- (t(w_obs) %*% w_obs - cc * t(z_obs) %*% z_obs >= -d)
  if (!condition) {
    print("Condition for passing F-screening not satisfied!")
    return(-Inf)
  }

  # Compute a term used in the noncentrality parameter for chi-squared sampling.
  term1 <- t(X[ , 1]) %*% P_diff %*% X[ , 1]
  W0 <- rnorm(B)
  Z_chi <- rchisq(B, df = n - p)

  # Define the likelihood function that depends on beta1.
  lik_fun <- function(theta) {
    beta1 <- theta[1]
    sigma_sq <- exp(theta[2])
    # Compute the mean of the projected beta1 effect.
    mean_w <- as.vector(t(U.1) %*% (X[ , 1] * beta1))
    numerator <- dnorm(w_obs, mean = mean_w, sd = sqrt(sigma_sq), log = TRUE) + sum(dnorm(z_obs, mean = 0, sd = sqrt(sigma_sq), log = TRUE))
    # denominator
    ncp_sqrt <- beta1 * sqrt(as.numeric(term1) / sigma_sq)
    W_chi <- (W0 + ncp_sqrt)^2
    select <- (W_chi - cc * Z_chi >= -d / sigma_sq)
    denom <- mean(select)
    denom <- pmax(denom, 1e-12)
    log_denom <- log(denom)
    if (!is.finite(denom) || denom <= 0) return(-Inf)
    denominator <- log(denom)
    # Return the likelihood value as the difference between the numerator and denominator.
    return(numerator - denominator)
  }
  return(lik_fun)
}


