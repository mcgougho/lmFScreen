#' Optimize Beta1 Parameter via Maximum Likelihood
#'
#' This function finds the maximum likelihood estimate (MLE) of `beta1` by optimizing
#' the likelihood function computed via Monte Carlo integration.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of response variables.
#' @param sigma_sq The variance parameter.
#' @param alpha_ov A parameter used in the likelihood calculation.
#' @param interval A numeric vector of length 2 specifying the range for `beta1` optimization (default: `c(-10, 10)`).
#' @param B The number of Monte Carlo samples (default: `1000000`).
#' @param seed A seed value for reproducibility (default: `12345`).
#'
#' @return A list containing:
#'   \item{beta1}{The optimal value of `beta1` that maximizes the likelihood.}
#'   \item{max_likelihood}{The maximum likelihood value computed at `beta1`.}
#'
#' @export
compute_MLE <- function(X, y, sigma_sq,  alpha_ov, interval = c(-10,10), B = 1000000, seed = 12345) {
  # Create the likelihood function based on the current data and parameters.
  # Note: compute_likelihood_function returns a function that computes the likelihood for a given beta1.
  lik_fun <- compute_likelihood_function(X, y, sigma_sq = sigma_sq, alpha_ov = alpha_ov, seed = seed, B = B)

  # Define a negative likelihood function to be minimized.
  neg_likelihood <- function(beta1) {
    -lik_fun(beta1)
  }

  # Use the optimize function to minimize the negative likelihood (thus maximizing the original likelihood).
  opt_result <- optimize(neg_likelihood, interval = interval)

  # Extract the best beta1 value from the optimization result.
  best_beta1 <- opt_result$minimum

  # Calculate the maximum likelihood value at the best beta1.
  best_value <- lik_fun(best_beta1)

  # Return the optimal beta1 and the corresponding maximum likelihood value.
  return(list(beta1 = best_beta1, max_likelihood = best_value))
}


#' Compute the Likelihood Function for Beta1 Optimization
#'
#' This function constructs a likelihood function (to be optimized) using the given data.
#' The likelihood is computed via a Monte Carlo integration method that involves generating
#' random samples from chi-squared distributions.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of response variables.
#' @param sigma_sq The variance parameter.
#' @param alpha_ov The significance level of the overall F-test.
#' @param B The number of Monte Carlo samples (default is 100000).
#' @param seed A seed value for reproducibility (default is 12345).
#'
#' @return A function that computes the (log) likelihood for a given `beta1`.
#'
#' @export
compute_likelihood_function <- function(X, y, sigma_sq, alpha_ov = 0.05, B = 1000000, seed = 12345) {
  n <- dim(X)[1]
  p <- dim(X)[2]
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

  # Define the likelihood function that depends on beta1.
  lik_fun <- function(beta1) {
    set.seed(seed)
    # Compute the mean of the projected beta1 effect.
    mean_w <- as.vector(t(U.1) %*% (X[ , 1] * beta1))
    # Compute the numerator as the log-density under the normal distribution.
    numerator <- dnorm(w_obs, mean = mean_w, sd = sqrt(sigma_sq), log = TRUE)
    # Generate Monte Carlo samples from a noncentral chi-squared distribution.
    W_chi <- rchisq(B, df = 1, ncp = beta1^2 * term1)
    # Generate samples from a central chi-squared distribution.
    Z_chi <- rchisq(B, df = n - p)
    # Evaluate the selection indicator based on the inequality condition.
    select <- (W_chi - cc * Z_chi >= -d / sigma_sq)
    # Compute the denominator by taking the logarithm of the mean of the selection indicator.
    denominator <- log(mean(select))
    # Return the likelihood value as the difference between the numerator and denominator.
    return(numerator - denominator)
  }
  return(lik_fun)
}
