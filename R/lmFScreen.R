#' lmFScreen: Selective Inference for Linear Regression
#'
#' This function fits a selective inference linear regression model using a specified formula and dataset.
#' It removes the intercept if present, centers the data, and then calls `Fscreen.fit` for model fitting.
#'
#' @param formula A formula specifying the model (e.g., `y ~ X`).
#' @param data An optional data frame containing the variables in the model.
#' @param alpha The significance level for hypothesis testing (default: 0.05).
#' @param alpha_ov The significance level for the overall F-test (default: 0.05).
#' @param sigma_sq The variance parameter (if NULL, estimated internally using debiased estimate).
#' @param B The number of Monte Carlo samples used for inference (default: 100000).
#' @param seed An optional seed for reproducibility (default: NULL).
#'
#' @return A list containing the selective and naive estimates, confidence intervals, and p-values for the model coefficients.
#'
#' @details
#' This function follows these steps:
#' 1. Converts the input formula into a model frame and extracts `X` and `y`.
#' 2. Centers `X` and `y` if an intercept is detected.
#' 3. Calls `Fscreen.fit` to compute selective inference results.
#'
#' @examples
#' data(mtcars)
#' result <- lmFScreen(mpg ~ wt + hp, data = mtcars)
#' summary(result)
#'
#' @export
lmFScreen <- function(formula, data, alpha = 0.05, alpha_ov = 0.05, sigma_sq = NULL, seed = NULL, compute_CI = TRUE, compute_est = TRUE, B = 100000) {
  # Handle missing data argument (use parent frame like lm)
  mf <- model.frame(formula, data = if (missing(data)) parent.frame() else data)

  # Check if the formula contains an intercept
  intercept <- attr(terms(formula), "intercept")

  # Generate the design matrix
  X <- model.matrix(formula, mf)

  # Extract the response variable (y)
  y <- model.response(mf)

  if(intercept==1) {
    X <- X[,-1, drop = FALSE]
    Xy_centered <- get_Xy_centered(X,y)
    X <- Xy_centered$X
    y <- Xy_centered$y
  }

  output <- lmFScreen.fit(X, y, alpha = alpha, alpha_ov = alpha_ov, test_cols = 1:ncol(X), sigma_sq = sigma_sq, seed = seed, compute_CI = compute_CI, compute_est = compute_est, B = B)

  return(output)
}


#' Fit linear model with F-screening Using Selective Inference
#'
#' This function performs selective inference for linear regression by computing maximum likelihood estimates,
#' confidence intervals, and p-values for each predictor variable.
#'
#' @param X A matrix of predictor variables.
#' @param y A vector of response variables.
#' @param alpha The significance level for hypothesis testing (default: 0.05).
#' @param alpha_ov The significance level for the overall F-test (default: 0.05).
#' @param test_cols Indices of predictors to test (default: all columns of `X`).
#' @param sigma_sq The variance parameter (if NULL, estimated internally).
#' @param B The number of Monte Carlo samples used for inference (default: 100000).
#' @param seed An optional seed for reproducibility (default: NULL).
#'
#' @return A list containing:
#'   \item{coefficients}{Selective estimates for each predictor.}
#'   \item{CIs}{Confidence intervals for each coefficient.}
#'   \item{pvalues}{Selective p-values for hypothesis testing.}
#'   \item{naive_coefs}{Naive coefficient estimates from standard OLS.}
#'
#' @details
#' This function follows these steps:
#' 1. If `sigma_sq` is not provided, estimates it using `get_variance_estimate`.
#' 2. Computes naive estimates using OLS regression.
#' 3. For each predictor, computes the selective MLE using `compute_MLE`.
#' 4. Obtains confidence intervals using `get_CI`.
#' 5. Computes p-values using `get_pselb`.
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10)
#' y <- rnorm(10)
#' result <- Fscreen.fit(X, y)
#' summary(result)
#'
#' @export
lmFScreen.fit <- function(X, y, alpha = 0.05, alpha_ov = 0.05, test_cols = 1:ncol(X), sigma_sq = NULL, seed = NULL, compute_CI = TRUE, compute_est = TRUE, B = 100000) {

  n <- dim(X)[1]
  p <- dim(X)[2]
  U <- svd(X)$u
  yPy <- as.numeric(sum((t(U) %*% y)^2))
  rss <- as.numeric(sum(y^2) - yPy)

  # Check if data passed F-screening
  overall_F_statistic <- (yPy / rss)
  cc <- p/(n-p) * qf(1 - alpha_ov, p, (n - p))
  if (overall_F_statistic < cc) {
    cat("Did not pass Fscreening!\n")
    cat("overall F statistic: ", overall_F_statistic * (n-p)/p, "\n" )
    return(NA)
  }

  if(is.null(seed)) {
    seed <- rnorm(1, 0, 10000)
    set.seed(seed)
  }

  # Compute debiased estimate if estimate of sigma is not given
  if(is.null(sigma_sq)) {
    sigma_sq <- get_variance_estimate(X,y,alpha_ov)
  }

  # Allocate empty vectors/matrices for coefficient estimates, confidence intervals, and pvalues
  beta <- pvalues <- rep(NA, length(test_cols))
  CIs  <- naive_CIs <- matrix(NA, nrow = length(test_cols), ncol=2)
  naive_estimates <- naive_pvalues <- rep(NA, length(test_cols))
  if (!is.null(colnames(X))) {
    names(beta) <- names(pvalues) <- names(naive_estimates) <- names(naive_pvalues) <- colnames(X)[test_cols]
    rownames(CIs) <- rownames(naive_CIs) <- colnames(X)[test_cols]
  } else {
    default_names <- paste0("X", test_cols)
    names(beta) <- names(pvalues) <- names(naive_estimates) <- names(naive_pvalues) <- default_names
    rownames(CIs) <- rownames(naive_CIs) <- default_names
  }

  # Get naive estimates of coefficients, confidence intervals, and pvalues using standard OLS output from lm
  naive_lm_coefs <- summary(lm(y ~ X + 0))$coefficients[test_cols, , drop = FALSE]
  X_orig <- X
  for(col in test_cols){

    # If test_cols is a strict subset of columns of X, need to re-index
    i <- which(test_cols == col)

    # get naive pvalues
    naive_pvalue <- naive_lm_coefs[i, 4]
    naive_pvalues[i] <- naive_pvalue

    # If compute_est == FALSE, we skip to getting pvalues
    if (compute_est){
      X <- X_orig[,c(col,setdiff(1:ncol(X),col))]

      # Save naive estimates
      naive_point_est <- naive_lm_coefs[i, 1]
      naive_estimates[i] <- naive_point_est
      naive_se_est <- naive_lm_coefs[i, 2]

      # If compute_CI == TRUE, save naive confidence intervals
      if (compute_CI){
        naive_CI <- c(
          naive_point_est - qnorm(1 - alpha / 2) * naive_se_est,
          naive_point_est + qnorm(1 - alpha / 2) * naive_se_est
        )
        naive_CIs[i, ] <- naive_CI
      }

      # Obtain the point estimate for beta1 using maximum likelihood.
      # interval uses naive estimates to get an interval to look for conditional MLE in
      interval <- c(naive_point_est - 10 * naive_se_est, naive_point_est + 10 * naive_se_est)
      point_est <- as.numeric(compute_MLE(X, y, sigma_sq = sigma_sq, interval = interval, seed = seed, alpha_ov = alpha_ov, B = B)[[1]])
      beta[i] <- point_est

      # Get selective confidence intervals if compute_CI == TRUE
      pselb <- get_pselb(X, y, sigma_sq = sigma_sq, yPy = yPy, rss = rss, alpha_ov = alpha_ov, seed = seed, B = B, verbose = FALSE)
      if (compute_CI){
        CI <- get_CI(pselb=pselb, point_est=point_est, naive_se_est=naive_se_est, alpha=alpha)
        CIs[i,1] <- CI[1]
        CIs[i,2] <- CI[2]
      }
    }

    # Get pvalue for test of H_0:\beta_1 = 0.
    B_max <- B * 10
    pselb <- get_pselb(X, y, sigma_sq = sigma_sq, yPy = yPy, rss = rss, alpha_ov = alpha_ov, seed = seed, B = B_max, verbose = FALSE)
    P_value <- max(pselb(0), 1/B_max) # return 1/B instead of zero
    pvalues[i] <- P_value
  }

  output <- list(
    "selective coefficients" = beta,
    "selective CIs" = CIs,
    "selective pvalues" = pvalues,
    "naive coefficients" = naive_estimates,
    "naive CIs" = naive_CIs,
    "naive pvalues" = naive_pvalues,
    "alpha" = alpha,
    "alpha_ov" = alpha_ov
  )
  class(output) <- "lmFScreen"
  return(output)
}




