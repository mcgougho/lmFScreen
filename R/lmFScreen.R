#' lmFScreen: Selective Inference for Linear Regression via F-screening
#'
#' Fits a selective inference linear regression model using an R formula interface.
#' If an intercept is present in the model, the data are projected to remove the intercept
#' before conducting inference.
#'
#' @param formula A formula specifying the linear model (e.g., y ~ x1 + x2).
#' @param data An optional data frame containing the variables in the model.
#' @param alpha Significance level for confidence intervals and hypothesis tests (default: 0.05).
#' @param alpha_ov Significance level for the overall F-test used for screening (default: 0.05).
#' @param sigma_sq Optional noise variance. If NULL, it is estimated using a corrected residual variance.
#' @param seed Optional seed for reproducibility.
#' @param compute_CI Logical; whether to compute selective confidence intervals (default: TRUE).
#' @param compute_est Logical; whether to compute selective point estimates (default: TRUE).
#' @param B Number of Monte Carlo samples used for selective inference (default: 100000).
#'
#' @return An object of class lmFScreen, which includes:
#'   - Selective coefficients, confidence intervals, and p-values
#'   - Naive (unadjusted) estimates, confidence intervals, and p-values
#'   - Model-level settings such as alpha and alpha_ov
#'
#' @details
#' This function performs the following steps:
#' 1. Converts the formula into a model matrix and response vector.
#' 2. Projects out the intercept if one is included in the formula.
#' 3. Calls lmFScreen.fit to compute selective inference results for all predictors.
#'
#' @examples
#' data(mtcars)
#' result <- lmFScreen(mpg ~ wt + hp, data = mtcars)
#' summary(result)
#' coef(result)
#' confint(result)
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


#' Fit a Linear Model with Selective Inference After F-screening
#'
#' Performs selective inference in linear regression by estimating selective
#' coefficients, confidence intervals, and p-values conditional on passing an overall F-test.
#'
#' @param X A numeric matrix of predictors.
#' @param y A numeric response vector.
#' @param alpha Significance level for confidence intervals and hypothesis tests (default: 0.05).
#' @param alpha_ov Significance level for the overall F-test used for screening (default: 0.05).
#' @param test_cols Indices of predictors to test (default: all columns of X).
#' @param sigma_sq Optional noise variance. If NULL, it is estimated using a corrected residual variance.
#' @param B Number of Monte Carlo samples used for selective inference (default: 100000).
#' @param seed Optional seed for reproducibility.
#' @param compute_CI Logical; whether to compute selective confidence intervals (default: TRUE).
#' @param compute_est Logical; whether to compute selective point estimates (default: TRUE).
#'
#' @return A list of class lmFScreen containing:
#'   - Selective coefficients, confidence intervals, and p-values
#'   - Naive (OLS) coefficients, confidence intervals, and p-values
#'   - Model settings: alpha and alpha_ov
#'
#' @details
#' The procedure follows these steps:
#' 1. Computes the overall F-statistic and checks whether the data passes the F-screening threshold.
#' 2. If passed, computes naive OLS estimates and standard errors.
#' 3. For each predictor:
#'    - Computes the selective MLE using Monte Carlo approximation
#'    - Constructs a selective confidence interval using root-finding
#'    - Computes a selective p-value using selective resampling
#'
#' The function skips selective estimation if compute_est is FALSE, and skips confidence interval construction if compute_CI is FALSE.
#'
#' @examples
#' data(mtcars)
#' X <- cbind(mtcars$wt, mtcars$hp)
#' y <- mtcars$mpg
#' svdP <- svd(rep(1,nrow(mtcars)), nu = nrow(mtcars))
#' tol <- nrow(mtcars) * max(svdP$d) * .Machine$double.eps
#' r <- sum(svdP$d > tol)
#' U_full <- svdP$u
#' U_perp <- U_full[, (r+1):ncol(U_full)]
#' X <- t(U_perp) %*% X
#' y <- t(U_perp) %*% y
#' result <- lmFScreen.fit(X,y)
#' summary(result)
#' coef(result)
#' confint(result)
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




