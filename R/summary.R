#' Summarize an lmFScreen Model
#'
#' Displays a structured summary of an lmFScreen model, including both selective
#' and naive estimates, confidence intervals, and p-values.
#'
#' The output includes:
#'
#' - A header indicating the number of predictors
#'
#' - A table of selective estimates with confidence intervals and p-values
#'
#' - A table of naive estimates with confidence intervals and p-values
#'
#' - Significance stars are added to p-values:
#'   * for p < 0.05
#'   ** for p < 0.01
#'   *** for p < 0.001
#'
#' @param object An object of class lmFScreen containing model results.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list containing:
#'
#'   - Selective: A data frame with selective estimates, confidence intervals, and p-values
#'
#'   - Standard: A data frame with standard estimates, confidence intervals, and p-values
#'
#' @examples
#' data(mtcars)
#' result <- lmFScreen(mpg ~ wt + hp, data = mtcars)
#' summary(result)
#'
#' @export
summary.lmFScreen <- function(object, ...) {
  if (!inherits(object, "lmFScreen")) {
    stop("summary.lmFScreen can only be used with objects of class 'lmFScreen'")
  }

  overall_Fstat <- object[["overall F-stat"]]

  # Extract selective values
  beta <- object[["selective coefficients"]]
  CIs <- object[["selective CIs"]]
  pvalues <- object[["selective pvalues"]]
  alpha <- object$alpha

  # Extract standard values
  naive_estimates <- object[["naive coefficients"]]
  naive_CIs <- object[["naive CIs"]]
  naive_pvalues <- object[["naive pvalues"]]

  # Function to add significance stars
  significance_stars <- function(p) {
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "")))
  }

  # Print summary
  cat("lmFScreen Model Summary", "\n")
  cat("--------------------------------------\n")
  cat(sprintf("Overall F-statistic:   %8.4f\n", overall_Fstat))
  cat("--------------------------------------\n\n")
  cat(sprintf("Number of post hoc tests: %d\n", length(beta)))
  cat("--------------------------------------\n\n")

  # Print Selective Estimates
  cat("Selective Estimates:\n")
  cat("Predictor       Estimate     Lower.CI    Upper.CI    P-value\n")
  cat("-------------------------------------------------------------\n")
  for (i in seq_along(beta)) {
    cat(sprintf(" %-12s %10.6f  %10.4f  %10.4f  %10.4f%s\n",
                names(beta)[i], beta[i], CIs[i, 1], CIs[i, 2], pvalues[i], significance_stars(pvalues[i])))
  }
  cat("\n")

  # Print Naive Estimates
  cat("Standard Estimates:\n")
  cat("Predictor       Estimate     Lower.CI    Upper.CI    P-value\n")
  cat("-------------------------------------------------------------\n")
  for (i in seq_along(beta)) {
    cat(sprintf(" %-12s %10.6f  %10.4f  %10.4f  %10.4f%s\n",
                names(naive_estimates)[i], naive_estimates[i], naive_CIs[i, 1], naive_CIs[i, 2], naive_pvalues[i], significance_stars(naive_pvalues[i])))
  }
  cat("\n")

  cat("\nSignificance levels: * < 0.05  ** < 0.01  *** < 0.001\n")

  invisible(NULL)
}




#' Extract Coefficients from an lmFScreen Model
#'
#' This S3 method prints and returns the selective and naive coefficient estimates
#' from an lmFScreen model.
#'
#' @param object An object of class lmFScreen.
#' @param ... Currently unused.
#'
#' @return Invisibly returns a data frame containing:
#'
#' - Predictor: Names of the predictor variables
#'
#' - Estimate (Selective): Selective coefficient estimates
#'
#' - Estimate (Naive): Naive coefficient estimates
#'
#'
#' @examples
#' data(mtcars)
#' mod <- lmFScreen(mpg ~ wt + hp, data = mtcars)
#' coef(mod)
#'
#'
#' @export
#' @method coef lmFScreen
coef.lmFScreen <- function(object, ...) {
  if (!inherits(object, "lmFScreen")) {
    stop("coef.lmFScreen can only be used with objects of class 'lmFScreen'")
  }

  selective_est <- object[["selective coefficients"]]
  naive_est <- object[["naive coefficients"]]

  out <- data.frame(
    Predictor = names(selective_est),
    Selective.Est = as.numeric(selective_est),
    Naive.Est = as.numeric(naive_est),
    row.names = NULL
  )

  cat("\nCoefficients from lmFScreen\n----------------------------\n")
  print(out)

  invisible(out)
}





#' Compute Confidence Intervals for an lmFScreen Model
#'
#' Prints and returns confidence intervals for both selective and naive estimates.
#' The confidence level is taken from the alpha value stored in the lmFScreen object.
#'
#' @param object An object of class lmFScreen containing model results.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list containing:
#'
#'   - Selective: A data frame of confidence intervals for selective estimates
#'
#'   - Standard: A data frame of confidence intervals for standard estimates
#'
#'
#' @examples
#' data(mtcars)
#' result <- lmFScreen(mpg ~ wt + hp, data = mtcars)
#' confint(result)
#'
#' @export
confint.lmFScreen <- function(object, ...) {
  if (!inherits(object, "lmFScreen")) {
    stop("confint.lmFScreen can only be used with objects of class 'lmFScreen'")
  }

  # Use default level from object if not provided
  level <- 1 - object$alpha


  alpha <- 1 - level

  # Extract confidence intervals
  CIs_selective <- object[["selective CIs"]]
  beta_names <- names(object[["selective coefficients"]])

  # Extract naive intervals (already precomputed)
  naive_CIs <- object[["naive CIs"]]

  # Print header
  cat("\nlmFScreen Model Confidence Intervals\n")
  cat("------------------------------------------------------\n")
  cat(sprintf("Confidence Level: %.2f%%\n", level * 100))
  cat(sprintf("Number of predictors: %d\n", length(beta_names)))
  cat("------------------------------------------------------\n\n")

  # Print Selective Intervals
  cat("Selective Confidence Intervals:\n")
  cat("Predictor       Lower.CI     Upper.CI\n")
  cat("---------------------------------------------\n")
  for (i in seq_along(beta_names)) {
    cat(sprintf(" %-14s %10.4f    %10.4f\n",
                beta_names[i], CIs_selective[i, 1], CIs_selective[i, 2]))
  }
  cat("\n")

  # Print Naive Intervals
  cat("Standard Confidence Intervals:\n")
  cat("Predictor       Lower.CI     Upper.CI\n")
  cat("---------------------------------------------\n")
  for (i in seq_along(beta_names)) {
    cat(sprintf(" %-14s %10.4f    %10.4f\n",
                beta_names[i], naive_CIs[i, 1], naive_CIs[i, 2]))
  }
  cat("\n")

  invisible(list(Selective = CIs_selective, Naive = naive_CIs))
}





#' Print an lmFScreen Fit
#' @description print coefficients from lmFScreen model
#' @param object An object of class lmFScreen containing model results.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list containing:
#'
#'   - Selective: A data frame with selective estimates
#'
#'   - Standard: A data frame with standard estimates
#'
#' @examples
#' data(mtcars)
#' result <- lmFScreen(mpg ~ wt + hp, data = mtcars)
#' print(result)
#'
#' @export
print.lmFScreen <- function(object, ...) {
  if (!inherits(object, "lmFScreen")) {
    stop("summary.lmFScreen can only be used with objects of class 'lmFScreen'")
  }

  overall_Fstat <- object[["overall F-stat"]]

  # Extract selective values
  beta <- object[["selective coefficients"]]

  # Print Selective Estimates
  cat("Selective Estimates:\n")
  cat("Predictor       Estimate \n")
  cat("---------------------------\n")
  for (i in seq_along(beta)) {
    cat(sprintf(" %-12s %10.6f  \n",
                names(beta)[i], beta[i]))
  }
  cat("\n")

  invisible(NULL)
}







