#' Summarize an `lmFScreen` Model
#'
#' This function provides a structured summary of an `lmFScreen` model, displaying
#' both **selective** and **naive** estimates, confidence intervals, and p-values.
#'
#' The output includes:
#' - A header summarizing the number of predictors.
#' - A table of **selective estimates** with confidence intervals and p-values.
#' - A table of **naive estimates** with confidence intervals and p-values.
#' - Significance stars appended to p-values:
#'   - `*` for p < 0.05
#'   - `**` for p < 0.01
#'   - `***` for p < 0.001
#'
#' @param object An object of class `"lmFScreen"`, containing model results.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{Selective}{A dataframe with selective estimates, CIs, and p-values.}
#'   \item{Naive}{A dataframe with naive estimates, CIs, and p-values.}
#' }
#'
#' @examples
#' # Assuming `fscreen_model` is an lmFScreen object:
#' summary(fscreen_model)
#'
#' @export
summary.lmFScreen <- function(object, ...) {
  if (!inherits(object, "lmFScreen")) {
    stop("summary.lmFScreen can only be used with objects of class 'lmFScreen'")
  }


  beta <- object[["selective coefficients"]]
  CIs <- object[["selective CIs"]]
  pvalues <- object[["selective pvalues"]]
  alpha <- object$alpha

  # Extract naive values
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
  cat("\n", "lmFScreen Model Summary", "\n")
  cat("--------------------------------------\n")
  cat(sprintf("Number of predictors: %d\n", length(beta)))
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
  cat("Naive Estimates:\n")
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

#' Extract Coefficients from an `lmFScreen` Model
#'
#' This function returns and prints the **selective** and **naive** coefficient estimates.
#'
#' @param object An object of class `"lmFScreen"`, containing model results.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns a dataframe containing:
#' \describe{
#'   \item{Predictor}{The predictor names.}
#'   \item{Estimate (Selective)}{The selective model coefficient estimates.}
#'   \item{Estimate (Naive)}{The naive model coefficient estimates.}
#' }
#'
#' @examples
#' coef(fscreen_model)
#'
#' @export
coef.lmFScreen <- function(object, ...) {
  if (!inherits(object, "lmFScreen")) {
    stop("coef.lmFScreen can only be used with objects of class 'lmFScreen'")
  }

  # Extract coefficients
  selective_est <- object[["selective coefficients"]]
  naive_est <- object[["naive coefficients"]]

  # Print header
  cat("\nlmFScreen Model Coefficients\n")
  cat("--------------------------------------------------\n")
  cat(sprintf("Number of predictors: %d\n", length(selective_est)))
  cat("--------------------------------------------------\n\n")

  # Print table header
  cat("Predictor       Selective.Est   Naive.Est\n")
  cat("---------------------------------------------\n")

  # Print each row formatted
  for (i in seq_along(selective_est)) {
    cat(sprintf(" %-14s %14.6f %12.6f\n",
                names(selective_est)[i],
                selective_est[i],
                naive_est[i]))
  }

  cat("\n")
  invisible(NULL)
}


#' Compute Confidence Intervals for an `lmFScreen` Model
#'
#' This function returns and prints confidence intervals for both **selective** and **naive** estimates.
#'
#' @param object An object of class `"lmFScreen"`, containing model results.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{Selective}{A dataframe with confidence intervals for selective estimates.}
#'   \item{Naive}{A dataframe with confidence intervals for naive estimates.}
#' }
#'
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
  cat("Naive Confidence Intervals:\n")
  cat("Predictor       Lower.CI     Upper.CI\n")
  cat("---------------------------------------------\n")
  for (i in seq_along(beta_names)) {
    cat(sprintf(" %-14s %10.4f    %10.4f\n",
                beta_names[i], naive_CIs[i, 1], naive_CIs[i, 2]))
  }
  cat("\n")

  invisible(list(Selective = CIs_selective, Naive = naive_CIs))
}
