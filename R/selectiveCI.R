#' Compute Selective Confidence Interval Using Root Finding with uniroot
#'
#' This function calculates the confidence interval (CI) for a selective estimate
#' by finding the values at which the selective p-value function `pselb(b)` equals
#' the significance level `alpha`.
#'
#' @param pselb A function that returns the selective p-value for a given `b`.
#' @param point_est The point estimate of the parameter.
#' @param naive_se_est The naive standard error estimate of the parameter.
#' @param alpha The significance level for the confidence interval (default: 0.05).
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#'
#' @details
#' The function finds the confidence interval bounds by solving for `b` where
#' `pselb(b) = alpha`. This is done using `uniroot()`, which finds the roots
#' of the function `p_diff(b) = pselb(b) - alpha`.
#'
#' Debugging statements are included to check:
#' - The p-values at key points in the search range.
#' - The computed confidence interval bounds.
#' - The p-values at the computed bounds.
#'
#'
#' @export
get_CI <- function(pselb, point_est, naive_se_est, alpha){

  pselb_deterministic <- memoise::memoise(pselb) # rerunning psel(b) with same b will give same output
  p_diff <- function(b) {
    p_val <- pselb_deterministic(b)
    return(p_val - alpha)
  }

  # Find the lower bound of the confidence interval by finding the root of p_diff in a range.
  r <- 10
  b_range_lower <- c(point_est - r * naive_se_est, point_est)
  if (p_diff(point_est) <= 0){
    print("Empty confidence interval!")
    return(c(Inf, -Inf))
  }

  while (p_diff(point_est - r * naive_se_est) >= 0 & r < 1000){
    cat("pval at lower (2): ", p_diff(point_est - r * naive_se_est), "\n")
    r <- r * 2
    if (r >= 1000){
      print("infinite CI!")
      return(c(-Inf, Inf))
    }
    b_range_lower <- c(point_est - r * naive_se_est, point_est)
  }

  root_lower <- uniroot(p_diff, interval = b_range_lower)$root

  # Find the upper bound of the confidence interval by finding the root of p_diff in a range.
  b_range_upper <- c(point_est, point_est + r * naive_se_est)

  while (p_diff(point_est + r * naive_se_est) >= 0 & r < 1000){
    cat("pval at upper: ", p_diff(point_est + r * naive_se_est), "\n")
    r <- r * 2
    if (r >= 1000){
      print("infinite CI!")
      return(c(-Inf, Inf))
    }
    b_range_upper <- c(point_est, point_est + r * naive_se_est)
  }
  root_upper <- uniroot(p_diff, interval = b_range_upper)$root

  # Debuggin
  #cat("Check P-values: ", pselb(point_est - 10 * naive_se_est), pselb(0.01982056), pselb(point_est), pselb(point_est + 10 * naive_se_est), "\n")
  #cat("selective interval: " , c(root_lower, root_upper), "\n")
  #cat("psel at upper CI: ", pselb(root_upper), "\n")
  #cat("psel at lower CI: ", pselb(root_lower), "\n")
  #cat("root lower: ", root_lower, "\n")
  #cat("root upper: ", root_upper, "\n")
  #cat("interval: " , c(root_lower, root_upper), "\n")

  # Combine the lower bound, point estimate, and upper bound into a confidence interval vector.
  CI <- c(root_lower, root_upper)
  return(CI)
}

