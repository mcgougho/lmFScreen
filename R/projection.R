#' Center X and y by Removing the Intercept
#'
#' This function centers the predictor matrix `X` and response vector `y` by
#' projecting out the intercept. It ensures that the transformed data does not
#' include a mean shift component.
#'
#' @param X A numeric matrix of predictors (dimensions: n x p).
#' @param y A numeric response vector of length n.
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{The transformed predictor matrix with the intercept projected out.}
#'   \item{y}{The transformed response vector with the intercept projected out.}
#' }
#'
#' @details This function uses singular value decomposition (SVD) to project out
#' the intercept and return a centered version of `X` and `y`.
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10)
#' y <- rnorm(10)
#' centered_data <- get_Xy_centered(X, y)
#' X_centered <- centered_data$X
#' y_centered <- centered_data$y
#'
#' @export
get_Xy_centered <- function(X, y) {
  n <- length(y)
  project_out <- cbind(rep(1, n))  # Intercept projection matrix
  svdP <- svd(project_out, nu = nrow(project_out))
  tol <- max(dim(project_out)) * max(svdP$d) * .Machine$double.eps
  r <- sum(svdP$d > tol)
  U_full <- svdP$u
  U_perp <- U_full[, (r+1):ncol(U_full)]
  X <- t(U_perp) %*% X
  y <- t(U_perp) %*% y
  return(list(X = X, y = y))
}

