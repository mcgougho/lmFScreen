#' Project Out the Intercept from X and y
#'
#' This function projects the predictor matrix `X` and response vector `y` onto the
#' orthogonal complement of the intercept (constant) vector. This removes any component
#' of the data explained by the intercept, effectively centering the data in the
#' sense of removing the mean direction.
#'
#' @param X A numeric matrix of predictors with dimensions (n x p).
#' @param y A numeric response vector of length n.
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{The predictor matrix after projection (dimensions: (n - 1) x p).}
#'   \item{y}{The response vector after projection (length: n - 1).}
#' }
#'
#' @details
#' This function uses the singular value decomposition (SVD) of the intercept vector
#' to construct an orthonormal basis for the space orthogonal to the intercept.
#' It then projects `X` and `y` onto that space, effectively removing any contribution
#' from the intercept.
#'
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 10)
#' y <- rnorm(10)
#' centered_data <- lmFScreen:::get_Xy_centered(X, y)
#' X_centered <- centered_data$X
#' y_centered <- centered_data$y
#'
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

