#' @title Deriving the second partial derivatives of the log likelihood function of the linear-binomial model in \code{aRD()} (Hessian matrix)
#' @description \code{aRD()} derives the second partial derivatives of the log likelihood function of the linear-binomial model.
#' @usage hess(theta, y, x)
#' @param theta A numeric vector containing the initial values of the model parameters.
#' @param y A numeric vector containing the dependent variable of the model.
#' @param x The model matrix.
#' @return A numeric matrix containing the second partial derivatives of the log likelihood function of the linear-binomial model (Hessian matrix).
#' @author Thomas Wolf, Julius Johannes Weise, Stefan Wagenpfeil
#' @export


hess <- function(theta, y, x) {
  p <- x %*% theta
  p[p >= 1] <- 1 - 1e-5
  p[p <= 0] <- 1e-5
  s <- ((y - p)^2) / (p^2 * (1 - p)^2)

  # Weight each row of x by sqrt(s)
  Xs <- sweep(x, 1, sqrt(s), "*")

  # Compute Hessian matrix as cross-product of weighted design matrix
  im <- crossprod(Xs)

  colnames(im) <- names(theta)
  rownames(im) <- names(theta)
  return(im)
}
