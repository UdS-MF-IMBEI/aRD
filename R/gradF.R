#' @title Deriving the first derivatives of the log likelihood function of the linear-binomial model in \code{aRD()}
#' @description \code{gradF()} derives the first derivatives of the log likelihood function of the linear-binomial model.
#' @usage gradF(theta, y, x)
#' @param theta A numeric vector containing the initial values of the model parameters.
#' @param y A numeric vector containing the dependent variable of the model.
#' @param x The model matrix.
#' @return A numeric vector containing the first derivatives of the log likelihood function of the linear-binomial model.
#' @author Thomas Wolf, Julius Johannes Weise, Stefan Wagenpfeil
#' @export

gradF <- function(theta, y, x) {
  p <- x %*% theta
  p[p >= 1] <- 1 - 1e-5
  p[p <= 0] <- 1e-5
  s <- (y - p) / (p * (1 - p))
  deriv1 <- as.vector(t(x) %*% s)
  return(-deriv1)
}
