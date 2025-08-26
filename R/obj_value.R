#' @title Log-Likelihood Function Value for the linear-binomial model in \code{aRD()}
#' @description
#' \code{obj_value()} computes the value of the log-likelihood function for the linear-binomial model
#' at a given parameter vector \code{theta}.
#' @usage obj_value(theta, y, x)
#' @param theta A numeric vector containing the initial values of the model parameters.
#' @param y A numeric vector containing the dependent variable of the model.
#' @param x The model matrix.
#' @return Numeric scalar representing the value of the log-likelihood function evaluated at \code{theta}.
#' @author Thomas Wolf, Julius Johannes Weise, Stefan Wagenpfeil
#' @export



obj_value <- function(theta, y, x) {
  p <- x %*% theta
  p[p >= 1] <- 1 - 1e-5
  p[p <= 0] <- 1e-5
  sum(y * log(p) + (1 - y) * log(1 - p))
}
