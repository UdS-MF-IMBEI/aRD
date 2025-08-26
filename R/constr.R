#' @title Setting the linear inequality constraints for \code{aRD()}
#' @description \code{constr()} sets the linear inequality constraints for \code{aRD()}.
#' @usage constr(x, version = 1)
#' @param x A model matrix.
#' @param version switch for constraints
#' @return A matrix containing the linear inequality constraints for \code{aRD()}.
#' @author Thomas Wolf, Julius Johannes Weise, Stefan Wagenpfeil
#' @export

constr <- function(x, version = 1) {
  x <- as.matrix(x)
  if (version == 0) {
    # Version 0: design matrix
    Amat <- unname(x)
  } else if (version == 1) {
    # Version 1: Generates constraints by creating all 2^m combinations of min and max values
    # for predictors (excluding intercept), forming linear inequalities that restrict coefficients
    # within the observed data range. Suitable for constraining predictions to realistic values. (default)
    colMax <- matrixStats::colMaxs(x)
    colMin <- matrixStats::colMins(x)
    const <- expand.grid(lapply(2:length(colMax), function(i) {c(colMin[i], colMax[i])}))
    Amat <- unname(as.matrix(cbind(rep(1, times = nrow(const)), const)))
  } else {
    stop("Invalid value for 'version'. Only 0 or 1 are allowed.")
  }
  return(Amat)
}

