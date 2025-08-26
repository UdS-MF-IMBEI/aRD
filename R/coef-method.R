#' @title Extracting the estimated model parameters of \code{aRD()}
#' @description For objects of class \code{"aRD"}, \code{coef()} extracts the estimated model parameters of \code{aRD()}.
#' @param object An object of class \code{"aRD"}.
#' @return A numeric vector containing the estimated model parameters.
#' @aliases coef,aRD-method
#' @author Thomas Wolf, Julius Johannes Weise, Stefan Wagenpfeil
#' @export

setMethod(f = "coef",
          signature = "aRD",
          definition = function(object) {
            return(object@coefficients)
          }
)
