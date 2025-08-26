#' @title Estimating confidence intervals of the estimated model parameters of \code{aRD()}
#' @description For objects of class \code{"aRD"}, \code{confint()} estimates confidence intervals of the estimated model parameters of \code{aRD()}.
#' @param object An object of class \code{"aRD"}.
#' @param parm A specification of which model parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all model parameters are considered.
#' @param level A numeric value that indicates the level of confidence.
#' @param method A character giving the estimation method of the confidence intervals (\code{"bca"} or \code{"wald"}).
#' @param R A positive integer giving the number of bootstrap replicates.
#' @param maxit A positive integer giving the maximum number of iterations.
#' @param verbose Logical; if TRUE, a detailed output is printed, including iteration logs,
#' objective values and solver messages (e.g., from \code{osqp}).
#' @param solver A character string specifying the solver to use. Options are "osqp" or "quadprog".
#' @param verbose Logical; if TRUE, a detailed output is printed, including iteration logs,
#' objective values and solver messages (e.g., from \code{osqp}).
#' @param polish_final Logical; if TRUE and solver is "osqp", performs a final polish step.
#' @param eps_abs Absolute tolerance for solver convergence (only used with "osqp").
#' @param eps_rel Relative tolerance for solver convergence (only used with "osqp").
#' @param conswitch Specifies how the constraint matrix is constructed:
#' \describe{
#'   \item{1 (default)}{Generates all possible combinations of minimum and maximum values for the predictors (excluding the intercept), resulting in \eqn{2^{m-1}} constraints. 
#'   This formulation constrains model predictions within the observed data range, making it suitable for both risk factor identification and prediction (prognosis).}
#'   \item{0}{Uses the raw design matrix \code{x} as the constraint matrix, resulting in \eqn{n} constraints. 
#'   This is primarily suitable for identifying risk factors, but not for prediction tasks, as predictions are not bounded to realistic ranges.}
#' }
#' @param use_nearPD Logical; if \code{TRUE}, the Hessian matrix is projected to the nearest positive definite matrix 
#'   using \code{Matrix::nearPD()} to ensure numerical stability, especially when using solvers that require a positive definite matrix 
#'   (e.g., \code{quadprog}). If \code{FALSE}, the raw Hessian is used directly, which is faster but may lead to numerical issues 
#'   if the matrix is not positive definite. Default is \code{TRUE}.
#' @details \code{confint} provides Wald (default) and bias-corrected accelerated bootstrap confidence intervals of the estimated model parameters of \code{aRD()}.
#' @return A matrix with columns giving the lower and upper confidence limits of each estimated model parameter.
#' @aliases confint,aRD-method
#' @author Thomas Wolf, Julius Johannes Weise, Stefan Wagenpfeil
#' @export

setMethod(f = "confint",
          signature = "aRD",
          definition = function(object, parm, level = .95, method = "wald", R = 1000L, maxit = NULL, verbose = NULL, 
                                solver = NULL, eps_abs = NULL, eps_rel = NULL, polish_final = NULL, conswitch = NULL, use_nearPD = NULL) {
            if (!is.numeric(level)) {
              stop("\"level\" must be a numeric value")
            }

            if (length(level) != 1L) {
              stop("single numeric value for \"level\" expected")
            }

            if (!is.character(method)) {
              stop("\"method\" must be a character string")
            }

            if (length(method) != 1L) {
              stop("single character string for \"method\" expected")
            }

            if (!(method %in% c("bca", "wald"))) {
              stop("\"method\" is misspecified. Currently available confidence interval estimation procedures are: \"bca\" and \"wald\"")
            }

            if (!is.integer(R)) {
              stop("\"R\" must be a positive integer")
            }

            if (length(R) != 1L) {
              stop("single positive integer for \"R\" expected")
            }

            if (R < 1000L) {
              stop("\"R\" must be a positive integer equal to or greater than 1000")
            }
            
            if (!is.null(maxit) && (!is.numeric(maxit) || length(maxit) != 1L || maxit <= 0)) {
              stop("\"maxit\" must be a single positive number or NULL")
            }
            
            if (!is.null(verbose) && (!is.logical(verbose) || length(verbose) != 1L)) {
              stop("\"verbose\" must be TRUE, FALSE, or NULL")
            }
            
            if (!is.null(solver) && (!is.character(solver) || !(solver %in% c("osqp", "quadprog")))) {
              stop("\"solver\" must be one of: \"osqp\", \"quadprog\", or NULL")
            }
            
            if (!is.null(eps_abs) && (!is.numeric(eps_abs) || length(eps_abs) != 1L || eps_abs <= 0 || eps_abs >= 1)) {
              stop("\"eps_abs\" must be a single numeric value between 0 and 1 (exclusive), or NULL")
            }
            
            if (!is.null(eps_rel) && (!is.numeric(eps_rel) || length(eps_rel) != 1L || eps_rel <= 0 || eps_rel >= 1)) {
              stop("\"eps_rel\" must be a single numeric value between 0 and 1 (exclusive), or NULL")
            }
            
            if (!is.null(polish_final) && (!is.logical(polish_final) || length(polish_final) != 1L)) {
              stop("\"polish_final\" must be TRUE, FALSE, or NULL")
            }
            
            if (!is.null(conswitch) && !conswitch %in% c(0, 1)) {
              stop("\"conswitch\" must be 0 or 1")
            }
            
            if (!is.null(use_nearPD) && (!is.logical(use_nearPD) || length(use_nearPD) != 1L)) {
              stop("\"use_nearPD\" must be TRUE or FALSE or NULL")
            }
            

            else {
              cf <- coef(object)
              pnames <- names(cf)
              if (missing(parm)) {
                parm <- pnames
              }
              else if (is.numeric(parm)) {
                parm <- pnames[parm]
              }
              alpha <- (1 - level) / 2
              p <- c(alpha, 1 - alpha)
              ci <- array(data = NA, dim = c(length(parm), 2L), dimnames = list(parm, paste(x = format(x = 100 * p, digits = 3, scientific = FALSE, trim = TRUE), "%", sep = "")))

              if (method == "bca") {
                call_list <- as.list(object@call)
                
                if (is.null(maxit)) {
                  maxit <- if (!is.null(call_list$maxit)) call_list$maxit else 80000L
                }
                
                if (is.null(verbose)) {
                  verbose <- if (!is.null(call_list$verbose)) eval(call_list$verbose, envir = parent.frame()) else FALSE
                }
                
                if (is.null(solver)) {
                  solver <- if (!is.null(call_list$solver)) call_list$solver else "quadprog"
                }
                
                if (is.null(eps_abs)) {
                  eps_abs <- if (!is.null(call_list$eps_abs)) call_list$eps_abs else 1e-5
                }
                
                if (is.null(eps_rel)) {
                  eps_rel <- if (!is.null(call_list$eps_rel)) call_list$eps_rel else 1e-5
                }
                
                if (is.null(polish_final)) {
                  polish_final <- if (!is.null(call_list$polish_final)) call_list$polish_final else TRUE
                }
                if (is.null(conswitch)) {
                  conswitch <- if (!is.null(call_list$conswitch)) call_list$conswitch else 1
                }
                
                if (is.null(use_nearPD)) {
                  use_nearPD <- if (!is.null(call_list$use_nearPD)) call_list$use_nearPD else TRUE
                }
                
                f <- function(data, indices) {
                  
                  dat <- data[indices, ]
                  
                  
                  fit <- aRD(formula = object@formula, data = dat, maxit = maxit, verbose = verbose, 
                             solver = solver, eps_abs = eps_abs, eps_rel = eps_rel, polish_final = polish_final, conswitch = conswitch, use_nearPD = use_nearPD)
                  return(coef(fit)[parm])
                }
                
                b <- boot::boot(data = object@data, statistic = f, R = R)
                ci[] <- matrix(unlist(lapply(1:ncol(b$t), function(i) {boot::boot.ci(b, conf = level, type = "bca", index = i)$bca[4:5]})), ncol = 2, byrow = TRUE)
              }

              if (method == "wald") {
                se <- sqrt(diag(solve(hess(cf, object@y, object@x))))[parm]
                ci[] <- cf[parm] + se %o% stats::qnorm(p)
              }
              return(ci)
            }
          }
)
