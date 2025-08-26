#' @title Bootstrapping for Risk Differences estimated by \code{aRD()}
#' @description \code{bootaRD()} applies nonparametric bootstrapping to an object of class \code{"aRD"} 
#' and computes bias-corrected accelerated confidence intervals (BCa) for the estimated risk differences.
#' 
#' @usage bootaRD(object, ci_level = 0.95, R = 1000L, maxit = NULL, verbose = NULL, 
#' solver = NULL, eps_abs = NULL, eps_rel = NULL, 
#' polish_final = NULL, conswitch = NULL, use_nearPD = NULL)
#' @param object An object of the class \code{"aRD"}.
#' @param ci_level A value between 0 and 1 indicating the confidence interval.
#' Provides bias-corrected accelerated bootstrap confidence intervals
#' of the original estimated model parameters of \code{aRD()}.
#' @param R A positive integer greater than or equal to 1000 giving the number of bootstrap replicates.
#' @param maxit A positive integer giving the maximum number of iterations in the \code{aRD()} algorithm.
#' If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#' @param verbose Logical; if TRUE, a detailed output is printed, including iteration logs,
#' objective values and solver messages (e.g., from \code{osqp}).
#' If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#' @param solver A character string specifying the solver to use. Options are "osqp" or "quadprog".
#' If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#' @param eps_abs Absolute tolerance for solver convergence (only used with "osqp").
#' If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#' @param eps_rel Relative tolerance for solver convergence (only used with "osqp").
#' If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#' @param polish_final Logical; if TRUE and solver is "osqp", performs a final polish step.
#' If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#' @param conswitch Specifies how the constraint matrix is constructed:
#' \describe{
#'   \item{1 (default)}{Generates all possible combinations of minimum and maximum values for the predictors (excluding the intercept), resulting in \eqn{2^{m-1}} constraints. 
#'   This formulation constrains model predictions within the observed data range, making it suitable for both risk factor identification and prediction (prognosis).}
#'   \item{0}{Uses the raw design matrix \code{x} as the constraint matrix, resulting in \eqn{n} constraints. 
#'   This is primarily suitable for identifying risk factors, but not for prediction tasks, as predictions are not bounded to realistic ranges.}
#' }
#' If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#' @param use_nearPD Logical; if \code{TRUE}, the Hessian matrix is projected to the nearest positive definite matrix 
#'   using \code{Matrix::nearPD()} to ensure numerical stability, especially when using solvers that require a positive definite matrix 
#'   (e.g., \code{quadprog}). If \code{FALSE}, the raw Hessian is used directly, which is faster but may lead to numerical issues 
#'   if the matrix is not positive definite. Default is \code{TRUE}.
#'   If \code{NULL} (the default), the value stored in the \code{aRD} object is passed internally to \code{bootaRD()}.
#'   
#' @return An object of class \code{"aRD_boot"}, which is a list containing:
#' \describe{
#'   \item{Call_aRD}{The original call to the \code{aRD()} function used to fit the model.}
#'   \item{Successful_Bootstraps}{The number of bootstrap replicates that were completed successfully.}
#'   \item{message}{A character string with a status message indicating how many bootstrap samples succeeded.}
#'   \item{Coefficients}{A matrix with the original estimated model parameters (Orig. Est.), the mean of the bootstrap estimates (Boot. Est.), 
#' the standard error of the bootstrap estimates (Boot. SE), the difference between the bootstrap mean and the original estimate (bias),
#' the Risk Difference (equal to the estimate; RD), and the bias-corrected accelerated confidence intervals at the specified level.}
#'   \item{Bootstrap_Object}{An object of class \code{"boot"} (from the \pkg{boot} package) containing the full bootstrap output, 
#'     including replicates and metadata. This can be used for further analyses or plotting.}
#' }
#' 
#' @importFrom boot boot boot.ci
#' @importFrom checkmate assert_class assert_integer assert_numeric
#' @importFrom aRD aRD
#' @importFrom stats complete.cases sd
#' 
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- rnorm(100, 50, 10)
#' y <- rbinom(100, 1, exp(-4 + x * 0.04))
#' library(aRD)
#' fit <- aRD(formula = y ~ x, data = data.frame(y = y, x = x))
#' result <- bootaRD(fit, ci_level = 0.90)
#' print(result)
#' }
#' 
#' @author Julius Johannes Weise, Thomas Wolf, Stefan Wagenpfeil
#' @export

#Definition of the new function bootaRD
bootaRD <- function(object, ci_level = 0.95, R = 1000L, maxit = NULL, verbose = NULL, 
                    solver = NULL, eps_abs = NULL, eps_rel = NULL, polish_final = NULL, conswitch = NULL, use_nearPD = NULL) {
  
  #Input validation
  checkmate::assert_class(object, "aRD")
  checkmate::assert_numeric(ci_level, lower = 0, upper = 1, any.missing = FALSE, len = 1)
  if (ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be strictly between 0 and 1 (exclusive)")
  }
  checkmate::assert_integer(R, lower = 1000, any.missing = FALSE, len = 1)
  
  if (!is.null(maxit)) {
    checkmate::assert_integer(maxit, lower = 200L, any.missing = FALSE, len = 1)
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
  
  
  #Determine number of coefficients
  coef_len <- length(object@coefficients)
  
  #Import aRD model information
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
  
  
  #Bootstrap help function
  f <- function(data, indices) {
    dat <- data[indices, ]
    model <- tryCatch(aRD(formula = object@formula, data = dat, maxit = maxit, verbose = verbose, 
                          solver = solver, eps_abs = eps_abs, eps_rel = eps_rel, 
                          polish_final = polish_final, conswitch = conswitch, use_nearPD = use_nearPD), error = function(e) NULL)
    if (is.null(model)) return(rep(NA, coef_len))
    return(coef(model))
  }
  
  #Perform bootstrap
  fit <- boot::boot(data = object@data, statistic = f, R = R)
  
  
  #Calculation of confidence intervals for each bootstrap result
  ci <- lapply(1:ncol(fit$t), function(i) {
    ci_result <- boot::boot.ci(fit, conf = ci_level, type = "bca", index = i)
    if (is.null(ci_result$bca)) {
      return(c(NA, NA))
    } else {
      return(ci_result$bca[4:5]) #Extract the lower (4) and upper (5) confidence interval
    }
  })
  
  #List of confidence intervals in a matrix format
  ci_matrix <- do.call(rbind, ci)
  
  #Extract the percentage level
  ci_level_lower <- (1 - ci_level) / 2 * 100  #Lower percentage level
  ci_level_upper <- 100 - ci_level_lower      #Upper percentage level
  
  #Dynamic naming of the confidence interval columns
  colnames(ci_matrix) <- c(
    paste(format(ci_level_lower, digits = 3, scientific = FALSE, trim = TRUE), "%", sep = ""),
    paste(format(ci_level_upper, digits = 3, scientific = FALSE, trim = TRUE), "%", sep = "")   
  )
  
  #Calculation of the bootstrap results
  bootstrap_replicates <- fit$t
  successful_bootstraps <- sum(complete.cases(bootstrap_replicates))
  
  #Calculation of mean values and standard errors
  original_Est <- fit$t0
  bootstrap_means <- colMeans(bootstrap_replicates, na.rm = TRUE)
  bias <- bootstrap_means - original_Est
  bootstrap_se <- apply(bootstrap_replicates, 2, sd, na.rm = TRUE)
  
  #Output success message
  message_text <- if (successful_bootstraps == R) {
    sprintf("%d of %d bootstraps completed successfully.", R, R)
  } else {
    sprintf("ATTENTION: Only %d out of %d bootstraps were completed successfully.
            This may be due to a restriction in one or more variables in the model.",
            successful_bootstraps, R)
  }
  
  #Create table
  coef.table <- cbind(as.matrix(original_Est), as.matrix(bootstrap_means), as.matrix(bootstrap_se), as.matrix(bias), as.matrix(bootstrap_means), ci_matrix)
  colnames(coef.table) <- c("Orig. Est.", "Boot. Est.", "Boot. SE", "Bias","Boot. RD", colnames(ci_matrix))
  
  
  result <- list(
    Call_aRD = object@call,
    Successful_Bootstraps = successful_bootstraps,
    message = message_text,
    Coefficients = coef.table,
    Bootstrap_Object = fit
  )
  
  class(result) <- "aRD_boot"
  
  return(result)
}