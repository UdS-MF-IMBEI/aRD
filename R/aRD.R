#' @title Adjusted Risk differences via specifically penalized likelihood
#' @description \code{aRD()} fits a linear-binomial model using a modified Newton-type algorithm (aRD algorithm) for solving the maximum likelihood estimation problem under linear box constraints.
#' @usage aRD(formula, data, maxit = 80000L, solver = "quadprog",
#' verbose = FALSE, polish_final = TRUE, eps_abs = 1e-5, 
#' eps_rel = 1e-5, conswitch = 1, use_nearPD = TRUE)
#' @param formula An object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param maxit A positive integer giving the maximum number of iterations.
#' @param solver A character string specifying the solver to use. Options are "osqp" () or "quadprog" ().
#' @param verbose Logical; if TRUE, a detailed output is printed, including iteration logs, armijo steps,
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
#' @importFrom osqp osqp osqpSettings
#' @importFrom Matrix Matrix nearPD
#' @importFrom quadprog solve.QP
#' @importFrom stats terms

#' @return An object of S4 class \code{"aRD"} containing the following slots:
#' \item{call}{An object of class \code{"call"}.}
#' \item{formula}{An object of class \code{"formula"}.}
#' \item{coefficients}{A numeric vector containing the estimated model parameters.}
#' \item{iter}{A positive integer indicating the number of iterations.}
#' \item{converged}{A logical constant that indicates whether the model has converged.}
#' \item{y}{A numerical vector containing the dependent variable of the model.}
#' \item{x}{The model matrix.}
#' \item{data}{A data frame containing the variables in the model.}
#' @references Wagenpfeil S, Schöpe J, Bekhit A (2025) Estimation of adjusted relative risks in log-binomial regression using the BSW algorithm. In: Mau J, Mukhin S, Wang G, Xu S (Eds.) Biokybernetika. DE GRUYTER, Berlin, Germany, pp. 665–676.
#' Wagenpfeil S (1991) Implementierung eines SQP-Verfahrens mit dem Algorithmus von Ritter und Best. Diplomarbeit, TUM, Munich, Germany
#' Stellato B, Banjac G, Goulart P, Boyd S, Bansal V (2024). *osqp: Quadratic Programming Solver using the 'OSQP' Library*. R package version 0.6.3.3. \doi{10.32614/CRAN.package.osqp}. Available at: \url{https://CRAN.R-project.org/package=osqp}.
#' Turlach BA, Weingessel A, Moler C (2019). *quadprog: Functions to Solve Quadratic Programming Problems*. R package version 1.5-8. \doi{10.32614/CRAN.package.quadprog}. Available at: \url{https://CRAN.R-project.org/package=quadprog}.
#' @author Thomas Wolf, Julius Johannes Weise, Stefan Wagenpfeil
#' @examples
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n, 50, 1)
#' y <- rbinom(n, 1, -1.5 + 0.04 * x)
#' fit <- aRD(formula = y ~ x, data = data.frame(y = y, x = x), 
#' solver = "quadprog", verbose = TRUE, maxit = 80000L, conswitch = 1)
#' summary(fit)
#' @export


aRD <- function(formula, data, maxit = 80000L, solver = "quadprog", verbose = FALSE, polish_final = TRUE, eps_abs = 1e-5, eps_rel = 1e-5, conswitch = 1, use_nearPD = TRUE) {
  call <- match.call()
  
  #Input validation
  solver <- match.arg(solver, choices = c("osqp", "quadprog"))
  if (!inherits(x = formula, what = "formula")) stop("\"formula\" must be of class \"formula\"")
  if (!is.data.frame(data)) stop("\"data\" must be a data.frame")
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit <= 0L || maxit != as.integer(maxit)) stop("single positive integer for \"maxit\" expected")
  if (!is.logical(verbose) || length(verbose) != 1L) stop("\"verbose\" must be TRUE or FALSE")
  if (!is.logical(polish_final) || length(polish_final) != 1L) stop("\"polish_final\" must be TRUE or FALSE")
  if (!is.null(conswitch) && !conswitch %in% c(0, 1)) {
    stop("\"conswitch\" must be 0 or 1")
  }
  if (!is.null(eps_abs) && (!is.numeric(eps_abs) || length(eps_abs) != 1L ||  eps_abs <= 0 || eps_abs >= 1)) stop("\"eps_abs\" must be a single numeric value between 0 and 1 (exclusive)")
  if (!is.null(eps_rel) && (!is.numeric(eps_rel) || length(eps_rel) != 1L ||  eps_rel <= 0 || eps_rel >= 1)) stop("\"eps_rel\" must be a single numeric value between 0 and 1 (exclusive)")
  if (!is.logical(use_nearPD) || length(use_nearPD) != 1L)
    stop("\"use_nearPD\" must be TRUE or FALSE")
  #Definitions
  data <- stats::model.frame(formula = formula, data = data)
  y <- unname(stats::model.matrix(stats::as.formula(paste("~", all.vars(formula)[1])), data = data)[,-1])
  x <- stats::model.matrix(object = formula, data = data)

  armijo_params <- list(c1 = 1e-4, rho = 0.5, min_alpha = 1e-14)
  theta <- c(mean(y), rep(0, times = ncol(x) - 1))
  Amat <- constr(x, version = conswitch)
  Amat_sparse <- Matrix(Amat, sparse = TRUE)
  l <- rep(0, nrow(Amat))
  u <- rep(1, nrow(Amat))
  converged <- FALSE
  iter <- 0

  #Initialize settings for "osqp"
  if (solver == "osqp") {
    
      settings <- osqpSettings(
      verbose = verbose,
      eps_abs = eps_abs,
      eps_rel = eps_rel,
      max_iter = 1000000L,
      polish = FALSE,
    )
  }
  #Adjust constraint matrix Amat for solver "quadprog"
  else if (solver == "quadprog") {
    Amat_neg <- -Amat
    Amat_t <- rbind(Amat, Amat_neg)
    Amat_Q <- t(Amat_t)
    b <- c(l, -u)
  }
  
  #Main algorithm
  while (!converged && iter < maxit) {
    iter <- iter + 1
    Dmat <- hess(theta, y, x)
    if (use_nearPD) {
      Dmat <- Matrix::nearPD(Dmat)$mat
    }
    Dmat_sparse <- Matrix(Dmat, sparse = TRUE)
    grad <- gradF(theta, y, x)
    dvec <- as.numeric(grad - Dmat %*% theta)

    #Solve quadratic problem
    if (solver == "osqp") {
      if (iter == 1) {
        model <- osqp(P = Dmat_sparse, q = dvec, A = Amat_sparse, l = l, u = u, settings)
      } else {
        model$Update(Px = Dmat_sparse@x, q = dvec)
      }
      model$WarmStart(x = theta)
      fit <- model$Solve()
      p <- fit$x - theta
    } else if (solver == "quadprog") {
      fit <- quadprog::solve.QP(Dmat = Dmat, dvec = -dvec, Amat = Amat_Q, bvec = b)
      p <- fit$solution - theta
    }

    # Armijo line search
    alpha <- 1
    armijo_steps <- 0
    f_old <- obj_value(theta, y, x)
    gTp <- sum(grad * p)
    f_new <- obj_value(theta + alpha * p, y, x)

    while (f_new < f_old + armijo_params$c1 * alpha * gTp) {
      alpha <- alpha * armijo_params$rho
      armijo_steps <- armijo_steps + 1
      if (alpha < armijo_params$min_alpha) stop("Armijo line search failed: step size too small")
      f_new <- obj_value(theta + alpha * p, y, x)
    }

    if (verbose) {
      cat(sprintf("Iter %d | Obj: %.6f -> %.6f | alpha: %.6f | Armijo steps: %d\n",
                  iter, f_old, f_new, alpha, armijo_steps))
    }

    theta_new <- theta + alpha * p
    converged <- all(abs(theta_new - theta) < 1e-4)
    #Additional solution polish after the last iteration
    if (converged && polish_final && solver == "osqp") {
      polish_settings <- osqpSettings(
        verbose = verbose,
        eps_abs = eps_abs,
        eps_rel = eps_rel,
        max_iter = 1000000L,
        polish = TRUE
      )
      model$UpdateSettings(polish_settings)
      model$WarmStart(x = theta)
      theta_polished <- model$Solve()$x
      theta_new <- theta_polished
     
    }

    theta <- theta_new
    names(theta) <- colnames(x)
  }

  if (iter == maxit && !converged) stop("Maximum number of iterations reached without convergence")

  #return object of class "aRD"
  return(methods::new(Class = "aRD", call = call, formula = formula,
                      coefficients = theta, iter = iter, converged = converged,
                      y = y, x = x, data = data))
}
