#' @title Variable Selection (Forward or Backward) for linear-Binomial Models
#' 
#' @description
#' Performs forward or backward variable selection based on Wald test p-values for models estimated using \code{aRD()}. 
#' In each step, a new model is fitted using \code{aRD()}, and variables are added or removed based on the significance level defined by \code{alpha}.
#'
#' @usage variable_selection_aRD(model, selection = c("backward", "forward"), 
#'                               alpha = 0.157, print_models = FALSE, maxit = NULL, 
#'                               verbose = NULL, solver = NULL, eps_abs = NULL, 
#'                               eps_rel = NULL, polish_final = NULL, conswitch = NULL, 
#'                               use_nearPD = NULL)
#'
#' @param model A model object from \code{aRD()} with full data and formula.
#' @param selection Character string, either \code{"backward"} or \code{"forward"}. Determines the direction of model selection. If not specified, backward elimination is performed by default.
#' @param alpha P-value threshold for variable inclusion (forward) or exclusion (backward). Defaults to 0.157, as recommended by Heinze, G., Wallisch, C., & Dunkler, D. (2018).
#' @param print_models Logical; whether to print each model during selection. Defaults to FALSE.
#' @param maxit Maximum number of iterations in the aRD() algorithm. If NULL, defaults to 200L or value from original model call.
#' @param verbose Logical; if TRUE, shows detailed solver output. If NULL, taken from original model or defaults to FALSE.
#' @param solver Solver to be used in aRD(), default is "osqp" unless specified in the original model.
#' @param eps_abs Numeric. Absolute tolerance used as a convergence criterion for the solver. If not specified, taken from original model or defaults to \code{1e-5}.
#' @param eps_rel Numeric. Relative tolerance used as a convergence criterion for the solver. If not specified, taken from original model or defaults to \code{1e-5}.
#' @param polish_final Logical; if TRUE and solver is "osqp", performs a final polish step.
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
#' @references Heinze, G., Wallisch, C., & Dunkler, D. (2018). Variable selection – A review and recommendations for the practicing statistician. Biometrical Journal, 60(3), 431–449.
#' @importFrom stats model.frame model.matrix as.formula
#' 
#' @return An object of class \code{"aRD_selection"}, which is a list containing:
#' \describe{
#'   \item{final_model}{An object of class \code{aRD} representing the final model selected through the variable selection process.}
#'   \item{model_list}{A list of intermediate \code{aRD} model objects fitted during each step of the selection.}
#'   \item{skipped_models}{A named list of models that failed to converge and were skipped during the selection. Each entry includes the attempted formula.}
#'   \item{final_formula}{The final model formula used in the last step.}
#'   \item{EPV}{Estimated events-per-variable (EPV) of the final model, used as a diagnostic for model stability.}
#'   \item{warnings}{Optional warning messages about convergence issues or model stability (e.g., low EPV or skipped variables).}
#' }
#' 
#' @author Julius Johannes Weise, Thomas Wolf, Stefan Wagenpfeil
#' @examples
#' set.seed(123)
#' x1 <- rnorm(500, 50, 10)
#' x2 <- rnorm(500, 30, 5)
#' x3 <- rnorm(500, 40, 8)
#' x4 <- rnorm(500, 60, 12)
#' logit <- (-4 + x1 * 0.04 + x3 * 0.04)
#' p <- 1 / (1 + exp(-logit))
#' y <- rbinom(500, 1, p)
#' df <- data.frame(y, x1, x2, x3, x4)
#' fit <- aRD(formula = y ~ x1 + x2 + x3 + x4, data = df)
#' result <- variable_selection_aRD(fit, selection = "forward", alpha = 0.1)
#' print(result)
#' 
#' @export

variable_selection_aRD <- function(model, selection = c("backward", "forward"),
                                   alpha = 0.157, print_models = FALSE, 
                                   maxit = NULL, verbose = NULL, solver = NULL, 
                                   eps_abs = NULL, eps_rel = NULL, polish_final = NULL, 
                                   conswitch = NULL, use_nearPD = NULL) {
  # Input validation
  selection <- match.arg(selection)
  
  if (!inherits(model, "aRD")) {
    stop("\"model\" must be an object of class \"aRD\" as returned by aRD()")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("\"alpha\" must be a single numeric value between 0 and 1 (exclusive)")
  }
  
  if (!is.logical(print_models) || length(print_models) != 1L) {
    stop("\"print_models\" must be TRUE or FALSE")
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
  
  
  data <- model@data
  full_formula <- model@formula
  
  terms_obj <- terms(full_formula, data = data)
  response <- as.character(attr(terms_obj, "variables")[[2]])
  all_predictors <- attr(terms_obj, "term.labels")
  
  call_list <- as.list(model@call)
  
  # Warnungsliste vorbereiten
  warn_text <- character()
  
  if (is.null(maxit)) {
    maxit <- if (!is.null(call_list$maxit)) call_list$maxit else 100000L
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
  
  model_list <- list()
  model_step <- 1
  
  skipped_models <- list()
  
  #Backward Selection
  if (selection == "backward") {
    
    predictors <- all_predictors
    
    while (length(predictors) > 0) {
      current_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
      current_model <- aRD(current_formula, data, maxit = maxit, verbose = verbose, 
                           solver = solver, eps_abs = eps_abs, eps_rel = eps_rel, 
                           polish_final = polish_final,conswitch = conswitch, use_nearPD = use_nearPD)
      
      model_list[[length(model_list) + 1]] <- current_model
      
      if (!current_model@converged) {
        warning("Model did not converge. Stopping selection.")
        break
      }
      
      if (print_models) {
        cat("\n--- Model step:", model_step, "---\n")
        try(summary(current_model), silent = TRUE)
      }
      
      model_step <- model_step + 1
      
      coefs <- coef(current_model)
      se <- sqrt(diag(solve(aRD::hess(coefs, current_model@y, current_model@x))))
      z <- coefs / se
      p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      p <- p[-1]  
      
      worst_var <- names(p)[which.max(p)]
      
      if (p[worst_var] > alpha) {
        if (length(predictors) > 1) {
          predictors <- setdiff(predictors, worst_var)
        } else {
          warn_text <- c(warn_text, "No variable remained below the p-value threshold.")
          break 
        }
      }
      else {
        break
      }
    }
    
    final_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
  } 
  
  # Forward Selection
  else if (selection == "forward") {
    
    selected_predictors <- c()
    remaining_predictors <- all_predictors
    
    repeat {
      p_values <- c()
      candidate_models <- list()
      
      for (var in remaining_predictors) {
        test_formula <- as.formula(paste(response, "~", paste(c(selected_predictors, var), collapse = " + ")))
        
        test_model <- tryCatch({
          aRD(test_formula, data, maxit = maxit, verbose = verbose,
              solver = solver, eps_abs = eps_abs, eps_rel = eps_rel,
              polish_final = polish_final, conswitch = conswitch, use_nearPD = use_nearPD)
        }, error = function(e) {
          return(NULL)
        })
        
        coefs <- coef(test_model)
        se <- tryCatch({
          sqrt(diag(solve(aRD::hess(coefs, test_model@y, test_model@x))))
        }, error = function(e) {
          return(NULL)
        })
        
        z <- coefs / se
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
        
        p_values[var] <- p[var]
        candidate_models[[var]] <- test_model
        
        if (is.null(test_model) || is.null(se) || !test_model@converged) {
          skipped_models[[paste0("Model step ", model_step, " variable ", var)]] <- list(
            formula = test_formula
          )
          next
        }
      }
      
      if (length(p_values) == 0) {
        if (length(remaining_predictors) == 0) {
          warn_text <- c(warn_text, "All variables have already been included.")
        } else {
          warn_text <- c(warn_text, "No further variables could be added due to failed convergence of all candidate models.")
        }
        break
      }
      
      best_var <- names(p_values)[which.min(p_values)]
      
      if (p_values[best_var] <= alpha) {
        selected_predictors <- c(selected_predictors, best_var)
        remaining_predictors <- setdiff(remaining_predictors, best_var)
        best_model <- candidate_models[[best_var]]
        model_list[[length(model_list) + 1]] <- best_model
        
        if (print_models) {
          cat("\n--- Model step:", model_step, "---\n")
          try(summary(best_model), silent = TRUE)
        }
        model_step <- model_step + 1
      } else {
        break
      }
    }
    
    # Warning if there are no predictors above the threshold
    if (length(selected_predictors) == 0) {
      warning("No variable met the significance criterion or all models failed to converge")
      return(invisible(NULL))
    }
    
    final_formula <- as.formula(paste(response, "~", paste(selected_predictors, collapse = " + ")))
  }
  
  final_model <- aRD(final_formula, data, maxit = maxit, verbose = verbose,
                     solver = solver, eps_abs = eps_abs, eps_rel = eps_rel,
                     polish_final = polish_final, conswitch = conswitch, use_nearPD = use_nearPD)
  
  
  # Stability Investigations
  y <- data[[response]]
  num_events <- sum(y == 1)
  num_predictors <- length(coef(final_model)) - 1
  epv <- num_events / max(1, num_predictors)
  
  
  # Warnung bei niedriger EPV
  if (epv < 25) {
    warn_text <- c(warn_text, "EPV < 25: Model may be unstable.")
  }
  
  # Warnung bei übersprungenen Modellen
  if (length(skipped_models) > 0) {
    skipped_names <- names(skipped_models)
    warn_text <- c(
      warn_text,
      paste0(length(skipped_models), " model(s) were skipped due to non-convergence."))
  }
  
  # Ergebnisobjekt vorbereiten
  out <- list(
    final_model = final_model,
    model_list = model_list,
    skipped_models = skipped_models,
    final_formula = final_formula,
    EPV = epv,
    warnings = if (length(warn_text) > 0) paste(warn_text, collapse = "\n") else NULL
  )
  
  class(out) <- "aRD_selection"
  invisible(out)
  
}