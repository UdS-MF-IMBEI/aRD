#' @method print aRD_selection
#' @export
print.aRD_selection <- function(x, ...) {
  cat("\nFINALES MODELL:\n\n")
  
  # Hier wird automatisch alles gedruckt wie gewünscht:
  summary(x$final_model, ...)
  
  if (length(x$skipped_models) > 0) {
    cat("\nThe following variables were skipped due to non-convergence or errors:\n")
    skipped_names <- names(x$skipped_models)
    cat(paste("- ", skipped_names, collapse = "\n"), "\n")
  }
  
  cat("\nStability Investigations:\n")
  cat("EPV (Events per Variable):", round(x$EPV, 2), "\n")
  
  if (!is.null(x$warnings) && any(nzchar(x$warnings))) {
    cat("\nWarnings:\n")
    cat(x$warnings, "\n")
  }
  
  
  invisible(x)
}