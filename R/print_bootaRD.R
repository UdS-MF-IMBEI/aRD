#' @method print aRD_boot
#' @export
print.aRD_boot <- function(x, ...) {
  cat("Call:\n")
  print(x$Call_aRD)
  
  cat("\nBootstrap Summary:\n")
  cat(x$message, "\n")
  
  cat("\nBootstrap Coefficients:\n")
  print(round(x$Coefficients, 4))
  
  invisible(x)
}