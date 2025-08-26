test_that("variable selection, no warning & expected number coefficients", {
  # Create data
  set.seed(123)
  x1 <- rnorm(1000, 50, 10)
  x2 <- rnorm(1000, 30, 5)
  x3 <- rnorm(1000, 20, 8)
  x4 <- rnorm(1000, 60, 12)
  x5 <- rnorm(1000, 55, 7)
  logit <- -4 + x1 * 0.04 + x2 * 0.05 + x3 * 0.02 + x4 * -0.03 + x5 * 0.01
  p <- 1 / (1 + exp(-logit)) 
  y <- rbinom(1000, 1, p)  
  data <- data.frame(y, x1, x2, x3, x4, x5)
  
  # Fitting the model
  fit <- aRD(formula = y ~ x1 + x2 + x3, data = data)
  
  # Tests
  result_back <- variable_selection_aRD(fit, alpha = 0.05, print_models = FALSE)
  expect_no_warning(result_back$final_model)
  expect_equal(length(coef(result_back$final_model)), 2)  # 1 predictor + intercept = 2
  expect_equal(length(result_back$model_list), 3)         # List contains three models
  
  result_for <- variable_selection_aRD(fit, alpha = 0.05, selection = "forward", print_models = FALSE)
  expect_no_warning(result_for$final_model)
  expect_equal(length(coef(result_for$final_model)), 2)   # 1 predictor + intercept = 2
  expect_equal(length(result_for$model_list), 1)          # List contains one model
})