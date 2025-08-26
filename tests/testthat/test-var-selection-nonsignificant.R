# Example Data
set.seed(123)
x1 <- rnorm(100, 50, 10)
x2 <- rnorm(100, 30, 5)
x3 <- rnorm(100, 40, 8)
y <- rbinom(100, 1, 0.50)

df <- data.frame(y, x1, x2, x3)

# Fitting the model
fit <- aRD(formula = y ~ x1 + x2 + x3, data = df)

# Test for backward and forward selection
test_that("Expect warning", {
  
  expect_no_warning({
    variable_selection_aRD(fit, alpha = 0.05, print_models = F)
  })
  
  expect_warning({
    variable_selection_aRD(fit, alpha = 0.05, selection = "forward", print_models = F)
  })
})