test_that("Test for correct input validation in bootaRD", {
  # Create an example of an aRD object
  set.seed(123)
  x <- rnorm(100, 50, 10)
  y <- rbinom(100, 1, exp(-4 + x * 0.04))
  
  fit <- aRD(formula = y ~ x, data = data.frame(y = y, x = x))
  
  # Test for invalid ci_level (less than 0 or greater than 1)
  expect_error(bootaRD(fit, ci_level = 0, R = 1000L), "ci_level")
  expect_error(bootaRD(fit, ci_level = 1, R = 1000L), "ci_level")
  
  # Test for invalid R (less than 1000)
  expect_error(bootaRD(fit, ci_level = 0.95, R = 500L), "Assertion on 'R' failed: Element 1 is not >= 1000.")
  expect_error(bootaRD(fit, ci_level = 0.95, R = 1000))
})
