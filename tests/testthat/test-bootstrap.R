test_that("Bootstrap process runs successfully", {
  
  # Create example data
  set.seed(123)
  x <- rnorm(100, 50, 10)
  y <- rbinom(100, 1, exp(-4 + x * 0.04))
  
  # Fit an aRD model
  fit <- aRD(formula = y ~ x, data = data.frame(y = y, x = x))
  
  # Run the bootstrap procedure
  result <- bootaRD(fit, ci_level = 0.95, R = 1000L)
  
  # Extract the bootstrap object
  bootobj <- result$Bootstrap_Object
  
  # Check that the original estimates were computed
  expect_true(is.numeric(bootobj$t0))  
  
  # Check that 1000 bootstrap replicates were performed
  expect_equal(length(bootobj$t[,1]), 1000)  
  
  # Check that at least 95% of bootstrap replicates are complete
  expect_gt(sum(complete.cases(bootobj$t)), 950)  
})
