test_that("Bootstrap process runs successfully with multiple variables", {
  
  # Create example data
  set.seed(123)
  
  # Generate 5 independent variables from normal distributions
  x1 <- rnorm(1000, 50, 10)
  x2 <- rnorm(1000, 30, 5)
  x3 <- rnorm(1000, 40, 8)
  x4 <- rnorm(1000, 60, 12)
  x5 <- rnorm(1000, 55, 7)
  
  # Compute the dependent variable y using a logistic transformation
  logit <- -4 + x1 * 0.04 + x2 * 0.05 + x3 * 0.02 + x4 * -0.03 + x5 * 0.01
  p <- 1 / (1 + exp(-logit))  
  y <- rbinom(1000, 1, p)  
  
  data <- data.frame(y, x1, x2, x3, x4, x5)
  
  # Fit the aRD model with a subset of predictors
  fit <- aRD(formula = y ~ x1 + x2 + x3, data = data)
  
  # Run the bootstrap procedure
  result <- bootaRD(fit, ci_level = 0.90, R = 1000L)
  
  # Extract the bootstrap object
  bootobj <- result$Bootstrap_Object
  
  # Checks 
  expect_true(is.numeric(bootobj$t0))  
  expect_equal(length(bootobj$t[,1]), 1000)  
  expect_gt(sum(complete.cases(bootobj$t)), 950) 
})