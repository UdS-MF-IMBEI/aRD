test_that("Contingency Table: aRD Provides Correct Estimate", {
  library(aRD)

  # Data
  df <- data.frame(
    y = rep(c(0, 1), each = 250),
    x = rep(c(0, 1, 0, 1), times = c(200, 50, 50, 200))
  )

  # Calculated values
  RD <- 200 / 250 - 50 / 250
  SE <- sqrt(
    (200 / 250) * (1 - 200 / 250) / 250 +
      (50 / 250) * (1 - 50 / 250) / 250
  )
  expected_z <- RD / SE
  expected_CI <- RD + SE * qnorm(c(0.025, 0.975))

  # Model
  fit <- aRD(y ~ x, df)
  out <- NULL
  invisible(capture.output({
    out <- summary(fit)
  }))


  # Tests
  expect_equal(unname(coef(fit)[2]), RD, tolerance = 1e-5)
  expect_equal(unname(out$std.err[2]), SE, tolerance = 1e-5)
  expect_equal(unname(out$z.value[2]), expected_z, tolerance = 1e-5)
  expect_equal(unname(confint(fit)[2, ]), expected_CI, tolerance = 1e-5)
})
