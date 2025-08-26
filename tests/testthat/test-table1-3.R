test_that("aRD model works correctly for Datensatz 1", {
  data <- data.frame(
    Y = rep(c(1, 0), each = 3),
    X = rep(c(-1, 0, 1), 2),
    Freq = c(10, 18, 5, 8, 9, 0)
  )
  expanded_data <- data[rep(1:nrow(data), data$Freq), 1:2]
  fit <- aRD(Y ~ X, data = expanded_data)
  out <- NULL
  invisible(capture.output({
    out <- summary(fit)
  }))
  expect_equal(unname(coef(fit)[2]), 0.2531358, tolerance = 1e-5)
  expect_equal(unname(out$std.err[2]), 0.13029234, tolerance = 1e-5)
})

test_that("aRD model works correctly for table 2", {
  data2 <- data.frame(
    Y = rep(c(1, 0), each = 3),
    X = rep(c(-1, 0, 1), 2),
    Freq = c(0, 0, 0, 17, 21, 12)
  )
  expanded_data2 <- data2[rep(1:nrow(data2), data2$Freq), 1:2]
  fit2 <- aRD(Y ~ X, data = expanded_data2)
  out2 <- NULL
  invisible(capture.output({
    out2 <- summary(fit2)
  }))

  expect_equal(unname(coef(fit2)[2]), -1.196369e-23, tolerance = 1e-10)
  expect_equal(unname(out2$std.err[2]), 0.1873162, tolerance = 1e-5)
})

test_that("aRD model works correctly for table 3", {
  data3 <- data.frame(
    Y = rep(c(1, 0), each = 3),
    X = rep(c(-1, 0, 1), times = 2),
    Freq = c(2, 14, 2, 2, 3, 17)
  )
  expanded_data3 <- data3[rep(1:nrow(data3), data3$Freq), 1:2]
  fit3 <- aRD(Y ~ X, data = expanded_data3)
  out3 <- NULL
  invisible(capture.output({
    out3 <- summary(fit3)
  }))

  expect_equal(unname(coef(fit3)[2]), -0.3270672, tolerance = 1e-5)
  expect_equal(unname(out3$std.err[2]), 0.0775550, tolerance = 1e-5)
})
