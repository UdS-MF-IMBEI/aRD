test_that("aRD converges and returns expected structure", {
  # Beispiel-Daten
  data <- data.frame(
    Y = rep(c(1, 0), each = 3),
    X = rep(c(-1, 0, 1), 2),
    Freq = c(10, 18, 5, 8, 9, 0)
  )

  expanded_data <- data[rep(1:nrow(data), data$Freq), 1:2]

  # Modell fitten
  fit <- aRD(formula = Y ~ X, data = expanded_data, solver = "osqp", verbose = FALSE)

  # Tests
  expect_s4_class(fit, "aRD")             # S4-Klasse korrekt
  expect_true(fit@converged)              # Modell konvergiert
  expect_length(fit@coefficients, 2)      # 2 Koeffizienten (Intercept + X)
})
