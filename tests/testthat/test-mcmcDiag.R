test_that("mcmcDiag correctly computes R-hat diagnostics", {
  # Example data
  n <- 1000
  data <- data.frame(
    X_1 = runif(n, 0, 1),
    Z_1 = rbinom(n, 1, 0.35),
    Z_2 = rbinom(n, 1, 0.5),
    trt = rbinom(n, 1, 0.5)
  )
  data$Y <- 2 * data$Z_1 + 2 * data$Z_1 * data$trt + rnorm(n, 0, 0.1)

  # Define candidate variables
  candsplinevars <- c("X_1")
  candbinaryvars <- paste0("Z_", 1:2)
  candinter <- c(candsplinevars, candbinaryvars)

  # Simulate trial results using rjMCMC function
  mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
  prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))
  trial_results <- rjMCMC(data, candsplinevars, candbinaryvars, candinter, mcmc_specs, prior_params)

  # Test mcmcDiag function
  diag_results <- mcmcDiag(trial_results)

  # Check that the result is a list and contains expected elements
  expect_type(diag_results, "list")
  expect_true("Rhat_trt_eff_posterior" %in% names(diag_results))
  expect_true("Rhat_inter" %in% names(diag_results))
  expect_true("Rhat_trt" %in% names(diag_results))

  # Check that R-hat values are within acceptable range (1 is perfect, typically < 1.1 is acceptable)
  expect_true(all(diag_results$Rhat_trt_eff_posterior > 0))
  expect_true(diag_results$Rhat_inter > 0)
  expect_true(diag_results$Rhat_trt > 0)
  expect_true(all(diag_results$Rhat_trt_eff_posterior < 1.1))
  expect_true(diag_results$Rhat_inter < 1.1)
  expect_true(diag_results$Rhat_trt < 1.1)

  # Check that the diagnostic data frames are correctly structured
  expect_true("trt_eff_posterior_plot_df" %in% names(diag_results))
  expect_true("inter_trt_param_plot_df" %in% names(diag_results))
  expect_s3_class(diag_results$trt_eff_posterior_plot_df, "data.frame")
  expect_s3_class(diag_results$inter_trt_param_plot_df, "data.frame")
})
