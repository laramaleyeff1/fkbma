# Test interimAnalysis when it successfully runs the rjMCMC function
test_that("interimAnalysis runs successfully with valid inputs", {
  # Mock rjMCMC to simulate trial results
  rjMCMC <- list(
    success = TRUE,
    trt_eff_posterior = matrix(rnorm(100), nrow = 10, ncol = 10),
    vars_prop_summ = c(a = 0.1, b = 0.05, c = 0.2),
    included_vars = c("a", "b")
  )
  # Example data setup
  data <- data.frame(
    X_1 = runif(100),
    X_2 = runif(100),
    Z_1 = rbinom(100, 1, 0.5),
    Z_2 = rbinom(100, 1, 0.5),
    trt = rbinom(100, 1, 0.5)
  )
  candsplinevars <- c("X_1", "X_2")
  candbinaryvars <- c("Z_1", "Z_2")
  candinter <- c("X_1", "Z_1")

  # Mock MCMC and trial specs
  mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
  prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))
  trial_specs <- list(alpha = 0.05, B_1 = 0.95, B_2 = 0.8, e_1 = 0, pi_var = 0.1, enrich = TRUE)

  # Call interimAnalysis
  result <- interimAnalysis(data, candsplinevars, candbinaryvars, candinter, last = FALSE,
                            mcmc_specs = mcmc_specs, prior_params = prior_params, trial_specs = trial_specs, mcmc_results = rjMCMC)

  # Check the results
  expect_type(result, "list")
  expect_true(result$success)
  expect_true("stop_efficacy" %in% names(result))
  expect_true("stop_futility" %in% names(result))
  expect_true("prop_eff" %in% names(result))
})

# Test interimAnalysis when the MCMC fails
test_that("interimAnalysis handles MCMC failure", {
  # Define an rjMCMC function that simulates failure
  rjMCMC <- list(success = FALSE)


  # Example data setup
  data <- data.frame(
    X_1 = runif(100),
    X_2 = runif(100),
    Z_1 = rbinom(100, 1, 0.5),
    Z_2 = rbinom(100, 1, 0.5),
    trt = rbinom(100, 1, 0.5)
  )
  candsplinevars <- c("X_1", "X_2")
  candbinaryvars <- c("Z_1", "Z_2")
  candinter <- c("X_1", "Z_1")

  # Mock MCMC and trial specs
  mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
  prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))
  trial_specs <- list(alpha = 0.05, B_1 = 0.95, B_2 = 0.8, e_1 = 0, pi_var = 0.1, enrich = TRUE)

  # Call interimAnalysis
  result <- interimAnalysis(data, candsplinevars, candbinaryvars, candinter, last = FALSE,
                            mcmc_specs = mcmc_specs, prior_params = prior_params, trial_specs = trial_specs,
                            mcmc_results = rjMCMC)

  # Expect the function to handle the failure and return success = FALSE
  expect_type(result, "list")
  expect_false(result$success)
})
