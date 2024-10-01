# Define a simple baseline and treatment pattern
baseline_pattern <- function(X_1, X_2, Z_1, Z_2) {
  return(X_1 + X_2 + Z_1 + Z_2)
}

trt_pattern <- function(X_1, X_2, Z_1, Z_2) {
  return(X_1 * Z_1)
}

# Test RunSimulations function
test_that("runSimulations works when argument lengths match the function definitions", {
  # Prepare the inputs for runSimulations
  sim_iters <- 1
  n_pool <- 5000
  n_test <- 1000

  candsplinevars <- c("X_1", "X_2")
  candbinaryvars <- paste0("Z_", 1:2)
  candinter <- c(candsplinevars, candbinaryvars)

  # Create the test and pool data
  data_test <- data.frame(
    X_1 = runif(n_test, 0, 1),
    X_2 = runif(n_test, 0, 1),
    Z_1 = rbinom(n_test, 1, 0.35),
    Z_2 = rbinom(n_test, 1, 0.5)
  )

  data_pool <- data.frame(
    X_1 = runif(n_pool, 0, 1),
    X_2 = runif(n_pool, 0, 1),
    Z_1 = rbinom(n_pool, 1, 0.35),
    Z_2 = rbinom(n_pool, 1, 0.5),
    trt = rbinom(n_pool, 1, 0.5)
  )

  # MCMC and prior specifications
  mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
  prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))

  # Trial specifications
  trial_specs <- list(alpha = 0.05, B_1 = 0.95, B_2 = 0.8, e_1 = 0, pi_var = 0.1, enrich = TRUE)

  # Call the runSimulations function, expect it to run without error
  result <- runSimulations(sim_iters, data_test, data_pool, baseline_pattern, trt_pattern,
                           candsplinevars, candbinaryvars, candinter, true_tailoring_vars = c("X_1"),
                           interim_n = c(100, 200, 300), sigma_eps = 0.1, mcmc_specs, prior_params, trial_specs, parallel = FALSE)

  # Check the structure of the result
  expect_type(result, "data.frame")
  expect_true("success" %in% names(result))
  expect_true("final_trial_size" %in% names(result))
})

test_that("runSimulations throws error when argument lengths do not match in trt_pattern", {
  # Define an invalid treatment pattern (should require more arguments)
  invalid_trt_pattern <- function(X_1, Z_1, Z_2) {
    return(X_1 * Z_1 * Z_2)
  }

  # Prepare the inputs for runSimulations
  sim_iters <- 10
  n_pool <- 5000
  n_test <- 1000

  candsplinevars <- c("X_1", "X_2")
  candbinaryvars <- paste0("Z_", 1:2)
  candinter <- c(candsplinevars, candbinaryvars)

  data_test <- data.frame(
    X_1 = runif(n_test, 0, 1),
    X_2 = runif(n_test, 0, 1),
    Z_1 = rbinom(n_test, 1, 0.35),
    Z_2 = rbinom(n_test, 1, 0.5)
  )

  data_pool <- data.frame(
    X_1 = runif(n_pool, 0, 1),
    X_2 = runif(n_pool, 0, 1),
    Z_1 = rbinom(n_pool, 1, 0.35),
    Z_2 = rbinom(n_pool, 1, 0.5),
    trt = rbinom(n_pool, 1, 0.5)
  )

  # MCMC and prior specifications
  mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
  prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))

  # Trial specifications
  trial_specs <- list(alpha = 0.05, B_1 = 0.95, B_2 = 0.8, e_1 = 0, pi_var = 0.1, enrich = TRUE)

  # Expect an error due to mismatch in number of arguments for invalid_trt_pattern
  expect_error(runSimulations(sim_iters, data_test, data_pool, baseline_pattern, invalid_trt_pattern,
                              candsplinevars, candbinaryvars, candinter, true_tailoring_vars = c("X_1"),
                              interim_n = c(100, 200, 300), sigma_eps = 0.1, mcmc_specs, prior_params, trial_specs, parallel = FALSE),
               "Error: pattern functions should take 4 arguments, but 3 were provided.")
})

test_that("runSimulations throws error when argument lengths do not match in baseline_pattern", {
  # Define an invalid baseline pattern (should require more arguments)
  invalid_baseline_pattern <- function(X_1, Z_1) {
    return(X_1 + Z_1)
  }

  # Prepare the inputs for runSimulations
  sim_iters <- 10
  n_pool <- 5000
  n_test <- 1000

  candsplinevars <- c("X_1", "X_2")
  candbinaryvars <- paste0("Z_", 1:2)
  candinter <- c(candsplinevars, candbinaryvars)

  data_test <- data.frame(
    X_1 = runif(n_test, 0, 1),
    X_2 = runif(n_test, 0, 1),
    Z_1 = rbinom(n_test, 1, 0.35),
    Z_2 = rbinom(n_test, 1, 0.5)
  )

  data_pool <- data.frame(
    X_1 = runif(n_pool, 0, 1),
    X_2 = runif(n_pool, 0, 1),
    Z_1 = rbinom(n_pool, 1, 0.35),
    Z_2 = rbinom(n_pool, 1, 0.5),
    trt = rbinom(n_pool, 1, 0.5)
  )

  # MCMC and prior specifications
  mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
  prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))

  # Trial specifications
  trial_specs <- list(alpha = 0.05, B_1 = 0.95, B_2 = 0.8, e_1 = 0, pi_var = 0.1, enrich = TRUE)

  # Expect an error due to mismatch in number of arguments for invalid_baseline_pattern
  expect_error(runSimulations(sim_iters, data_test, data_pool, invalid_baseline_pattern, trt_pattern,
                              candsplinevars, candbinaryvars, candinter, true_tailoring_vars = c("X_1"),
                              interim_n = c(100, 200, 300), sigma_eps = 0.1, mcmc_specs, prior_params, trial_specs, parallel = FALSE),
               "Error: pattern functions should take 4 arguments, but 2 were provided.")
})
