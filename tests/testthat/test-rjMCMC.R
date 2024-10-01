test_that("rjMCMC runs successfully with valid inputs", {
  # Example dataset
  set.seed(123)
  n <- 1000
  data <- data.frame(
    X_1 = runif(n, 0, 1),
    Z_1 = rbinom(n, 1, 0.35),
    Z_2 = rbinom(n, 1, 0.5),
    Z_3 = rbinom(n, 1, 0.65),
    Z_4 = rbinom(n, 1, 0.2),
    Z_5 = rbinom(n, 1, 0.35),
    trt = rbinom(n, 1, 0.5)
  )
  data$Y <- 2 * data$Z_1 + 2 * data$Z_1 * data$trt + rnorm(n, 0, 0.1)

  candsplinevars <- c("X_1")
  candbinaryvars <- paste0("Z_", 1:5)
  candinter <- c(candsplinevars, candbinaryvars)

  mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
  prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))

  # Test rjMCMC function
  result <- rjMCMC(data, candsplinevars, candbinaryvars, candinter, mcmc_specs, prior_params)

  # Check that the result is a list and contains expected elements
  expect_type(result, "list")
  expect_true("success" %in% names(result))
  expect_true("trt_eff_posterior" %in% names(result))
  expect_true("vars_prop_summ" %in% names(result))

  # Check for convergence
  expect_true(result$success)

  # Check that output dimensions are correct
  expect_equal(dim(result$trt_eff_posterior), c(n, mcmc_specs$B))
  expect_equal(length(result$vars_prop_summ),
               length(candsplinevars) + length(candbinaryvars) + length(candinter))
})

test_that("rjMCMC throws an error with incorrect inputs", {
  # Invalid dataset (missing 'trt' column)
  invalid_data <- data.frame(
    X_1 = runif(100),
    Y = rnorm(100)
  )

  expect_error(
    rjMCMC(invalid_data, c("X_1"), NULL, NULL),
    "data must contain columns Y and trt"
  )

  # Invalid candidate variables not in data
  data <- data.frame(
    X_1 = runif(100),
    trt = rbinom(100, 1, 0.5),
    Y = rnorm(100)
  )

  expect_error(
    rjMCMC(data, c("X_2"), NULL, NULL),
    "X_2 are not columns of data"
  )
})
