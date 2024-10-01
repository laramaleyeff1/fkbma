
library(doParallel)
library(plyr)
library(purrr)

sim_iters = 10

# n_pool: Size of the population pool from which we sample during the trial
n_pool = 5000
n_test = 10000

# Names of the candidate continuous and binary variables
candsplinevars = c("X_1")
candbinaryvars = paste0("Z_", 1:5)
candinter = c(candsplinevars, candbinaryvars)

# interim_n: vector with length equal to the number of interim analyses to be
# performed, denoting how many individuals are added to the sample at each one.
# Here, the total sample size is 500 with one interim analysis at 300 individuals.
interim_n = c(300,200)

baseline_pattern = function(x1,z1,z2,z3,z4,z5) {
  return(2*z1)
}
trt_pattern = function(x1,z1,z2,z3,z4,z5) {
  return(2*z1)
}

# Generate external dataset to assess accuracy of treatment recommendations
# at the end of each trial
data_test =  data.frame(X_1 = runif(n_test,0,1),
                        Z_1 = rbinom(n_test,1,0.35),
                        Z_2 = rbinom(n_test,1,0.5),
                        Z_3 = rbinom(n_test,1,0.65),
                        Z_4 = rbinom(n_test,1,0.2),
                        Z_5 = rbinom(n_test,1,0.35))

data_pool = data.frame(X_1 = runif(n_pool,0,1),
                       Z_1 = rbinom(n_pool,1,0.35),
                       Z_2 = rbinom(n_pool,1,0.5),
                       Z_3 = rbinom(n_pool,1,0.65),
                       Z_4 = rbinom(n_pool,1,0.2),
                       Z_5 = rbinom(n_pool,1,0.35),
                       trt = rbinom(n_pool,1,0.5))

trial_specs = list(alpha = 0.05,
                   B_1 = 0.95,
                   B_2 = 0.8,
                   e_1 = 0,
                   b_1 = 0,
                   b_2 = 0,
                   pi_var = 0.1,
                   enrich = TRUE
                   )
mcmc_specs = list(B = 1000,
                  burnin = 1000,
                  thin = 1,
                  chains = 2,
                  sigma_v = 0.1,
                  bma = TRUE)
prior_params = list(lambda_1 = 0.1,
                    lambda_2 = 1,
                    a_0 = 0.01,
                    b_0 = 0.01,
                    degree = 3,
                    k_max = 9,
                    w = 1,
                    sigma_B = sqrt(20))
true_tailoring_vars = c("X_1")

#' Check Function Arguments for Pattern Matching
#'
#' This internal function checks whether the number of arguments a function expects
#' matches the number of provided variables. If there is a mismatch, the function
#' will stop and return an error message. Otherwise, it returns `TRUE` indicating
#' that the check passed.
#'
#' @noRd
#' @param pattern A function whose arguments are to be checked (e.g., treatment pattern).
#' @param variables A vector of variables to check against the expected arguments of the pattern.
#'
#' @return A logical value: `TRUE` if the number of arguments matches, otherwise stops execution with an error message.
#'
#' @keywords internal
check_fn_args <- function(pattern, variables) {
  # Get the number of arguments that trt_pattern expects
  num_args_pattern <- length(formals(pattern))

  # Check if the number of arguments matches the length of candinter
  if (num_args_pattern != length(variables)) {
    stop(paste("Error: trt_pattern expects", num_args_pattern,
               "arguments, but", length(variables), "were provided."))
  } else {
    return(TRUE)  # If the check passes, return TRUE
  }
}
#' Run Simulations for Adaptive Trial Designs
#'
#' @description
#' This function runs simulations for adaptive clinical trial designs using Bayesian methods,
#' where participants are enrolled and analyzed based on predefined patterns. The simulations
#' allow for multiple iterations to evaluate the performance of the trial design, using
#' free-knot B-splines, Bayesian model averaging, and sequential interim analyses.
#'
#' @param sim_iters Integer specifying the number of simulations to run.
#' @param data_test A data frame representing the external testing dataset, which includes the
#' variables defined in `candsplinevars` and `candbinaryvars`, along with a true treatment outcome (`true_trt`).
#' @param data_pool A data frame representing the population data available for enrollment in the trial,
#' which includes variables for outcome (`Y`), treatment group (`trt`), and the candidate spline and binary variables.
#' @param baseline_pattern A function describing the baseline outcome pattern for participants based on covariates.
#' @param trt_pattern A function describing the treatment effect pattern for participants.
#' @param candsplinevars A vector of continuous candidate predictive variables (e.g., age, lab values).
#' @param candbinaryvars A vector of binary candidate predictive variables (e.g., gender, disease presence).
#' @param candinter A vector specifying which of `candsplinevars` and `candbinaryvars` are tailoring variables.
#' @param true_tailoring_vars A vector containing the true tailoring variables used in the simulation.
#' @param interim_n A vector specifying the number of participants to enroll at each interim analysis.
#' @param sigma_eps Numeric value representing the standard deviation of noise added to the outcome.
#' @param mcmc_specs A list of MCMC specifications, including:
#'                   - `B`: Number of posterior samples.
#'                   - `burnin`: Number of burn-in samples.
#'                   - `thin`: The thinning parameter.
#'                   - `chains`: The number of MCMC chains.
#'                   - `sigma_v`: Proposal variance for "jump" terms.
#'                   - `bma`: Boolean indicating whether to include Bayesian model averaging.
#' @param prior_params A list of prior distribution parameters, including:
#'                     - `lambda_1`: Prior for the number of terms in the model.
#'                     - `lambda_2`: Prior for the number of knots in each spline.
#'                     - `a_0`, `b_0`: Parameters for inverse gamma prior on individual-level variance.
#'                     - `degree`: Degree of B-splines.
#'                     - `k_max`: Maximum number of knots for each spline term.
#'                     - `w`: Window to propose knot location changes.
#'                     - `sigma_B`: Prior Normal variance for model coefficients.
#' @param trial_specs A list of trial specifications, including:
#'                    - `alpha`: Significance level for hypothesis testing.
#'                    - `B_1`, `B_2`: Threshold probabilities for stopping the trial for efficacy and futility.
#'                    - `e_1`: Minimum posterior quantile value for enrollment eligibility.
#'                    - `pi_var`: Minimum proportion of participants showing efficacy.
#'                    - `enrich`: Boolean indicating if adaptive enrichment is used.
#' @param parallel Boolean indicating whether to run the simulations in parallel.
#'
#' @return A data frame summarizing the simulated trial results, including:
#' \describe{
#'   \item{success}{Boolean indicating whether the trial was successful.}
#'   \item{final_stop_efficacy}{Boolean if the trial stopped for efficacy at the final interim analysis.}
#'   \item{final_stop_futility}{Boolean if the trial stopped for futility at the final interim analysis.}
#'   \item{prop_eff_final}{Proportion of participants in the effective subgroup at the final analysis.}
#'   \item{final_trial_size}{Total number of participants enrolled by the trial's conclusion.}
#'   \item{effect_and_subgroup}{Boolean indicating efficacy specific to a subgroup.}
#'   \item{effect_overall}{Boolean indicating overall efficacy across all participants.}
#'   \item{subgroup_spec}{Boolean if subgroup-specific analyses were performed.}
#'   \item{selected_markers}{String listing the tailoring variables selected by the model.}
#'   \item{correct_marker}{Boolean indicating if the selected tailoring variables match the true tailoring variables.}
#'   \item{include_correct_marker}{Boolean if the selected tailoring variables include the true tailoring variables.}
#'   \item{accuracy}{Proportion of correct treatment decisions based on the final trial results.}
#'   \item{mean_subgroup_ate}{Mean treatment effect in the effective subgroup.}
#'   \item{mean_subgroup_ate_gr_cutoff}{Proportion of subgroup with treatment effect greater than a specified cutoff.}
#'   \item{mean_subgroup_ate_l_cutoff}{Proportion of subgroup with treatment effect less than a specified cutoff.}
#'   \item{time}{Time taken for each iteration.}
#' }
#' @details
#' This function evaluates the performance of an adaptive clinical trial by running multiple
#' iterations of the trial. It uses a free-knot B-spline model to fit the data, with the option
#' of Bayesian model averaging. The trial is enriched over time based on interim analyses and
#' participants are selectively enrolled based on whether they are predicted to benefit from the
#' treatment.
#'
#' @examples
#' sim_iters <- 10
#' n_pool <- 5000
#' n_test <- 10000
#'
#' # Data setup
#' candsplinevars <- c("X_1", "X_2")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' data_test <- data.frame(X_1 = runif(n_test, 0, 1),
#'                         X_2 = runif(n_test, 0, 1),
#'                         Z_1 = rbinom(n_test, 1, 0.35),
#'                         Z_2 = rbinom(n_test, 1, 0.5),
#'                         Z_3 = rbinom(n_test, 1, 0.65),
#'                         Z_4 = rbinom(n_test, 1, 0.2),
#'                         Z_5 = rbinom(n_test, 1, 0.35))
#'
#' data_pool <- data.frame(X_1 = runif(n_pool, 0, 1),
#'                         X_2 = runif(n_pool, 0, 1),
#'                         Z_1 = rbinom(n_pool, 1, 0.35),
#'                         Z_2 = rbinom(n_pool, 1, 0.5),
#'                         Z_3 = rbinom(n_pool, 1, 0.65),
#'                         Z_4 = rbinom(n_pool, 1, 0.2),
#'                         Z_5 = rbinom(n_pool, 1, 0.35),
#'                         trt = rbinom(n_pool, 1, 0.5))
#'
#' # Generate outcomes
#' baseline_pattern <- function(X_1, X_2, Z_1, Z_2, Z_3, Z_4, Z_5) {
#'   2 * Z_1 + 2 * Z_1 * X_1
#' }
#' trt_pattern <- function(X_1, X_2, Z_1, Z_2, Z_3, Z_4, Z_5) {
#'   1 * Z_1 + 0.5 * X_1 * Z_1
#' }
#'
#' # MCMC and prior specifications
#' mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
#' prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3,
#'                      k_max = 9, w = 1, sigma_B = sqrt(20))
#'
#' # Trial specifications
#' trial_specs <- list(alpha = 0.05, B_1 = 0.95, B_2 = 0.8, e_1 = 0, b_1 = 0, b_2 = 0,
#'                     pi_var = 0.1, enrich = TRUE)
#'
#' # Run simulation
#' results <- runSimulations(sim_iters, data_test, data_pool, baseline_pattern, trt_pattern,
#'                           candsplinevars, candbinaryvars, candinter,
#'                           true_tailoring_vars = c("X_1"), interim_n = c(100, 200, 300),
#'                           sigma_eps = 0.1, mcmc_specs, prior_params, trial_specs, parallel = FALSE)
#'
#' @export
runSimulations <- function(sim_iters,
                           data_test,
                           data_pool,
                           baseline_pattern,
                           trt_pattern,
                           candsplinevars,
                           candbinaryvars,
                           candinter,
                           true_tailoring_vars,
                           interim_n,
                           sigma_eps,
                           mcmc_specs,
                           prior_params,
                           trial_specs,
                           parallel) {
  check_fn_args(trt_pattern,length(candinter))
  check_fn_args(baseline_pattern,length(candsplinevars) + length(candbinaryvars))
  data_test$trt = 1
  data_test$truth = pmap(data_test[,candinter], trt_pattern)
  data_test$true_trt = as.numeric(data_test$truth > 0)

  if (parallel) {
    ncores = Sys.getenv("SLURM_CPUS_PER_TASK")
    registerDoParallel(cores=ncores)
  }
  results = foreach(i=1:sim_iters, .combine='rbind.fill') %dopar% {
    data_pool$Y = pmap(data_pool[,c(candsplinevars,candbinaryvars)],baseline_pattern) +
      pmap(data_pool[,candinter],trt_pattern)*data_pool$trt + rnorm(length(data_pool),0,sigma_eps)

    start_time_fk_bma = Sys.time()
    results_fk_bma = runTrial(data_pool,
                              data_test,
                              candsplinevars,
                              candbinaryvars,
                              candinter,
                              true_tailoring_vars,
                              interim_n,
                              mcmc_specs,
                              prior_params,
                              trial_specs
    )
    end_time_fk_bma = Sys.time()
    results_fk_bma$time = difftime(end_time_fk_bma, start_time_fk_bma, units = "mins")
    colnames(results_fk_bma) = paste0("fk_bma",colnames(results_fk_bma))
  }

  for (name in names(trial_specs)) {
    results[[name]] <- trial_specs[[name]]
  }

  # Add elements from mcmc_specs to results
  for (name in names(mcmc_specs)) {
    results[[name]] <- mcmc_specs[[name]]
  }

  # Add elements from prior_params to results
  for (name in names(prior_params)) {
    results[[name]] <- prior_params[[name]]
  }
  results$n_test = n_test
  results$n_pool = n_pool
  results$interim_n = paste(interim_n,collapse=",")
  results$sim_iters = sim_iters
  results$candbinaryvars = paste(candbinaryvars,collapse=",")
  results$candsplinevars = paste(candsplinevars,collapse=",")
  return(results)
}
