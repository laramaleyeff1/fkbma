library(doParallel)
library(plyr)
library(purrr)

#' Helper function to execute a simulated Bayesian adaptive design using FK-BMA
#'
#' Conducts an adaptive clinical trial w/ interim analyses. This function enrolls participants
#' based on interim results, dynamically adjusts trial enrollment criteria using enrichment methods,
#' and can stop the trial early for efficacy or futility. It returns performance metrics such as accuracy
#' and subgroup treatment effects.
#'
#' @noRd
#' @param data_pool A data frame representing the population with the following columns:
#' \describe{
#'   \item{trt}{The group indicator (1 = experimental, 0 = control).}
#'   \item{Y}{The outcome variable.}
#'   \item{candsplinevars}{Columns named according to continuous candidate variables.}
#'   \item{candbinaryvars}{Columns named according to binary candidate variables.}
#' }
#' @param data_test A large, external testing data frame for assessing external validity. Must include:
#' \describe{
#'   \item{candsplinevars}{Continuous predictive variables (same as data_pool).}
#'   \item{candbinaryvars}{Binary predictive variables (same as data_pool).}
#'   \item{true_trt}{Indicator if treatment is truly effective (1 = yes, 0 = no).}
#' }
#' @param candsplinevars A vector of continuous candidate predictive variables (e.g., age, lab values).
#' @param candbinaryvars A vector of binary candidate predictive variables (e.g., gender, disease presence).
#' @param candinter A vector of candidate tailoring variables (subset of candsplinevars and candbinaryvars).
#' @param true_tailoring_vars A vector of true tailoring variables used for performance metric calculation.
#' @param interim_n A vector specifying the number of participants enrolled at each interim analysis.
#' @param mcmc_specs A list of MCMC specifications (optional), e.g., number of iterations, burn-in, thinning, etc.
#' @param prior_params A list of prior distribution parameters for the Bayesian model (optional).
#' @param trial_specs A list of trial specifications, e.g., thresholds for efficacy, futility, and enrichment (optional).
#'
#' @return A data frame summarizing the trial results, including:
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
#' }
#' @keywords internal
#'
#' @details The function adjusts trial enrollment dynamically based on interim analyses, enabling enrichment
#' by enrolling participants expected to benefit more from the treatment. At each interim analysis, the
#' trial can stop early based on efficacy or futility thresholds. The function tracks trial performance metrics,
#' such as the accuracy of treatment decisions and subgroup-specific results.
#'
#' @examples
#' candbinaryvars = paste0("Z_",1:5)
#' candsplinevars = c("X_1")
#' candinter = c(candsplinevars, candbinaryvars)
#'
#' n = 1000
#' data_pool = data.frame(X_1 = runif(n, 0, 1),
#'                        Z_1 = rbinom(n, 1, 0.35),
#'                        Z_2 = rbinom(n, 1, 0.5),
#'                        Z_3 = rbinom(n, 1, 0.65),
#'                        Z_4 = rbinom(n, 1, 0.2),
#'                        Z_5 = rbinom(n, 1, 0.35),
#'                        trt = rbinom(n, 1, 0.5))
#'
#' data_pool$Y = 2 * data_pool$Z_1 + 2 * data_pool$Z_1 * data_pool$trt + rnorm(n, 0, 0.1)
#'
#' data_test = data_pool  # For simplicity, use same structure as data_pool for testing.
#' data_test$true_trt = ifelse(data_test$Z_1 > 0, 1, 0)
#'
#' trial_specs = list(alpha = 0.05,
#'                    B_1 = 0.95,
#'                    B_2 = 0.8,
#'                    e_1 = 0,
#'                    b_1 = 0,
#'                    b_2 = 0,
#'                    pi_var = 0.1,
#'                    enrich = TRUE)
#'
#' mcmc_specs = list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
#'
#' prior_params = list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))
#'
#' trial_results = runTrial(data_pool, data_test, candsplinevars, candbinaryvars, candinter,
#'                          c("Z_1"), c(100, 200, 300), mcmc_specs, prior_params, trial_specs)
runTrial <- function(data_pool,
                     data_test,
                     candsplinevars,
                     candbinaryvars,
                     candinter,
                     true_tailoring_vars,
                     interim_n,
                     mcmc_specs = NULL,
                     prior_params = NULL,
                     trial_specs = NULL
) {
  subgroup_spec = FALSE
  for (interim in 1:length(interim_n)) {
    max_interim = interim
    print(paste("interim: ", interim))
    # If we are adaptively enriching trial population, after the first interim analysis
    # we restrict enrollment to individuals expected to benefit from treatment, i.e.
    # the effective subspace criterion from previous interim analysis
    if (interim > 1 & trial_specs$enrich) {
      elibible_results <- eligibleEnroll(data_pool,
                                         candsplinevars,
                                         candbinaryvars,
                                         candinter,
                                         interim_analysis_results$trial_results,
                                         trial_specs$alpha,
                                         trial_specs$e_1)
      data_indx = sample(elibible_results$eligible_idx,interim_n[interim])
      data = rbind(data, data_pool[data_indx,])
    }

    if (interim == 1 | !trial_specs$enrich) {
      data_indx = sample(nrow(data_pool), interim_n[interim])
      data = data_pool[data_indx,]
    }

    assign(paste0("true_prop_eff_",interim),mean(data$truth>0))

    data_pool = data_pool[-data_indx,]

    interim_analysis_results = interimAnalysis(data,
                                               candsplinevars,
                                               candbinaryvars,
                                               candinter,
                                               last = (interim == length(interim_n)),
                                               mcmc_specs,
                                               prior_params,
                                               trial_specs)

    if (!interim_analysis_results$success) {
      return(data.frame(success = FALSE))
    }

    assign(paste0("prop_eff_",interim),interim_analysis_results$prop_eff)
    assign(paste0("stop_efficacy_",interim),interim_analysis_results$stop_efficacy)
    assign(paste0("stop_futility_",interim),interim_analysis_results$stop_futility)

    # If the proportion effective is ever less than one,
    # we are performing a subgroup-specific analysis.

    # For example, if prop_eff = 0.5 in interim analysis 1,
    # then we restrict the entry criteria and prop_eff may be
    # equal to 1 in the second cohort. An overall analysis is
    # performed if prop_eff=1 for all interim analyses
    if (interim_analysis_results$prop_eff < 1) {
      subgroup_spec = TRUE
    }

    # If we have reached the last interim analysis, or the trial is being
    # stopped early for efficacy or futility: compute external accuracy and break
    # out of the outer for loop
    if (interim == length(interim_n) |
        interim_analysis_results$stop_efficacy |
        interim_analysis_results$stop_futility) {
      final_trial_size = sum(interim_n[1:interim])
      # Assess accuracy based on large, external dataset
      quantile_test = getPosteriorQuantiles(data_test,
                                            candsplinevars,
                                            candbinaryvars,
                                            candinter,
                                            interim_analysis_results$trial_results,
                                            trial_specs$alpha
      )
      # Here, "trt_decision" refers to the treatment decisions that would be made based on the observed
      # trial data
      data_test$trt_decision = as.numeric(quantile_test > trial_specs$e_1)
      accuracy = sum(as.numeric(data_test$trt_decision == data_test$true_trt))/nrow(data_test)
      break
    }
  }

  # correct_marker tells us if the selecting tailoring variables are identical to the true
  # tailoring variables
  correct_marker = identical(true_tailoring_vars,interim_analysis_results$included_vars)
  # identical() doesn't handle length 0 vectors, so we adjust for this
  if (length(true_tailoring_vars) == 0) {
    correct_marker = (length(interim_analysis_results$included_vars)==0)
  }

  # check if the selected tailoring variables contain the true tailoring variables; i.e.
  # = TRUE even if we have selected extra
  include_correct_marker = all(true_tailoring_vars %in% interim_analysis_results$included_vars)

  returned = data.frame(success = TRUE,
                        final_stop_efficacy = interim_analysis_results$stop_efficacy,
                        final_stop_futility = interim_analysis_results$stop_futility,
                        prop_eff_final = interim_analysis_results$prop_eff,
                        final_trial_size = final_trial_size,
                        effect_and_subgroup = (interim_analysis_results$stop_efficacy & subgroup_spec),
                        effect_overall = (interim_analysis_results$stop_efficacy & !subgroup_spec),
                        subgroup_spec = subgroup_spec,
                        selected_markers = paste(interim_analysis_results$included_vars,collapse=","),
                        correct_marker = correct_marker,
                        include_correct_marker = include_correct_marker,
                        accuracy = accuracy
  )

  for (i in 1:max_interim) {
    returned[[paste0("stop_efficacy_",i)]] = eval(parse(text=paste0("stop_efficacy_",i)))
    returned[[paste0("stop_futility_",i)]] = eval(parse(text=paste0("stop_futility_",i)))
    returned[[paste0("prop_eff_",i)]] = eval(parse(text=paste0("prop_eff_",i)))
    returned[[paste0("true_prop_eff_",i)]] = eval(parse(text=paste0("true_prop_eff_",i)))
  }

  return(returned)
}

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
    stop(paste("Error: pattern functions should take", length(variables),
               "arguments, but", num_args_pattern, "were provided."))
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
#' The input variables must be in the order c(candsplinevars, candbinaryvars).
#' @param trt_pattern A function describing the treatment effect pattern for participants.
#' The input variables must be in the order candinter.
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
  check_fn_args(trt_pattern,candinter)
  check_fn_args(baseline_pattern,c(candsplinevars,candbinaryvars))
  default_trial_specs <- list(alpha = 0.05,
                              B_1 = 0.95,
                              B_2 = 0.8,
                              e_1 = 0,
                              b_1 = 0,
                              b_2 = 0,
                              pi_var = 0.1,
                              enrich = TRUE)

  # If trial_specs is NULL, use default_trial_specs
  if (is.null(default_trial_specs)) {
    trial_specs <- default_trial_specs
  }

  # Check if each parameter is missing, if so, set it to the default value
  for (param in names(default_trial_specs)) {
    if (is.null(trial_specs[[param]])) {
      trial_specs[[param]] <- default_trial_specs[[param]]
      print(paste("Setting", param, "to default of", default_trial_specs[[param]]))
    }
  }
  data_test$trt = 1
  data_test$truth = unlist(pmap(data_test[,candinter], trt_pattern))
  data_test$true_trt = as.numeric(data_test$truth > 0)

  if (parallel) {
    ncores = Sys.getenv("SLURM_CPUS_PER_TASK")
    registerDoParallel(cores=ncores)
  }
  results = foreach(i=1:sim_iters, .combine='rbind.fill') %dopar% {
    data_pool$Y = unlist(pmap(data_pool[,c(candsplinevars,candbinaryvars)],baseline_pattern)) +
      unlist(pmap(data_pool[,candinter],trt_pattern))*data_pool$trt + rnorm(nrow(data_pool),0,sigma_eps)

    data_pool$truth = unlist(pmap(data_pool[,candinter], trt_pattern))

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
    results_fk_bma
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
  results$n_test = nrow(data_test)
  results$n_pool = nrow(data_pool)
  results$interim_n = paste(interim_n,collapse=",")
  results$sim_iters = sim_iters
  results$candbinaryvars = paste(candbinaryvars,collapse=",")
  results$candsplinevars = paste(candsplinevars,collapse=",")
  return(results)
}
