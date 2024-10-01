######################################################################################################
#           Code to run a Bayesian adaptive design with interim decision rules described in
#                                 Section 2.4 of Maleyeff et al. (2024).
#                                   Contact: laramaleyeff@gmail.com
#                                       Last updated: October 2024
######################################################################################################

library(matrixStats)
library(splines)
library(dplyr)

#' Performs Interim Analysis for an Adaptive Clinical Trial
#'
#' This function performs interim analysis for an adaptive clinical trial using MCMC results.
#' It determines whether the trial should stop for efficacy or futility based on pre-specified
#' thresholds and trial results generated from the rjMCMC procedure.
#'
#' @param data A data frame containing the trial observations, including the outcome, treatment,
#' and candidate spline and binary variables.
#' @param candsplinevars A vector of continuous candidate predictive variables (e.g., age, lab values).
#' @param candbinaryvars A vector of binary candidate predictive variables (e.g., gender, disease presence).
#' @param candinter A vector of candidate tailoring variables, which is a subset of \code{candsplinevars} and \code{candbinaryvars}.
#' @param last Boolean indicating whether this is the last interim analysis. If \code{TRUE}, futility is not considered.
#' @param mcmc_specs List of MCMC specifications (optional), containing:
#' \itemize{
#'   \item \code{B}: Number of posterior samples (default: 2000).
#'   \item \code{burnin}: Number of burn-in samples (default: 10000).
#'   \item \code{thin}: Thinning parameter (default: 5).
#'   \item \code{chains}: Number of chains (default: 1).
#'   \item \code{sigma_v}: Proposal variance for "jump" terms (default: 0.1).
#'   \item \code{bma}: Boolean indicating whether to include Bayesian model averaging (default: TRUE).
#' }
#' @param prior_params List of prior distribution parameters (optional), including:
#' \itemize{
#'   \item \code{lambda_1}: Prior parameter for the number of terms in the model (default: 0.1).
#'   \item \code{lambda_2}: Prior parameter for the number of knots in each spline (default: 1).
#'   \item \code{a_0}: Shape parameter for inverse gamma prior on variance (default: 0.01).
#'   \item \code{b_0}: Rate parameter for inverse gamma prior on variance (default: 0.01).
#'   \item \code{degree}: Degree of B-splines (default: 3).
#'   \item \code{k_max}: Maximum number of knots for each spline term (default: 9).
#'   \item \code{w}: Window to propose knot location changes (default: 1).
#'   \item \code{sigma_B}: Prior variance for model coefficients (default: \code{sqrt(20)}).
#' }
#' @param trial_specs List of trial specifications (optional), including thresholds for efficacy, futility,
#' and the proportion of effective participants (default values are used if not provided). It includes:
#' \itemize{
#'   \item \code{alpha}: Significance level for hypothesis testing (default: 0.05).
#'   \item \code{B_1}: Threshold probability for stopping the trial for efficacy (default: 0.95).
#'   \item \code{B_2}: Threshold probability for stopping the trial for futility (default: 0.8).
#'   \item \code{e_1}: Minimum posterior quantile value required for enrollment eligibility (default: 0).
#'   \item \code{b_1}: Minimum average treatment effect to stop the trial for efficacy (default: 0).
#'   \item \code{b_2}: Maximum average treatment effect under which the trial can stop for futility (default: 0).
#'   \item \code{pi_var}: Minimum proportion of participants required to demonstrate efficacy for the trial to continue (default: 0.1).
#'   \item \code{enrich}: Boolean indicating whether the trial is enriched for more responsive participants (default: TRUE).
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{success}: Boolean indicating whether the MCMC procedure was successful.
#'   \item \code{included_vars}: Vector of tailoring variables included in the effective subgroup.
#'   \item \code{stop_efficacy}: Boolean indicating whether the trial should stop for efficacy.
#'   \item \code{stop_futility}: Boolean indicating whether the trial should stop for futility.
#'   \item \code{prop_eff}: Proportion of participants in the effective subgroup.
#'   \item \code{trial_results}: Results from the MCMC procedure if successful.
#' }
#'
#' @details This function runs the MCMC procedure using \code{rjMCMC} until convergence to generate
#' posterior distributions of treatment effects. The posterior quantiles are compared to a predefined
#' threshold, and participants whose quantiles exceed this threshold form the "effective subgroup."
#' Based on predefined thresholds for efficacy and futility, the function decides whether to stop the
#' trial early. If it is the last interim analysis, futility is not considered.
#' @examples
#' candbinaryvars = paste0("Z_",1:5)
#' candsplinevars = c("X_1")
#' candinter = c(candsplinevars, candbinaryvars)
#'
#' n = 1000
#' data = data.frame(X_1 = runif(n,0,1),
#'                    Z_1 = rbinom(n,1,0.35),
#'                    Z_2 = rbinom(n,1,0.5),
#'                    Z_3 = rbinom(n,1,0.65),
#'                    Z_4 = rbinom(n,1,0.2),
#'                    Z_5 = rbinom(n,1,0.35),
#'                    trt = rbinom(n,1,0.5))
#'
#' data$Y = 2 * data$Z_1 + 2 * data$Z_1 * data$trt + rnorm(n, 0, 0.1)
#'
#' trial_specs = list(alpha = 0.05,
#'                    B_1 = 0.95,
#'                    B_2 = 0.8,
#'                    e_1 = 0,
#'                    b_1 = 0,
#'                    b_2 = 0,
#'                    pi_var = 0.1,
#'                    enrich = TRUE)
#' mcmc_specs = list(B = 1000,
#'                   burnin = 1000,
#'                   thin = 1,
#'                   chains = 2,
#'                   sigma_v = 0.1,
#'                   bma = TRUE)
#' prior_params = list(lambda_1 = 0.1,
#'                     lambda_2 = 1,
#'                     a_0 = 0.01,
#'                     b_0 = 0.01,
#'                     degree = 3,
#'                     k_max = 9,
#'                     w = 1,
#'                     sigma_B = sqrt(20))
#'
#' trial_results = rjMCMC(data, candsplinevars, candbinaryvars, candinter, mcmc_specs, prior_params)
#'
#' interim_results = interimAnalysis(data, candsplinevars, candbinaryvars, candinter, last = FALSE,
#'                                   mcmc_specs, prior_params, trial_specs)
#'
#' @export
interimAnalysis <- function(data,
                            candsplinevars,
                            candbinaryvars,
                            candinter,
                            last,
                            mcmc_specs = NULL,
                            prior_params = NULL,
                            trial_specs = NULL) {

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

  for (i in 1:5) {
    trial_results = rjMCMC(data,
                            candsplinevars,
                            candbinaryvars,
                            candinter,
                            mcmc_specs,
                            prior_params)
    if (trial_results$success) {
      print("rjMCMC procedure was successful!")
      break
    }
  }

  if (trial_results$success) {
    included_vars <- sub("_[^_]+$", "", grep("inter",
                        names(which(trial_results$vars_prop_summ > trial_specs$pi_var)), value = TRUE))

    trt_eff_posterior = trial_results$trt_eff_posterior
    quantiles = rowQuantiles(trt_eff_posterior,probs = c(trial_specs$alpha))

    data_subgroup = data[which(quantiles>trial_specs$e_1),]
    trt_eff_posterior_subgroup = trt_eff_posterior[which(quantiles>trial_specs$e_1),]

    trt_eff_per_person = rowMeans(trt_eff_posterior)

    prop_eff = nrow(data_subgroup)/nrow(data)
    print(paste("Prevalence of effective subspace:", prop_eff))
    print(paste("Included vars: ", paste(included_vars,collapse=",")))

    stop_efficacy = FALSE
    stop_futility = FALSE

    if (prop_eff > trial_specs$pi_var) {
      # Compute the average treatment effect for the sensitive subgroup
      combined_subgroup_ate = colMeans(trt_eff_posterior_subgroup)

      if (mean(combined_subgroup_ate > trial_specs$b_1) > trial_specs$B_1) {
        stop_efficacy = TRUE
      }

      if (((mean(combined_subgroup_ate < trial_specs$b_2) > trial_specs$B_2)) & !last & !stop_efficacy) {
        stop_futility = TRUE
      }
    } else {
      stop_futility = TRUE
    }

    if (stop_futility) {
      print("Stop trial for futility")
    }

    if (stop_efficacy) {
      print("Stop trial for efficacy")
    }

    return(list(success = TRUE,
                included_vars = included_vars,
                stop_efficacy = stop_efficacy,
                stop_futility = stop_futility,
                prop_eff = prop_eff,
                trial_results = trial_results
               )
    )
  } else {
    return(list(success=FALSE))
  }


}

#' Determines Eligibility of New Participants for Enrollment in a Trial
#'
#' This function determines the eligibility of new participants for enrollment in a clinical trial
#' based on their posterior quantiles relative to a threshold. Participants whose posterior quantiles
#' exceed the specified enrollment threshold are considered eligible.
#'
#' @param data_new A data frame containing the new participant data, including variables specified by
#' \code{candsplinevars} and \code{candbinaryvars}.
#' @param candsplinevars A vector of names for continuous candidate predictive variables (e.g., age, lab values).
#' @param candbinaryvars A vector of names for binary candidate predictive variables (e.g., gender, disease presence).
#' @param candinter A vector of candidate tailoring variables, typically a subset of \code{candsplinevars} and \code{candbinaryvars}.
#' @param trial_results A list containing trial results, typically generated from MCMC, including the posterior distributions
#' for treatment effects.
#' @param alpha The significance or confidence level for computing posterior quantiles.
#' @param e_1 The enrollment threshold. Participants whose posterior quantiles exceed this threshold will be considered
#' eligible for enrollment.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{eligible_data_new}: A data frame containing the data of participants eligible for enrollment.
#'   \item \code{eligible_idx}: A vector of indices representing eligible participants in \code{data_new}.
#' }
#'
#' @details This function calculates the posterior quantiles of new participants using their data and the trial results.
#' These quantiles are compared to the enrollment threshold (\code{e_1}). Participants whose quantiles exceed \code{e_1}
#' are considered eligible for enrollment in the trial.
#'
#' @examples
#' # Example usage similar to interim analysis example
#' candbinaryvars = paste0("Z_",1:5)
#' candsplinevars = c("X_1")
#' candinter = c(candsplinevars, candbinaryvars)
#'
#' n = 1000
#' data = data.frame(X_1 = runif(n,0,1),
#'                    Z_1 = rbinom(n,1,0.35),
#'                    Z_2 = rbinom(n,1,0.5),
#'                    Z_3 = rbinom(n,1,0.65),
#'                    Z_4 = rbinom(n,1,0.2),
#'                    Z_5 = rbinom(n,1,0.35),
#'                    trt = rbinom(n,1,0.5))
#'
#' data$Y = 2 * data$Z_1 + 2 * data$Z_1 * data$trt + rnorm(n, 0, 0.1)
#'
#' n = 500
#' data_new = data.frame(X_1 = runif(n,0,1),
#'                       Z_1 = rbinom(n,1,0.35),
#'                       Z_2 = rbinom(n,1,0.5),
#'                       Z_3 = rbinom(n,1,0.65),
#'                       Z_4 = rbinom(n,1,0.2),
#'                       Z_5 = rbinom(n,1,0.35),
#'                       trt = rbinom(n,1,0.5))
#'
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
#' mcmc_specs = list(B = 1000,
#'                   burnin = 1000,
#'                   thin = 1,
#'                   chains = 2,
#'                   sigma_v = 0.1,
#'                   bma = TRUE)
#'
#' prior_params = list(lambda_1 = 0.1,
#'                     lambda_2 = 1,
#'                     a_0 = 0.01,
#'                     b_0 = 0.01,
#'                     degree = 3,
#'                     k_max = 9,
#'                     w = 1,
#'                     sigma_B = sqrt(20))
#'
#' trial_results = rjMCMC(data, candsplinevars, candbinaryvars, candinter, mcmc_specs, prior_params)
#'
#' eligible_results = eligibleEnroll(data_new, candsplinevars, candbinaryvars, candinter, trial_results, alpha = 0.05, e_1 = 0)
#'
#' @export
eligibleEnroll <- function(data_new,
                           candsplinevars,
                           candbinaryvars,
                           candinter,
                           trial_results,
                           alpha,
                           e_1) {
  quantile_interim <- getPosteriorQuantiles(data_new,
                                            candsplinevars,
                                            candbinaryvars,
                                            candinter,
                                            trial_results,
                                            alpha)
  eligible_idx = which(quantile_interim>e_1)
  eligible_data_new = data_new[eligible_idx,]
  return(list(eligible_data_new = eligible_data_new,
         eligible_idx = eligible_idx))
}

#' Conduct an Adaptive Clinical Trial with Interim Analyses
#'
#' Conducts an adaptive clinical trial using interim analyses. This function enrolls participants
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
#'   \item{mean_subgroup_ate}{Mean treatment effect in the effective subgroup.}
#'   \item{mean_subgroup_ate_gr_cutoff}{Proportion of subgroup with treatment effect greater than a specified cutoff.}
#'   \item{mean_subgroup_ate_l_cutoff}{Proportion of subgroup with treatment effect less than a specified cutoff.}
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

    if (interim == 1 | !enrich) {
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
                        accuracy = accuracy,
                        mean_subgroup_ate = interim_analysis_results$mean_subgroup_ate,
                        mean_subgroup_ate_gr_cutoff = interim_analysis_results$mean_subgroup_ate_gr_cutoff,
                        mean_subgroup_ate_l_cutoff = interim_analysis_results$mean_subgroup_ate_l_cutoff
  )

  for (i in 1:max_interim) {
    returned[[paste0("stop_efficacy_",i)]] = eval(parse(text=paste0("stop_efficacy_",i)))
    returned[[paste0("stop_futility_",i)]] = eval(parse(text=paste0("stop_futility_",i)))
    returned[[paste0("prop_eff_",i)]] = eval(parse(text=paste0("prop_eff_",i)))
    returned[[paste0("true_prop_eff_",i)]] = eval(parse(text=paste0("true_prop_eff_",i)))
  }

  return(returned)
}
