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
#' @param mcmc_results (Optional) Precomputed MCMC results to use instead of running the MCMC process within the function. Default is `NULL`.
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
                            trial_specs = NULL,
                            mcmc_results = NULL
                            ) {

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

  if (is.null(mcmc_results)) {
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
  } else {
    trial_results = mcmc_results
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
