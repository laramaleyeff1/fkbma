# Function:       getPosteriorQuantiles
# Author:         Lara Maleyeff
# Description:    Find the alpha-quantile of the posterior treatment effect distribution for
#                 each individual (row). It uses the model that was fit in the most recent interim analysis (original data) to
#                 find the alpha-quantiles of external data (used for either the adaptive enrichment steps or
#                 accuracy calculations)
# Parameters:     data              A data frame with one individual per row and information on candsplineinter and candbinaryinter
#                 candsplinevars    Vector with names of continuous variables (not used)
#                 candbinaryvars    Vector with names of binary variables (not used)
#                 trial_results     Results from most recent interim analysis, a data frame containing:
#                                   - included_vars: selected tailoring variables
#                                   - candsplineinter: candidate spline tailoring variables
#                                   - candbinaryinter: candidate binary tailoring variables
#                                   - trt_param: posterior distribution of treatment effect
#                                   - data_fit: original data used to fit the model, used for interpolation
#                                   - splines_fitted: posterior distribution of fitted values for each spline
#                                   - binary_param: posterior distribution of binary coefficients
#                 alpha             Cutoff for the effective subspace
#
# Returns:        alpha-row quantile for each individuals to then be compared with e_1. If
#                 the alpha-quantile of the treatment effect is > e_1 then that individual's
#                 variable combination is in the effective subspace
getPosteriorQuantiles <- function(data,
                                  candsplinevars,
                                  candbinaryvars,
                                  candinter,
                                  trial_results,
                                  alpha
) {

  candsplineinter = intersect(candsplinevars,candinter)
  candbinaryinter = intersect(candbinaryvars,candinter)

  combined_posterior = matrix(1,nrow=nrow(data)) %*% trial_results$inter_trt_param[,"trt"]

  if (length(candsplineinter) > 0) {
    for (m in 1:length(candsplineinter)) {
      interpolated_splines_fitted <- sapply(1:nrow(trial_results$splines_fitted[[candsplineinter[m]]]), function(j) {
        approx(x = trial_results$data_fit[,candsplineinter[m]],
               y = trial_results$splines_fitted[[candsplineinter[m]]][j,], xout = data[,candsplineinter[m]],rule=2)$y
      })
      combined_posterior = combined_posterior + interpolated_splines_fitted
    }
  }
  if (length(candbinaryinter) > 0) {
    combined_posterior = combined_posterior +
      as.matrix(data[,candbinaryinter]) %*% t(trial_results$binary_param[,paste0(candbinaryinter,"_inter")])
  }

  return(rowQuantiles(combined_posterior,probs = c(alpha)))
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
