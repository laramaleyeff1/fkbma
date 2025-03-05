#' Summarize Results from Reversible Jump MCMC (rjMCMC)
#'
#' This function provides a detailed summary of the results from the `rjMCMC` procedure, including
#' model information, parameter estimates, posterior inclusion probabilities, convergence diagnostics,
#' and plots for spline terms. The function also prints the model formula with `fbs()` notation for spline terms,
#' indicating the use of free-knot B-splines.
#'
#' @param object An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#' @param digits Number of digits in summary output (default = 3)
#' @param level Credible interval level (default = 0.95)
#' @param pip_cutoff Posterior inclusion probability cutoff for reporting effective sample size
#' and R-squared (default = 0.10)
#' @param ... Additional arguments to be passed to other methods or functions.
#'
#' @return Prints the following summary information:
#' \describe{
#'   \item{Model Formula}{The model formula with spline terms wrapped in `fbs()`, indicating free-knot B-splines, and interaction terms appropriately formatted.}
#'   \item{Convergence Diagnostics}{Reports any convergence issues based on Geweke diagnostics.}
#'   \item{MCMC Sampler Arguments}{Displays MCMC sampler arguments, including the number of posterior samples, burn-in, thinning, and chains.}
#'   \item{Parameter Estimates}{Posterior mean, standard error, 95% credible intervals, effective sample size (ESS), Gelman-Rubin statistic (Rhat), and posterior inclusion probabilities (PIP) for binary parameters, exposure intercept, and exposure effect.}
#'   \item{Gaussian Family Parameters}{Posterior summary for the residual standard error (sigma).}
#'   \item{Posterior Inclusion Probabilities for Splines}{Prints the posterior inclusion probabilities for spline terms.}
#'   \item{Plots for Fitted Exposure Effects}{Plots the mean and 95% credible intervals for each spline term vs fitted exposure effects.}
#' }
#'
#' @details The function produces detailed summaries similar to those from `brms`, including
#' diagnostics, estimates, posterior inclusion probabilities, and spline effects. The spline terms
#' are wrapped in `fbs()` notation, indicating the use of free-knot B-splines in the model. If the sampler
#' did not converge, a warning is issued. The function also allows the user to view diagnostic plots for fitted
#' exposure effects.
#'
#' @importFrom stats quantile sd
#' @importFrom rstan Rhat
#' @importFrom coda effectiveSize mcmc mcmc.list
#'
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' summary(results)
#' }
#'
#' @export
summary.rjMCMC <- function(object, digits = 3, level = 0.95, pip_cutoff = 0.1,...) {
  results <- object
  factor_var = results$factor_var
  outcome = results$outcome
  # Print the formula and data information, with fbs() around spline variables
  cat("Model Information:\n")
  formula_parts <- c()

  # Handle spline variables, including those in candinter
  if (length(results$candsplinevars) > 0) {
    spline_terms <- paste0("fbs(", results$candsplinevars, ")")
    formula_parts <- c(formula_parts, spline_terms)
  }

  # Add binary variables
  if (length(results$candbinaryvars) > 0) {
    formula_parts <- c(formula_parts, results$candbinaryvars)
  }

  # Handle interaction terms (some of which may be splines)
  if (length(results$candinter) > 0) {
    interaction_terms <- sapply(results$candinter, function(var) {
      if (var %in% results$candsplinevars) {
        return(paste0("fbs(", var, "):",factor_var))
      } else {
        return(paste0(var, ":" ,factor_var))
      }
    })
    formula_parts <- c(formula_parts, factor_var, interaction_terms)
  } else {
    formula_parts <- c(formula_parts, factor_var)
  }

  # Build the formula string
  formula_string <- paste(outcome,"~", paste(formula_parts, collapse = " + "))
  cat("Formula: ", formula_string, "\n")
  cat("Note: fbs() indicates a free-knot B-spline.\n")
  cat("Data: ", deparse(substitute(object$data_fit)), "\n")
  cat("Number of observations: ", nrow(results$data_fit), "\n")

  # Extract sampler arguments
  mcmc_specs <- results$mcmc_specs
  cat("MCMC Sampler Arguments:\n")
  cat("  - iter: ", mcmc_specs$iter, "\n")
  cat("  - warmup: ", mcmc_specs$warmup, "\n")
  cat("  - thin: ", mcmc_specs$thin, "\n")
  cat("  - chains: ", mcmc_specs$chains, "\n")

  B = ceiling((mcmc_specs$iter - mcmc_specs$warmup) / mcmc_specs$thin)*mcmc_specs$chains
  # Create a summary for binary_param, fixed_param, and sigma_sq
  cat("\nParameter Estimates:\n")

  # Convert matrices to data.frame to ensure column names are preserved
  combined_params <- data.frame(intercept = results$fixed_param[, "intercept"],
                                factor_var = results$fixed_param[, factor_var],
                                results$binary_param, check.names = FALSE)
  names(combined_params)[2] = factor_var

  # Posterior inclusion probabilities for binary and spline parameters
  vars_prop_summ <- results$vars_prop_summ

  # Separate the posterior inclusion probabilities for intercept and trt
  intercept_pip <- 1
  trt_pip <- 1
  alpha = 1-level
  # Extract PIPs for other parameters and match them to the corresponding parameter names
  other_pips <- vars_prop_summ[names(vars_prop_summ) %in% colnames(combined_params)]
  other_pips <- other_pips[order(-other_pips)]  # Sort parameters by posterior inclusion probabilities (descending)

  # Rearrange the combined_params so that intercept and trt are first, followed by other parameters in order of PIPs
  param_order <- c("intercept", factor_var, names(other_pips))
  combined_params <- combined_params[, param_order]

  # Add PIP column to the parameter summary table
  all_pips <- c(intercept = intercept_pip, trt = trt_pip, other_pips)
  summary_table <- function(param_mat, pips) {
    mean_val <- apply(param_mat, 2, mean)
    sd_val <- apply(param_mat, 2, stats::sd)
    lower_ci <- apply(param_mat, 2, function(x) stats::quantile(x, alpha/2))
    upper_ci <- apply(param_mat, 2, function(x) stats::quantile(x, 1-alpha/2))

    rhat <- apply(param_mat, 2, function(x) {
      chain_matrix <- matrix(x, nrow = B, ncol = mcmc_specs$chains, byrow = FALSE)
      rstan::Rhat(chain_matrix)
    })


    ess <- apply(param_mat, 2, function(x) {
      chain_matrix <- matrix(x, nrow = B, ncol = mcmc_specs$chains, byrow = FALSE)
      chain_list <- split(chain_matrix, col(chain_matrix))
      # Convert each chain into an 'mcmc' object
      mcmc_chains <- lapply(chain_list, function(chain) coda::mcmc(matrix(chain, ncol = 1)))
      # Combine the individual 'mcmc' objects into an 'mcmc.list'
      mcmc_list <- coda::mcmc.list(mcmc_chains)
      coda::effectiveSize(mcmc_list)
    })

    df = data.frame(
      Estimate = round(mean_val,digits),
      `Est.Error` = round(sd_val,digits),
      `l-95% CI` = round(lower_ci,digits),
      `u-95% CI` = round(upper_ci, digits),
      Eff.Sample = round(ess, digits) ,
      Rhat = round(rhat,digits),
      PIP = round(pips,digits),
      check.names = FALSE
    )
    df$Rhat[df$PIP < pip_cutoff] <- NA
    df$Eff.Sample[df$PIP < pip_cutoff] <- NA

    colnames(df)[3:4] <- c(paste0("l-",level*100,"% CI"),
                           paste0("u-",level*100,"% CI"))
    df
  }


  cat("\nNon-spline Parameters:\n")
  combined_summary <- summary_table(combined_params, pips = all_pips)
  print(combined_summary)
  cat("\nPIP = posterior inclusion probability\n")


  cat("\nGaussian Family Specific Parameters:\n")
  sigma_sq_summary <- summary_table(sqrt(results$sigma_sq), pips = 1)
  rownames(sigma_sq_summary) <- c("sigma")
  print(sigma_sq_summary)

  # Separate Posterior Inclusion Probabilities for `candsplinevars`
  cat("\nPosterior Inclusion Probabilities for Splines:\n")
  spline_pips <- vars_prop_summ[names(vars_prop_summ) %in% c(results$candsplinevars,
                                                             paste0(results$candsplinevars,paste0(":",factor_var)))]
  print(spline_pips)
}

#' Print a summary of results from from Reversible Jump MCMC (rjMCMC)
#'
#' This function provides a detailed summary of the results from the `rjMCMC` procedure, including
#' model information, parameter estimates, posterior inclusion probabilities, convergence diagnostics,
#' and plots for spline terms. The function also prints the model formula with `fbs()` notation for spline terms,
#' indicating the use of free-knot B-splines.
#'
#' @param x An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#' @param ... Additional arguments to be passed to other methods or functions.
#'
#' @return Prints the following summary information:
#' \describe{
#'   \item{Model Formula}{The model formula with spline terms wrapped in `fbs()`, indicating free-knot B-splines, and interaction terms appropriately formatted.}
#'   \item{Convergence Diagnostics}{Reports any convergence issues based on Geweke diagnostics.}
#'   \item{MCMC Sampler Arguments}{Displays MCMC sampler arguments, including the number of posterior samples, burn-in, thinning, and chains.}
#'   \item{Parameter Estimates}{Posterior mean, standard error, 95% credible intervals, effective sample size (ESS), Gelman-Rubin statistic (Rhat), and posterior inclusion probabilities (PIP) for binary parameters, exposure intercept, and exposure effect.}
#'   \item{Gaussian Family Parameters}{Posterior summary for the residual standard error (sigma).}
#'   \item{Posterior Inclusion Probabilities for Splines}{Prints the posterior inclusion probabilities for spline terms.}
#'   \item{Plots for Fitted Exposure Effects}{Plots the mean and 95% credible intervals for each spline term vs fitted exposure effects.}
#' }
#'
#' @details The function produces detailed summaries similar to those from `brms`, including
#' diagnostics, estimates, posterior inclusion probabilities, and spline effects. The spline terms
#' are wrapped in `fbs()` notation, indicating the use of free-knot B-splines in the model. If the sampler
#' did not converge, a warning is issued. The function also allows the user to view diagnostic plots for fitted
#' exposure effects.
#'
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'            outcome = "Y", factor_var = "trt")
#' print(results)
#' }
#' @export
print.rjMCMC <- function(x,...) {
  results <- x
  summary(results)
}

#' Predict from Reversible Jump MCMC (rjMCMC) Model Results
#'
#' This function generates posterior predictions from an rjMCMC model based on the provided data.
#' It combines the fixed effects, spline terms, and binary parameters, and includes residual variance
#' in the predictions.
#'
#' @param object An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#' @param newdata A data frame for which predictions are to be made. If NA, the original fitted data is used.
#' @param ... Additional arguments to be passed to other methods or functions.
#'
#' @return A matrix of predicted values.
#'
#' @importFrom stats rnorm approx
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'            outcome = "Y", factor_var = "trt")
#' predict(results)
#' }
#' @export
predict.rjMCMC <- function(object, newdata = NULL, ...) {
  results <- object
  factor_var = results$factor_var
  if (is.null(newdata)) {
    newdata = results$data_fit
  }
  fitted_posterior =  results$fixed_param %*% t(cbind(1,newdata[,factor_var]))

  all_splines = names(results$splines_fitted)
  all_binary = colnames(results$binary_param)

  # Check if newdata contains required columns
  required_columns <- c(results$candbinaryvars, results$candsplinevars, factor_var)
  missing_columns <- setdiff(required_columns, colnames(newdata))

  if (length(missing_columns) > 0) {
    stop(paste("Error: The following columns are missing in newdata:", paste(missing_columns, collapse = ", ")))
  }

  if (length(all_splines) > 0) {
    for (m in 1:length(all_splines)) {
      if (grepl(paste0(":",factor_var), all_splines[m])) {
        interpolated_splines_fitted <- sapply(1:nrow(results$splines_fitted[[all_splines[m]]]), function(j) {
          stats::approx(x = results$data_fit[,gsub(paste0(":",factor_var,"$"), "", all_splines[m])]*results$data_fit[,factor_var],
                 y = results$splines_fitted[[all_splines[m]]][j,]*results$data_fit[,factor_var], xout = newdata[,gsub(paste0(":",factor_var,"$"), "", all_splines[m])]*newdata[,factor_var],rule=2,ties="mean")$y
        })
      } else {
        interpolated_splines_fitted <- sapply(1:nrow(results$splines_fitted[[all_splines[m]]]), function(j) {
          stats::approx(x = results$data_fit[,all_splines[m]],
                 y = results$splines_fitted[[all_splines[m]]][j,], xout = newdata[,all_splines[m]],rule=2,ties="mean")$y
        })
      }
      fitted_posterior = fitted_posterior + t(interpolated_splines_fitted)
    }
  }

  interaction_terms <- colnames(results$binary_param)[grepl(paste0(":",factor_var,"$"), colnames(results$binary_param))]

  # Create interaction columns in newdata based on the matching regular terms
  for (interaction_term in interaction_terms) {
    # Get the base name before paste0(":",factor_var)
    base_name <- gsub(paste0(":",factor_var,"$"), "", interaction_term)

    # If the base term exists in newdata, create the interaction term
    if (base_name %in% colnames(newdata)) {
      newdata[[interaction_term]] <- newdata[[base_name]] * newdata[,factor_var]
    }
  }
  if (length(all_binary) > 0) {
    fitted_posterior = fitted_posterior +
      results$binary_param %*% t(as.matrix(newdata[,all_binary]))
  }

  residual_se <- sqrt(results$sigma_sq)
  predictive_posterior <- fitted_posterior + matrix(stats::rnorm(length(fitted_posterior), mean = 0, sd = residual_se), nrow = nrow(fitted_posterior), ncol = ncol(fitted_posterior))

  return(predictive_posterior)
}

#' Fitted values from Reversible Jump MCMC (rjMCMC) Model Results
#'
#' This function generates posterior fitted values from an rjMCMC model based on the provided data.
#' It combines the fixed effects, spline terms, and binary parameters.
#'
#' @param object An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#' @param newdata A data frame for which fitted values are to be computed. If NA, the original fitted data is used.
#' @param ... Additional arguments to be passed to other methods or functions.
#'
#' @return A matrix of fitted values.
#' @importFrom stats approx
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'            outcome = "Y", factor_var = "trt")
#'
#' newdata = data.frame(Z_1 = 1, Z_2 = 1, Z_3 = 1, Z_4 = 1, Z_5 = 1,
#'                    trt = 1, X_1 = seq(0,1,by=0.01))
#' fitted(results)
#' fitted(results,newdata)
#' }
#' @export
fitted.rjMCMC <- function(object, newdata = NULL,...) {
  results <- object
  factor_var = results$factor_var
  if (is.null(newdata)) {
    newdata = results$data_fit
  }

  all_splines = names(results$splines_fitted)
  all_binary = colnames(results$binary_param)

  # Check if newdata contains required columns
  required_columns <- c(results$candbinaryvars, results$candsplinevars, factor_var)
  missing_columns <- setdiff(required_columns, colnames(newdata))


  if (length(missing_columns) > 0) {
    stop(paste("Error: The following columns are missing in newdata:", paste(missing_columns, collapse = ", ")))
  }

  fitted_posterior =  results$fixed_param %*% t(cbind(1,newdata[,factor_var]))

  if (length(all_splines) > 0) {
    for (m in 1:length(all_splines)) {
      if (grepl(paste0(":",factor_var), all_splines[m])) {
        interpolated_splines_fitted <- sapply(1:nrow(results$splines_fitted[[all_splines[m]]]), function(j) {
          stats::approx(x = results$data_fit[,gsub(paste0(":",factor_var,"$"), "", all_splines[m])]*results$data_fit[,factor_var],
                 y = results$splines_fitted[[all_splines[m]]][j,]*results$data_fit[,factor_var], xout = newdata[,gsub(paste0(":",factor_var,"$"), "", all_splines[m])]*newdata[,factor_var],rule=2,ties="mean")$y
        })
      } else {
        interpolated_splines_fitted <- sapply(1:nrow(results$splines_fitted[[all_splines[m]]]), function(j) {
          stats::approx(x = results$data_fit[,all_splines[m]],
                 y = results$splines_fitted[[all_splines[m]]][j,], xout = newdata[,all_splines[m]],rule=2,ties="mean")$y
        })
      }
      fitted_posterior = fitted_posterior + t(interpolated_splines_fitted)
    }
  }

  interaction_terms <- colnames(results$binary_param)[grepl(paste0(":",factor_var,"$"), colnames(results$binary_param))]

  # Create interaction columns in newdata based on the matching regular terms
  for (interaction_term in interaction_terms) {
    # Get the base name before paste0(":",factor_var)
    base_name <- gsub(paste0(":",factor_var,"$"), "", interaction_term)

    # If the base term exists in newdata, create the interaction term
    if (base_name %in% colnames(newdata)) {
      newdata[[interaction_term]] <- newdata[[base_name]] * newdata[,factor_var]
    }
  }
  if (length(all_binary) > 0) {
    fitted_posterior = fitted_posterior +
      results$binary_param %*% t(as.matrix(newdata[,all_binary]))
  }

  return(fitted_posterior)
}

#' Fitted exposure effect values from Reversible Jump MCMC (rjMCMC) Model Results
#'
#' This function generates posterior fitted exposure effects from an rjMCMC model based on the provided data.
#' It combines the fixed effects, spline terms, and binary parameters.
#'
#' @param results An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#' @param newdata A data frame for which fitted values are to be computed. If NA, the original fitted data is used.
#'
#' @return A matrix of fitted values.
#' @importFrom stats approx
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' newdata = data.frame(Z_1 = 1, Z_2 = 1, Z_3 = 1, Z_4 = 1, Z_5 = 1,
#'                    trt = 1, X_1 = seq(0,1,by=0.01))
#' fittedExposureEff(results)
#' fittedExposureEff(results,newdata)
#' }
#' @export
fittedExposureEff <- function(results, newdata = NULL) {
  factor_var = results$factor_var
  if (is.null(newdata)) {
    newdata = results$data_fit
  }

  # Check if newdata contains required columns
  required_columns <- c(results$candinter)
  missing_columns <- setdiff(required_columns, colnames(newdata))


  if (length(missing_columns) > 0) {
    stop(paste("Error: The following columns are missing in newdata:", paste(missing_columns, collapse = ", ")))
  }

  fitted_posterior =  results$fixed_param[,factor_var] %*% matrix(1,ncol=nrow(newdata))
  candsplineinter = results$candsplineinter
  if (length(candsplineinter) > 0) {
    for (m in 1:length(candsplineinter)) {
      interpolated_splines_fitted <- sapply(1:nrow(results$splines_fitted[[paste0(candsplineinter[m],paste0(":",factor_var))]]), function(j) {
        stats::approx(x = results$data_fit[,candsplineinter[m]],
               y = results$splines_fitted[[paste0(candsplineinter[m],paste0(":",factor_var))]][j,], xout = newdata[,candsplineinter[m]],rule=2,ties="mean")$y
      })
      fitted_posterior = fitted_posterior + t(interpolated_splines_fitted)
    }
  }

  candbinaryinter <- colnames(results$binary_param)[grepl(paste0(":",factor_var,"$"), colnames(results$binary_param))]

  if (length(candbinaryinter) > 0) {
    fitted_posterior = fitted_posterior +
      results$binary_param[,candbinaryinter] %*% t(as.matrix(newdata[,gsub(paste0(":",factor_var,"$"), "",candbinaryinter)]))
  }

  return(fitted_posterior)
}


#' Predict Exposure Effect
#'
#' This function predicts the exposure effect for new data based the Reversible Jump MCMC (rjMCMC) results.
#'
#' @param results An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#' @param newdata A data frame for which predicted values are to be computed If NA, the original fitted data is used.
#'
#' @return A matrix of predictive posterior samples for the exposure effect, where each row corresponds to a posterior sample
#'   and each column corresponds to an observation in \code{newdata}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Checks if the required columns in \code{results$candinter} are present in \code{newdata}.
#'   \item Computes the fitted posterior exposure effect based on main exposure effects, spline interactions, and binary interactions.
#'   \item Adds noise to the fitted posterior using the residual variance \code{results$sigma_sq} to generate predictive posterior samples.
#' }
#'
#' Spline interactions are handled by interpolating the spline coefficients for the values in \code{newdata}.
#' @importFrom stats rnorm approx
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' newdata = data.frame(Z_1 = 1, Z_2 = 1, Z_3 = 1, Z_4 = 1, Z_5 = 1,
#'                    trt = 1, X_1 = seq(0,1,by=0.01))
#' predictExposureEff(results)
#' predictExposureEff(results,newdata)
#' }
#' @export
predictExposureEff <- function(results, newdata = NULL) {
  factor_var = results$factor_var
  if (is.null(newdata)) {
    newdata = results$data_fit
  }

  # Check if newdata contains required columns
  required_columns <- c(results$candinter)
  missing_columns <- setdiff(required_columns, colnames(newdata))


  if (length(missing_columns) > 0) {
    stop(paste("Error: The following columns are missing in newdata:", paste(missing_columns, collapse = ", ")))
  }

  fitted_posterior =  results$fixed_param[,factor_var] %*% matrix(1,ncol=nrow(newdata))
  candsplineinter = results$candsplineinter
  if (length(candsplineinter) > 0) {
    for (m in 1:length(candsplineinter)) {
      interpolated_splines_fitted <- sapply(1:nrow(results$splines_fitted[[paste0(candsplineinter[m],paste0(":",factor_var))]]), function(j) {
        stats::approx(x = results$data_fit[,candsplineinter[m]],
               y = results$splines_fitted[[paste0(candsplineinter[m],paste0(":",factor_var))]][j,], xout = newdata[,candsplineinter[m]],rule=2,ties="mean")$y
      })
      fitted_posterior = fitted_posterior + t(interpolated_splines_fitted)
    }
  }

  candbinaryinter <- colnames(results$binary_param)[grepl(paste0(":",factor_var,"$"), colnames(results$binary_param))]

  if (length(candbinaryinter) > 0) {
    fitted_posterior = fitted_posterior +
      results$binary_param[,candbinaryinter] %*% t(as.matrix(newdata[,gsub(paste0(":",factor_var,"$"), "",candbinaryinter)]))
  }

  residual_se <- sqrt(results$sigma_sq)
  predictive_posterior <- fitted_posterior + matrix(stats::rnorm(length(fitted_posterior), mean = 0, sd = residual_se), nrow = nrow(fitted_posterior), ncol = ncol(fitted_posterior))

  return(predictive_posterior)
}

#' Compute Posterior Inclusion Probabilities (PIPs) for rjMCMC Results
#'
#' This function returns the posterior inclusion probabilities (PIPs) for all variables,
#' including the intercept and exposure, based on the results from an rjMCMC model.
#'
#' @param results An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#'
#' @return A numeric vector with the PIPs for the intercept, exposure, and other variables.
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' pip(results)
#' }
#' @export
pip <- function(results) {
  factor_var = results$factor_var
  pip <- c(intercept = 1, factor_var = 1, results$vars_prop_summ)
  names(pip)[2] = factor_var
  return(pip)
}

#' Extract Posterior Mean Coefficients from rjMCMC Results
#'
#' This function extracts the posterior means of the intercept, exposure parameters,
#' and binary parameters from the results of an rjMCMC model.
#'
#' @param object An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#' @param ... Additional arguments to be passed to other methods or functions.
#'
#' @return A numeric vector containing the posterior mean of the intercept, exposure, and binary parameters.
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' coef(results)
#' }
#' @export
coef.rjMCMC <- function(object, ...) {
  results <- object
  c(apply(results$fixed_param, 2, mean),
          apply(results$binary_param, 2, mean))
}

#' Credible Intervals for rjMCMC Results
#'
#' This function calculates the posterior mean and credible intervals for parameters
#' from the rjMCMC results, including both intercept/exposure parameters and binary parameters.
#' The credible intervals are computed based on the specified confidence level.
#'
#' @param results An object of class rjMCMC containing the output from the `rjMCMC` procedure, including posterior samples.
#' @param level The level for the credible intervals (default is 0.95).
#'
#' @return A data frame with estimates, lower, and upper bounds of the credible intervals.
#' @importFrom stats quantile
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' credint(results)
#' }
#' @export
credint <- function(results, level = 0.95) {
  alpha = 1-level
  df = data.frame(Estimate = c(apply(results$fixed_param, 2, mean),
                          apply(results$binary_param, 2, mean)),
             Lower = c(apply(results$fixed_param, 2, stats::quantile,prob=alpha/2),
                       apply(results$binary_param, 2, stats::quantile, prob=1-alpha/2)),
             Upper = c(apply(results$fixed_param, 2,stats::quantile,prob=alpha/2),
                      apply(results$binary_param, 2, stats::quantile, prob=1-alpha/2)))
  colnames(df) = c("Estimate", paste0("Q",alpha/2*100), paste0("Q",(1-alpha/2)*100))
  return(df)
}

#' Generate Histogram and Trace Plots for MCMC Samples
#'
#' This internal function creates histogram or trace plots for MCMC samples from an `rjMCMC` model,
#' displaying the posterior distributions or sampling trajectories of specified variables.
#'
#' @param results A fitted model object from \code{rjMCMC}.
#' @param variables A vector of variable names to include in the plot, representing the estimand or covariates of interest.
#' @param sample_type Character string specifying the sample type: \code{"estimand"}, \code{"fitted"}, or \code{"predictive"}.
#'        When set to \code{"estimand"}, the function plots individual parameters from the posterior.
#'        For \code{"fitted"} or \code{"predictive"}, it uses posterior samples for one individual’s
#'        fitted or predictive values.
#' @param plot_type Character string specifying the plot type: \code{"hist"} for histogram or \code{"trace"} for trace plots.
#'        \code{"hist"} shows the distribution of parameter values, while \code{"trace"} displays MCMC sampling trajectories.
#' @param aux_vars A list of auxiliary variables and their values. Used to generate data for \code{fitted} or \code{predictive} values.
#'        Each element name should correspond to a model variable.
#' @param facet_by A vector of variable names to facet by in the plot. Automatically set to binary model variables if \code{NULL}.
#'
#' @details
#' - **Sample and Plot Compatibility**:
#'     - For \code{sample_type = "estimand"}, only \code{plot_type = "hist"} or "trace" is compatible, as it represents the posterior
#'       distribution or MCMC trajectory of individual parameters.
#'     - For \code{sample_type = "fitted"} or \code{"predictive"}, this function plots either the main effect
#'       (no interaction with exposure) or exposure effect for an individual using auxiliary variable values.
#' - **Effect Types**:
#'     - For \code{sample_type = "fitted"} or \code{"predictive"}, \code{"outcome"} generates plots without exposure interaction,
#'       and \code{"exposure_effect"} plots the interaction effect with exposure.
#'
#' @note This is an internal function not intended for direct use by package users.
#'
#' @return A \code{ggplot2} object displaying histograms or trace plots of the MCMC samples.
#' @importFrom ggplot2 ggplot geom_histogram geom_line labs facet_wrap theme_minimal ggtitle aes
#' @importFrom stats quantile predict
#' @importFrom rlang .data
#' @keywords internal
plotHistTrace <- function(results,
                          variables = NULL,
                          sample_type = "fitted",
                          effect_type = "exposure_effect",
                          plot_type = "cred",
                          aux_vars = list(),
                          facet_by = NULL) {

  chains <- results$mcmc_specs$chains
  num_persons <- ncol(results$trt_eff_posterior)
  B_per_chain <- nrow(results$trt_eff_posterior) / chains
  # Generate histogram and trace plots for specified covariate values
  hist_trace_plots <- list()

  if (sample_type == "estimand") {
    binary_vars <- intersect(c(variables,facet_by), colnames(results$binary_param))
    fixed_vars <- intersect(c(variables,facet_by), colnames(results$fixed_param))
    sigma_vars <- intersect(c(variables,facet_by), "sigma")

    # Prepare data for fixed_param plots (only if inter_vars exist)
    if (length(fixed_vars) > 0) {
      fixed_param_plot_df <- expand.grid(
        iterations = seq_len(B_per_chain),
        chains = as.factor(seq_len(chains)),
        param = fixed_vars
      )
      fixed_param_plot_df$value <- as.vector(results$fixed_param[, fixed_vars])
    } else {
      fixed_param_plot_df <- NULL
    }

    # Prepare data for binary_param plots (only if binary_vars exist)
    if (length(binary_vars) > 0) {
      binary_param_plot_df <- expand.grid(
        iterations = seq_len(B_per_chain),
        chains = as.factor(seq_len(chains)),
        param = binary_vars
      )
      binary_param_plot_df$value <- as.vector(results$binary_param[, binary_vars])
    } else {
      binary_param_plot_df <- NULL
    }

    if (length(sigma_vars) > 0) {
      sigma_param_plot_df <- expand.grid(
        iterations = seq_len(B_per_chain),
        chains = as.factor(seq_len(chains)),
        param = sigma_vars
      )
      sigma_param_plot_df$value <- sqrt(results$sigma_sq)
    } else {
      sigma_param_plot_df <- NULL
    }

    data_list <- list(fixed_param_plot_df, binary_param_plot_df, sigma_param_plot_df)
    plot_df <- do.call(rbind, Filter(Negate(is.null), data_list))

    if (is.null(plot_df) || nrow(plot_df) == 0) {
      stop("No data available for plotting.")
    }
  } else {
    newdata <- aux_vars

    # Adjust `facet_by` variables
    for (var in c(facet_by,variables)) {
      if (var %in% results$candbinaryvars) {
        # Set binary variables to c(0, 1)
        newdata[[var]] <- c(0, 1)
      } else if (var %in% colnames(results$data_fit)) {
        # Set continuous variables to two quantiles of X_1 in results$data_fit
        quantiles <- stats::quantile(results$data_fit[[var]], probs = c(0.25, 0.75))
        newdata[[var]] <- quantiles
      }
    }

    expanded_data <- expand.grid(newdata)
    mcmc_samples_all <- if (sample_type == "fitted") {
      if (effect_type == "outcome") fitted(results, expanded_data) else fittedExposureEff(results, expanded_data)
    } else {
      if (effect_type == "outcome") predict(results, expanded_data) else predictExposureEff(results, expanded_data)
    }

    label_text = rep(NA,nrow(expanded_data))
    for (i in 1:nrow(expanded_data)) {
      label_text[i] <- paste(sapply(names(expanded_data), function(name) {
        paste0(name, "=", round(expanded_data[i, name],2))
      }), collapse = ", ")
    }


    plot_df = expand.grid(
      iterations = seq_len(B_per_chain),
      chains = as.factor(seq_len(chains)),
      param = label_text
    )

    plot_df$value = as.vector(mcmc_samples_all)
    if (is.null(plot_df) || nrow(plot_df) == 0) {
      stop("No data available for plotting.")
    }
  }



  if (plot_type == "hist") {
    title <- paste0(
      "Histogram of ", sample_type, " MCMC samples for ",
      if (effect_type == "outcome") "main effect" else if (sample_type == "estimand") "interaction effect" else "exposure effect"
    )
    p_hist <- plot_df %>%
      ggplot2::ggplot(ggplot2::aes(x=.data$value, fill = .data$chains)) +
      ggplot2::geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
      ggplot2::labs(x = "Value", y = "Count", fill = "Chain") +
      ggplot2::facet_wrap(~param, scales = "free",ncol=1) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(title)
    print(p_hist)
  } else if (plot_type == "trace") {
    title <- paste0(
      "Trace plot of ", sample_type, " MCMC samples for ",
      if (effect_type == "outcome") "main effect" else if (sample_type == "estimand") "interaction effect" else "exposure effect"
    )
    p_trace <- plot_df %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$iterations, y = .data$value, col = .data$chains)) +
      ggplot2::geom_line(alpha = 0.7, linetype = "dashed") +
      ggplot2::labs(x = "Iteration", y = "Value", col = "Chain") +
      ggplot2::facet_wrap(~param, scales = "free_y",ncol=1) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(title)
    print(p_trace)
  }
}

#' Plotting function for \code{rjMCMC} results
#'
#' This function generates plots for model results from \code{rjMCMC} based on specified sample type, effect type, and plot type.
#' The function is flexible for various combinations of \code{sample_type}, \code{effect_type}, and \code{plot_type}, as outlined below.
#'
#' @param x A fitted model object from \code{rjMCMC}.
#' @param ... Additional arguments to be passed to other methods or functions.
#' @param variables A vector of variable names to include in the plot. Automatically set to continuous variables if \code{NULL}.
#' @param sample_type Character string specifying the type of sample: "fitted", "predictive", or "estimand".
#'                    "fitted" and "predictive" are compatible with \code{plot_type = "cred"}.
#'                    "estimand" is compatible with \code{plot_type = "hist"} or "trace" (only used for individual parameter trajectories).
#' @param effect_type Character string indicating the effect type: "exposure_effect" or "outcome".
#'                    For "exposure_effect", the function plots the fitted or predictive effect of exposure; for "outcome", it plots
#'                    outcome values without interaction with exposure.
#' @param plot_type Character string specifying the plot type: "cred" for credible interval plots, or "hist"/"trace" for histogram or trace plots
#'                  of individual parameters (only for \code{sample_type = "estimand"}).
#' @param level Numeric value for the credible interval level (default is 0.95).
#' @param aux_vars A list of auxiliary variables and their fixed values. Each element name must match a model variable.
#' @param facet_by A vector of variable names to facet by in the plot. Automatically set to binary model variables if \code{NULL}.
#' @param pip_cutoff Numeric threshold for the posterior inclusion probability (PIP) of model variables to include in the plot.
#'
#' @details
#' - **Sample and Plot Compatibility**:
#'     - For \code{sample_type = "estimand"}, only \code{plot_type = "hist"} or "trace" is allowed, as these are designed to visualize the posterior distribution
#'       or MCMC trajectory of individual parameters. Parameters like \code{intercept}, \code{trt}, and \code{sigma} are agnostic to \code{effect_type} as they
#'       do not interact with exposure.
#'     - \code{plot_type = "cred"} is designed for use with \code{sample_type = "fitted"} or "predictive" and shows credible intervals for the outcome (y-axis)
#'       across biomarker values (x-axis) by covariate pattern. \code{effect_type} controls whether the exposure effect or main effect is displayed.
#' - **Effect Types**:
#'     - \code{outcome} plots either the fitted or predictive values without exposure interaction, allowing for exposure (\code{trt}) values to be specified.
#'     - \code{exposure_effect} plots the interaction of the exposure effect across different covariate patterns.
#'
#' @return A \code{ggplot2} object or a grid of plots.
#'
#' @importFrom ggplot2 ggplot geom_point geom_errorbar xlab ylab ggtitle geom_ribbon facet_wrap aes
#' @importFrom dplyr group_by summarize select across all_of %>% n_distinct mutate
#' @importFrom cowplot plot_grid
#' @importFrom stats as.formula predict fitted
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' plot(results, sample_type = "fitted", effect_type = "exposure_effect", plot_type = "cred")
#' plot(results, sample_type = "estimand", plot_type = "hist")
#' }
#' @export
plot.rjMCMC <- function(x,
                        ...,
                        variables = NULL,
                        sample_type = "fitted",
                        effect_type = "exposure_effect",
                        plot_type = "cred",
                        level = 0.95,
                        aux_vars = list(),
                        facet_by = NULL,
                        pip_cutoff = 0.1) {
  results <- x
  factor_var = results$factor_var
  # Check if any elements of `facet_by` are in `variables`
  common_elements <- intersect(facet_by, variables)
  # Throw an error if there are common elements
  if (length(common_elements) > 0) {
    stop(paste("Error: The following variables are in both `facet_by` and `variables`, which is not allowed:",
               paste(common_elements, collapse = ", ")))
  }

  # Check compatibility of `sample_type` and `plot_type`
  if (sample_type == "estimand" && plot_type == "cred") {
    # Automatically set `plot_type` to "trace" and provide an informative message
    message("`plot_type` set to 'trace' for use with `sample_type = estimand`. Only 'trace' and 'hist' plot types are available for `estimand` sample_type.")
    plot_type <- "trace"
  }

  # Validate `aux_vars`
  if (!is.list(aux_vars)) {
    stop("aux_vars must be a list where each element corresponds to a model variable.")
  }

  # List of all model variables
  model_vars <- c(results$candsplinevars, results$candbinaryvars)
  invalid_aux_vars <- setdiff(names(aux_vars), c(factor_var,model_vars))
  if (length(invalid_aux_vars) > 0) {
    stop("The following variables in aux_vars are not part of the model variables: ", paste(invalid_aux_vars, collapse = ", "))
  }

  # Validate `plot_type`
  valid_plot_types <- c("cred", "hist", "trace")
  if (!plot_type %in% valid_plot_types) {
    stop(paste("Invalid sample_type:", plot_type,
               "- must be one of", paste(valid_plot_types, collapse = ", ")))
  }

  # Validate `sample_type`
  valid_sample_types <- c("fitted", "predictive", "estimand")
  if (!sample_type %in% valid_sample_types) {
    stop(paste("Invalid sample_type:", sample_type,
               "- must be one of", paste(valid_sample_types, collapse = ", ")))
  }

  # Validate `effect_type`
  valid_effect_types <- c("exposure_effect", "outcome")
  if (!effect_type %in% valid_effect_types) {
    stop(paste("Invalid effect_type:", effect_type,
               "- must be one of", paste(valid_effect_types, collapse = ", ")))
  }

  # Exclude intercept and trt, which always have pip = 1
  pip_cov = pip(results)[-c(1,2)]
  pip_cov = pip_cov[pip_cov > pip_cutoff]
  # (1) Get variables without paste0(":",factor_var) for pip_cov_main
  pip_cov_main <- names(pip_cov)[!grepl(paste0(":",factor_var), names(pip_cov))]

  # (2) Get variables with paste0(":",factor_var) and remove paste0(":",factor_var) for pip_cov_inter
  pip_cov_inter <- gsub(paste0(":",factor_var), "", names(pip_cov)[grepl(paste0(":",factor_var), names(pip_cov))])

  # Set `variables` automatically if it is NULL
  if (is.null(variables)) {
    message("Automatically setting variables to be continuous model variables with pip > pip_cutoff.")
    variables <- intersect(if (effect_type == "exposure_effect") pip_cov_inter else pip_cov_main, results$candsplinevars)
    # Set `facet_by` automatically if it is NULL, only in the case when variables is also NULL
    # (reserves the right for users to specify facet_by = NULL generally)
    if (is.null(facet_by)) {
      message("Automatically setting facet_by to be binary model variables with pip > pip_cutoff.")
      facet_by <- intersect(if (effect_type == "exposure_effect") pip_cov_inter else pip_cov_main, results$candbinaryvars)
    }
  }

  # Identify missing variables not in aux_vars, variables, or facet_by
  provided_vars <- c(names(aux_vars), variables, facet_by)
  missing_vars <- if (effect_type == "outcome") setdiff(c(model_vars,factor_var), provided_vars) else setdiff(model_vars, provided_vars)

  # Set missing variables to 0 and notify user
  if (length(missing_vars) > 0 & sample_type != "estimand") {
    message("The following variables were not provided and will be set to 0: ", paste(missing_vars, collapse = ", "))
    aux_vars[missing_vars] <- 0
  }

  # Error handling if variables is empty
  if (length(variables) == 0) {
    stop("Error: variables is not provided and there are no continuous variables with pip > pip_cutoff to plot.")
  }

  if (plot_type == "cred") {
    # Initialize a list to store plots
    plots <- list()

    # Loop over each variables and generate a plot
    for (var in variables) {

      # Generate `newdata` from `aux_vars`
      newdata <- aux_vars

      if (length(variables) > 0) {
        # Set auxiliary variables to 0 except for the current var and facet_by variables
        within_plot_aux_vars <- setdiff(variables, var)
        vars_to_set <- setdiff(within_plot_aux_vars, names(aux_vars))

        # Set those variables to 0 in `newdata` and print a message if any were set
        if (length(vars_to_set) > 0) {
          newdata[vars_to_set] <- lapply(vars_to_set, function(x) 0)
          message_text <- paste0("The following variables will be set to 0 in the ", var, " plot: ",
                                 paste(vars_to_set, collapse = ", "))
          warning(message_text)
        }
      }

      # Adjust `var` if it is continuous
      if (var %in% colnames(results$data_fit)) {
        if (var %in% results$candbinaryvars) {
          newdata[[var]] <- c(0, 1)
        } else {
          newdata[[var]] <- seq(min(results$data_fit[[var]]), max(results$data_fit[[var]]), length.out = 100)
        }
      }

      # Adjust `facet_by` variables
      for (facet_var in facet_by) {
        if (facet_var %in% names(newdata)) {
          # Skip if already provided in aux_vars
          next
        } else if (facet_var %in% results$candbinaryvars) {
          # Set binary variables to c(0, 1)
          newdata[[facet_var]] <- c(0, 1)
        } else if (facet_var %in% colnames(results$data_fit)) {
          # Set continuous variables to three quantiles of X_1 in results$data_fit
          quantiles <- stats::quantile(results$data_fit[[facet_var]], probs = c(0.25, 0.5, 0.75))
          newdata[[facet_var]] <- quantiles
        }
      }

      # Create newdata as an expand.grid from the list `newdata`
      expanded_data <- expand.grid(newdata)

      # Ensure other variables are identical across rows
      if (!is.null(facet_by)) {
        vars_to_check <- setdiff(colnames(expanded_data), c(var, facet_by))
        inconsistent_rows <- expanded_data %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(facet_by))) %>%
          dplyr::summarize(dplyr::across(dplyr::all_of(vars_to_check), ~ dplyr::n_distinct(.) > 1), .groups = "drop") %>%
          dplyr::select(-all_of(facet_by)) %>%
          rowSums() %>%
          sum()

        if (inconsistent_rows > 0) {
          stop("Within each level of facet_by, all other variables must be constant except for variables.")
        }
      }

      if (effect_type == "outcome") {
        expanded_data <- expanded_data %>%
          mutate(across(all_of(results$candinter),
                        ~ . * !!sym(facet_var),
                        .names = "{.col}:{facet_var}"))
      }



      # Get mcmc_estimates based on sample_type and effect_type
      mcmc_estimates <- if (sample_type == "fitted") {
        if (effect_type == "outcome") fitted(results, expanded_data) else fittedExposureEff(results, expanded_data)
      } else {
        if (effect_type == "outcome") predict(results, expanded_data) else predictExposureEff(results, expanded_data)
      }

      fixed_values <- sapply(newdata, function(col) unique(col))
      fixed_values <- fixed_values[names(fixed_values)!=var]
      if (effect_type == "exposure_effect") {
        fixed_values <- fixed_values[names(fixed_values) != factor_var]
      }

      if (!is.null(facet_by)) {
        fixed_values <- fixed_values[!(names(fixed_values) %in% facet_by)]
      }


      fixed_values_text <- paste(names(fixed_values), "=", fixed_values, collapse = ", ")

      # Set the plot title
      title <- if (effect_type == "outcome") {
        effect_type_text <- if (sample_type == "fitted") "Fitted Outcome" else "Predictive Outcome"
        paste(effect_type_text, "of", var, "for", fixed_values_text)
      } else {
        effect_type_text <- if (sample_type == "fitted") "Fitted Exposure Effect" else "Predictive Exposure Effect"
        paste(effect_type_text, "of", var, "for", fixed_values_text)
      }

      mean <- apply(mcmc_estimates, 2, mean)
      alpha <- 1 - level
      lower_ci <- apply(mcmc_estimates, 2, function(x) stats::quantile(x, alpha / 2))
      upper_ci <- apply(mcmc_estimates, 2, function(x) stats::quantile(x, 1 - alpha / 2))

      plot_data <- data.frame(
        x = expanded_data[[var]],
        mean = mean,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        expanded_data[facet_by]
      )

      ylab <- ifelse(effect_type == "outcome",
                  paste0("Outcome (", level * 100, "% Credible Interval)"),
                  paste0("Effect (", level * 100, "% Credible Interval)"))

      # Generate plot for each var and add it to the plots list
      p <- if (var %in% results$candbinaryvars) {
        ggplot2:: ggplot(plot_data, ggplot2::aes(x = as.factor(x), y = mean)) +
          ggplot2::geom_point() +
          ggplot2::geom_errorbar(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
          ggplot2::xlab(var) +
          ggplot2::ylab(ylab) +
          ggplot2::ggtitle(title)
      } else {
        ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = mean)) +
          ggplot2::geom_line() +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
          ggplot2::xlab(var) +
          ggplot2::ylab(ylab) +
          ggplot2::ggtitle(title)
      }

      if (!is.null(facet_by)) {
        p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", paste(facet_by, collapse = " + "))), labeller = "label_both")
      }

      plots[[var]] <- p
    }

    # Combine all plots in a grid
    cowplot::plot_grid(plotlist = plots)
  } else {
    plotHistTrace(results,
                  variables,
                  sample_type,
                  effect_type,
                  plot_type,
                  aux_vars,
                  facet_by)
  }
}


#' Get Effective Subspace
#'
#' This function identifies the "effective subspace" where exposure is effective
#' based on posterior inference results from the FK-BMA model. It analyzes interaction terms
#' between exposure and covariates, allowing for both binary and continuous variables.
#'
#' @param results A fitted model object from \code{rjMCMC}.
#' @param newdata Optional. A new dataset for evaluating the effective subspace.
#'   If \code{NULL}, the function uses \code{results$data_fit}.
#' @param alpha Numeric. The alpha level used for computing quantiles. Default is \code{0.05}.
#' @param pip_cutoff Numeric. The minimum Posterior Inclusion Probability (PIP)
#'   threshold for selecting covariates. Default is \code{0.1}.
#'
#' @return A list with the following components:
#'   \item{quantiles}{A vector of quantile values for the exposure effect in the new dataset.}
#'   \item{is_effective_subspace}{A logical vector indicating whether the exposure effect is
#'   positive in the effective subspace.}
#'
#' @details
#' - The function computes the posterior exposure effect for each observation in the
#'   dataset using the \code{fittedExposureEff} function and evaluates its quantiles at the
#'   specified \code{alpha} level.
#' - Binary variables with high posterior inclusion probabilities (PIP) are used to define
#'   subgroups, and the corresponding effective subspaces for a continuous variable
#'   are identified by checking where the exposure effect quantiles are strictly positive.
#' - If the number of binary variables is \code{<= 3} and there is exactly one continuous variable,
#'   the function describes the effective subspace in terms of disjoint intervals.
#' - For more complex cases, a warning is issued suggesting alternative methods such as
#'   Bayesian regression trees for interpretation.
#'
#' @importFrom matrixStats colQuantiles
#' @importFrom stats quantile
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' getEffectiveSubspace(results)
#' }
#' @export
getEffectiveSubspace <- function(results,
                                 newdata = NULL,
                                 alpha = 0.05,
                                 pip_cutoff = 0.1) {
  factor_var = results$factor_var
  if (is.null(newdata)) {
    newdata = results$data_fit
  }
  mcmc_estimates = fittedExposureEff(results, newdata)
  quantiles = matrixStats::colQuantiles(mcmc_estimates,probs=alpha)
  newdata$quantiles = quantiles
  pip_cov = pip(results)[-c(1,2)]
  pip_cov = pip_cov[pip_cov > pip_cutoff]

  # (2) Get variables with paste0(":",factor_var) and remove paste0(":",factor_var) for pip_cov_inter
  pip_cov_inter <- gsub(paste0(":",factor_var), "", names(pip_cov)[grepl(paste0(":",factor_var), names(pip_cov))])
  binary_vars = intersect(pip_cov_inter,results$candbinaryvars)
  continuous_var = intersect(pip_cov_inter,results$candsplinevars)
  all_vars <- c(binary_vars, continuous_var)
  # Initialize descriptions
  descriptions <- list()

  # Create all possible combinations of binary variable levels
  binary_levels <- expand.grid(lapply(binary_vars, function(var) unique(newdata[[var]])))
  colnames(binary_levels) <- binary_vars

  if (length(binary_vars) <= 3 && length(continuous_var) == 1) {
    message("Effective subspace descriptions:\n")
    for (i in 1:nrow(binary_levels)) {
      # Subset data for the current binary combination
      subset_data <- newdata
      for (var in binary_vars) {
        subset_data <- subset_data[subset_data[[var]] == binary_levels[i, var], , drop = FALSE]
      }

      if (nrow(subset_data) > 0) {
        # Extract continuous variable and corresponding quantiles
        continuous_data <- subset_data[[continuous_var]]
        continuous_quantiles <- subset_data$quantiles
        sorted_indices <- order(continuous_data)
        sorted_data <- continuous_data[sorted_indices]
        sorted_quantiles <- continuous_quantiles[sorted_indices]

        # Check if the quantiles never cross zero
        if (all(sorted_quantiles > 0)) {
          # All quantiles > 0: entire range is effective
          binary_description <- paste0(
            paste0(binary_vars, " = ", binary_levels[i, ], collapse = ", "),
            ": ", continuous_var, " in [", round(min(sorted_data), 2), ", ", round(max(sorted_data), 2), "]"
          )
          descriptions[[length(descriptions) + 1]] <- binary_description
        } else if (all(sorted_quantiles <= 0)) {
          # All quantiles <= 0: this binary combination is not in the effective subspace
          next
        } else {
          # Quantiles cross zero: find disjoint intervals where quantiles > 0
          transitions <- c(1,which(diff(sign(sorted_quantiles)) != 0),length(sorted_quantiles))
          positive_intervals <- list()
          for (j in 1:(length(transitions)-1)) {
            if (sorted_quantiles[transitions[j]+1]>0) {
              positive_intervals <- c(
                positive_intervals,
                list(c(sorted_data[transitions[j]], sorted_data[transitions[j+1]]))
              )
            }
          }

          # Format the intervals
          interval_strings <- sapply(positive_intervals, function(interval) {
            paste0("[", round(interval[1], 2), ", ", round(interval[2], 2), "]")
          })
          binary_description <- paste0(
            paste0(binary_vars, " = ", binary_levels[i, ], collapse = ", "),
            ": ", continuous_var, " in ", paste(interval_strings, collapse = " or ")
          )
          descriptions[[length(descriptions) + 1]] <- binary_description
        }
      }
    }
    message(paste(descriptions, collapse = ", "))
  } else {
    warning(
      "The decision rules may be too complex to interpret. Consider using other methods, such as a Bayesian regression tree, to identify effective subspaces."
    )
  }


  return(list(quantiles = quantiles, is_effecitve_subspace = quantiles > 0))
}


#' Calculate and Print Rhat Diagnostics for rjMCMC Results
#'
#' This function calculates the Rhat diagnostic for convergence based on the posterior samples
#' of individual exposure effects, intercept, and main exposure effect from an rjMCMC model.
#' It prints the median, minimum, and maximum Rhat values for the exposure effects, as well as
#' the Rhat for the intercept and exposure effect.
#'
#' @param results An object of class rjMCMC containing the output from the `rjMCMC` procedure, which includes:
#' \describe{
#'   \item{fixed_param}{Matrix of posterior samples for exposure intercept and main effect.}
#'   \item{binary_param}{Matrix of posterior samples for binary variable parameters.}
#'   \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'   \item{vars_prop_summ}{Posterior inclusion probabilities for candidate variables.}
#'   \item{splines_fitted}{List of matrices containing fitted values for spline terms across iterations.}
#'   \item{data_fit}{Original dataset used in the `rjMCMC` procedure.}
#'   \item{candsplineinter}{Names of continuous candidate predictive spline variables.}
#'   \item{candsplinevars}{Names of continuous candidate spline variables.}
#'   \item{candbinaryvars}{Names of binary candidate variables.}
#'   \item{candinter}{Names of interaction terms, which can include spline variables.}
#'   \item{mcmc_specs}{MCMC sampler specifications, including the number of iterations, burn-in, thinning, and chains.}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{Rhat_trt_eff_posterior}{Vector of R-hat values for each individual's exposure effect.}
#'   \item{Rhat_inter}{R-hat value for the intercept parameter.}
#'   \item{Rhat_factor_var}{R-hat value for the main effect of exposure.}
#' }
#'
#' @details
#' This function calculates R-hat statistics to assess MCMC convergence for both exposure effects and model parameters.
#' Diagnostic plots are generated to visually inspect the chains across iterations.
#'
#' @importFrom rstan Rhat
#' @importFrom stats median
#' @examples
#' \donttest{
#' # Example dataset
#' data("simulated_data")
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:5)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
#'                   outcome = "Y", factor_var = "trt")
#' rhats(results)
#' }
#' @export
rhats <- function(results) {
  factor_var = results$factor_var
  chains = results$mcmc_specs$chains
  num_persons = ncol(results$trt_eff_posterior)
  B_per_chain = nrow(results$trt_eff_posterior)/chains
  Rhat_trt_eff_posterior <- numeric(ncol(results$trt_eff_posterior))
  # Loop over each column (parameter)
  for (col in seq_len(num_persons)) {
    # Split the column into num_chains matrices of dimensions (num_iterations x num_chains)
    chain_matrix <- matrix(results$trt_eff_posterior[, col], nrow = B_per_chain, ncol = chains, byrow = FALSE)

    # Apply the Rhat function to the matrix
    Rhat_trt_eff_posterior[col] <- rstan::Rhat(chain_matrix)
  }

  message(paste0("Rhat for individual exposure effects: (Median [Range]) ",
               round(stats::median(Rhat_trt_eff_posterior),4), " [", round(min(Rhat_trt_eff_posterior),4),",", round(max(Rhat_trt_eff_posterior),4),"]"))

  chain_matrix_inter <- matrix(results$fixed_param[, 1], nrow = B_per_chain, ncol = chains, byrow = FALSE)
  Rhat_inter <- rstan::Rhat(chain_matrix_inter)
  message(paste0("Rhat for intercept: ", round(Rhat_inter,4) ))

  chain_matrix_trt <- matrix(results$fixed_param[, 2], nrow = B_per_chain, ncol = chains, byrow = FALSE)
  Rhat_factor_var <- rstan::Rhat(chain_matrix_trt)
  message(paste0("Rhat for main effect of ",factor_var,": ", round(Rhat_factor_var,4) ))

  return(list(
    Rhat_trt_eff_posterior = Rhat_trt_eff_posterior,
    Rhat_inter = Rhat_inter,
    Rhat_factor_var = Rhat_factor_var
  ))
}
