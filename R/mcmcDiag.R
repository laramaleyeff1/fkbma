library(ggplot2)
library(rstan)
library(dplyr)
#' MCMC Diagnostics for Treatment Effect Posterior Distributions
#'
#' Perform MCMC diagnostics for treatment effect posterior distributions, intercept, and main treatment effect
#' using R-hat diagnostics. This function computes the R-hat values to assess the convergence of MCMC chains for
#' each individual's treatment effect and model parameters. It also generates diagnostic plots showing the
#' evolution of the posterior distributions across iterations for selected individuals and parameters.
#'
#' @param trial_results A list containing the output from `rjMCMC`, including:
#' \describe{
#'   \item{trt_eff_posterior}{Matrix (rows = iterations, columns = individuals) of posterior treatment effects.}
#'   \item{inter_trt_param}{Matrix of posterior estimates for intercept and main treatment effect.}
#' }
#' @param chains Integer indicating the number of chains used during MCMC sampling.
#' @param ids (Optional) Vector specifying the indices of individuals for plotting treatment effects. If not provided,
#' 8 individuals will be randomly selected.
#'
#' @return A list containing:
#' \describe{
#'   \item{Rhat_trt_eff_posterior}{Vector of R-hat values for each individual's treatment effect.}
#'   \item{Rhat_inter}{R-hat value for the intercept parameter.}
#'   \item{Rhat_trt}{R-hat value for the main effect of treatment.}
#'   \item{trt_eff_posterior_plot_df}{Data frame for the treatment effect plot.}
#'   \item{inter_trt_param_plot_df}{Data frame for the intercept and treatment effect plot.}
#' }
#'
#' @details
#' This function calculates R-hat statistics to assess MCMC convergence for both treatment effects and model parameters.
#' Diagnostic plots are generated to visually inspect the chains across iterations.
#'
#' @export
#'
#' @examples
#' # Example using a dataset and rjMCMC results
#' n <- 1000
#' data <- data.frame(
#'   X_1 = runif(n, 0, 1),
#'   Z_1 = rbinom(n, 1, 0.35),
#'   Z_2 = rbinom(n, 1, 0.5),
#'   trt = rbinom(n, 1, 0.5)
#' )
#' data$Y <- 2 * data$Z_1 + 2 * data$Z_1 * data$trt + rnorm(n, 0, 0.1)
#'
#' candsplinevars <- c("X_1")
#' candbinaryvars <- paste0("Z_", 1:2)
#' candinter <- c(candsplinevars, candbinaryvars)
#'
#' # Simulating trial results using rjMCMC function
#' mcmc_specs <- list(B = 1000, burnin = 1000, thin = 1, chains = 2, sigma_v = 0.1, bma = TRUE)
#' prior_params <- list(lambda_1 = 0.1, lambda_2 = 1, a_0 = 0.01, b_0 = 0.01, degree = 3, k_max = 9, w = 1, sigma_B = sqrt(20))
#'
#' trial_results <- rjMCMC(data, candsplinevars, candbinaryvars, candinter, mcmc_specs, prior_params)
#'
#' # Performing MCMC diagnostics
#' mcmcDiag(trial_results)
mcmcDiag <- function(trial_results, ids = NULL) {
  chains = trial_results$mcmc_specs$chains
  num_persons = ncol(trial_results$trt_eff_posterior)
  B_per_chain = nrow(trial_results$trt_eff_posterior)/chains
  Rhat_trt_eff_posterior <- numeric(ncol(trial_results$trt_eff_posterior))
  # Loop over each column (parameter)
  for (col in seq_len(num_persons)) {
    # Split the column into num_chains matrices of dimensions (num_iterations x num_chains)
    chain_matrix <- matrix(trial_results$trt_eff_posterior[, col], nrow = B_per_chain, ncol = chains, byrow = FALSE)

    # Apply the Rhat function to the matrix (requires the coda package)
    Rhat_trt_eff_posterior[col] <- Rhat(chain_matrix)  # Assuming coda::Rhat
  }

  print(paste0("Rhat for individual trtment effects: (Median [Range]) ",
               round(median(Rhat_trt_eff_posterior),4), " [", round(min(Rhat_trt_eff_posterior),4),",", round(max(Rhat_trt_eff_posterior),4),"]"))

  chain_matrix_inter <- matrix(trial_results$inter_trt_param[, 1], nrow = B_per_chain, ncol = chains, byrow = FALSE)
  Rhat_inter <- Rhat(chain_matrix_inter)
  print(paste0("Rhat for intercept: ", round(Rhat_inter,4) ))

  chain_matrix_trt <- matrix(trial_results$inter_trt_param[, 2], nrow = B_per_chain, ncol = chains, byrow = FALSE)
  Rhat_trt <- Rhat(chain_matrix_trt)
  print(paste0("Rhat for main effect of treatment: ", round(Rhat_trt,4) ))

  if (is.null(ids)) {
    ids <- sample(num_persons,8)
  }
  trt_eff_posterior_plot_df <- expand.grid(
    iterations = seq_len(B_per_chain),
    chains = as.factor(seq_len(chains)),
    person = paste0("Person ", ids)
  )

  trt_eff_posterior_plot_df$value <- as.vector(trial_results$trt_eff_posterior[,ids])
  p1 <- trt_eff_posterior_plot_df %>%
    ggplot(aes(x=iterations,y=value,col=chains)) +
    geom_line(alpha=0.7,linetype="dashed") +
    labs(x="Iteration", y="Treatment effect",col="") +
    facet_wrap(~person,scales="free_y")

  print(p1)
  inter_trt_param_plot_df <- expand.grid(
    iterations = seq_len(B_per_chain),
    chains = as.factor(seq_len(chains)),
    param = c("Intercept", "Main effect of treatment")
  )

  inter_trt_param_plot_df$value <- as.vector(trial_results$inter_trt_param)
  readline(prompt = "Press [Enter] to see the next plot...")
  p2 <- inter_trt_param_plot_df %>%
    ggplot(aes(x=iterations,y=value,col=chains)) +
    geom_line(alpha=0.7,linetype="dashed") +
    labs(x="Iteration", y="Treatment effect",col="") +
    facet_wrap(~param,scales="free_y")
  print(p2)

  return(list(
    Rhat_trt_eff_posterior = Rhat_trt_eff_posterior,
    Rhat_inter = Rhat_inter,
    Rhat_trt = Rhat_trt,
    trt_eff_posterior_plot_df = trt_eff_posterior_plot_df,
    inter_trt_param_plot_df = inter_trt_param_plot_df
  ))
}
