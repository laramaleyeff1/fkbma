% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_maleyeff.R
\name{rjMCMC}
\alias{rjMCMC}
\title{Run Reversible Jump MCMC (rjMCMC) Procedure}
\usage{
rjMCMC(
  data,
  candsplinevars,
  candbinaryvars,
  candinter,
  mcmc_specs = NULL,
  prior_params = NULL
)
}
\arguments{
\item{data}{A data frame containing the observations, including the following columns:
\describe{
  \item{trt}{Group indicator (=1 for experimental group; =0 for control group).}
  \item{Y}{Continuous-valued outcome.}
  \item{candsplinevars}{All candidate spline variables as described in `candsplinevars`.}
  \item{candbinaryvars}{All candidate binary variables as described in `candbinaryvars`.}
}}

\item{candsplinevars}{A vector of names for continuous predictive candidate variables (default = NULL).}

\item{candbinaryvars}{A vector of names for binary predictive candidate variables (default = NULL).}

\item{candinter}{A vector indicating which of the candidate variables are tailoring (default = NULL).}

\item{mcmc_specs}{A list containing:
\describe{
  \item{B}{Number of posterior samples (default = 2000).}
  \item{burnin}{Number of burn-in samples (default = 10000).}
  \item{thin}{Thinning parameter (default = 5).}
  \item{chain}{Number of chains (default = 1).}
  \item{sigma_v}{Proposal variance for "jump" terms (default = 0.1).}
  \item{bma}{Boolean indicating whether to include Bayesian model averaging step (default = TRUE).}
}}

\item{prior_params}{A list containing prior parameters:
\describe{
  \item{lambda_1}{Prior parameter for the number of terms in the model (default = 0.1).}
  \item{lambda_2}{Prior parameter for the number of knots in each spline (default = 1).}
  \item{a_0}{Shape parameter for inverse gamma prior on individual-level variance (default = 0.01).}
  \item{b_0}{Rate parameter for inverse gamma prior on individual-level variance (default = 0.01).}
  \item{degree}{Degree of B-splines (default = 3).}
  \item{k_max}{Maximum number of knots for each spline term (default = 9).}
  \item{w}{Window for proposing knot location changes (default = 1).}
  \item{sigma_B}{Prior normal variance for model coefficients (default = sqrt(20)).}
}}
}
\value{
If the procedure is successful, a list containing:
\describe{
  \item{success}{Boolean indicating whether the procedure was successful based on Geweke convergence.}
  \item{accept_var}{Matrix of whether proposed variable addition/removal was accepted for each iteration.}
  \item{accept_add_knot}{Matrix of whether proposed knot addition was accepted for each iteration.}
  \item{accept_remove_knot}{Matrix of whether proposed knot removal was accepted for each iteration.}
  \item{accept_move_knot}{Matrix of whether proposed knot position change was accepted for each iteration.}
  \item{trt_eff_posterior}{Posterior distribution of treatment effects for each individual.}
  \item{splines_fitted}{Fitted values for interaction spline terms, used for prediction in new data.}
  \item{binary_param}{Posterior distribution of binary variable parameters.}
  \item{inter_trt_param}{Posterior distribution of treatment intercept and main effect.}
  \item{sigma_sq}{Posterior distribution of model standard deviation.}
  \item{k}{Posterior distribution of the number of knots for each spline term.}
  \item{vars_prop_summ}{Posterior inclusion probability for each term.}
}
}
\description{
This internal function performs a Reversible Jump MCMC (rjMCMC) procedure to generate the posterior distribution for one chain,
using Bayesian model averaging and free-knot B-splines.
}
\examples{
# Example dataset
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

results <- rjMCMC(data, candsplinevars, candbinaryvars, candinter, mcmc_specs, prior_params)
}
