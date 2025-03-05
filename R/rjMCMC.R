######################################################################################################
#                       Implementation of proposed rjMCMC procedure of Maleyeff et al. (2024)
#                                   Contact: laramaleyeff@gmail.com
######################################################################################################

#' Run Reversible Jump MCMC (rjMCMC) Procedure
#'
#' This function performs a Reversible Jump MCMC (rjMCMC) procedure to generate the posterior distribution,
#' using Bayesian model averaging and free-knot B-splines.
#'
#' @param data A data frame containing the observations, including the following columns:
#' \describe{
#'   \item{factor_var}{Exposure group binary indicator.}
#'   \item{outcome}{Continuous-valued outcome.}
#'   \item{candsplinevars}{All candidate spline variables as described in `candsplinevars`.}
#'   \item{candbinaryvars}{All candidate binary variables as described in `candbinaryvars`.}
#' }
#' @param candsplinevars A vector of names for continuous predictive candidate variables (default = NULL).
#' @param candbinaryvars A vector of names for binary predictive candidate variables (default = NULL).
#' @param candinter A vector indicating which of the candidate variables are tailoring (default = NULL).
#' @param outcome A string indicating the name of the continuous outcome variable in data.
#' @param factor_var A string indicating the name of the binary exposure variable in data.
#' @param mcmc_specs A list containing:
#' \describe{
#'   \item{iter}{Number of total iterations per chain (including warmup; default = 4000).}
#'   \item{warmup}{Number of warmup/burn-in samples per chain (default = 2000).}
#'   \item{thin}{Thinning parameter (default = 1).}
#'   \item{chain}{Number of chains (default = 4).}
#'   \item{sigma_v}{Proposal variance for "jump" terms (default = 0.1).}
#'   \item{bma}{Boolean indicating whether to include Bayesian model averaging step (default = TRUE).}
#' }
#' @param prior_params A list containing prior parameters:
#' \describe{
#'   \item{lambda_1}{Prior parameter for the number of terms in the model (default = 0.1).}
#'   \item{lambda_2}{Prior parameter for the number of knots in each spline (default = 1).}
#'   \item{a_0}{Shape parameter for inverse gamma prior on individual-level variance (default = 0.01).}
#'   \item{b_0}{Rate parameter for inverse gamma prior on individual-level variance (default = 0.01).}
#'   \item{degree}{Degree of B-splines (default = 3).}
#'   \item{k_max}{Maximum number of knots for each spline term (default = 9).}
#'   \item{w}{Window for proposing knot location changes (default = 1).}
#'   \item{sigma_B}{Prior normal variance for model coefficients (default = sqrt(20)).}
#' }
#'
#' @return An rjMCMC object with the following components:
#'   \describe{
#'     \item{accept_var}{Matrix indicating the acceptance of variable inclusion/removal for each iteration.}
#'     \item{accept_add_knot}{Matrix indicating acceptance of knot addition for each spline term across iterations.}
#'     \item{accept_remove_knot}{Matrix indicating acceptance of knot removal for each spline term across iterations.}
#'     \item{accept_move_knot}{Matrix indicating acceptance of knot movement for each spline term across iterations.}
#'     \item{splines_fitted}{List of matrices, one per spline interaction term, containing fitted spline values across iterations.}
#'     \item{binary_param}{Matrix containing posterior samples of binary variable parameters.}
#'     \item{inter_trt_param}{Matrix containing posterior samples of the treatment intercept and main effect.}
#'     \item{sigma_sq}{Matrix of posterior samples for the residual variance (sigma squared).}
#'     \item{vars_prop}{Matrix indicating the inclusion of variables across iterations (1 for included, 0 for excluded).}
#'     \item{vars_prop_summ}{Posterior inclusion probabilities for all candidate variables (spline and binary).}
#'     \item{k}{Matrix indicating the number of knots for each spline term across iterations.}
#'     \item{trt_eff_posterior}{Matrix of posterior treatment effect estimates, including spline effects.}
#'     \item{data_fit}{The original dataset passed to the function.}
#'     \item{candsplineinter}{A character vector indicating the spline interaction terms.}
#'     \item{candsplinevars}{A character vector of candidate spline variables.}
#'     \item{candbinaryvars}{A character vector of candidate binary variables.}
#'     \item{candinter}{A character vector of interaction terms with treatment (can include splines and binary variables).}
#'     \item{mcmc_specs}{The MCMC specifications used in the procedure.}
#'     \item{prior_params}{The prior parameters used in the procedure.}
#'   }
#'
#' @references
#'   Maleyeff, L., Golchi, S., Moodie, E. E. M., & Hudson, M. (2024) "An adaptive
#'   enrichment design using Bayesian model averaging for selection and threshold-identification
#'   of predictive variables" <doi:10.1093/biomtc/ujae141>
#'
#' @importFrom stats glm sigma quantile rgamma rnorm
#' @importFrom splines bs
#' @importFrom MASS mvrnorm
#' @importFrom utils setTxtProgressBar setTxtProgressBar
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
#' }
#' @export
rjMCMC <- function(data,
                   candsplinevars,
                   candbinaryvars,
                   candinter,
                   outcome,
                   factor_var,
                   mcmc_specs = NULL,
                   prior_params = NULL) {

  check_binary_columns_result = check_binary_columns(data,candbinaryvars)
  if (!(all(check_binary_columns_result$binary_check))) {
    stop(paste0(paste(check_binary_columns_result$non_binary_columns,collapse=", " ), " of candbinaryvars are not binary"))
  }

  default_mcmc_specs <- list(
    iter = 4000,
    warmup = NULL,  # Placeholder for later assignment
    thin = 1,
    chains = 4,
    sigma_v = 0.1,
    bma = TRUE
  )

  # If mcmc_specs is NULL, use default_mcmc_specs
  if (is.null(mcmc_specs)) {
    mcmc_specs <- default_mcmc_specs
  }

  # Set warmup to iter / 2 if it's NULL
  if (is.null(mcmc_specs$warmup)) {
    mcmc_specs$warmup <- mcmc_specs$iter / 2
  }

  # Check if each parameter is missing, if so, set it to the default value
  for (param in names(default_mcmc_specs)) {
    if (is.null(mcmc_specs[[param]])) {
      mcmc_specs[[param]] <- default_mcmc_specs[[param]]
      message(paste("Setting", param, "to default of", default_mcmc_specs[[param]]))
    }
  }

  if(mcmc_specs$warmup > mcmc_specs$iter) {
    stop("Error: warmup should not be greater than iter")
  }

  default_prior_params <- list(
    lambda_1 = 0.1,
    lambda_2 = 1,
    a_0 = 0.01,
    b_0 = 0.01,
    degree = 3,
    k_max = 9,
    w = 1,
    sigma_B = sqrt(20)
  )

  # If prior_params is NULL, use default_prior_params
  if (is.null(prior_params)) {
    prior_params <- default_prior_params
  }

  # Check if each parameter is missing, if so, set it to the default value
  for (param in names(default_prior_params)) {
    if (is.null(prior_params[[param]])) {
      prior_params[[param]] <- default_prior_params[[param]]
      message(paste("Setting", param, "to default of", default_prior_params[[param]]))
    }
  }

  if (!all(c(outcome, factor_var) %in% colnames(data))) {
    stop(paste("data must contain columns", outcome, "and", factor_var))
  }

  data$Y = data[,outcome]
  data$trt = data[,factor_var]

  if (!all(c(candsplinevars,candbinaryvars) %in% colnames(data))) {
    stop(paste0(paste(setdiff(c(candsplinevars,candbinaryvars),colnames(data)),collapse=", "), " are not columns of data - candsplinevars and candbinaryvars must correspond to columns in data"))
  }

  if(!all(candinter %in% c(candsplinevars,candbinaryvars))) {
    stop("candinter must be a subset of c(candsplinevars,candbinaryvars)")
  }

  if (!(all(data$trt %in% c(0, 1)) && length(unique(data$trt)) == 2)) {
    stop(paste(factor_var, "must be binary"))
  }

  B_per_chain = ceiling((mcmc_specs$iter - mcmc_specs$warmup) / mcmc_specs$thin)
  B = B_per_chain * mcmc_specs$chains
  burnin = mcmc_specs$warmup
  thin = mcmc_specs$thin
  bma = mcmc_specs$bma
  sigma_v = mcmc_specs$sigma_v

  degree = prior_params$degree
  k_max = prior_params$k_max
  lambda_1 = prior_params$lambda_1
  lambda_2 = prior_params$lambda_2
  a_0 = prior_params$a_0
  b_0 = prior_params$b_0
  sigma_B = prior_params$sigma_B
  w = prior_params$w

  iterations = B_per_chain * thin + burnin

  candsplinevars_ext = c(paste0(candsplinevars,"_main", recycle0 = T),
                         paste0(intersect(candsplinevars, candinter),"_inter", recycle0 = T))
  candbinaryvars_ext = c(paste0(candbinaryvars,"_main", recycle0 = T),
                         paste0(intersect(candbinaryvars, candinter),"_inter", recycle0 = T))
  n_cand_vars = length(candsplinevars_ext) + length(candbinaryvars_ext)
  candvars_ext = c(paste0(candsplinevars,"_main", recycle0 = T),
                   paste0(candbinaryvars,"_main", recycle0 = T),
                   paste0(intersect(candsplinevars, candinter),"_inter", recycle0 = T),
                   paste0(intersect(candbinaryvars, candinter),"_inter", recycle0 = T))
  candsplineinter = intersect(candsplinevars, candinter)
  candmain = c(candsplinevars,candbinaryvars)

  inter_trt_param_all = matrix(nrow = B, ncol = 2, dimnames = list(NULL,c("intercept",factor_var)))
  binary_param_all = matrix(nrow = B, ncol = length(candbinaryvars_ext), dimnames = list(NULL,candbinaryvars_ext))
  sigma_sq_all = array(dim = c(B, 1))
  splines_fitted_all = list()
  if (length(candsplinevars_ext) > 0) {
    for (l in 1:length(candsplinevars_ext)) {
      splines_fitted_all[[candsplinevars_ext[l]]] <- matrix(0, nrow = B, ncol = nrow(data))
    }
  }
  vars_prop_all = matrix(nrow = B, ncol = n_cand_vars, dimnames = list(NULL,candvars_ext))
  # Store k per spline
  k_all = matrix(nrow = B, ncol = length(candsplinevars_ext), dimnames = list(NULL, candsplinevars_ext))

  # Initialize acceptance matrices
  accept_var_all = matrix(nrow = B, ncol = 2, dimnames = list(NULL, c("add var", "remove var")))
  accept_add_knot_all = matrix(nrow = B, ncol = length(candsplinevars_ext), dimnames = list(NULL, candsplinevars_ext))
  accept_remove_knot_all = matrix(nrow = B, ncol = length(candsplinevars_ext), dimnames = list(NULL, candsplinevars_ext))
  accept_move_knot_all = matrix(nrow = B, ncol = length(candsplinevars_ext), dimnames = list(NULL, candsplinevars_ext))


  for (chain in 1:mcmc_specs$chains) {
    pb <- utils::txtProgressBar(min = 1, max = iterations, style = 3, width = 50, char = "=")
    message(paste0("Chain ", chain, ": "))
    inter_trt_param = array(dim = c(iterations,2))
    sigma_sq = array(dim = c(iterations,1))
    accept_var = matrix(nrow = iterations, ncol = 2)
    accept_add_knot = matrix(nrow = iterations, ncol = length(candsplinevars_ext))
    accept_remove_knot = matrix(nrow = iterations, ncol = length(candsplinevars_ext))
    accept_move_knot = matrix(nrow = iterations, ncol = length(candsplinevars_ext))

    k = matrix(nrow = iterations, ncol = length(candsplinevars_ext))

    vars_prop = matrix(nrow = iterations, ncol = n_cand_vars)
    vars_prop[1,] <- rep(1,n_cand_vars)

    spline_param = list()
    splines_fitted = list()

    curspline_ext = candsplinevars_ext
    curbinary_ext = candbinaryvars_ext
    curinter = candinter
    curmain = candmain
    curvars_ext = candvars_ext

    if (length(candbinaryvars_ext) > 0) {
      binary_param = array(dim = c(iterations,length(candbinaryvars_ext)))
      colnames(binary_param) = candbinaryvars_ext
      if (length(intersect(candbinaryvars, candinter))== 0) {
        binary_mod_mat = data[candbinaryvars]
      } else {
        binary_mod_mat = cbind(data[candbinaryvars],
                               data[intersect(candbinaryvars, candinter)]*data$trt)
      }

      colnames(binary_mod_mat) = candbinaryvars_ext
    } else {
      binary_param = as.matrix(rep(0,iterations))
      binary_mod_mat = 0
    }

    if (length(candsplinevars_ext) > 0) {
      knotscand = list()
      knotscur_idx = vector("list", length = length(candsplinevars_ext))
      names(knotscur_idx) = candsplinevars_ext

      for (i in 1:length(candsplinevars_ext)) {
        knotscand_i = stats::quantile(data[[sub("_[^_]+$", "", candsplinevars_ext[i])]],
                               seq(0,1,length.out=k_max+2))[-c(1,k_max+2)]

        knotscand[[candsplinevars_ext[i]]] = knotscand_i
      }

      spline_mod_mat = list()
      spline_mod_mat_raw = list()
      for (i in 1:length(candsplinevars_ext)) {
        mod_mat_i = splines::bs(data[[sub("_[^_]+$", "", candsplinevars_ext[i])]],
                       degree = degree,
                       knots = knotscand[[candsplinevars_ext[i]]][knotscur_idx[[candsplinevars_ext[i]]]],
                       intercept = F)
        spline_mod_mat_raw[[candsplinevars_ext[i]]] = mod_mat_i

        if (grepl("_inter", candsplinevars_ext[i])) {
          mod_mat_i = mod_mat_i*data$trt
        }

        spline_mod_mat[[candsplinevars_ext[i]]] = mod_mat_i
      }
      k[1,] <- unlist(lapply(knotscur_idx,length))
    } else {
      spline_mod_mat = 0
    }

    if (length(candbinaryvars_ext) == 0 & length(candsplinevars_ext) == 0) {
      mod_start <- stats::glm(Y ~ trt, data=data)
    } else if (length(candbinaryvars_ext) == 0) {
      mod_start <- stats::glm(Y ~ trt +
                         do.call(cbind, spline_mod_mat),data=data)
    } else if (length(candsplinevars_ext) == 0) {
      mod_start <- stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat), data=data)
      binary_param[1,] = coef(mod_start)[-c(1,2)]
    } else {
      mod_start <- stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat) +
                         do.call(cbind, spline_mod_mat),data=data)
      binary_param[1,] = coef(mod_start)[3:(length(candbinaryvars_ext)+2)]
    }


    inter_trt_param[1,] = c(coef(mod_start)[1], coef(mod_start)[2])
    sigma_sq[1] = stats::sigma(mod_start)^2
    if (length(candsplinevars_ext) > 0) {
      ncoef_perx = k[1,] + degree
      coefs = split(coef(mod_start)[-c(1:(length(candbinaryvars_ext)+2))],rep(1:(length(ncoef_perx)),ncoef_perx))

      spline_ols_param = list()
      for (i in 1:length(candsplinevars_ext)) {
        spline_ols_param[[candsplinevars_ext[i]]] = unlist(coefs[i])
      }

      spline_param[[1]] = spline_ols_param

      for (l in 1:length(candsplinevars_ext)) {
        splines_fitted[[candsplinevars_ext[l]]] = matrix(nrow=iterations,
                                                      ncol=nrow(data))
        splines_fitted[[candsplinevars_ext[l]]][1,] = splinesFitted(0,
                                                                    spline_param[[1]][[candsplinevars_ext[l]]],
                                                                    spline_mod_mat_raw[[candsplinevars_ext[l]]])

      }
    }


    for (i in 1:(iterations-1)) {
      utils::setTxtProgressBar(pb, i)
      if (length(candsplinevars_ext) > 0) {
        spline_param[[i+1]] = spline_param[[i]]
        k[i+1,] = k[i,]
        for (j in 1:length(candsplinevars_ext)) {
          j_name = candsplinevars_ext[j]
          if (j_name %in% curspline_ext) {
            knotscand_x = knotscand[[j_name]]
            knotscur_idx_x = knotscur_idx[[j_name]]
            knotscur_x = knotscand_x[knotscur_idx_x]


            # move knot
            move_knot = moveKnot(data$Y, data, degree, w, j_name, j, knotscand_x, knotscur_x,
                                 knotscur_idx_x, binary_param[i,],
                                 binary_mod_mat, spline_mod_mat, spline_mod_mat_raw,
                                 spline_ols_param[[j_name]],
                                 inter_trt_param[i,], spline_param[[i+1]], sigma_sq[i])

            knotscur_idx_x = move_knot$knotscur_idx_x
            spline_ols_param[[j_name]] = move_knot$spline_ols_param_x
            spline_mod_mat[[j_name]] = move_knot$spline_mod_mat_x
            spline_mod_mat_raw[[j_name]] = move_knot$spline_mod_mat_raw_x
            knotscur_x = move_knot$knotscur_x
            accept_move_knot[i,] = move_knot$accept

            v <- stats::rnorm(1,0,sigma_v)
            u_1 <- runif(1)

            # Add a knot
            if (u_1 < 0.5) {
              add_remove_knot = addKnot(data$Y, data, degree, v, sigma_v, sigma_B, lambda_2, j_name, j,
                                        k[i,j], k_max, candsplinevars,
                                        candbinaryvars, knotscand_x, knotscur_x,
                                        knotscur_idx_x,binary_param[i,],
                                        binary_mod_mat, spline_mod_mat, spline_mod_mat_raw,
                                        spline_ols_param[[j_name]],
                                        inter_trt_param[i,], spline_param[[i+1]], sigma_sq[i])
              accept_add_knot[i,] = add_remove_knot$accept
            }

            # remove knot
            if (u_1 >= 0.5) {
              add_remove_knot = removeKnot(data$Y, data, degree, v, sigma_v, sigma_B, lambda_2, j_name, j,
                                           k[i,j], k_max, candsplinevars,
                                           candbinaryvars, knotscand_x, knotscur_x,
                                           knotscur_idx_x, binary_param[i,],
                                           binary_mod_mat, spline_mod_mat, spline_mod_mat_raw,
                                           spline_ols_param[[j_name]],
                                           inter_trt_param[i,], spline_param[[i+1]], sigma_sq[i])
              accept_remove_knot[i,] = add_remove_knot$accept
            }

            # Update spline parameters
            spline_param[[i+1]][[j_name]] = add_remove_knot$spline_param_x

            # Update k parameters
            k[i+1,j] = add_remove_knot$k

            # Update knot indices
            knotscur_idx_x = add_remove_knot$knotscur_idx_x

            # Update knot values
            knotscur_x = add_remove_knot$knotscur_x

            # Update current ols spline parameters (used to scale coefficients)
            spline_ols_param[[j_name]] = add_remove_knot$spline_ols_param_x

            # Update spline model matrix
            spline_mod_mat[[j_name]] = add_remove_knot$spline_mod_mat_x
            spline_mod_mat_raw[[j_name]] = add_remove_knot$spline_mod_mat_raw_x

            # Update all splines
            curspline_idx = which(candsplinevars_ext %in% curvars_ext)
            curspline_idx = curspline_idx[curspline_idx!=j]
            if (is.list(binary_mod_mat)) {
              X_rest = cbind(1,data$trt,do.call(cbind, spline_mod_mat[curspline_idx]),as.matrix(binary_mod_mat))
              beta_rest = c(inter_trt_param[i,],do.call(c, spline_param[[i+1]][curspline_idx]),binary_param[i,])
            } else {
              X_rest = cbind(1,data$trt,do.call(cbind, spline_mod_mat[curspline_idx]))
              beta_rest = c(inter_trt_param[i,],do.call(c, spline_param[[i+1]][curspline_idx]))
            }

            XtX_block_spline = t(spline_mod_mat[[j_name]]) %*% spline_mod_mat[[j_name]]

            Sigma_block_inv = XtX_block_spline / sigma_sq[i] + 1/sigma_B^2
            Sigma_block = solve(Sigma_block_inv)
            mean_block = Sigma_block %*% (t(spline_mod_mat[[j_name]]) %*% (data$Y - X_rest %*% beta_rest)) / sigma_sq[i]

            # Sample from the conditional posterior
            spline_param[[i+1]][[j_name]] <- MASS::mvrnorm(1, mu = mean_block, Sigma = Sigma_block)
            knotscur_idx[[j_name]] = knotscur_idx_x
          }
        }
      } else {
        spline_param[[i+1]] = 0
      }

      if (bma) {
        internotinmodel = c(setdiff(candinter,curinter))
        eligibletoremove = c(paste_(setdiff(curmain, curinter), "_main"),
                             paste_(curinter, "_inter"))

        eligibletoadd = c(paste_(c(setdiff(candmain,curmain)),"_main"),
                          paste_(internotinmodel[internotinmodel %in% curmain],"_inter"))

        n_cur_vars = length(curmain) + length(curinter)
        u_2 = runif(1)
        if (u_2 < 0.5) {
          add_remove_var = addVar(data$Y,
                                  data$trt,
                                  n_cand_vars,
                                  n_cur_vars,
                                  sigma_v,
                                  sigma_B,
                                  lambda_1,
                                  eligibletoremove,
                                  eligibletoadd,
                                  curvars_ext,
                                  curmain,
                                  curinter,
                                  curbinary_ext,
                                  curspline_ext,
                                  candbinaryvars_ext,
                                  inter_trt_param[i,],
                                  binary_param[i,],
                                  binary_mod_mat,
                                  spline_param[[i+1]],
                                  spline_mod_mat,
                                  sigma_sq[i]
                                  )
          accept_var[i,1] = add_remove_var$accept
        }

        if (u_2 >= 0.5) {
          add_remove_var = removeVar(data$Y,
                                     data$trt,
                                     n_cand_vars,
                                     n_cur_vars,
                                     sigma_v,
                                     sigma_B,
                                     lambda_1,
                                     eligibletoremove,
                                     eligibletoadd,
                                     curvars_ext,
                                     curmain,
                                     curinter,
                                     curbinary_ext,
                                     curspline_ext,
                                     candbinaryvars_ext,
                                     inter_trt_param[i,],
                                     binary_param[i,],
                                     binary_mod_mat,
                                     spline_param[[i+1]],
                                     spline_mod_mat,
                                     sigma_sq[i]
                                     )
          accept_var[i,2] = add_remove_var$accept

        }

        spline_param[[i+1]] = add_remove_var$spline_param_i
        inter_trt_param[i+1,] = add_remove_var$inter_trt_param_i
        binary_param[i+1,] = add_remove_var$binary_param_i
        curmain = add_remove_var$curmain
        curinter = add_remove_var$curinter
        curvars_ext = add_remove_var$curvars_ext
        curbinary_ext = add_remove_var$curbinary_ext
        curspline_ext = add_remove_var$curspline_ext
        vars_prop[i+1,] = as.numeric(candvars_ext %in% curvars_ext)
      } else {
        inter_trt_param[i+1,] = inter_trt_param[i,]
        binary_param[i+1,] =  binary_param[i,]
        spline_param[[i+1]] = spline_param[[i+1]]
        vars_prop[i+1,] = vars_prop[i,]
      }



      for (j in 1:2) {
        if (j==1) {
          values = matrix(rep(1,nrow(data)),ncol=1)
        } else if (j==2) {
          values = matrix(data$trt,ncol=1)
        }
        # V_beta_j = 1/ (X_j^t X_j / sigma2 + 1 / sigma2_beta_j)
        V_beta_j <- 1 / ((t(values) %*% values) / sigma_sq[i] + 1 / sigma_B^2)
        # beta_j_hat = V_beta_j * X_j^t (Y-X_[-j]hat(beta)_[-j]) / sigma2[i]
        beta_j_hat <- V_beta_j * t(values) %*% (data$Y - computeFittedValues(data$trt,
                                                                             inter_trt_param[i+1,],
                                                                             binary_param[i+1,],
                                                                             binary_mod_mat,
                                                                             spline_param[[i+1]],
                                                                             spline_mod_mat) + values %*% inter_trt_param[i+1, j]) / sigma_sq[i]
        inter_trt_param[i+1, j] <- stats::rnorm(1, beta_j_hat, sqrt(V_beta_j))
      }

      if (length(curbinary_ext) > 0) {
        for (j in 1:ncol(binary_param)) {
          if (candbinaryvars_ext[j] %in% curbinary_ext) {
            values = matrix(binary_mod_mat[,j],ncol=1)
            # V_beta_j = 1/ (X_j^t X_j / sigma2 + 1 / sigma2_beta_j)
            V_beta_j <- 1 / ((t(values) %*% values) / sigma_sq[i] + 1 / sigma_B^2)
            # beta_j_hat = V_beta_j * X_j^t (Y-X_[-j]hat(beta)_[-j]) / sigma2[i]
            beta_j_hat <- V_beta_j * t(values) %*% (data$Y - computeFittedValues(data$trt,
                                                                                 inter_trt_param[i+1,],
                                                                                 binary_param[i+1,],
                                                                                 binary_mod_mat,
                                                                                 spline_param[[i+1]],
                                                                                 spline_mod_mat) + values %*% binary_param[i+1,j]) / sigma_sq[i]
            binary_param[i+1, j] <- stats::rnorm(1, beta_j_hat, sqrt(V_beta_j))
          }
        }
      }

      # Update sigma^2
      resids = data$Y - computeFittedValues(data$trt,
                                           inter_trt_param[i+1,],
                                           binary_param[i+1,],
                                           binary_mod_mat,
                                           spline_param[[i+1]],
                                           spline_mod_mat

      )
      sigma_sq[i+1] <- 1/stats::rgamma(1, shape=nrow(data)/2+a_0, rate=sum(resids^2)/2+b_0)

      #
      # Update treatment effect, use the "raw" model.matrix which contains the
      # spline matrix for the interaction effect, before being multiplied by the
      # treatment indicator
      if (length(candsplinevars_ext) > 0) {
        for (l in 1:length(candsplinevars_ext)) {
          splines_fitted[[candsplinevars_ext[l]]][i+1,] = splinesFitted(0,
                                                                   spline_param[[i+1]][[candsplinevars_ext[l]]],
                                                                   spline_mod_mat_raw[[candsplinevars_ext[l]]])

        }
      }
    }

    close(pb)

    ### Post-processing: Remove burn-in and apply thinning
    index_to_keep <- seq(from = burnin + 1, to = iterations, by = thin)

    ### Append results of each chain (thinned and post-burnin) to global variables
    inter_trt_param_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- inter_trt_param[index_to_keep, ]
    sigma_sq_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- sigma_sq[index_to_keep, ]
    vars_prop_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- vars_prop[index_to_keep, ]

    accept_var_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- accept_var[index_to_keep, ]
    accept_add_knot_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- accept_add_knot[index_to_keep, ]
    accept_remove_knot_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- accept_remove_knot[index_to_keep, ]
    accept_move_knot_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- accept_move_knot[index_to_keep, ]

    binary_param_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- binary_param[index_to_keep, ]

    if (length(candsplinevars_ext) > 0) {
      for (l in 1:length(candsplinevars_ext)) {
        spline_name <- candsplinevars_ext[l]
        splines_fitted_all[[spline_name]][((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- splines_fitted[[spline_name]][index_to_keep, ]
      }
    }

    k_all[((chain - 1) * B_per_chain + 1):(chain * B_per_chain), ] <- k[index_to_keep, ]

  }
############################################ MCMC diagnostics ###############################################
  rename_vars <- function(var_name) {
    var_name <- gsub("_main$", "", var_name)  # Remove "_main"
    var_name <- gsub("_inter$", paste0(":",factor_var), var_name)  # Replace "_inter" with ":trt"
    return(var_name)
  }
  colnames(binary_param_all) <- unname(sapply(colnames(binary_param_all), rename_vars))
  colnames(accept_remove_knot_all) <- unname(sapply(colnames(accept_remove_knot_all), rename_vars))
  colnames(accept_add_knot_all) <- unname(sapply(colnames(accept_add_knot_all), rename_vars))
  colnames(accept_move_knot_all) <- unname(sapply(colnames(accept_move_knot_all), rename_vars))
  colnames(vars_prop_all) <- unname(sapply(colnames(vars_prop_all), rename_vars))
  colnames(k_all) <- unname(sapply(colnames(k_all), rename_vars))

  names(splines_fitted_all) <- sapply(names(splines_fitted_all), rename_vars)

  trt_eff_posterior = t(matrix(1,nrow=nrow(data)) %*% t(matrix(inter_trt_param_all[,2])))
  if (length(candsplineinter) > 0) {
    for (m in 1:length(candsplineinter)) {
      trt_eff_posterior = trt_eff_posterior + splines_fitted_all[[paste0(candsplineinter[m],paste0(":",factor_var))]]
    }
  }

  candbinaryinter = intersect(candbinaryvars,candinter)
  if (length(candbinaryinter) > 0) {
    trt_eff_posterior = trt_eff_posterior +
      binary_param_all[,paste0(candbinaryinter,paste0(":",factor_var))] %*% t(as.matrix(data[,candbinaryinter]))
  }


  vars_prop_summ = colMeans(vars_prop_all)
  names(vars_prop_summ) = colnames(vars_prop_all)

  obj <- list(
    accept_var = accept_var_all,
    accept_add_knot = accept_add_knot_all,
    accept_remove_knot = accept_remove_knot_all,
    accept_move_knot = accept_move_knot_all,
    splines_fitted = splines_fitted_all,
    binary_param = binary_param_all,
    fixed_param = inter_trt_param_all,
    sigma_sq = sigma_sq_all,
    vars_prop = vars_prop_all,
    vars_prop_summ = vars_prop_summ,
    k = k_all,
    trt_eff_posterior = trt_eff_posterior,
    data_fit = data,
    candsplineinter = candsplineinter,
    candsplinevars = candsplinevars,
    candbinaryvars = candbinaryvars,
    candinter = candinter,
    outcome = outcome,
    factor_var = factor_var,
    mcmc_specs = mcmc_specs,
    prior_params = prior_params)
  class(obj) <- c("rjMCMC", class(obj))
  return(obj)
}

