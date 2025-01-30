# Helper Functions for FK-BMA

#' Check if Candidate Binary Variables Are Binary
#'
#' This internal function verifies whether the candidate binary variables in the given data are truly binary
#' (i.e., take values 0 or 1).
#'
#' @noRd
#' @param data A data frame containing the variables.
#' @param candbinaryvars A vector of candidate binary variable names.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{binary_check}{A named logical vector indicating whether each candidate binary variable is binary.}
#'   \item{non_binary_columns}{A vector of variable names that are not binary.}
#' }
#' @importFrom stats setNames
#' @keywords internal
check_binary_columns <- function(data, candbinaryvars) {
  # Initialize a named logical vector to store the results
  result <- stats::setNames(logical(length(candbinaryvars)), candbinaryvars)
  non_binary_cols <- c()

  # Check each column
  for (col in candbinaryvars) {
    if (col %in% colnames(data)) {
      is_binary <- all(data[[col]] %in% c(0, 1)) && length(unique(data[[col]])) == 2
      result[col] <- is_binary
      if (!is_binary) {
        non_binary_cols <- c(non_binary_cols, col)
      }
    } else {
      result[col] <- FALSE
      non_binary_cols <- c(non_binary_cols, col)
    }
  }
  return(list(binary_check = result, non_binary_columns = non_binary_cols))
}

#' Sample One Element from a Vector
#'
#' This internal function samples one element from a vector. If the vector has only one element, it returns that element.
#' @noRd
#' @param x A vector to sample from.
#'
#' @return A sampled element or the single element if the vector has only one element.
#' @keywords internal
sample_ <- function(x) {
  if (length(x) == 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}

#' Compute Fitted Values for the Model
#'
#' This internal function computes the fitted values of a regression model based on the treatment effect, binary parameters,
#' and spline parameters.
#'
#' @noRd
#' @param trt A vector indicating treatment status.
#' @param inter_trt_param_ A vector of treatment parameters (intercept and treatment effect).
#' @param binary_param_ A vector of binary variable coefficients.
#' @param binary_mod_mat_ A matrix of binary covariates.
#' @param spline_param_minus_j A list of spline parameters excluding the j-th variable.
#' @param spline_mat_minus_j A list of spline matrices excluding the j-th variable.
#'
#' @return A vector of fitted values.
#' @keywords internal
computeFittedValues <- function(trt,
                                inter_trt_param_,
                             binary_param_,
                             binary_mod_mat_,
                             spline_param_minus_j,
                             spline_mat_minus_j) {

  binary_eff = 0
  if (is.list(binary_mod_mat_)) {
    binary_eff = do.call(cbind,binary_mod_mat_) %*% binary_param_
  }

  spline_eff = 0
  if (is.list(spline_mat_minus_j)) {
    spline_eff = do.call(cbind, spline_mat_minus_j) %*% unlist(spline_param_minus_j)
  }

  return(cbind(1,trt) %*% inter_trt_param_ + binary_eff + spline_eff)
}

#' Compute Fitted Splines
#'
#' This internal function computes the fitted spline values given the treatment and spline matrix.
#'
#' @noRd
#' @param trt A vector of treatment indicators.
#' @param spline A vector of spline coefficients.
#' @param mat A matrix for the spline basis functions.
#'
#' @return A vector of fitted spline values.
#' @keywords internal
splinesFitted <- function(trt, spline, mat) {
  return(trt + mat %*% spline)
}

#' Append Suffix to Variable Names
#'
#' This internal function appends a suffix to the given variable names.
#'
#' @noRd
#' @param set A vector of variable names.
#' @param suffix A string to append to the variable names.
#'
#' @return A vector of variable names with the appended suffix.
#' @keywords internal
paste_ <- function(set, suffix) {
  if (length(set) == 0) {
    return(NULL)
  } else {
    return(paste0(set,suffix))
  }

}

#' Custom Log-Likelihood Calculation
#'
#' This internal function calculates the log-likelihood of the data given the model parameters.
#'
#' @noRd
#' @param Y The outcome variable.
#' @param trt The treatment variable.
#' @param inter_trt_param_ A vector of treatment effect parameters.
#' @param binary_param_ A vector of binary variable parameters.
#' @param binary_mod_mat_ A matrix of binary covariates.
#' @param spline_param_ A list of spline parameters.
#' @param spline_mat_ A list of spline matrices.
#' @param var The variance of the model.
#'
#' @return The log-likelihood of the data given the model.
#' @importFrom stats dnorm
#' @keywords internal
logLikelihoodCustom <- function(Y,
                                trt,
                                inter_trt_param_,
                                binary_param_,
                                binary_mod_mat_,
                                spline_param_,
                                spline_mat_,
                                var) {
  pred = computeFittedValues(trt,
                             inter_trt_param_,
                             binary_param_,
                             binary_mod_mat_,
                             spline_param_,
                             spline_mat_)

  singlelikelihoods = stats::dnorm(Y, mean = pred, sd = sqrt(var), log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

#' Move a Spline Knot
#'
#' This internal function moves a knot in a spline and evaluates whether the move improves the model fit.
#'
#' @noRd
#' @param Y The outcome variable.
#' @param data A data frame containing the variables.
#' @param degree The degree of the spline.
#' @param w The window size for proposing knot changes.
#' @param ... Additional arguments related to the current spline and binary parameters.
#'
#' @return A list containing the updated knot indices, spline coefficients, and acceptance indicator.
#' @importFrom splines bs
#' @importFrom stats glm
#' @keywords internal
moveKnot <- function(Y,
                     data,
                     degree,
                     w,
                     varcur,
                     j,
                     knotscand_x,
                     knotscur_x,
                     knotscur_idx_x,
                     binary_param,
                     binary_mod_mat,
                     spline_mod_mat,
                     spline_mod_mat_raw,
                     spline_ols_param_x,
                     inter_trt_param_i,
                     spline_param_i,
                     sigma_sq_i) {
  spline_mod_mat_x = spline_mod_mat[[varcur]]
  spline_mod_mat_raw_x = spline_mod_mat_raw[[varcur]]

  if (length(knotscur_idx_x) > 0) {
    knot_tomove = knotscand_x[sample_(knotscur_idx_x)]
    knots_avail = setdiff(knotscand_x, knotscur_x)
    window_propose = knots_avail[knots_avail<=(knot_tomove+w) &
                                   knots_avail>=(knot_tomove-w)]
  } else {
    window_propose = NULL
  }

  if (length(window_propose) > 0) {
    knot_newvalue = sample_(window_propose)

    knots_propose = sort(setdiff(union(knotscur_x, knot_newvalue),knot_tomove))
    knotscur_idx_propose = sort(setdiff(union(which(knotscand_x == knot_newvalue),
                                              knotscur_idx_x),
                                        which(knotscand_x == knot_tomove)))

    knots_avail_reverse = setdiff(knotscand_x, knots_propose)
    window_reverse = knots_avail_reverse[knots_avail_reverse<=(knot_newvalue+w) &
                                           knots_avail_reverse>=(knot_newvalue-w)]

    offset = computeFittedValues(data$trt,
                                 inter_trt_param_i,
                              binary_param,
                              binary_mod_mat,
                              spline_param_i[-j],
                              spline_mod_mat[-j])

    new_spline = splines::bs(data[[sub("_[^_]+$", "", varcur)]],
                    degree = degree,
                    knots = knots_propose)

    new_spline_raw = new_spline
    if (length(grep("inter",varcur))>0) {
      new_spline = new_spline*data$trt
    }

    mod_propose <- stats::glm(Y ~ -1 + new_spline,
                       data = data,
                       offset = offset)

    spline_mod_mat_propose = spline_mod_mat
    spline_mod_mat_propose[[varcur]] = new_spline

    log_prob = logLikelihoodCustom(Y,
                                   data$trt,
                                   inter_trt_param_i,
                                   binary_param,
                                   binary_mod_mat,
                                   spline_param_i,
                                   spline_mod_mat_propose,
                                   sigma_sq_i) -
      logLikelihoodCustom(Y,
                          data$trt,
                          inter_trt_param_i,
                          binary_param,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sigma_sq_i) +
      log(length(window_propose)) -
      log(length(window_reverse))


    gamma <- stats::runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of which knots
      knotscur_idx_x = knotscur_idx_propose
      knotscur_x = knots_propose
      # keep track of OLS model coefs
      spline_ols_param_x = stats::coef(mod_propose)
      # keep track of model matrices
      spline_mod_mat_x = new_spline
      spline_mod_mat_raw_x = new_spline_raw

      accept = 1
    } else {
      accept = 0
    }
  } else {
    accept = 0
  }

  return(list(knotscur_idx_x = knotscur_idx_x,
              knotscur_x = knotscur_x,
              spline_ols_param_x = spline_ols_param_x,
              spline_mod_mat_x = spline_mod_mat_x,
              spline_mod_mat_raw_x = spline_mod_mat_raw_x,
              accept = accept))
}


#' Add a Spline Knot
#'
#' This internal function adds a knot to a spline and evaluates whether the addition improves the model fit.
#'
#' @noRd
#' @param Y The outcome variable.
#' @param data A data frame containing the variables.
#' @param degree The degree of the spline.
#' @param ... Additional arguments related to the current spline and binary parameters.
#'
#' @return A list containing the updated knot indices, spline coefficients, and acceptance indicator.
#' @importFrom splines bs
#' @importFrom stats glm coef runif
#' @keywords internal
addKnot <- function(Y,
                    data,
                    degree,
                    v,
                    sigma_v,
                    sigma_B,
                    lambda_2,
                    varcur,
                    j,
                    k_ij,
                    k_max,
                    candsplinevars,
                    candbinaryvars,
                    knotscand_x,
                    knotscur_x,
                    knotscur_idx_x,
                    binary_param,
                    binary_mod_mat,
                    spline_mod_mat,
                    spline_mod_mat_raw,
                    spline_ols_param_x,
                    inter_trt_param_i,
                    spline_param_i,
                    sigma_sq_i) {
  spline_mod_mat_x = spline_mod_mat[[varcur]]
  spline_mod_mat_raw_x = spline_mod_mat_raw[[varcur]]
  if (k_ij < k_max) {
    index_propose = sample_(setdiff(1:k_max,knotscur_idx_x))
    knot_propose = knotscand_x[index_propose]
    knots_propose = sort(union(knotscur_x, knot_propose))
    k_interval_propose = which(knots_propose == knot_propose)
    knotscur_idx_propose = sort(union(index_propose, knotscur_idx_x))
    k_cur_propose = k_ij + 1


    offset = computeFittedValues(data$trt,
                                 inter_trt_param_i,
                              binary_param,
                              binary_mod_mat,
                              spline_param_i[-j],
                              spline_mod_mat[-j])


    new_spline = splines::bs(data[[sub("_[^_]+$", "", varcur)]],
                    degree = degree,
                    knots = knots_propose)

    new_spline_raw = new_spline
    if (length(grep("inter",varcur))>0) {
      new_spline = new_spline*data$trt
    }

    mod_propose <- stats::glm(Y ~ -1 + new_spline,
                       data = data,
                       offset = offset)

    spline_propose = stats::coef(mod_propose)
    spline_propose[-(k_interval_propose+1)] = spline_propose[-(k_interval_propose+1)] +
      spline_param_i[[j]] -
      spline_ols_param_x

    spline_propose[(k_interval_propose+1)] = spline_propose[(k_interval_propose+1)] + v

    splines_propose = spline_param_i
    splines_propose[[varcur]] = spline_propose

    spline_mod_mat_propose = spline_mod_mat
    spline_mod_mat_propose[[varcur]] = new_spline

    log_prob = logLikelihoodCustom(Y,
                                   data$trt,
                                   inter_trt_param_i,
                                   binary_param,
                                   binary_mod_mat,
                                   splines_propose,
                                   spline_mod_mat_propose,
                                   sigma_sq_i) -
      logLikelihoodCustom(Y,
                          data$trt,
                          inter_trt_param_i,
                          binary_param,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sigma_sq_i) +
      log((lambda_2*sigma_v)/(k_cur_propose*sigma_B)) -
      (v^2/2*sigma_v^2) + (1/(2*sigma_B^2))*(t(spline_param_i[[varcur]]) %*% spline_param_i[[varcur]] -
                                               t(spline_propose) %*% (spline_propose))

    gamma <- stats::runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coeffients
      spline_param_x = spline_propose
      # keep track of number of knots
      k = k_cur_propose
      # keep track of which knots
      knotscur_idx_x = knotscur_idx_propose
      # keep track of OLS model coeffieints
      spline_ols_param_x = stats::coef(mod_propose)
      # keep track of model matrices
      spline_mod_mat_x = new_spline
      spline_mod_mat_raw_x = new_spline_raw
      accept = 1

    } else {
      spline_param_x = spline_param_i[[varcur]]
      k = k_ij
      accept = 0
    }
  } else {
    spline_param_x = spline_param_i[[varcur]]
    k = k_ij
    accept = 0
  }

  return(list(spline_param_x = spline_param_x,
              k = k,
              knotscur_idx_x = knotscur_idx_x,
              knotscur_x = knotscur_x,
              spline_ols_param_x = spline_ols_param_x,
              spline_mod_mat_x = spline_mod_mat_x,
              spline_mod_mat_raw_x = spline_mod_mat_raw_x,
              accept = accept))

}

#' Remove a Spline Knot
#'
#' This internal function removes a knot from a spline and evaluates whether the removal improves the model fit.
#'
#' @noRd
#' @param Y The outcome variable.
#' @param data A data frame containing the variables.
#' @param degree The degree of the spline.
#' @param ... Additional arguments related to the current spline and binary parameters.
#'
#' @return A list containing the updated knot indices, spline coefficients, and acceptance indicator.
#' @importFrom splines bs
#' @importFrom stats glm coef runif
#' @keywords internal
removeKnot <- function(Y,
                       data,
                       degree,
                    v,
                    sigma_v,
                    sigma_B,
                    lambda_2,
                    varcur,
                    j,
                    k_ij,
                    k_max,
                    candsplinevars,
                    candbinaryvars,
                    knotscand_x,
                    knotscur_x,
                    knotscur_idx_x,
                    binary_param,
                    binary_mod_mat,
                    spline_mod_mat,
                    spline_mod_mat_raw,
                    spline_ols_param_x,
                    inter_trt_param_i,
                    spline_param_i,
                    sigma_sq_i) {
  spline_mod_mat_x = spline_mod_mat[[varcur]]
  spline_mod_mat_raw_x = spline_mod_mat_raw[[varcur]]

  if (k_ij > 0) {
    index_propose = sample_(knotscur_idx_x)
    knot_propose = knotscand_x[index_propose]
    knots_propose = setdiff(knotscur_x, knot_propose)

    k_interval_propose = which(knotscur_x == knot_propose)
    knotscur_idx_propose = setdiff(knotscur_idx_x, index_propose)
    k_cur_propose = k_ij - 1

    offset = computeFittedValues(data$trt,
                                 inter_trt_param_i,
                              binary_param,
                              binary_mod_mat,
                              spline_param_i[-j],
                              spline_mod_mat[-j])

    new_spline = splines::bs(data[[sub("_[^_]+$", "", varcur)]],
                    degree = degree,
                    knots = knots_propose)

    new_spline_raw = new_spline
    if (length(grep("inter",varcur))>0) {
      new_spline = new_spline*data$trt
    }

    mod_propose <- stats::glm(Y ~ -1 + new_spline,
                       data = data,
                       offset = offset)
    spline_propose = stats::coef(mod_propose)
    spline_propose = spline_propose +
      spline_param_i[[j]][-(k_interval_propose+1)] -
      spline_ols_param_x[-(k_interval_propose+1)]

    splines_propose = spline_param_i
    splines_propose[[varcur]] = spline_propose

    spline_mod_mat_propose = spline_mod_mat
    spline_mod_mat_propose[[varcur]] = new_spline

    log_prob = logLikelihoodCustom(Y,
                                   data$trt,
                                   inter_trt_param_i,
                                   binary_param,
                                   binary_mod_mat,
                                   splines_propose,
                                   spline_mod_mat_propose,
                                   sigma_sq_i) -
      logLikelihoodCustom(Y,
                          data$trt,
                          inter_trt_param_i,
                          binary_param,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sigma_sq_i) +
      log(sigma_B*k_ij/(lambda_2*sigma_v)) -
      (v^2/2*sigma_v^2) + (1/(2*sigma_B^2)) *
      (t(spline_param_i[[varcur]]) %*% spline_param_i[[varcur]] -
         t(spline_propose) %*% (spline_propose))

    gamma <- stats::runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coefs
      spline_param_x = spline_propose
      # keep track of number of knots
      k = k_cur_propose
      # keep track of which knots
      knotscur_idx_x = knotscur_idx_propose
      # keep track of OLS model coefs
      spline_ols_param_x = stats::coef(mod_propose)
      # keep track of model matrices
      spline_mod_mat_x = new_spline
      spline_mod_mat_raw_x = new_spline_raw
      accept = 1

    } else {
      spline_param_x = spline_param_i[[varcur]]
      k = k_ij
      accept = 0
    }
  } else {
    spline_param_x = spline_param_i[[varcur]]
    k = k_ij
    accept = 0
  }

  return(list(spline_param_x = spline_param_x,
              k = k,
              knotscur_idx_x = knotscur_idx_x,
              knotscur_x = knotscur_x,
              spline_ols_param_x = spline_ols_param_x,
              spline_mod_mat_x = spline_mod_mat_x,
              spline_mod_mat_raw_x = spline_mod_mat_raw_x,
              accept = accept))

}

#' Add a New Variable to the Model
#'
#' This internal function adds a new variable (either binary or spline) to the model and evaluates whether the addition improves the model fit.
#'
#' @noRd
#' @param Y The outcome variable.
#' @param trt The treatment variable.
#' @param ... Additional arguments related to the current binary and spline parameters.
#'
#' @return A list containing the updated variable parameters and acceptance indicator.
#' @importFrom stats dnorm glm coef runif
#' @keywords internal
addVar <- function(Y,
                   trt,
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
                   inter_trt_param_i,
                   binary_param_i,
                   binary_mod_mat,
                   spline_param_i,
                   spline_mod_mat,
                   sigma_sq_i) {
  if (length(eligibletoadd) > 0) {
    to_add = sample_(eligibletoadd)
    spline_param_propose = spline_param_i
    binary_param_propose = binary_param_i
    inter_trt_param_propose = inter_trt_param_i

    curbinary_propose = curbinary_ext
    curspline_propose = curspline_ext

    if (to_add %in% candbinaryvars_ext) {
      curbinary_propose =  c(curbinary_ext, to_add)
      index_add = which(to_add == curbinary_propose)

      dim_x = 1
      v <- rnorm(dim_x,0,sigma_v)
      binary_param_propose[to_add] = v
      if (length(curspline_ext) == 0) {
        mod_propose = stats::glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_propose]))

        if (length(curbinary_ext) == 0) {
          mod_prev = stats::glm(Y ~ trt)
        } else {
          mod_prev = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]))
          binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
            stats::coef(mod_prev)[-c(1,2)] +
            stats::coef(mod_propose)[setdiff(3:(length(curbinary_propose)+2),2+index_add)]
        }
      } else {
        mod_propose = stats::glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_propose]) +
                           do.call(cbind, spline_mod_mat[curspline_ext]))

        if (length(curbinary_ext) == 0) {
          mod_prev = stats::glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_ext]))
        } else {
          mod_prev = stats::glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_ext]) +
                             do.call(cbind, spline_mod_mat[curspline_ext]))
          binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
            stats::coef(mod_prev)[3:(length(curbinary_ext)+2)] +
            stats::coef(mod_propose)[setdiff(3:(length(curbinary_ext)+2),2+index_add)]
        }

        spline_ols_prev = convertSplineCoefs(stats::coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                             spline_mod_mat,
                                             curspline_ext)
        spline_ols_propose = convertSplineCoefs(stats::coef(mod_propose)[-c(1:(length(curbinary_propose)+2))],
                                                spline_mod_mat,
                                                curspline_ext)

        for (m in 1:length(curspline_ext)) {
          spline_param_propose[[curspline_ext[m]]] =
            spline_param_i[[curspline_ext[m]]] -
            spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
        }
      }

      binary_param_propose[to_add] = v + stats::coef(mod_propose)[2+index_add]
    } else {
      curspline_propose =  c(curspline_ext, to_add)
      dim_x = length(spline_param_i[[to_add]])
      index_add = which(to_add == curspline_propose)

      v <- rnorm(dim_x,0,sigma_v)
      if (length(curbinary_ext) == 0) {
        mod_propose = stats::glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_propose]))

        spline_ols_propose = convertSplineCoefs(stats::coef(mod_propose)[-c(1,2)],
                                             spline_mod_mat,
                                             curspline_propose)

        if (length(curspline_ext) == 0) {
          mod_prev = stats::glm(Y ~ trt)
        } else {
          mod_prev = stats::glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_ext]))
          spline_ols_prev = convertSplineCoefs(stats::coef(mod_prev)[-c(1,2)],
                                                  spline_mod_mat,
                                                  curspline_ext)

          for (m in 1:length(curspline_ext)) {
            spline_param_propose[[curspline_ext[m]]] =
              spline_param_i[[curspline_ext[m]]] -
              spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
          }
        }

      } else {
        mod_propose = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]) +
                           do.call(cbind, spline_mod_mat[curspline_propose]))

        spline_ols_propose = convertSplineCoefs(stats::coef(mod_propose)[-c(1:(length(curbinary_ext)+2))],
                                             spline_mod_mat,
                                             curspline_propose)

        if (length(curspline_ext) == 0) {
          mod_prev = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]))
        } else {
          mod_prev = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]) +
                          do.call(cbind, spline_mod_mat[curspline_ext]))
          spline_ols_prev = convertSplineCoefs(stats::coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                                  spline_mod_mat,
                                                  curspline_ext)

          for (m in 1:length(curspline_ext)) {
            spline_param_propose[[curspline_ext[m]]] =
              spline_param_i[[curspline_ext[m]]] -
              spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
          }
        }
        binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
          stats::coef(mod_prev)[3:(length(curbinary_ext)+2)] +
          stats::coef(mod_propose)[3:(length(curbinary_ext)+2)]
      }
      spline_param_propose[[to_add]] = v + spline_ols_propose[[to_add]]
    }

    inter_trt_param_propose = inter_trt_param_i - stats::coef(mod_prev)[c(1,2)] + stats::coef(mod_propose)[c(1,2)]
    log_prob = logLikelihoodCustom(Y,
                                   trt,
                                   inter_trt_param_propose,
                                   binary_param_propose,
                                   binary_mod_mat,
                                   spline_param_propose,
                                   spline_mod_mat,
                                   sigma_sq_i) -
      logLikelihoodCustom(Y,
                          trt,
                          inter_trt_param_i,
                          binary_param_i,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sigma_sq_i) +
      sum(stats::dnorm(c(inter_trt_param_propose, binary_param_propose, unlist(spline_param_propose)), mean = 0, sd = sigma_B, log = T)) -
      sum(stats::dnorm(c(inter_trt_param_i, binary_param_i, unlist(spline_param_i)), mean = 0,
                sd = sigma_B, log = T)) +
      log(lambda_1/(n_cand_vars-n_cur_vars)) +
      log(length(eligibletoadd)) - log(length(eligibletoremove)+1) -
      sum(stats::dnorm(v,0,sigma_v,log=T))

    curmain_propose = curmain
    curinter_propose = curinter
    if (length(grep("main", to_add)>0)) {
      curmain_propose = c(sub("_[^_]+$", "", to_add), curmain)
    } else {
      curinter_propose = c(sub("_[^_]+$", "", to_add), curinter)
    }

    gamma <- stats::runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coeffients
      spline_param_i = spline_param_propose
      binary_param_i = binary_param_propose
      inter_trt_param_i = inter_trt_param_propose
      curmain = curmain_propose
      curinter = curinter_propose
      curvars_ext = c(curbinary_propose, curspline_propose)
      curbinary_ext = curbinary_propose
      curspline_ext = curspline_propose
      accept = 1
    } else {
      accept = 0
    }

  } else {
    accept = 0
  }
  return(list(spline_param_i = spline_param_i,
              inter_trt_param_i = inter_trt_param_i,
              binary_param_i = binary_param_i,
              curmain = curmain,
              curinter = curinter,
              curvars_ext = curvars_ext,
              curbinary_ext = curbinary_ext,
              curspline_ext = curspline_ext,
              accept = accept))
}


#' Convert Spline Coefficients
#'
#' This internal function converts the OLS coefficients of splines to the format required by the model.
#'
#' @noRd
#' @param spline_coefs The spline coefficients from the OLS model.
#' @param spline_mat The matrix of spline basis functions.
#' @param vars The names of the spline variables.
#'
#' @return A list of converted spline coefficients.
#' @keywords internal
convertSplineCoefs <- function(spline_coefs,
                               spline_mat,
                               vars) {
  ncoef_perx = unlist(lapply(spline_mat[vars],ncol))

  coef_prop = split(spline_coefs,
                    rep(1:(length(ncoef_perx)),ncoef_perx))
  spline_ols = list()
  for (m in 1:length(vars)) {
    spline_ols[[vars[m]]] = unlist(coef_prop[m])
  }
  return(spline_ols)
}

#' Remove a Variable from the Model
#'
#' This internal function removes a variable (either binary or spline) from the model and evaluates whether the removal improves the model fit.
#'
#' @noRd
#' @param Y The outcome variable.
#' @param trt The treatment variable.
#' @param ... Additional arguments related to the current binary and spline parameters.
#'
#' @return A list containing the updated variable parameters and acceptance indicator.
#' @importFrom stats dnorm glm runif coef
#' @keywords internal
removeVar <- function(Y,
                      trt,
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
                     inter_trt_param_i,
                     binary_param_i,
                     binary_mod_mat,
                     spline_param_i,
                     spline_mod_mat,
                     sigma_sq_i) {
  if (length(eligibletoremove) > 0) {
    to_remove = sample_(eligibletoremove)
    spline_param_propose = spline_param_i

    binary_param_propose = binary_param_i
    inter_trt_param_propose = inter_trt_param_i
    curbinary_propose = curbinary_ext
    curspline_propose = curspline_ext

    if (to_remove %in% candbinaryvars_ext) {
      curbinary_propose =  setdiff(curbinary_ext, to_remove)
      dim_x = 1
      v <- rnorm(dim_x,0,sigma_v)
      binary_param_propose[to_remove] = 0

      if (length(curspline_ext) == 0) {
        binary_mod_mat_prev = binary_mod_mat[curbinary_ext]
        index_remove = which(to_remove == curbinary_ext)
        mod_prev = stats::glm(Y ~ trt + do.call(cbind, binary_mod_mat_prev))

        if (length(curbinary_propose) == 0) {
          mod_propose = stats::glm(Y ~ trt)
        } else {
          mod_propose = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_propose]))
          binary_param_propose[curbinary_propose] = binary_param_i[curbinary_propose] -
            stats::coef(mod_prev)[setdiff(3:(length(curbinary_ext)+2),2+index_remove)] +
            stats::coef(mod_propose)[-c(1,2)]
        }
      } else {
        mod_prev = stats::glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_ext]) + do.call(cbind, spline_mod_mat[curspline_ext]))

        if (length(curbinary_propose) == 0) {
          mod_propose = stats::glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_ext]))
        } else {
          mod_propose = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_propose]) +
                             do.call(cbind, spline_mod_mat[curspline_ext]))
          index_remove = which(to_remove == curbinary_ext)
          binary_param_propose[curbinary_propose] = binary_param_i[curbinary_propose] -
            stats::coef(mod_prev)[setdiff(3:(length(curbinary_ext)+2),2+index_remove)] +
            stats::coef(mod_propose)[3:(length(curbinary_propose)+2)]
        }

        spline_ols_prev = convertSplineCoefs(stats::coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                         spline_mod_mat,
                                         curspline_ext)
        spline_ols_propose = convertSplineCoefs(stats::coef(mod_propose)[-c(1:(length(curbinary_propose)+2))],
                                         spline_mod_mat,
                                         curspline_ext)

        for (m in 1:length(curspline_ext)) {
          spline_param_propose[[curspline_ext[m]]] =
            spline_param_i[[curspline_ext[m]]] -
            spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
        }
      }

      binary_param_propose[to_remove] = 0
    } else {
      dim_x = length(spline_param_i[[to_remove]])
      curspline_propose =  setdiff(curspline_ext, to_remove)

      v <- rnorm(dim_x,0,sigma_v)

      if (length(curbinary_ext) == 0) {
        mod_mat_prev = do.call(cbind, spline_mod_mat[curspline_ext])
        mod_prev = stats::glm(Y ~ trt + mod_mat_prev)

        spline_ols_prev = convertSplineCoefs(stats::coef(mod_prev)[-c(1,2)],
                                             spline_mod_mat,
                                             curspline_ext)

        if (length(curspline_propose) == 0) {
          mod_propose = stats::glm(Y ~ trt)
        } else {
          mod_mat_propose = do.call(cbind, spline_mod_mat[curspline_propose])
          mod_propose = stats::glm(Y ~ trt + mod_mat_propose)
          spline_ols_propose = convertSplineCoefs(stats::coef(mod_propose)[-c(1,2)],
                                                  spline_mod_mat,
                                                  curspline_propose)

          for (m in 1:length(curspline_propose)) {
            spline_param_propose[[curspline_propose[m]]] =
              spline_param_i[[curspline_propose[m]]] -
              spline_ols_prev[[curspline_propose[m]]] + spline_ols_propose[[curspline_propose[m]]]
          }
        }

      } else {
        mod_mat_prev = do.call(cbind, spline_mod_mat[curspline_ext])

        mod_prev = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]) + mod_mat_prev)

        spline_ols_prev = convertSplineCoefs(stats::coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                             spline_mod_mat,
                                             curspline_ext)

        if (length(curspline_propose) == 0) {
          mod_propose = stats::glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]))
        } else {
          mod_propose = stats::glm(Y ~ trt +  do.call(cbind,binary_mod_mat[curbinary_ext]) +
                             do.call(cbind, spline_mod_mat[curspline_propose]))
          spline_ols_propose = convertSplineCoefs(stats::coef(mod_propose)[-c(1:(length(curbinary_ext)+2))],
                                                  spline_mod_mat,
                                                  curspline_propose)
          for (m in 1:length(curspline_propose)) {
            spline_param_propose[[curspline_propose[m]]] =
              spline_param_i[[curspline_propose[m]]] -
              spline_ols_prev[[curspline_propose[m]]] + spline_ols_propose[[curspline_propose[m]]]
          }


        }
        binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
          stats::coef(mod_prev)[3:(length(curbinary_ext)+2)] +
          stats::coef(mod_propose)[3:(length(curbinary_ext)+2)]
      }
      spline_param_propose[[to_remove]] = rep(0,dim_x)
    }

    inter_trt_param_propose = inter_trt_param_i - stats::coef(mod_prev)[c(1,2)] + stats::coef(mod_propose)[c(1,2)]
    #inter_trt_param_propose = inter_trt_param_i
    log_prob = logLikelihoodCustom(Y,
                                   trt,
                                   inter_trt_param_propose,
                                   binary_param_propose,
                                   binary_mod_mat,
                                   spline_param_propose,
                                   spline_mod_mat,
                                   sigma_sq_i) -
      logLikelihoodCustom(Y,
                          trt,
                          inter_trt_param_i,
                          binary_param_i,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sigma_sq_i) +
      sum(stats::dnorm(c(inter_trt_param_propose, binary_param_propose, unlist(spline_param_propose)), mean = 0, sd = sigma_B, log = T)) -
      sum(stats::dnorm(c(inter_trt_param_i, binary_param_i, unlist(spline_param_i)), mean = 0,
                sd = sigma_B, log = T)) +
      log((n_cand_vars-n_cur_vars + 1)/lambda_1) +
      log(length(eligibletoremove)) - log(length(eligibletoadd)+1) +
      sum(stats::dnorm(v,0,sigma_v,log=T))

    curmain_propose = curmain
    curinter_propose = curinter
    if (length(grep("main", to_remove)>0)) {
      curmain_propose = setdiff(curmain, sub("_[^_]+$", "", to_remove))
    } else {
      curinter_propose =  setdiff(curinter, sub("_[^_]+$", "", to_remove))
    }

    gamma <- stats::runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coefficients
      spline_param_i = spline_param_propose
      binary_param_i = binary_param_propose
      inter_trt_param_i = inter_trt_param_propose
      curmain = curmain_propose
      curinter = curinter_propose
      curvars_ext = c(curbinary_propose, curspline_propose)
      curbinary_ext = curbinary_propose
      curspline_ext = curspline_propose
      accept = 1
    } else {
      accept = 0
    }

  } else {
    accept = 0
  }
  return(list(spline_param_i = spline_param_i,
              inter_trt_param_i = inter_trt_param_i,
              binary_param_i = binary_param_i,
              curmain = curmain,
              curinter = curinter,
              curvars_ext = curvars_ext,
              curbinary_ext = curbinary_ext,
              curspline_ext = curspline_ext,
              accept = accept))
}

