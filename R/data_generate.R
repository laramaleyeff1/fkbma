#' A synthetic dataset with continuous and binary covariates and a binary treatment variable (`trt`).
#'
#' This dataset was generated using the formula:
#' \eqn{Y = 2 * Z_1 + 2 * X_1 + 2 * Z_1 * trt + \cos(X_1 \cdot 2 \pi) * trt + \epsilon},
#' where \eqn{\epsilon \sim N(0, 0.1)}.
#'
#' @docType data
#' @name simulated_data
#' @rdname simulated_data
#' @usage data(simulated_data)
#' @format A data frame with 1000 rows and 8 variables:
#' \describe{
#'   \item{X_1}{A continuous variable sampled from \code{U(0, 1)}.}
#'   \item{Z_1}{A binary variable sampled from \code{Bernoulli(0.35)}.}
#'   \item{Z_2}{A binary variable sampled from \code{Bernoulli(0.5)}.}
#'   \item{Z_3}{A binary variable sampled from \code{Bernoulli(0.65)}.}
#'   \item{Z_4}{A binary variable sampled from \code{Bernoulli(0.2)}.}
#'   \item{Z_5}{A binary variable sampled from \code{Bernoulli(0.35)}.}
#'   \item{trt}{A binary treatment variable sampled from \code{Bernoulli(0.5)}.}
#'   \item{Y}{The outcome variable calculated using the formula above.}
#' }
#' @examples
#' data(simulated_data)
#' head(simulated_data)
#' @keywords datasets
"simulated_data"

#' @noRd
#' @keywords internal
.generateData <- function(seed = 1234L) {
  set.seed(seed)
  n = 1000
  simulated_data = data.frame(X_1 = stats::runif(n,0,1),
                    Z_1 = stats::rbinom(n,1,0.35),
                    Z_2 = stats::rbinom(n,1,0.5),
                    Z_3 = stats::rbinom(n,1,0.65),
                    Z_4 = stats::rbinom(n,1,0.2),
                    Z_5 = stats::rbinom(n,1,0.35),
                    trt = stats::rbinom(n,1,0.5))

  simulated_data$Y = 2 * simulated_data$Z_1 + 2* simulated_data$X_1 + 2 * simulated_data$Z_1 * simulated_data$trt + cos(simulated_data$X_1*2*pi) * simulated_data$trt + stats::rnorm(n, 0, 0.1)
  return(simulated_data)
}


