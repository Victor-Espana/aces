#' @title Estimate Coefficients for Adaptive Constrained Enveloping Splines (ACES)
#'
#' @description
#'
#' This function estimates a vector of coefficients for the Adaptive Constrained Enveloping Splines (ACES) model. These coefficients enforce the desired shape properties, such as monotonicity and concavity to fit the data within the ACES framework.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier.
#'
#' @param B
#' A \code{matrix} of linear basis functions derived from input variables.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_scores
#' A \code{matrix} containing DEA-VRS efficiency scores, calculated using an output-oriented radial model. For models with multiple outputs, each column corresponds to the scores for one specific output.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity.
#'
#' @return
#'
#' A \code{vector} containing the estimated coefficients for the ACES model.

estimate_coefficients <- function (
    model_type,
    B,
    y_obs,
    dea_scores,
    it_list,
    Bp_list,
    shape
    ) {

  if (model_type == "envelopment") {

    coefs <- estimate_coefficients_envelopment (
      B = B,
      y_obs = y_obs,
      dea_scores = dea_scores,
      it_list = it_list,
      Bp_list = Bp_list,
      shape = shape
    )

  } else {

    coefs <- estimate_coefficients_stochastic (
      B = B,
      y_obs = y_obs,
      it_list = it_list,
      Bp_list = Bp_list,
      shape = shape
    )

  }

  return(coefs)

}

#' @title Estimate Coefficients in Adaptive Constrained Enveloping Splines: Envelopment Version
#'
#' @description
#'
#' This function estimates a vector of coefficients for the Adaptive Constrained Enveloping Splines (ACES) model in its envelopment version. It solves a mathematical programming problem to determine the coefficients that best fit the observed data while respecting certain shape constraints such as monotonicity and/or concavity. The model aims to envelop the data points, particularly when estimating production functions or efficiency frontiers.
#'
#' @param B
#' A \code{matrix} of linear basis functions derived from input variables.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_scores
#' A \code{matrix} containing DEA-VRS efficiency scores, calculated using an output-oriented radial model. For models with multiple outputs, each column corresponds to the scores for one specific output.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity.
#'
#' @return
#'
#' A \code{vector} containing the estimated coefficients for the ACES model.

estimate_coefficients_envelopment <- function (
    B,
    y_obs,
    dea_scores,
    it_list,
    Bp_list,
    shape
    ) {

  # monotonicity
  mono <- shape[["mono"]]

  # concavity
  conc <- shape[["conc"]]

  # dea weights
  deaw <- as.matrix(1 / dea_scores)

  # sample size
  N <- nrow(B)

  # number of basis functions (inputs)
  p <- ncol(B)

  # number of outputs
  nY <- ncol(y_obs)

  # number of inputs
  nX <- length(Bp_list)

  # structure of coefficients
  coefs <- matrix(NA, nrow = p, ncol = nY)

  for (out in 1:nY) {

    # select an output
    y_ind <- y_obs[, out]

    # ==================================================== #
    # vars: c(coef_0, coef_1, ..., coef_P, e_1, ... , e_n) #
    # ==================================================== #

    objVal <- c(rep(0, p), deaw[, out])

    # ==================================================== #
    # A: envelopment + concavity + monotonicity            #
    # ==================================================== #

    # envelopment: y_hat - e = y
    EMat <- cbind(B, diag(rep(- 1, N), N))

    # matrix of coefficients
    Amat <- rbind(EMat)

    # generate non-decreasing monotonicity matrix
    if (mono) {

      MMat <- monotonocity_matrix (
        it_list = it_list,
        Bp_list = Bp_list,
        N = N
      )

      # add MMat to Amat
      Amat <- rbind(Amat, MMat)

    }

    # generate concavity matrix
    if (conc) {

      CMat <- concavity_matrix(
        Bp_list = Bp_list,
        N = N
      )

      # add CMat to Amat
      Amat <- rbind(Amat, CMat)

    }

    # ==================================================== #
    # b: envelopment + concavity + monotonicity            #
    # ==================================================== #

    bvec <- c(y_ind, rep(0, nrow(Amat) - N))

    # ==================================================== #
    # directions of inequalities                           #
    # ==================================================== #

    dirs <- c(rep("==", N), rep(">=", nrow(Amat) - N))

    # ==================================================== #
    # lower and upper bounds                               #
    # ==================================================== #

    bnds <- list(lower = list(ind = 1:p, val = rep(- Inf, p)))

    # ==================================================== #
    # solution of the optimization problem                 #
    # ==================================================== #

    sols <- Rglpk_solve_LP (
      obj = objVal,
      mat = Amat,
      dir = dirs,
      rhs = bvec,
      bounds = bnds
    )

    if (sols[["status"]] == 0) {

      coefs[, out] <- sols$solution[1:p]

    } else {

      coefs[, out] <- rep(0, p)

    }

  }

  return(coefs)

}

#' @title Estimate Coefficients in Adaptive Constrained Enveloping Splines Fitting: Stochastic Version
#'
#' @description
#'
#' This function estimates a vector of coefficients for the Adaptive Constrained Enveloping Splines (ACES) model in its stochastic version. It solves a quadratic programming problem to determine the coefficients that best fit the observed data while respecting certain shape constraints such as monotonicity and/or concavity.
#'
#' @param B
#' A \code{matrix} of linear basis functions derived from input variables.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity.
#'
#' @importFrom quadprog solve.QP
#'
#' @return
#'
#' A \code{vector} containing the estimated coefficients for the ACES model.

estimate_coefficients_stochastic <- function (
    B,
    y_obs,
    it_list,
    Bp_list,
    shape
    ) {

  # monotonicity
  mono <- shape[["mono"]]

  # concavity
  conc <- shape[["conc"]]

  # sample size
  N <- nrow(B)

  # number of basis functions
  p <- ncol(B)

  # number of outputs
  nY <- ncol(y_obs)

  # number of inputs
  nX <- length(Bp_list)

  # structure of coefficients
  coefs <- matrix(NA, nrow = p, ncol = nY)

  for (out in 1:nY) {

    # select an output
    y_ind <- y_obs[, out]

    # ==================================================== #
    # vars: c(coef_0, coef_1, ..., coef_P, e_1, ... , e_n) #
    # ==================================================== #

    # ================= #
    # D: Quadratic part #
    # ================= #

    # identity matrix
    Dmat <- diag(rep(1, p + N))

    # near 0 for alpha
    Dmat[1:p, 1:p] <- diag(rep(1e-12, p), p)
    # diag(Dmat[(p + 1):(N + p), (p + 1):(N + p)]) <- 1 / y_ind

    # ============== #
    # d: Linear part #
    # ============== #

    dvec <- rep(0, p + N)

    # ==================================================== #
    # A: fitting + concavity + monotonicity                #
    # ==================================================== #

    # fitting matrix: y_hat - e = y
    FMat <- cbind(B, diag(rep(- 1, N), N))
    Amat <- FMat

    # generate monotonicity matrix
    if (mono) {

      MMat <- monotonocity_matrix (
        it_list = it_list,
        Bp_list = Bp_list,
        N = N
      )

      Amat <- rbind(Amat, MMat)

    }

    # generate concavity matrix
    if (conc) {

      CMat <- concavity_matrix (
        Bp_list = Bp_list,
        N = N
      )

      Amat <- rbind(Amat, CMat)

    }

    # ==================================================== #
    # b: fitting + concavity + monotonicity                #
    # ==================================================== #

    if (ppto_val) {
      bvec <- c(y_ind, 1, rep(0, nrow(Amat) - N - 1))

    } else {
      bvec <- c(y_ind, rep(0, nrow(Amat) - N))

    }

    # ==================================================== #
    # solution of the optimization problem                 #
    # ==================================================== #

    sols <- solve.QP (
      Dmat = Dmat, dvec = dvec,
      Amat = t(Amat), bvec = bvec,
      meq = N
    )

    coefs[, out] <- sols$solution[1:p]

  }

  return(coefs)
}
