#' @title Estimate Coefficients for Adaptive Constrained Enveloping Splines (ACES) Fitting
#'
#' @description
#'
#' This function solves a Mathematical Programming problem to obtain a vector of coefficients that impose the required shaped on the Adaptive Constrained Enveloping Splines estimator.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function estimates.
#'
#' @param B
#' A \code{matrix} of linear basis functions.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_eff
#' An indicator vector with ones for efficient DMUs and zeros for inefficient DMUs.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @return
#'
#' A \code{vector} of estimated coefficients.

estimate_coefficients <- function (
    model_type,
    B,
    y_obs,
    dea_eff,
    it_list,
    Bp_list,
    shape
    ) {

  if (model_type == "env") {

    coefs <- estim_coefs_env (
      B = B,
      y_obs = y_obs,
      dea_eff = dea_eff,
      it_list = it_list,
      Bp_list = Bp_list,
      shape = shape
    )

  } else {

    coefs <- estim_coefs_sto (
      B = B,
      y_obs = y_obs,
      it_list = it_list,
      Bp_list = Bp_list,
      shape = shape
    )

  }

  return(coefs)
}

#' @title Estimate Coefficients in Adaptive Constrained Enveloping Splines Fitting: Envelopment Version
#'
#' @description
#'
#' This function solves a Linear Programming problem to obtain a vector of coefficients that impose the required shaped on the Adaptive Constrained Enveloping Splines estimator in the envelopment version.
#'
#' @param B
#' A \code{matrix} of linear basis functions.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_eff
#' An indicator vector with ones for efficient DMUs and zeros for inefficient DMUs.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @return
#'
#' A \code{vector} of estimated coefficients.

estim_coefs_env <- function (
    B,
    y_obs,
    dea_eff,
    it_list,
    Bp_list,
    shape
    ) {

  # monotonicity
  mono = shape[["mono"]]

  # concavity
  conc = shape[["conc"]]

  # pass through the origin
  ptto = shape[["ptto"]]

  # ===================== #
  # called from smoothing #
  # ===================== #

  if (ptto != FALSE && (is.null(it_list) && is.null(Bp_list))) {

    # f(0) = 0 vector of coefficients
    origin_vec <- c(B[nrow(B), ], rep(0, nrow(B) - 1))

    # B matrix
    B <- B[1:(nrow(B) - 1), ]

    # y_obs
    y_obs <- as.matrix(y_obs[1:(length(y_obs) - 1)])

  }

  # sample size
  N <- nrow(B)

  # number of efficient DMUs
  N_eff <- length(dea_eff)

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

    objVal <- c(rep(0, p), rep(1, N))

    # ==================================================== #
    # A: envelopment + concavity + monotonicity            #
    # ==================================================== #

    # envelopment: y_hat + e = y
    EMat <- cbind(B, diag(rep(- 1, N), N))
    EMat <- EMat[dea_eff, ]

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

    bvec <- c(y_ind[dea_eff], rep(0, nrow(Amat) - N_eff))

    # ==================================================== #
    # directions of inequalities                           #
    # ==================================================== #

    dirs <- c(rep("==", N_eff), rep(">=", nrow(Amat) - N_eff))

    # ==================================================== #
    # Prediction: f(0) = 0                                 #
    # ==================================================== #

    if (ptto != FALSE) {

      # f(0) = 0
      if (!is.null(it_list) && !is.null(Bp_list)) {

        origin_vec <- predict_origin (
          it_list = it_list,
          Bp_list = Bp_list,
          n_bf = p,
          N = N
        )

        # add origin_vec to Amat
        Amat <- rbind(Amat, origin_vec)

        # update right-hand terms
        bvec <- c(bvec, if (ptto == "0") 0 else 1)

        # update vector of directions
        dirs <- c(rep("==", N_eff), rep(">=", nrow(Amat) - N_eff - 1), "==")

      } else {

        # add f(0) = 0 to Amat
        Amat <- rbind(Amat, origin_vec)

        # update right-hand terms
        bvec <- c(bvec, if (ptto == "0") 0 else 1)

        # add f(0) = 0 to dirs
        dirs <- c(rep("==", N_eff), rep(">=", nrow(Amat) - N_eff - 1), "==")
      }

    }

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

    coefs[, out] <- sols$solution[1:p]

  }

  return(coefs)
}

#' @title Estimate Coefficients in Adaptive Constrained Enveloping Splines Fitting: Stochastic Version
#'
#' @description
#'
#' This function solves a Quadratic Programming problem to obtain a vector of coefficients that impose the required shaped on the Adaptive Constrained Enveloping Splines estimator in the stochastic version
#'
#' @param B
#' A \code{matrix} of linear basis functions.
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
#' @param mono
#' A \code{logical} value indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#'
#' @param conc
#' A \code{logical} value indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @param ptto
#' A \code{logical} value indicating whether the estimator should satisfy f(0) = 0.
#'
#' @importFrom quadprog solve.QP
#'
#' @return
#'
#' A \code{vector} of estimated coefficients.

estim_coefs_sto <- function (
    B,
    y_obs,
    it_list,
    Bp_list,
    mono,
    conc,
    ptto
    ) {

  # called from smoothing
  if (ptto != FALSE && (is.null(it_list) && is.null(Bp_list))) {

    # f(0) = 0 vector of coefficients
    origin_vec <- c(B[nrow(B), ], rep(0, nrow(B) - 1))

    # B matrix
    B <- B[1:(nrow(B) - 1), ]

    # y_obs
    y_obs <- as.matrix(y_obs[1:(length(y_obs) - 1)])

  }

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

    # ==================================================== #
    # Prediction: f(0) = 0                                 #
    # ==================================================== #

    if (ptto != FALSE) {

      # f(0) = 0
      if (!is.null(it_list) && !is.null(Bp_list)) {

        origin_vec <- predict_origin (
          it_list = it_list,
          Bp_list = Bp_list,
          n_bf = p,
          N = N
        )

        # add origin_vec to Amat
        Amat <- rbind(Amat, origin_vec)

      } else {

        # add f(0) = 0 to Amat
        Amat <- rbind(Amat, origin_vec)

      }
    }

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

    if (ptto == "1") {
      bvec <- c(y_ind, 1, rep(0, nrow(Amat) - N - 1))

    } else {
      bvec <- c(y_ind, rep(0, nrow(Amat) - N))

    }

    # ==================================================== #
    # solution of the optimization problem                 #
    # ==================================================== #

    if (ptto != FALSE) {

      sols <- solve.QP (
        Dmat = Dmat, dvec = dvec,
        Amat = t(Amat), bvec = bvec,
        meq = N + 1
      )

    } else {

      sols <- solve.QP (
        Dmat = Dmat, dvec = dvec,
        Amat = t(Amat), bvec = bvec,
        meq = N
      )

    }

    coefs[, out] <- sols$solution[1:p]
  }

  return(coefs)
}
