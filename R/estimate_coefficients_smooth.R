#' @title Estimate Coefficients during Non-decreasing and Concave Smoothing Procedure
#'
#' @description
#'
#' This function solves a Linear Programming problem to obtain a set of coefficients that impose concavity and non-decreasing monotonicity in the Smoothing Adaptive Constrained Enveloping Splines estimator.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function will estimate.
#'
#' @param B
#' A \code{matrix} of cubic or quintic basis functions.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_eff
#' An indicator vector with 1s for efficient DMUs and 0s for inefficient DMUs.
#'
#' @param n_pair
#' An \code{integer} specifying the number of coefficients associated with the paired basis functions.
#'
#' @param n_lsub
#'  An \code{integer} specifying the number of coefficients associated with the unpaired left-side basis functions.
#'
#' @param monotonicity
#' A \code{logical} value indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#'
#' @param concavity
#' A \code{logical} value indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @param origin
#' A \code{logical} value indicating whether the estimator should satisfy f(0) = 0.
#'
#' @return
#'
#' A \code{vector} with the estimated coefficients.

estimate_coefficients_smoothed <- function (
    model_type,
    B,
    y_obs,
    dea_eff,
    n_pair,
    n_lsub,
    monotonicity,
    concavity,
    origin
    ) {

  if (model_type == "env") {
    coefs <- ecs_envelopment (
      B = B,
      y_obs = y_obs,
      dea_eff = dea_eff,
      n_pair = n_pair,
      n_lsub = n_lsub,
      monotonicity = monotonicity,
      concavity = concavity,
      origin = origin
    )

  } else {
    coefs <- ecs_stochastic (
      B = B,
      y_obs = y_obs,
      n_pair = n_pair,
      n_lsub = n_lsub,
      monotonicity = monotonicity,
      concavity = concavity,
      origin = origin
    )
  }

  return(coefs)
}

#' @title Estimate Coefficients in Smoothing Adaptive Constrained Enveloping Splines Fitting: Envelopment Version and Additive Error
#'
#' @description
#'
#' This function solves a Linear Programming problem to obtain a vector of coefficients that impose the required shaped on the Smoothing Adaptive Constrained Enveloping Splines estimator in the envelopment version, assuming an additive error structure.
#'
#' @param B
#' A \code{matrix} of cubic or quintic basis functions.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_eff
#' An indicator vector with 1s for efficient DMUs and 0s for inefficient DMUs.
#'
#' @param n_pair
#' An \code{integer} specifying the number of coefficients associated with the paired basis functions.
#'
#' @param n_lsub
#'  An \code{integer} specifying the number of coefficients associated with the unpaired left-side basis functions.
#'
#' @param monotonicity
#' A \code{logical} value indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#'
#' @param concavity
#' A \code{logical} value indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @param origin
#' A \code{logical} value indicating whether the estimator should satisfy f(0) = 0.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return
#'
#' A \code{vector} of estimated coefficients.

ecs_envelopment <- function (
    B,
    y_obs,
    dea_eff,
    n_pair,
    n_lsub,
    monotonicity,
    concavity,
    origin
    ) {

  if (origin != "FALSE") {
    # f(0) = 0 vector of coefficients
    Ovec <- c(B[nrow(B), ], rep(0, nrow(B) - 1))

    # B matrix
    B <- B[1:(nrow(B) - 1), ]

    # y_obs
    y_obs <- as.matrix(y_obs[1:(length(y_obs) - 1)])
  }

  # sample size
  N <- nrow(B)

  # number of efficient DMUs
  N_eff <- length(dea_eff)

  # number of basis functions: 1 + n_pair + n_lsub
  p <- ncol(B)

  # number of pairs
  pairs <- n_pair / 2

  # number of outputs
  nY <- ncol(y_obs)

  # initialize coefficients structure
  coefs <- matrix(NA, nrow = p, ncol = nY)

  for (out in 1:nY) {

    # vars: c(tau_0, alpha_1, beta_2, ... , alpha_(H-1), beta_H, w_1, ..., w_R, e_1, ... , e_n)
    objVal <- c(rep(0, p), rep(1, N))

    # select a variable
    y_ind <- y_obs[, out]

    # ========================================= #
    # A: envelopment + concavity + monotonicity #
    # ========================================= #

    # equality constraints: y_hat - e = y
    EMat <- cbind(B, diag(rep(- 1, N), N))
    EMat <- EMat[dea_eff, ]

    # matrix of coefficients
    Amat <- rbind(EMat)

    # concavity: paired basis functions
    if (concavity) {
      CMat <- matrix(0, nrow = pairs, ncol = p + N)

      if (pairs >= 1) {
        for (i in 1:pairs) {
          CMat[i, ] <- c (
            0,
            rep(c(0, 0), i - 1),
            - 1,
            - 1,
            rep(c(0, 0), pairs - i),
            rep(0, n_lsub + N)
            )
        }
      }

      # add concavity matrix to the envelopment matrix
      Amat <- rbind(Amat, CMat)
    }

    # non-decreasing monotonicity: paired basis functions
    if (monotonicity) {
      MMat <- cbind (
        rep(0, n_pair),
        diag(c(1, - 1), n_pair),
        matrix(0, n_pair, n_lsub),
        matrix(0, n_pair, N)
        )

      # add non-decreasing monotonicity matrix to the envelopment matrix
      Amat <- rbind(Amat, MMat)
    }

    # concavity & increasing monotonicity: left-side unpaired basis functions
    if (monotonicity | concavity) {
      MCMat <- cbind (
        matrix(0, n_lsub, 1 + n_pair),
        diag(- 1, n_lsub),
        matrix(0, n_lsub, N)
        )

      # add concavity & increasing monotonicity [LSBF] to envelopment matrix
      Amat <- rbind(Amat, MCMat)
    }

    # ==================================================== #
    # b: envelopment + concavity + monotonicity            #
    # ==================================================== #

    bvec <- c(y_ind[dea_eff], rep(0, nrow(Amat) - N_eff))

    # ==================================================== #
    # Directions of inequalities                           #
    # ==================================================== #

    dirs <- c(rep("==", N_eff), rep(">=", nrow(Amat) - N_eff))

    # ==================================================== #
    # Prediction: f(0) = 0                                 #
    # ==================================================== #

    if (origin != "FALSE") {
      # add f(0) = 0 to Amat
      Amat <- rbind(Amat, Ovec)

      # update right-hand terms
      if (origin == "0") {
        bvec <- c(bvec, 0)

      } else {
        bvec <- c(bvec, 1)
      }

      # add f(0) = 0 to dirs
      dirs <- c(rep("==", N_eff), rep(">=", nrow(Amat) - N_eff - 1), "==")
    }

    # ==================================================== #
    # Lower and upper bounds                               #
    # ==================================================== #

    bnds <- list (
      lower = list (
        ind = 1L:p,
        val = rep(- Inf, p)
        )
      )

    # ==================================================== #
    # Solution of the optimization problem                 #
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

#' @title Estimate Coefficients in Smoothing Adaptive Constrained Enveloping Splines Fitting: Stochastic Version and Additive Error
#'
#' @description
#'
#' This function solves a Linear Programming problem to obtain a vector of coefficients that impose the required shaped on the Smoothing Adaptive Constrained Enveloping Splines estimator in the stochastic version, assuming an additive error structure.
#'
#' @param B
#' A \code{matrix} of cubic or quintic basis functions.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param n_pair
#' An \code{integer} specifying the number of coefficients associated with the paired basis functions.
#'
#' @param n_lsub
#'  An \code{integer} specifying the number of coefficients associated with the unpaired left-side basis functions.
#'
#' @param monotonicity
#' A \code{logical} value indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#'
#' @param concavity
#' A \code{logical} value indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @param origin
#' A \code{logical} value indicating whether the estimator should satisfy f(0) = 0.
#'
#' @importFrom quadprog solve.QP
#'
#' @return
#'
#' A \code{vector} of estimated coefficients.

ecs_stochastic <- function (
    B,
    y_obs,
    dea_eff,
    n_pair,
    n_lsub,
    monotonicity,
    concavity,
    origin
    ) {

  if (origin != "FALSE") {
    # f(0) = 0 vector of coefficients
    Ovec <- c(B[nrow(B), ], rep(0, nrow(B) - 1))

    # B matrix
    B <- B[1:(nrow(B) - 1), ]

    # y_obs
    y_obs <- as.matrix(y_obs[1:(length(y_obs) - 1)])
  }

  # sample size
  N <- nrow(B)

  # number of basis functions: 1 + n_pair + n_lsub
  p <- ncol(B)

  # number of pairs
  pairs <- n_pair / 2

  # number of outputs
  nY <- ncol(y_obs)

  # initialize coefficients structure
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

    # ============== #
    # d: Linear part #
    # ============== #

    dvec <- rep(0, p + N)

    # ========================================= #
    # A: envelopment + concavity + monotonicity #
    # ========================================= #

    # fitting matrix: y_hat - e = y
    FMat <- cbind(B, diag(rep(- 1, N), N))
    Amat <- FMat

    # f(0) = 0
    if (origin != "FALSE") {
      Amat <- rbind(Amat, Ovec)
    }

    # concavity: paired basis functions
    if (concavity) {
      CMat <- matrix(0, nrow = pairs, ncol = p + N)

      if (pairs >= 1) {
        for (i in 1:pairs) {
          CMat[i, ] <- c (
            0,
            rep(c(0, 0), i - 1),
            - 1,
            - 1,
            rep(c(0, 0), pairs - i),
            rep(0, n_lsub + N)
          )
        }
      }

      # add concavity matrix to the envelopment matrix
      Amat <- rbind(Amat, CMat)
    }

    # non-decreasing monotonicity: paired basis functions
    if (monotonicity) {
      MMat <- cbind (
        rep(0, n_pair),
        diag(c(1, - 1), n_pair),
        matrix(0, n_pair, n_lsub),
        matrix(0, n_pair, N)
      )

      # add non-decreasing monotonicity matrix to the envelopment matrix
      Amat <- rbind(Amat, MMat)
    }

    # concavity & increasing monotonicity: left-side unpaired basis functions
    if (monotonicity | concavity) {
      MCMat <- cbind (
        matrix(0, n_lsub, 1 + n_pair),
        diag(- 1, n_lsub),
        matrix(0, n_lsub, N)
      )

      # add concavity & increasing monotonicity [LSBF] to envelopment matrix
      Amat <- rbind(Amat, MCMat)
    }

    # ==================================================== #
    # b: envelopment + concavity + monotonicity            #
    # ==================================================== #

    if (origin == "1") {
      bvec <- c(y_ind, 1, rep(0, nrow(Amat) - N - 1))

    } else {
      bvec <- c(y_ind, rep(0, nrow(Amat) - N))

    }

    # ==================================================== #
    # Solution of the optimization problem                 #
    # ==================================================== #

    if (origin != "FALSE") {
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
