#' @title Estimate Smoothed ACES Coefficients
#'
#' @description
#'
#' Solves the reduced linear program for a cubic- or quintic-smoothed ACES model.
#'
#' Error variables are removed analytically so that the optimization contains
#' only the paired and unpaired smooth-basis coefficients. Envelopment and the
#' requested monotonicity or concavity restrictions are imposed in the reduced
#' problem.
#'
#' @param B
#' Matrix of evaluated cubic or quintic basis functions.
#'
#' @param y_obs
#' Matrix of observed outputs.
#'
#' @param dea_scores
#' A matrix of output-oriented DEA-VRS scores, with one column per output.
#'
#' @param n_pair
#' Number of coefficients for paired basis functions.
#'
#' @param n_lsub
#' Number of coefficients for unpaired left-side basis functions.
#'
#' @param shape
#' A list with logical elements \code{mono} and \code{conc}.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return
#'
#' A numeric vector of estimated coefficients.

estimate_coefficients_smoothed <- function(
  B,
  y_obs,
  dea_scores,
  n_pair,
  n_lsub,
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

  # number of basis functions: 1 + n_pair + n_lsub
  p <- ncol(B)

  # number of pairs
  pairs <- n_pair / 2

  # number of outputs
  nY <- ncol(y_obs)

  # initialize coefficients structure
  coefs <- matrix(NA, nrow = p, ncol = nY)

  for (out in 1:nY) {
    # vars: c(tau_0, alpha_1, beta_2, ... , alpha_(H-1), beta_H, w_1, ..., w_R)

    # objective function: (w' B)
    objVal <- as.vector(crossprod(B, deaw[, out]))

    # select a variable
    y_ind <- y_obs[, out]

    # envelopment: y_hat >= y (equivalent to e >= 0)
    Amat <- B

    # generate concavity matrix for paired basis functions
    if (conc) {
      CMat <- matrix(0, nrow = pairs, ncol = p)

      if (pairs >= 1) {
        for (i in 1:pairs) {
          CMat[i, ] <- c(
            0,
            rep(c(0, 0), i - 1),
            -1,
            -1,
            rep(c(0, 0), pairs - i),
            rep(0, n_lsub)
          )
        }
      }

      Amat <- rbind(Amat, CMat)
    }

    # generate non-decreasing monotonicity matrix for paired basis functions
    if (mono) {
      MMat <- cbind(
        rep(0, n_pair),
        diag(c(1, -1), n_pair),
        matrix(0, n_pair, n_lsub)
      )

      Amat <- rbind(Amat, MMat)
    }

    # generate concavity & increasing monotonicity matrix for LF-U-BF
    if (mono | conc) {
      MCMat <- cbind(
        matrix(0, n_lsub, 1 + n_pair),
        diag(-1, n_lsub)
      )

      Amat <- rbind(Amat, MCMat)
    }

    # right-hand side
    bvec <- c(y_ind, rep(0, nrow(Amat) - N))

    # directions of inequalities
    dirs <- rep(">=", nrow(Amat))

    # lower and upper bounds
    bnds <- list(
      lower = list(
        ind = 1L:p,
        val = rep(-Inf, p)
      )
    )

    # solution of the optimization problem
    sols <- Rglpk_solve_LP(
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
