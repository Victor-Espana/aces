#' @title Estimate Coefficients during Forward and Backward Algorithms.
#'
#' @description This function solves a Linear Programming problem to obtain a set of coefficients (that impose concavity and monotonicity) in the Adaptive Constrained Enveloping Splines estimator.
#'
#' @param B A \code{matrix} of linear basis functions.
#' @param y A \code{matrix} of output data.
#' @param it_list A \code{list} containing the set of intervals by input.
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return A \code{vector} of estimated coefficients.

estimate_coefficients <- function (
    B, y, it_list, Bp_list, monotonicity, concavity, x0_y0
    ) {

  # Sample size
  N  <- nrow(B)
  # Number of basis functions
  p  <- ncol(B)
  # Number of outputs
  nY <- ncol(y)
  # Number of inputs
  nX <- length(Bp_list)

  coefs <- matrix(NA, nrow = p, ncol = nY)

  for (out in 1:nY) {

    # Select an output
    y_ind <- y[, out]

    # ==================================================== #
    # vars: c(coef_0, coef_1, ..., coef_P, e_1, ... , e_n) #
    # ==================================================== #
    objVal <- c(rep(0, p), rep(1, N))

    # ==================================================== #
    # A: envelopment + concavity + monotonicity            #
    # ==================================================== #

    # Envelopment: y_hat - e = y
    EMat <- cbind(B, diag(rep(- 1, N), N))

    # Matrix of coefficients
    Amat <- rbind(EMat)

    if (x0_y0) {
      # f(0, 0) = 0
      const_0_0 <- predict_0_0 (
        it_list = it_list,
        Bp_list = Bp_list,
        n_bf = p
      )
      const_0_0 <- c(const_0_0, rep(0, N))

      # Add cons0 to Amat
      Amat <- rbind(Amat, const_0_0)
    }

    # Increasing monotonicity
    if (monotonicity) {
      MMat <- monotonocity_matrix (
        it_list = it_list,
        Bp_list = Bp_list,
        N = N
        )

      # Add MMat to Amat
      Amat <- rbind(Amat, MMat)
    }

    # Concavity
    if (concavity) {
      CMat <- concavity_matrix(
        Bp_list = Bp_list,
        N = N
        )

      # Add CMat to Amat
      Amat <- rbind(Amat, CMat)
    }

    # ==================================================== #
    # b: envelopment + concavity + monotonicity            #
    # ==================================================== #
    bvec <- c(y_ind, rep(0, nrow(Amat) - N))

    # ==================================================== #
    # Directions of inequalities                           #
    # ==================================================== #
    if (x0_y0) {
      dirs <- c(rep("==", N + 1), rep(">=", nrow(Amat) - N - 1))
    } else {
      dirs <- c(rep("==", N), rep(">=", nrow(Amat) - N))
    }

    # ==================================================== #
    # Lower and upper bounds                               #
    # ==================================================== #
    bnds <- list(lower = list(ind = 1:p, val = rep(- Inf, p)))

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

    # if (sols$status == 1) {
      # sols <- Rglpk_solve_LP(
        # objVal, Amat, dirs, bvec,
        # bnds, control = list(presolve = T)
      # )
    # }

    coefs[, out] <- sols$solution[1:p]
  }

  return(coefs)
}

#' @title  Build Constraint to Predict f(0) = 0
#'
#' @description This function constructs the constraint to enforce f(0) = 0.
#'
#' @param it_list A \code{list} containing the set of intervals by input.
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#' @param n_bf Number of basis functions (including the intercept term).
#'
#' @return A \code{vector} to enforce f(0) = 0.

predict_0_0 <- function (
    it_list, Bp_list, n_bf
) {

  # =========== #
  # A: f(0) = 0 #
  # =========== #

  # number of inputs
  nX <- length(Bp_list)

  # basis activated at first interval (all betas for each variable)
  bf_it0 <- vector("list", nX)

  # Bp indexes for the first interval
  for (v in 1:nX) {
    n_int <- length(it_list[[v]])
    if (n_int == 0 || is.null(it_list[[v]][[1]])) next

    bf_it0[[v]] <- it_list[[v]][[1]][["Bp"]][["Bp_xi"]]
  }

  # constraint to ensure f(0) = 0
  const_0_0 <- c(1)

  for (v in 1:nX) {

    # number of paired basis functions
    n_pbf <- length(Bp_list[[v]][["paired"]])
    # number of right-side basis functions
    n_rbf <- length(Bp_list[[v]][["right"]])
    # number of left-side basis functions
    n_lbf <- length(Bp_list[[v]][["left"]])

    # basis functions for the v-th variable
    const_0_0_v <- rep(0, 2 * n_pbf + n_rbf + n_lbf)

    for (side in c("paired", "right", "left")) {
      if (is.null(Bp_list[[v]][[side]])) next
      for (j in 1:length(bf_it0[[v]])) {
        if (is.null(bf_it0[[v]][j])) next
        for (b in 1:length(Bp_list[[v]][[side]])) {
          if (bf_it0[[v]][j] %in% Bp_list[[v]][[side]][[b]][["Bp_xi"]]) {

            const_0_0_v[bf_it0[[v]][j]] <- Bp_list[[v]][[side]][[b]][["t"]]
          }
        }
      }
    }

    # add to const_0_0
    const_0_0 <- c(const_0_0, const_0_0_v)
  }

  return(const_0_0)
}

#' @title Build Constraint Matrix for Enforcing Non-decreasing Monotonicity
#'
#' @description This function constructs the constraint matrix (often denoted as matrix "A") to enforce non-decreasing monotonicity.
#'
#' @param it_list A \code{list} containing the set of intervals by input.
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#' @param N Sample size.
#'
#' @importFrom Matrix bdiag
#'
#' @return A \code{matrix} to enforce increasing monotonicity constraint on the estimated function.

monotonocity_matrix <- function (
    it_list, Bp_list, N
    ) {

  # ========================================= #
  # A: envelopment + concavity + MONOTONICITY #
  # ========================================= #

  # Number of inputs
  nX <- length(Bp_list)

  # Initialize monotonicity matrix
  MMat <- vector("list", nX)

  for (v in 1:nX) {

    # Number of paired basis functions
    n_pbf <- length(Bp_list[[v]][["paired"]])
    # Number of right-side basis functions
    n_rbf <- length(Bp_list[[v]][["right"]])
    # Number of left-side basis functions
    n_lbf <- length(Bp_list[[v]][["left"]])
    # Number of intervals by variable
    n_int <- length(it_list[[v]])

    if (n_pbf + n_rbf + n_lbf == 0) next

    MMat_xi <- matrix(0, nrow = n_int, ncol = 2 * n_pbf + n_rbf + n_lbf)

    for (it in 1:n_int) {
      Bp_xi <- it_list[[v]][[it]][["Bp"]][, "Bp_xi"]
      # constrained value
      con_val <- c()

      for (side in it_list[[v]][[it]][["Bp"]][, "status"]) {

        if (side == "right") {
          con_val <- c(con_val, 1)
        } else {
          con_val <- c(con_val, - 1)
        }
      }
      MMat_xi[it, Bp_xi] <- con_val
    }
    MMat[[v]] <- MMat_xi
  }

  # Drop NULL matrix
  MMat[sapply(MMat, is.null)] <- NULL

  # monotonicity matrix
  MMat <- as.matrix(Matrix::bdiag(MMat))
  MMat <- cbind(rep(0, nrow(MMat)), MMat, matrix(0, nrow = nrow(MMat), ncol = N))

  return(MMat)
}

#' @title  Build Constraint Matrix for Enforcing Concavity
#'
#' @description This function constructs the constraint matrix (often denoted as matrix "A") to enforce concavity.
#'
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#' @param N Sample size.
#'
#' @importFrom Matrix bdiag
#'
#' @return A \code{matrix} to enforce concavity constraint on the estimated function.

concavity_matrix <- function (
    Bp_list, N
    ) {

  # ========================================= #
  # A: envelopment + CONCAVITY + monotonicity #
  # ========================================= #

  # Number of inputs
  nX <- length(Bp_list)

  # Initialize concavity matrix
  CMat <- vector("list", nX)

  for (v in 1:nX) {

    # Number of paired basis functions
    n_pbf <- length(Bp_list[[v]][["paired"]])
    # Number of right-side basis functions
    n_rbf <- length(Bp_list[[v]][["right"]])
    # Number of left-side basis functions
    n_lbf <- length(Bp_list[[v]][["left"]])

    if (n_pbf > 0) {
      CMat1 <- matrix(0, nrow = n_pbf, ncol = 2 * n_pbf)
      for (i in 1:nrow(CMat1)) {
        CMat1[i, ] <- c(rep(c(0, 0), i - 1), - 1, - 1, rep(c(0, 0), n_pbf - i))
      }
      CMat[[v]][[1]] <- CMat1
    }

    if (n_rbf + n_lbf > 0) {
      CMat2 <- diag(- 1, n_rbf + n_lbf, n_rbf + n_lbf)
      CMat[[v]][[2]] <- CMat2
    }

    if (!is.null(CMat[[v]])) {
      CMat[[v]][sapply(CMat[[v]], is.null)] <- NULL
      CMat[[v]] <- as.matrix(Matrix::bdiag(CMat[[v]]))
    }
  }

  # Drop NULL matrix
  CMat[sapply(CMat, is.null)] <- NULL

  # Concavity matrix
  CMat <- as.matrix(Matrix::bdiag(CMat))
  CMat <- cbind(rep(0, nrow(CMat)), CMat, matrix(0, nrow = nrow(CMat), ncol = N))

  return(CMat)
}

#' @title Estimate Coefficients during Shape-Restricted Smoothing Procedure.
#'
#' @description This function solves a Linear Programming problem to obtain a set of coefficients (that impose concavity and monotonicity) in the Smoothing Adaptive Constrained Enveloping Splines estimator.
#'
#' @param B A \code{matrix} of cubic / quintic basis functions.
#' @param y The output \code{vector} in data.
#' @param h Number of coefficients associated with the paired basis functions.
#' @param r Number of coefficients associated with the unpaired left-side basis functions.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return A \code{vector} with the estimated coefficients.

estimate_coefficients_smoothed <- function (
    B, y, h, r, monotonicity, concavity, x0_y0
    ) {

  # sample size
  N <- nrow(B)
  # number of basis functions: 1 + h + r
  p <- ncol(B)
  # number of pairs
  pairs <- h / 2
  # number of outputs
  nY <- ncol(y)

  # initialize coefficients
  coefs <- matrix(NA, nrow = p, ncol = nY)

  for (out in 1:nY) {

    # vars: c(tau_0, alpha_1, beta_2, ... , alpha_(H-1), beta_H, w_1, ..., w_R, e_1, ... , e_n)
    objVal <- c(rep(0, p), rep(1, N))

    # select a variable
    y_ind <- y[, out]

    # ========================================= #
    # A: envelopment + concavity + monotonicity #
    # ========================================= #

    # Equality constraints: y_hat - e = y
    Emat <- cbind(B, diag(rep(- 1, N), N))

    Amat <- rbind(Emat)

    # Concavity [paired basis functions]
    if (concavity) {
      CMat <- matrix(0, nrow = pairs, ncol = p + N)

      if (pairs >= 1) {
        for (i in 1:pairs) {
          CMat[i, ] <- c(0, rep(c(0, 0), i - 1), - 1, - 1, rep(c(0, 0), pairs - i), rep(0, r + N))
        }
      }

      # add concavity matrix to envelopment matrix
      Amat <- rbind(Amat, CMat)
    }

    if (monotonicity) {
      # increasing monotonicity [paired basis functions]
      MMat <- cbind(rep(0, h), diag(c(1, - 1), h), matrix(0, h, r), matrix(0, h, N))

      # add monotonicity matrix to envelopment matrix
      Amat <- rbind(Amat, MMat)
    }

    if (monotonicity | concavity) {
      # concavity & increasing monotonicity [left-side basis functions]
      MCMat <- cbind(matrix(0, r, 1 + h), diag(- 1, r), matrix(0, r, N))

      # add concavity & increasing monotonicity [LSBF] to envelopment matrix
      Amat <- rbind(Amat, MCMat)
    }

    # ==================================================== #
    # b: envelopment + concavity + monotonicity            #
    # ==================================================== #
    bvec <- c(y_ind, rep(0, nrow(Amat) - N))

    # ==================================================== #
    # Directions of inequalities                           #
    # ==================================================== #
    dirs <- c(rep("==", N), rep(">=", nrow(Amat) - N))

    if (x0_y0) {
      # f(0) = 0
      const_0_0 <- c(1, rep(c(0, 1), pairs), rep(1, r), rep(0, N))

      Amat <- rbind(Amat, const_0_0)
      bvec <- c(bvec, 0)
      dirs <- c(rep("==", N), rep(">=", nrow(Amat) - N - 1), "==")
    }

    # ==================================================== #
    # Lower and upper bounds                               #
    # ==================================================== #
    bnds <- list(lower = list(ind = 1L:p, val = rep(- Inf, p)))

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


