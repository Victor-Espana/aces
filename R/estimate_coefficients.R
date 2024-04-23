#' @title Build Constraint to Predict f(0) = 0
#'
#' @description
#'
#' This function constructs the constraint to ensure that the function prediction at zero equals zero, i.e., f(0) = 0.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param n_bf
#' An \code{integer} specifying the number of basis functions, including the intercept term.
#'
#' @param N
#' Sample size.
#'
#' @return
#'
#' A \code{vector} of coefficients that enforces the constraint f(0) = 0.

predict_origin <- function (
    it_list,
    Bp_list,
    n_bf,
    N
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
  Ovec <- c(1)

  for (v in 1:nX) {

    # number of paired basis functions
    n_pbf <- length(Bp_list[[v]][["paired"]])

    # number of right-side basis functions
    n_rbf <- length(Bp_list[[v]][["right"]])

    # number of left-side basis functions
    n_lbf <- length(Bp_list[[v]][["left"]])

    # basis functions for the v-th variable
    Ovec_v <- rep(0, 2 * n_pbf + n_rbf + n_lbf)

    for (side in c("paired", "right", "left")) {
      if (is.null(Bp_list[[v]][[side]])) next

      for (j in 1:length(bf_it0[[v]])) {
        if (is.null(bf_it0[[v]][j])) next

        for (b in 1:length(Bp_list[[v]][[side]])) {
          if (bf_it0[[v]][j] %in% Bp_list[[v]][[side]][[b]][["Bp_xi"]]) {
            Ovec_v[bf_it0[[v]][j]] <- Bp_list[[v]][[side]][[b]][["t"]]
          }
        }
      }
    }

    # add to Ovec
    Ovec <- c(Ovec, Ovec_v)
  }

  # add variables "e"
  Ovec <- c(Ovec, rep(0, N))

  return(Ovec)
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

  # number of inputs
  nX <- length(Bp_list)

  # initialize concavity matrix
  CMat <- vector("list", nX)

  for (v in 1:nX) {

    # number of paired basis functions
    n_pbf <- length(Bp_list[[v]][["paired"]])
    # number of right-side basis functions
    n_rbf <- length(Bp_list[[v]][["right"]])
    # number of left-side basis functions
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
