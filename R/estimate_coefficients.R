#' @title Build Constraint to Predict f(0) = 0
#'
#' @description
#'
#' This function constructs the constraint to ensure that the function prediction at zero equals zero, i.e., f(0) = 0. In the context of the model, this means that the model prediction for the output variable should be zero when all input variables are zero.
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
#' A \code{vector} of coefficients that enforces the constraint f(0) = 0. This vector can be used to modify the model coefficients to satisfy the constraint.

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

  for (v in 1:nX) {

    if (length(it_list[[v]]) == 0 || is.null(it_list[[v]][[1]])) next

    bf_it0[[v]] <- it_list[[v]][[1]][["Bp"]][["Bp_xi"]]

  }

  # constraint to ensure f(0) = 0
  origin_vec <- c(1)

  for (v in 1:nX) {

    # verify bf_it0[[v]] is not NULL
    if (length(it_list[[v]]) == 0) next

    # basis for the v-th variable
    Bp_v <- Bp_list[[v]]

    # number of paired basis functions
    n_pbf <- length(Bp_v[["paired"]])

    # number of right-side basis functions
    n_rbf <- length(Bp_v[["right"]])

    # number of left-side basis functions
    n_lbf <- length(Bp_v[["left"]])

    # basis functions for the v-th variable
    origin_vec_v <- rep(0, 2 * n_pbf + n_rbf + n_lbf)

    # exist an interval but not a left-side basis function
    if (!is.null(bf_it0[[v]])) {

      for (side in c("paired", "right", "left")) {

        Bp_v_side <- Bp_v[[side]]

        if (!is.null(Bp_v_side)) {

          for (b in 1:length(Bp_v_side)) {

            Bp_xi <- Bp_v_side[[b]][["Bp_xi"]]
            t_val <- Bp_v_side[[b]][["t"]]

            for (j in 1:length(bf_it0[[v]])) {

              if (is.null(bf_it0[[v]][j])) next

              if (bf_it0[[v]][j] %in% Bp_xi) {

                origin_vec_v[bf_it0[[v]][j]] <- t_val

              }
            }
          }
        }
      }

    }

    # add to vector of constraints
    origin_vec <- c(origin_vec, origin_vec_v)

  }

  # add variables of error to satisfy size format
  origin_vec <- c(origin_vec, rep(0, N))

  return(origin_vec)

}

#' @title Build Constraint Matrix for Enforcing Non-decreasing Monotonicity
#'
#' @description
#' This function constructs the constraint matrix (often denoted as matrix "A") to enforce non-decreasing monotonicity in an estimated function. The constraint matrix enforces concavity by ensuring that the coefficients of the basis functions satisfy certain conditions.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param N
#' Sample size.
#'
#' @importFrom Matrix bdiag
#'
#' @return
#' A \code{matrix} representing the constraint matrix to enforce non-decreasing monotonicity on the estimated function.

monotonocity_matrix <- function (
    it_list,
    Bp_list,
    N
    ) {

  # ========================================= #
  # A: envelopment + concavity + MONOTONICITY #
  # ========================================= #

  # number of inputs
  nX <- length(Bp_list)

  # initialize monotonicity matrix
  MMat <- vector("list", nX)

  for (v in 1:nX) {

    # number of paired basis functions
    n_pbf <- length(Bp_list[[v]][["paired"]])

    # number of right-side basis functions
    n_rbf <- length(Bp_list[[v]][["right"]])

    # number of left-side basis functions
    n_lbf <- length(Bp_list[[v]][["left"]])

    # number of intervals by variable
    n_int <- length(it_list[[v]])

    if (n_pbf + n_rbf + n_lbf == 0) next

    # monotonicity matrix for the v-th input
    MMat_xi <- matrix(0, nrow = n_int, ncol = 2 * n_pbf + n_rbf + n_lbf)

    for (it in 1:n_int) {

      Bp_xi <- it_list[[v]][[it]][["Bp"]][, "Bp_xi"]

      # create constraint
      con_val <- c()

      for (side in it_list[[v]][[it]][["Bp"]][, "status"]) {
        con_val <- c(con_val, ifelse(side == "right", 1, -1))
      }

      MMat_xi[it, Bp_xi] <- con_val

    }

    MMat[[v]] <- MMat_xi

  }

  # drop NULL matrix
  MMat[sapply(MMat, is.null)] <- NULL

  # monotonicity matrix
  MMat <- as.matrix(Matrix::bdiag(MMat))

  MMat <- cbind (
    rep(0, nrow(MMat)),
    MMat,
    matrix(0, nrow = nrow(MMat), ncol = N)
    )

  return(MMat)

}

#' @title  Build Constraint Matrix for Enforcing Concavity
#'
#' @description
#' This function constructs the constraint matrix (often denoted as matrix "A") to enforce concavity in the estimated function. The constraint matrix enforces concavity by ensuring that the coefficients of the basis functions satisfy certain conditions.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param N
#' Sample size.
#'
#' @importFrom Matrix bdiag
#'
#' @return
#' A \code{matrix} representing the constraint matrix to enforce concavity on the estimated function.

concavity_matrix <- function (
    Bp_list,
    N
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

  # drop NULL matrix
  CMat[sapply(CMat, is.null)] <- NULL

  # concavity matrix
  CMat <- as.matrix(Matrix::bdiag(CMat))
  CMat <- cbind (
    rep(0, nrow(CMat)),
    CMat,
    matrix(0, nrow = nrow(CMat), ncol = N)
    )

  return(CMat)

}
