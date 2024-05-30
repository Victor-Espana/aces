#' @title Add a New Pair of Basis Functions in Adaptive Constrained Enveloping Splines
#'
#' @description
#'
#' This function adds the pair of basis functions that results in the largest reduction of the lack-of-fit criterion.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param xi_degree
#' A \code{matrix} indicating the degree of each input variable.
#'
#' @param dea_eff
#' An indicator vector with ones for efficient DMUs and zeros for inefficient DMUs.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function estimates.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.
#'
#' @param forward_model
#' A \code{list} containing the current set of basis functions (\code{BF_set}) and the B matrix (\code{B}).
#'
#' @param Bp_list
#' A \code{list} containing the current set of basis functions for each input variable.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param kn_list
#' A \code{list} containing the current set of selected knots for each input variable.
#'
#' @param kn_grid
#' A \code{list} with the grid of knots to perform ACES.
#'
#' @param L
#' A \code{integer} value specifying the minimum number of observations between two adjacent knots.
#'
#' @param Le
#' A \code{integer} value specifying the minimum number of observations before the first and after the final knot.
#'
#' @param err_min
#' A \code{numeric} value specifying the minimum error obtained by the forward algorithm in the current iteration.
#'
#' @param hd_cost
#' A \code{numeric} value specifying the minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.
#'
#' @return
#'
#' An updated \code{list} containing the matrix of basis functions (\code{B}), an updated \code{list} with information of the basis functions (\code{BF_set}), an updated \code{list} of selected knots (\code{knots_list}) and the minimum error obtained (\code{err_min}).

add_basis_function <- function (
    data,
    x,
    y,
    xi_degree,
    dea_eff,
    model_type,
    metric,
    forward_model,
    Bp_list,
    shape,
    kn_list,
    kn_grid,
    L,
    Le,
    err_min,
    hd_cost
    ) {

  # sample size
  N <- nrow(data)

  # number of inputs / outputs as inputs
  nX <- ncol(data) - length(y)

  # set of basis functions
  bf_set <- forward_model[["bf_set"]]

  # number of basis functions
  nbf <- length(bf_set)

  # matrix of basis functions
  B <- forward_model[["B"]]

  # signal to indicate improvement in the fitting
  improvement <- FALSE

  # basis function for the expansion: it is always the first in the additive model
  bf <- bf_set[[1]]

  for (xi in x) {

    # ================================ #
    # Create the set of eligible knots #
    # ================================ #

    knots <- set_knots (
      data = data,
      nX = nX,
      var_idx = xi,
      minspan = L,
      endspan = Le,
      kn_list = kn_list,
      bf = bf,
      kn_grid = kn_grid
      )

    # skip to the next step if there are not eligible knots
    if (is.null(knots)) next

    # shuffle the knots randomly
    knots <- sample(knots)

    for (i in 1:length(knots)) {

      Bp_list_aux <- Bp_list

      # ================== #
      # Create two new bfs #
      # ================== #

      new_pair <- create_linear_basis (
        data = data,
        var_idx = xi,
        knot = knots[i],
        B = B
        )

      # ============== #
      # Update Bp_list #
      # ============== #

      # number of paired basis functions for the input "xi"
      nbf_xi <- length(Bp_list_aux[[xi]][["paired"]])

      # add the 2 new basis functions to the Bp_list
      new_bf <- list (
        "id" = c(nbf + 1, nbf + 2),
        "Bp" = cbind(new_pair[[1]], new_pair[[2]]),
        "t"  = knots[i]
        )

      Bp_list_aux[[xi]][["paired"]][[nbf_xi + 1]] <- new_bf

      # sort basis function by variable, side and knot
      for (v in 1:nX) {
        for (side in c("paired", "right", "left")) {
          if (!is.null(Bp_list_aux[[v]][[side]])) {
            arrangement <- order(sapply(Bp_list_aux[[v]][[side]], "[[", "t"))
            Bp_list_aux[[v]][[side]] <- Bp_list_aux[[v]][[side]][arrangement]
          }
        }
      }

      # column indexes of the BFs in the B matrix by variable
      # (intercept column is not considered)
      for (v in 1:nX) {
        Bp_xi <- 1
        for (side in c("paired", "right", "left")) {
          if (!is.null(Bp_list_aux[[v]][[side]])) {
            for (l in 1:length(Bp_list_aux[[v]][[side]])) {
              if (side == "paired") {
                Bp_list_aux[[v]][[side]][[l]][["Bp_xi"]] <- c(Bp_xi, Bp_xi + 1)
                Bp_xi <- Bp_xi + 2
              } else {
                Bp_list_aux[[v]][[side]][[l]][["Bp_xi"]] <- Bp_xi
                Bp_xi <- Bp_xi + 1
              }
            }
          }
        }
      }

      # ============== #
      # Update it_list #
      # ============== #

      it_list <- set_intervals (
        Bp_list = Bp_list_aux
        )

      # =============== #
      # Update B matrix #
      # =============== #

      new_B <- matrix(1, nrow = N)
      for (v in 1:nX) {
        for (side in c("paired", "right", "left")) {
          if (!is.null(Bp_list_aux[[v]][[side]])) {
            for (l in 1:length(Bp_list_aux[[v]][[side]])) {
              bf_jp <- Bp_list_aux[[v]][[side]][[l]][["Bp"]]
              new_B <- cbind(new_B, bf_jp)
            }
          }
        }
      }

      # ===================== #
      # Estimate coefficients #
      # ===================== #

      coefs <- estimate_coefficients (
        model_type = model_type,
        B = new_B,
        y_obs = data[, y, drop = F],
        dea_eff = dea_eff,
        it_list = it_list,
        Bp_list = Bp_list_aux,
        shape = shape
        )

      # predictions
      y_hat <- matrix(NA, nrow = N, ncol = length(y))

      for (out in 1:length(y)) {
        y_hat[, out] <- new_B %*% coefs[, out]
      }

      # lack-of-fit measure
      err <- err_metric (
        y_obs = data[, y, drop = F],
        y_hat = y_hat[drop = F],
        metric = metric
        )

      # variable
      err[2] <- xi_degree[2, xi]

      # compute gcv
      nknots_xi <- sapply(kn_list, function(x) {
        if (is.null(x)) {
          return(0)

        } else {
          return(length(x))

          }
        }
      )

      GCV <- compute_gcv (
        y_obs = data[, y, drop = F],
        y_hat = y_hat,
        metric = metric,
        n_bf = ncol(new_B),
        d = 1,
        knots = sum(nknots_xi) + 1,
        xi_degree = NULL
      )

      # compute GRSq: Milborrow (2014)
      GRSq <- 1 - GCV / bf_set[[1]][["GCV"]]

      # check if the best basis function should be checked
      add <- FALSE

      if (is.na(err_min[2]) || err_min[2] == 1) {

        if (xi_degree[2, xi] == 1) {

          add <- err[1] < err_min[1]

        } else if (xi_degree[2, xi] > 1) {

          add <- (err[1] - err_min[1]) / err_min[1] < - hd_cost

        }

      } else {

        add <- err[1] < err_min[1]

      }

      if (add) {

        # model has improved
        improvement <- TRUE

        # minimum error
        err_min <- err

        # best B matrix
        best_B <- new_B

        # best Bp_list
        best_Bp_list <- Bp_list_aux

        # index
        tindex <- which(kn_grid[[xi]] == knots[i])

        # new pair of basis functions
        bf1 <- bf2 <-  bf

        # id
        bf1[["id"]] <- nbf + 1
        bf2[["id"]] <- nbf + 2

        # status
        bf1[["status"]] <- bf2[["status"]] <- "paired"

        # side
        bf1[["side"]] <- "R"
        bf2[["side"]] <- "L"

        # Bp
        bf1[['Bp']] <- new_pair[[1]]
        bf2[['Bp']] <- new_pair[[2]]

        # xi
        if (all(bf[['xi']] == - 1)) {
          bf1[['xi']] <- bf2[['xi']] <- c(xi)

        } else {
          bf1[['xi']] <- bf2[['xi']] <- c(bf[['xi']], xi)
        }

        # knot
        bf1[["t"]] <- bf2[["t"]] <- knots[i]

        # R
        bf1[['R']] <- bf2[['R']] <- err_min[1]

        # GCV
        bf1[['GCV']] <- bf2[['GCV']] <- GCV

        # GRSq
        bf1[['GRSq']] <- bf2[['GRSq']] <- GRSq

        # coefficients
        bf1[['coefs']] <- bf2[['coefs']] <- coefs

      }
    }
  }

  if (improvement) {

    # append new basis functions
    bf_set <- append(bf_set, list(bf1))
    bf_set <- append(bf_set, list(bf2))

    # knots for each variable
    len_knt <- length(bf1[["xi"]])
    var_pos <- bf1[["xi"]][len_knt]

    kn_list[[var_pos]] <- append (
      kn_list[[bf1[['xi']][len_knt]]],
      list(list(t = bf1[["t"]], index = tindex))
      )

    return(list(best_B, bf_set, kn_list, best_Bp_list, err_min))

  } else {

    return(improvement)

  }

}

#' @title Generate the Set of Eligible Knots
#'
#' @description
#' This function generates a vector of knots that can be used to create a new pair of basis functions during the forward algorithm.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param nX
#' Number of inputs.
#'
#' @param var_idx
#' Index of the input variable that creates the new pair of basis functions.
#'
#' @param minspan
#' Minimum number of observations between two adjacent knots.
#'
#' @param endspan
#' Minimum number of observations before the first and after the final knot.
#'
#' @param kn_list
#' A \code{list} containing the set of selected knots.
#'
#' @param bf
#' A \code{list} with the basis function for the expansion of the model.
#'
#' @param kn_grid
#' Grid of virtual knots to perform ACES.
#'
#' @importFrom dplyr %>%
#'
#' @return
#' Numeric vector with the knots values for the variable.

set_knots <- function (
    data,
    nX,
    var_idx,
    minspan,
    endspan,
    kn_list,
    bf,
    kn_grid
    ) {

  # sample size
  N <- nrow(data)

  # minimum number of observations between two adjacent knots
  sp1 <- ifelse(length(minspan) > 1, minspan[var_idx], minspan)

  # minimum number of observations before the first and after the final knot.
  sp2 <- ifelse(length(endspan) > 1, endspan[var_idx], endspan)

  # observations in the space of the basis function.
  nobs_bf <- bf[['Bp']] != 0

  # minimum span for Friedman's approach
  if (sp1 == - 1) {
    sp1 <- floor(min(N * 0.1, - log2(- (1 / (nX * sum(nobs_bf))) * log(1 - 0.05)) / 2.5))
  }

  # boolean: possible knots
  kn_idx <- rep(TRUE, length(kn_grid[[var_idx]]))

  # Set to FALSE the first and last "Le" observations
  if (sp2 != 0) {
    set_sp2_idx <- c(1:sp2, (length(kn_idx) - sp2 + 1):length(kn_idx))

  } else {
    set_sp2_idx <- c()

  }

  # Set to FALSE adjacent "L" observations
  if (!is.null(kn_list[[var_idx]])) {

    # used knots values
    used_kn_val <- sapply(kn_list[[var_idx]], "[[", "t")

    # used knots indexes
    used_kn_idx <- c()

    for (kn in 1:length(used_kn_val)) {
      used_kn_idx <- c (
        used_kn_idx,
        which(kn_grid[[var_idx]][order(kn_grid[[var_idx]])] == used_kn_val[kn])
      )
    }

    set_sp1_idx <- c()

    for (idx in 1:length(used_kn_idx)) {
      set_sp1_idx <- append (
        set_sp1_idx,
        (used_kn_idx[idx] - sp1):(used_kn_idx[idx] + sp1)
        )
    }

    # indexes must be greater than 0 and lower than N
    set_sp1_idx <- unique(sort(set_sp1_idx[set_sp1_idx > 0 & set_sp1_idx <= sum(nobs_bf)]))

  } else {

    set_sp1_idx <- c()

  }

  # set of observations to FALSE based on minspan and endspan
  kn_idx[sort(unique(c(set_sp1_idx, set_sp2_idx)))] <- FALSE

  # define knots
  if (all(kn_idx == FALSE)) {

    knots <- NULL

  } else {

    # minimum value of observations in the basis function
    minimum <- min(data[nobs_bf, var_idx])

    # maximum value of observations in the basis function
    maximum <- max(data[nobs_bf, var_idx])

    # exclude observations according to the minspan and the endspan
    varknots <- kn_grid[[var_idx]][order(kn_grid[[var_idx]])][kn_idx]

    # knots between minimum and maximum
    knots <- varknots[varknots >= minimum & varknots <= maximum]

  }

  if (length(knots) <= 1) knots <- NULL

  return(knots)

}

#' @title Generate a New Pair of Linear Basis Functions
#'
#' @description
#' This function generates two new linear basis functions from a variable and a knot.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param var_idx
#' Index of the input variable that creates the new pair of basis functions.
#'
#' @param knot
#' Knot for creating the new pair of basis functions.
#'
#' @param B
#' A \code{matrix} of basis functions on which the new pair of functions is added.
#'
#' @return
#' A \code{list} with the new pair of basis functions.

create_linear_basis <- function (
    data,
    var_idx,
    knot,
    B
    ) {

  # create (xi-t)+ and (t-xi)+
  hinge1 <- pmax(0, data[, var_idx] - knot)
  hinge2 <- pmax(0, knot - data[, var_idx])

  # two new basis functions
  bf1 <- B[, 1] * hinge1
  bf2 <- B[, 1] * hinge2

  return(list(bf1, bf2))

}
