#' @title Add a New Pair of Basis Functions in Adaptive Constrained Enveloping Splines
#'
#' @description This function adds the pair of basis functions that results in the largest reduction of the lack-of-fit criterion.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param x Column indexes of input variables in \code{data}.
#' @param xi_degree Matrix indicating the degree of each input variable.
#' @param y Column indexes of output variables in \code{data}.
#' @param metric Lack-of-fit criterion to evaluate the model performance.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#' @param forward_model A \code{list} containing the current set of basis functions (\code{BF_set}) and the B matrix (\code{B}).
#' @param kn_list A \code{list} containing the current set of selected knots for each input variable.
#' @param Bp_list A \code{list} containing the current set of basis functions for each input variable.
#' @param L Minimum number of observations between two adjacent knots.
#' @param Le Minimum number of observations before the first and after the final knot.
#' @param kn_grid Grid of knots to perform ACES.
#' @param err_min Minimum error obtained by the forward algorithm in the current iteration.
#' @param hd_cost Minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.
#'
#' @return An undated \code{list} containing the matrix of basis functions (\code{B}), an undated \code{list} with information of the basis functions (\code{BF_set}), an updated \code{list} of selected knots (\code{knots_list}) and the minimum error obtained (\code{err_min}).

add_basis_function <- function (
    data, x, xi_degree, y, metric, monotonicity, concavity, x0_y0, forward_model,
    kn_list, Bp_list, L, Le, kn_grid, err_min, hd_cost
    ) {

  # Samples in data
  N <- nrow(data)

  # Number of inputs / outputs as inputs
  nX <- length(x)

  # Set of basis functions
  bf_set <- forward_model[["bf_set"]]

  # Number of basis functions
  nbf <- length(bf_set)

  # Matrix of basis functions
  B <- forward_model[["B"]]

  # Signal of improvement
  improvement <- FALSE

  # Basis function for the expansion: additive model (always 1)
  bf <- bf_set[[1]]

  for (xi in 1:nX) {

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

    # skip to the next step if there are no knots
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

      # Number of paired basis functions for the input "xi"
      nbf_xi <- length(Bp_list_aux[[xi]][["paired"]])

      # Add the 2 new basis functions to the Bp_list
      new_bf <- list (
        "id" = c(nbf + 1, nbf + 2),
        "Bp" = cbind(new_pair[[1]], new_pair[[2]]),
        "t"  = knots[i]
        )

      Bp_list_aux[[xi]][["paired"]][[nbf_xi + 1]] <- new_bf

      # Sort basis function by variable | side | knot
      for (v in 1:nX) {
        for (side in c("paired", "right", "left")) {
          if (!is.null(Bp_list_aux[[v]][[side]])) {
            arrangement <- order(sapply(Bp_list_aux[[v]][[side]], "[[", "t"))
            Bp_list_aux[[v]][[side]] <- Bp_list_aux[[v]][[side]][arrangement]
          }
        }
      }

      # Column indexes of the basis function in B matrix by variable (The intercept column is not considered)
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
        B = new_B,
        y = data[, y, drop = F],
        it_list = it_list,
        Bp_list = Bp_list_aux,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0
        )

      # Predictions
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
        B = new_B,
        d = 1,
        knots = sum(nknots_xi),
        xi_degree = NULL
      )

      # compute GRSq
      GRSq <- 1 - GCV / bf_set[[1]][["GCV"]]

      # add the new basis function
      if (is.na(err_min[2]) || err_min[2] == 1) {
        if (xi_degree[2, xi] == 1 && err[1] < err_min[1]) {
          add <- TRUE
        } else if (xi_degree[2, xi] > 1 && (err[1] - err_min[1]) / err_min[1] < - hd_cost) {
          add <- TRUE
        } else {
          add <- FALSE
        }
      } else {
        if (err[1] < err_min[1]) {
          add <- TRUE
        } else {
          add <- FALSE
        }
      }

      if (add) {

        # Model has improved
        improvement <- TRUE
        # Minimum error
        err_min <- err
        # Best B matrix
        best_B <- new_B
        # Best Bp_list
        best_Bp_list <- Bp_list_aux

        # index
        tindex <- which(kn_grid[[xi]] == knots[i])

        # New pair of basis functions
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

        # t
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
#' @description This function generates a vector of knots that can be used to create a new pair of basis functions.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param nX Number of inputs.
#' @param var_idx Index of the input variable that creates the new pair of basis functions.
#' @param minspan Minimum number of observations between two adjacent knots.
#' @param endspan Minimum number of observations before the first and after the final knot.
#' @param kn_list A \code{list} containing the set of selected knots.
#' @param bf A \code{list} with the basis function for the expansion of the model.
#' @param kn_grid Grid of virtual knots to perform ACES.
#'
#' @importFrom dplyr %>%
#'
#' @return Numeric vector with the knots values for the variable.

set_knots <- function (
    data, nX, var_idx, minspan, endspan, kn_list, bf, kn_grid
    ) {

  # Sample size
  N <- nrow(data)

  # Minimum number of observations between two adjacent knots
  if (length(minspan) > 1) {
    sp1 <- minspan[var_idx]
  } else {
    sp1 <- minspan
  }

  # Minimum number of observations before the first and after the final knot.
  if (length(endspan) > 1) {
    sp2 <- endspan[var_idx]
  } else {
    sp2 <- endspan
  }

  # Observations in the space of the basis function.
  nobs_bf <- bf[['Bp']] != 0

  # Minimum span for Friedman's approach
  if (sp1 == - 1) {
    sp1 <- floor(min(N * 0.1, - log2(- (1 / (nX * sum(nobs_bf))) * log(1 - 0.05)) / 2.5))
  }

  # Boolean: possible knots
  kn_idx <- rep(TRUE, length(kn_grid[[var_idx]]))

  # Set to FALSE the first and last "Le" observations
  if (sp2 != 0) {
    set_sp2_idx <- c(1:sp2, (length(kn_idx) - sp2 + 1):length(kn_idx))
  } else {
    set_sp2_idx <- c()
  }

  # Set to FALSE adjacent "L" observations
  if (!is.null(kn_list[[var_idx]])) {
    # Used knots values
    used_kn_val <- sapply(kn_list[[var_idx]], "[[", "t")
    # Used knots indexes
    used_kn_idx <- c()

    for (kn in 1:length(used_kn_val)) {
      used_kn_idx <- c(
        used_kn_idx,
        which(kn_grid[[var_idx]][order(kn_grid[[var_idx]])] == used_kn_val[kn])
      )
    }

    set_sp1_idx <- c()

    for (idx in 1:length(used_kn_idx)) {
      set_sp1_idx <- append(
        set_sp1_idx,
        (used_kn_idx[idx] - sp1):(used_kn_idx[idx] + sp1)
        )
    }

    # Indexes must be greater than 0 and lower than N
    set_sp1_idx <- unique(sort(set_sp1_idx[set_sp1_idx > 0 & set_sp1_idx <= sum(nobs_bf)]))

  } else {
    set_sp1_idx <- c()
  }

  # Set of observations to FALSE based on minspand and endspan
  kn_idx[sort(unique(c(set_sp1_idx, set_sp2_idx)))] <- FALSE

  # Knots
  if (all(kn_idx == FALSE)) {
    knots <- NULL

  } else {
    # Minimum value of observations in BF
    minimum <- min(data[nobs_bf, var_idx])
    # Maximum value of observations in BF
    maximum <- max(data[nobs_bf, var_idx])

    # Exclude observations according to minspan and endspan
    varknots <- kn_grid[[var_idx]][order(kn_grid[[var_idx]])][kn_idx]

    # Knots between minimum and maximum
    knots <- varknots[varknots >= minimum & varknots <= maximum]
  }

  if (length(knots) <= 1) knots <- NULL

  return(knots)
}

#' @title Generate a new pair of Linear Basis Functions
#'
#' @description This function generates two new linear basis functions from a variable and a knot.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param var_idx Index of the input variable that creates the new pair of basis functions.
#' @param knot Knot for creating the new pair of basis functions.
#' @param B A \code{matrix} of basis functions on which the new pair of functions is added.
#'
#' @return A \code{list} with the new pair of basis functions.

create_linear_basis <- function (
    data, var_idx, knot, B
    ) {

  # Create (xi-t)+ and (t-xi)+
  hinge1 <- pmax(0, data[, var_idx] - knot)
  hinge2 <- pmax(0, knot - data[, var_idx])

  # two new basis functions
  bf1 <- B[, 1] * hinge1
  bf2 <- B[, 1] * hinge2

  return(list(bf1, bf2))
}


#' @title Generate the Set of Intervals
#'
#' @description This function determines the boundaries for the intervals and identifies what basis functions activates on each interval.
#'
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#'
#' @return A \code{list} with the boundaries of the intervals and the IDs of the activated basis functions within each interval.

set_intervals <- function (
    Bp_list
    ) {

  # Number of inputs
  nX <- length(Bp_list)

  # Set of intervals by variable.
  it_list <- vector("list", nX)

  for (v in 1:nX) {

    t_vector <- c(0, Inf)
    for (side in c("paired", "right", "left")) {
      if (!is.null(Bp_list[[v]][[side]])) {
        t_vector <- c(t_vector, sapply(Bp_list[[v]][[side]], "[[", "t"))
      }
    }

    # Sort knots vector in ascending order
    t_vector <- sort(t_vector)

    for (j in 1:(length(t_vector) - 1)) {
      # Lower bound in the interval
      lb <- t_vector[[j]]
      # Upper bound in the interval
      ub <- t_vector[[j + 1]]

      # Add lower and upper bounds
      it_list[[v]][[j]] <- list("Lb" = lb, "Ub" = ub)

      # Columns index of the basis function by variable
      Bp_xi <- c(); status <- c()
      for (side in c("paired", "right", "left")) {
        if (!is.null(Bp_list[[v]][[side]])) {
          for (l in 1:length(Bp_list[[v]][[side]])) {
            if (side == "paired") {
              if (Bp_list[[v]][[side]][[l]][["t"]] <= lb) {
                Bp_xi <- c(Bp_xi, Bp_list[[v]][[side]][[l]][["Bp_xi"]][[1]])
                status  <- c(status, "right")

              } else if (Bp_list[[v]][[side]][[l]][["t"]] >= ub) {
                Bp_xi <- c(Bp_xi, Bp_list[[v]][[side]][[l]][["Bp_xi"]][[2]])
                status  <- c(status, "left")
              }

            } else if (side == "right") {
              if (Bp_list[[v]][[side]][[l]][["t"]] <= lb) {
                Bp_xi <- c(Bp_xi, Bp_list[[v]][[side]][[l]][["Bp_xi"]])
                status  <- c(status, "right")
              }

            } else {
              if (Bp_list[[v]][[side]][[l]][["t"]] >= ub) {
                Bp_xi <- c(Bp_xi, Bp_list[[v]][[side]][[l]][["Bp_xi"]])
                status  <- c(status, "left")
              }
            }
          }
        }
      }
      if (is.null(Bp_xi)) {
        it_list[[v]][[j]] <- NULL
      } else {
        it_list[[v]][[j]][["Bp"]] <- data.frame(Bp_xi, status)
      }
    }
  }

  return(it_list)
}
