#' @title Add a New Pair of Basis Functions in Adaptive Constrained Enveloping Splines
#'
#' @description
#'
#' This function adds a pair of basis functions to the Adaptive Constrained Enveloping Splines (ACES) model, selecting the pair that leads to the largest reduction in the lack-of-fit criterion. The function updates the model by incorporating these basis functions, adjusting the knots, and refining the error.
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
#' @param z
#' Column indexes of contextual variables in \code{data}.
#'
#' @param xi_degree
#' A \code{matrix} indicating the degree of each input variable.
#'
#' @param compl_cost
#' A \code{numeric} value specifying the minimum percentage of improvement over the best 1 degree basis function to justify adding a higher-degree basis function.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier.
#'
#' @param dea_scores
#' A \code{matrix} containing DEA-VRS efficiency scores, calculated using an output-oriented radial model. For models with multiple outputs, each column corresponds to the scores for one specific output.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.
#'
#' @param forward_model
#' A \code{list} containing the current state of the forward model, including the matrix of basis functions (\code{B}) and the set of selected basis functions (\code{BF_set}).
#'
#' @param Bp_list
#' A \code{list} containing the basis functions for each input variable, with detailed information about their structure and placement in the model.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param kn_list
#' A \code{list} containing the current set of selected knots for each input variable.
#'
#' @param kn_grid
#' A \code{list} containing the available grid of knots used to construct the basis functions for each input variable.
#'
#' @param span
#' A \code{numeric} vector specifying the minimum number of observations between two adjacent knots (L) and the minimum number of observations before the first and after the final knot (Le).
#'
#' @param err_min
#' A \code{numeric} value specifying the minimum error obtained by the forward algorithm in the current iteration.
#'
#' @param var_imp
#' A \code{matrix} tracking the best residual reduction for each variable (columns) in each iteration (rows).
#'
#' @param quick_aces
#' A \code{logical} indicating whether to use the fast version of ACES.
#'
#' @return
#'
#' An updated \code{list} containing:
#' \itemize{
#'   \item \code{B}: The updated matrix of basis functions.
#'   \item \code{BF_set}: The updated set of basis functions included in the model.
#'   \item \code{kn_list}: The updated list of selected knots for each variable.
#'   \item \code{err_min}: The updated minimum error achieved in the current iteration.
#'   \item \code{var_imp}; The updated matrix of variable importance
#' }

add_basis_function <- function (
    data,
    x,
    y,
    z,
    xi_degree,
    compl_cost,
    model_type,
    dea_scores,
    metric,
    forward_model,
    Bp_list,
    shape,
    kn_list,
    kn_grid,
    span,
    err_min,
    var_imp,
    quick_aces
    ) {

  # number of inputs
  nX <- ncol(data) - length(y) - length(z)

  # number of contextual variables
  nZ <- ncol(data) - length(x) - length(y)

  # set of basis functions
  bf_set <- forward_model[["bf_set"]]

  # number of basis functions
  nbf <- length(bf_set)

  # signal to indicate improvement in the fitting
  improvement <- FALSE

  # initial error
  err_ini <- err_min

  # first basis function for the expansion (additive model)
  bf <- bf_set[[1]]

  # compute weighted mean reduction in error for each variable
  if (quick_aces && nrow(var_imp) > 1) {

    iters <- nrow(var_imp) - 1

    # weighted mean reduction for in error each variable: exponential decay
    weights <- 0.95 ^ (iters:1)
    weighted_reduction <- colSums(var_imp[1:iters, , drop = FALSE] * weights) / sum(weights)

    # normalize importance
    normalized_importance <- weighted_reduction / max(weighted_reduction)

  }

  for (xi in c(x, z)) {

    # ================================ #
    # Create the set of eligible knots #
    # ================================ #

    if (quick_aces) {

      # remove variables based on Spearman and Kendall correlation
      if (var_imp[1, xi] == - 1) next

      # get knots: do not apply minimum and end span
      knots <- kn_grid[[xi]]

      # reduce number of knots based on importance
      if (nrow(var_imp) > 1) {

        # threshold for important variables
        threshold <- sort(normalized_importance, decreasing = TRUE)[2] / 2

        # set a deactivation threshold
        is_active <- normalized_importance[xi] >= threshold

        if (is_active) {

          # proportion of knots to be selected
          kn_prop <- normalized_importance[xi]

          # sample of knots based on importance
          kn_indx <- sample(length(knots), size = round(kn_prop * length(knots)))

          # reduced set of knots
          knots <- knots[kn_indx]

        } else {

          next

        }
      }

    } else {

      knots <- set_knots (
        data = data,
        nX = nX + nZ,
        var_idx = xi,
        minspan = span[1],
        endspan = span[2],
        kn_list = kn_list,
        bf = bf,
        kn_grid = kn_grid
      )

    }

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
        knot = knots[i]
        )

      # ============== #
      # Update Bp_list #
      # ============== #

      # number of paired basis functions for the input "xi"
      nbf_xi <- length(Bp_list_aux[[xi]][["paired"]])

      # add 2 new basis functions to the Bp_list
      new_bf <- list (
        "id" = c(nbf + 1, nbf + 2),
        "Bp" = cbind(new_pair[[1]], new_pair[[2]]),
        "t"  = knots[i]
        )

      Bp_list_aux[[xi]][["paired"]][[nbf_xi + 1]] <- new_bf

      # sort basis function by variable, side and knot
      Bp_list_aux[[xi]][["paired"]] <- Bp_list_aux[[xi]][["paired"]][
        order(sapply(Bp_list_aux[[xi]][["paired"]], "[[", "t"))
      ]

      # column indexes of BFs in B matrix by variable (ignore intercept)
      for (v in x) {

        Bp_xi <- 1

        if (!is.null(Bp_list_aux[[v]][["paired"]])) {

          for (l in 1:length(Bp_list_aux[[v]][["paired"]])) {

            Bp_list_aux[[v]][["paired"]][[l]][["Bp_xi"]] <- c(Bp_xi, Bp_xi + 1)
            Bp_xi <- Bp_xi + 2

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

      new_B <- matrix(1, nrow = nrow(data))

      for (v in x) {

        if (!is.null(Bp_list_aux[[v]][["paired"]])) {

          for (l in 1:length(Bp_list_aux[[v]][["paired"]])) {

            bf_jp <- Bp_list_aux[[v]][["paired"]][[l]][["Bp"]]
            new_B <- cbind(new_B, bf_jp)

          }
        }
      }

      if (!is.null(z)) {

        # =============== #
        # Update Z matrix #
        # =============== #

        new_Z <- matrix(1, nrow = nrow(data))

        for (v in z) {

          if (!is.null(Bp_list_aux[[v]][["paired"]])) {

            for (l in 1:length(Bp_list_aux[[v]][["paired"]])) {

              bf_jp <- Bp_list_aux[[v]][["paired"]][[l]][["Bp"]]
              new_Z <- cbind(new_Z, bf_jp)

            }
          }
        }

        if (ncol(new_Z) > 1) {

          new_Z <- new_Z[, 2:ncol(new_Z)]

        } else {

          new_Z <- NULL

        }
      }

      # ===================== #
      # Estimate coefficients #
      # ===================== #

      coefs <- estimate_coefficients (
        model_type = model_type,
        B = new_B,
        Z = new_Z,
        y_obs = data[, y, drop = F],
        dea_scores = dea_scores,
        it_list = it_list,
        Bp_list = Bp_list_aux,
        shape = shape
        )

      # predictions
      y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

      for (out in 1:length(y)) {
        y_hat[, out] <- new_B %*% coefs[, out]
      }

      # lack-of-fit measure
      err <- err_metric (
        y_obs = data[, y, drop = F],
        y_hat = y_hat[drop = F],
        metric = metric,
        weight = 1 / dea_scores
        )

      # variable
      err[2] <- xi_degree[2, xi]

      # update variable importance
      reduction <- round(1 - err[1] / err_ini[1], 2)
      best_reduction <- var_imp[nrow(var_imp), xi]

      if (reduction > best_reduction) {

        var_imp[nrow(var_imp), xi] <- reduction

      }

      # check if the best basis function should be checked
      add <- FALSE

      if (is.na(err_min[2]) || err_min[2] == 1) {

        if (xi_degree[2, xi] == 1) {

          add <- err[1] < err_min[1]

        } else if (xi_degree[2, xi] > 1) {

          required_reduction <- err_min[1] * (1 - compl_cost)

          add <- err[1] < required_reduction

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

        # knot index
        tindex <- which(kn_grid[[xi]] == knots[i])

        # new pair of basis functions
        bf1 <- bf2 <- bf

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
        bf1[['xi']] <- bf2[['xi']] <- ifelse(all(bf[['xi']] == - 1), xi, c(bf[['xi']], xi))

        # knot
        bf1[["t"]] <- bf2[["t"]] <- knots[i]

        # R
        bf1[['R']] <- bf2[['R']] <- err_min[1]

        # coefficients
        bf1[['coefs']] <- bf2[['coefs']] <- coefs

      }
    }
  }

  if (improvement) {

    # append new basis functions
    bf_set <- append(bf_set, list(bf1))
    bf_set <- append(bf_set, list(bf2))

    # update the list of knots for each variable
    len_knt <- length(bf1[["xi"]])
    var_pos <- bf1[["xi"]][len_knt]

    kn_list[[var_pos]] <- append (
      kn_list[[bf1[['xi']][len_knt]]],
      list(list(t = bf1[["t"]], index = tindex))
      )

    return(list(best_B, bf_set, kn_list, best_Bp_list, err_min, var_imp))

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

#' @title Generate a Pair of Linear Basis Functions
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
#' @return
#' A \code{list} with the new pair of basis functions.

create_linear_basis <- function (
    data,
    var_idx,
    knot
    ) {

  # create (xi-t)+ and (t-xi)+
  hinge1 <- pmax(0, data[, var_idx] - knot)
  hinge2 <- pmax(0, knot - data[, var_idx])

  return(list(hinge1, hinge2))

}
