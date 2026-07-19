#' @title Add a Pair of ACES Basis Functions
#'
#' @description
#'
#' Adds the eligible pair of basis functions that gives the largest reduction in
#' lack of fit. Candidate knots and variables must satisfy the span, interaction,
#' and shape restrictions. In Quick ACES mode, the candidate set is reduced
#' before the best pair is selected.
#'
#' @param data
#' A matrix containing the prepared model variables.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param xi_degree
#' A matrix that records the degree of each input in every expanded variable.
#'
#' @param inter_cost
#' Minimum relative improvement over the best first-degree basis function
#' required to add a higher-degree basis function.
#'
#' @param dea_scores
#' A matrix of output-oriented DEA-VRS scores, with one column per output.
#'
#' @param fdh_scores
#' A matrix of output-oriented FDH scores, with one column per output.
#'
#' @param metric
#' Character string specifying the lack-of-fit measure.
#'
#' @param forward_model
#' Current forward model, including \code{B} and \code{BF_set}.
#'
#' @param Bp_list
#' Current basis functions, grouped by input.
#'
#' @param shape
#' A list with logical elements \code{mono} and \code{conc}.
#'
#' @param kn_list
#' Selected knots, grouped by input.
#'
#' @param kn_grid
#' Available knot candidates, grouped by input.
#'
#' @param span
#' A numeric vector containing \code{minspan} and \code{endspan}.
#'
#' @param err_min
#' Best error found in the current iteration.
#'
#' @param var_imp
#' Matrix whose rows contain the best relative error reduction found for each
#' prepared input in each forward iteration. Quick ACES uses its recent history
#' to reduce the next candidate search.
#'
#' @param quick_keep
#' Logical vector indicating which prepared inputs passed the Quick ACES
#' correlation screen.
#'
#' @param quick_aces
#' If \code{TRUE}, use Quick ACES to reduce the candidate search.
#'
#' @return
#'
#' An updated \code{list} containing:
#' \itemize{
#'   \item{\code{B}: Updated basis matrix.}
#'   \item{\code{BF_set}: Updated set of selected basis functions.}
#'   \item{\code{kn_list}: Updated knots.}
#'   \item{\code{err_min}: Updated best error.}
#'   \item{\code{var_imp}: Updated error reductions by iteration and input.}
#' }

add_basis_function <- function(
  data,
  x,
  y,
  xi_degree,
  inter_cost,
  dea_scores,
  fdh_scores,
  metric,
  forward_model,
  Bp_list,
  shape,
  kn_list,
  kn_grid,
  span,
  err_min,
  var_imp,
  quick_keep,
  quick_aces
) {
  # number of inputs
  nX <- ncol(data) - length(y)

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

    # weighted mean reduction for error in each variable: exponential decay
    weights <- 0.95^(iters:1)
    reduction_history <- var_imp[seq_len(iters), , drop = FALSE]
    reduction_history[is.na(reduction_history)] <- 0
    weighted_reduction <- colSums(reduction_history * weights) / sum(weights)

    # normalize importance; if no input improved, keep the full search active
    max_reduction <- suppressWarnings(max(weighted_reduction, na.rm = TRUE))

    if (is.finite(max_reduction) && max_reduction > 0) {
      normalized_importance <- pmax(0, weighted_reduction / max_reduction)
    } else {
      normalized_importance <- rep(1, length(weighted_reduction))
    }

    normalized_importance[!quick_keep] <- 0
  }

  for (xi in x) {
    # ================================ #
    # Create the set of eligible knots #
    # ================================ #

    if (quick_aces) {
      # remove variables based on Spearman and Kendall correlation
      if (!quick_keep[xi]) next

      # get knots: do not apply minimum and end span
      knots <- kn_grid[[xi]]

      # reduce number of knots based on importance
      if (nrow(var_imp) > 1) {
        # threshold for important variables
        if (length(x) >= 2) {
          threshold <- sort(normalized_importance, decreasing = TRUE)[2] / 2
        } else {
          threshold <- sort(normalized_importance, decreasing = TRUE)[1] / 2
        }

        # set a deactivation threshold
        is_active <- normalized_importance[xi] >= threshold

        if (is_active) {
          # proportion of knots to be selected
          kn_prop <- normalized_importance[xi]

          # sample of knots based on importance
          sample_size <- max(1L, round(kn_prop * length(knots)))
          kn_indx <- sample(length(knots), size = sample_size)

          # reduced set of knots
          knots <- knots[kn_indx]
        } else {
          next
        }
      }
    } else {
      knots <- set_knots(
        data = data,
        nX = nX,
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
    knots <- sort(knots)

    for (i in 1:length(knots)) {
      Bp_list_aux <- Bp_list

      # ================== #
      # Create two new bfs #
      # ================== #

      new_pair <- create_linear_basis(
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
      new_bf <- list(
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

      it_list <- set_intervals(
        Bp_list = Bp_list_aux
      )

      # =============== #
      # Update B matrix #
      # =============== #

      # compute insertion position using ORIGINAL Bp_list (before adding new knot)
      col_offset <- 1L
      for (v in x[x < xi]) {
        if (!is.null(Bp_list[[v]][["paired"]])) {
          col_offset <- col_offset + 2L * length(Bp_list[[v]][["paired"]])
        }
      }

      if (!is.null(Bp_list[[xi]][["paired"]])) {
        existing_knots_xi <- sapply(Bp_list[[xi]][["paired"]], "[[", "t")
        knots_before <- sum(existing_knots_xi < knots[i])
      } else {
        knots_before <- 0L
      }

      col_insert <- col_offset + knots_before * 2L

      existing_B <- forward_model[["B"]]

      if (col_insert < ncol(existing_B)) {
        new_B <- cbind(
          existing_B[, 1:col_insert, drop = FALSE],
          new_pair[[1]],
          new_pair[[2]],
          existing_B[, (col_insert + 1):ncol(existing_B), drop = FALSE]
        )
      } else {
        new_B <- cbind(existing_B, new_pair[[1]], new_pair[[2]])

      }

      # ===================== #
      # Estimate coefficients #
      # ===================== #

      coefs <- estimate_coefficients(
        B = new_B,
        y_obs = data[, y, drop = F],
        dea_scores = dea_scores,
        fdh_scores = fdh_scores,
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
      err <- err_metric(
        y_obs = data[, y, drop = F],
        y_hat = y_hat[drop = F],
        metric = metric,
        weight = 1 / dea_scores
      )

      # variable
      err[2] <- xi_degree[2, xi]

      # Keep the best counterfactual reduction for this variable in the current
      # iteration, whether or not this variable supplies the selected split.
      reduction <- if (
        is.finite(err_ini[1]) && err_ini[1] > 0 && is.finite(err[1])
      ) {
        max(0, 1 - err[1] / err_ini[1])
      } else {
        NA_real_
      }

      best_reduction <- var_imp[nrow(var_imp), xi]

      if (!is.na(reduction) &&
          (is.na(best_reduction) || reduction > best_reduction)) {
        var_imp[nrow(var_imp), xi] <- reduction
      }

      # check if the best basis function should be checked
      add <- FALSE

      if (is.na(err_min[2]) || err_min[2] == 1) {
        if (xi_degree[2, xi] == 1) {
          add <- err[1] < err_min[1]
        } else if (xi_degree[2, xi] > 1) {
          required_reduction <- err_min[1] * (1 - inter_cost)

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
        bf1[["Bp"]] <- new_pair[[1]]
        bf2[["Bp"]] <- new_pair[[2]]

        # xi
        bf1[["xi"]] <- bf2[["xi"]] <- ifelse(all(bf[["xi"]] == -1), xi, c(bf[["xi"]], xi))

        # knot
        bf1[["t"]] <- bf2[["t"]] <- knots[i]

        # R
        bf1[["LOF"]] <- bf2[["LOF"]] <- err_min[1]

        # coefficients
        bf1[["coefs"]] <- bf2[["coefs"]] <- coefs
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

    kn_list[[var_pos]] <- append(
      kn_list[[bf1[["xi"]][len_knt]]],
      list(list(t = bf1[["t"]], index = tindex))
    )

    return(list(best_B, bf_set, kn_list, best_Bp_list, err_min, var_imp))
  } else {
    return(improvement)
  }
}

#' @title Find Eligible Knots
#'
#' @description
#' Finds knot candidates that satisfy the spacing rules during the forward step.
#'
#' @param data
#' A matrix containing the prepared model variables.
#'
#' @param nX
#' Number of inputs.
#'
#' @param var_idx
#' Index of the input used to create the new basis functions.
#'
#' @param minspan
#' Minimum number of observations between two adjacent knots.
#'
#' @param endspan
#' Minimum number of observations before the first knot and after the last knot.
#'
#' @param kn_list
#' Selected knots, grouped by input.
#'
#' @param bf
#' Parent basis function to expand.
#'
#' @param kn_grid
#' Available knot candidates, grouped by input.
#'
#' @return
#' A numeric vector of eligible knots.

set_knots <- function(
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
  nobs_bf <- bf[["Bp"]] != 0

  # minimum span for Friedman's approach
  if (sp1 == -1) {
    sp1 <- floor(min(N * 0.1, -log2(-(1 / (nX * sum(nobs_bf))) * log(1 - 0.05)) / 2.5))
  }

  # boolean: possible knots
  kn_idx <- rep(TRUE, length(kn_grid[[var_idx]]))

  # Set to FALSE the first and last "Le" observations
  if (sp2 != 0) {
    set_sp2_idx <- c(1:sp2, (length(kn_idx) - sp2 + 1):length(kn_idx))
  } else {
    set_sp2_idx <- c()
  }

  # sorted grid of knots for the variable
  sorted_grid <- kn_grid[[var_idx]][order(kn_grid[[var_idx]])]

  # Set to FALSE adjacent "L" observations
  if (!is.null(kn_list[[var_idx]])) {
    # used knots values
    used_kn_val <- sapply(kn_list[[var_idx]], "[[", "t")

    # used knots indexes
    used_kn_idx <- which(sorted_grid %in% used_kn_val)

    # indexes within "L" observations of a used knot
    set_sp1_idx <- as.vector(outer(used_kn_idx, seq(-sp1, sp1), "+"))

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
    varknots <- sorted_grid[kn_idx]

    # knots between minimum and maximum
    knots <- varknots[varknots >= minimum & varknots <= maximum]
  }

  if (length(knots) <= 1) knots <- NULL

  return(knots)
}

#' @title Create Linear Basis Functions
#'
#' @description
#' Creates the two hinge functions on either side of a knot.
#'
#' @param data
#' A matrix containing the prepared model variables.
#'
#' @param var_idx
#' Index of the input used to create the basis functions.
#'
#' @param knot
#' Knot location.
#'
#' @return
#' A list containing the right and left hinge vectors.

create_linear_basis <- function(
  data,
  var_idx,
  knot
) {
  # create (xi-t)+ and (t-xi)+
  hinge1 <- pmax(0, data[, var_idx] - knot)
  hinge2 <- pmax(0, knot - data[, var_idx])

  return(list(hinge1, hinge2))
}
