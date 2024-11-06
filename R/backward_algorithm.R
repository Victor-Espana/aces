#' @title Backward Algorithm for Adaptive Constrained Enveloping Splines
#'
#' @description
#'
#' This function implements the Backward Algorithm in Adaptive Constrained Enveloping Splines. It creates a portfolio of sub-models by iteratively removing basis functions one by one.
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
#' A \code{list} containing the information relative to the forward step of the ACES model.
#'
#' @param Bp_list
#' A \code{list} containing the current set of basis functions for each input variable.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param kn_penalty
#' An \code{integer} specifying the Generalized Cross Validation (GCV) penalty per knot.
#'
#' @return
#'
#' A \code{list} containing a portfolio of Adaptive Constrained Enveloping Splines sub-models.

aces_pruning <- function (
    data,
    x,
    y,
    z,
    xi_degree,
    model_type,
    dea_scores,
    metric,
    forward_model,
    Bp_list,
    shape,
    kn_penalty
    ) {

  # number of inputs
  nX <- length(x)

  # number of terms to drop
  terms <- length(forward_model[["basis"]]) - 1

  # knots sorted by variable and value
  knots <- update_knot_list (
    Bp_list = Bp_list
  )

  # B matrix
  B_mat <- forward_model[["Bmatx"]]

  # vector of coefficients from the forward step of ACES
  coefs <- forward_model[["basis"]][[terms + 1]][["coefs"]]

  # measures of performance
  LOF <- forward_model[["basis"]][[terms + 1]][["LOF"]]

  GCV <- compute_gcv (
    y_obs = data[, y, drop = F],
    LOF = LOF,
    n_bf = ncol(B_mat),
    kn_penalty = kn_penalty,
    kn_list = knots,
    xi_degree = xi_degree
  )

  # set of ACES sub-models with p terms
  models <- list (
    list (
      "id" = length(forward_model[["basis"]]),
      "B" = B_mat,
      "LOF" = LOF,
      "GCV" = GCV,
      "t" = knots,
      "coefs" = coefs,
      "Bp_list" = Bp_list
      )
    )

  # ======= #
  # Pruning #
  # ======= #

  # B matrix
  B <- forward_model[["Bmatx"]]

  # set of basis functions
  basis <- forward_model[["basis"]]

  while (terms > 0) {

    # error
    err <- Inf

    # basis functions to be removed excluding:
    # - B1(X) = 1
    # - basis function when appearing alone in an interval: excepting right-hand extreme interval
    bf_to_drop <- basis_to_drop (
      Bp_list = Bp_list,
      conc = shape[["conc"]],
      nX = nX
      )

    if (nrow(bf_to_drop) == 0) break

    # mapping between basis function ID and matrix column index by variable
    id_Bpxi <- mapping_basis (
      Bp_list = Bp_list,
      nX = nX
      )

    # add id column to basis function to drop
    bf_to_drop <- merge(bf_to_drop, id_Bpxi, by = c("xi", "Bp_xi"))

    # shuffle rows
    bf_to_drop <- bf_to_drop[sample(1:nrow(bf_to_drop)), ]

    for (bf in 1:nrow(bf_to_drop)) {

      # ======================= #
      # Update the set of basis #
      # ======================= #

      new_basis <- basis

      # id to drop the basis function from new_basis
      id <- bf_to_drop[bf, "id"]

      # select the basis function to be dropped
      dropTerm <- sapply(new_basis, function(x) x[["id"]] == id)

      # number of 0's in the basis function to be dropped (in case of same GCV)
      basis_0s <- sum(new_basis[dropTerm][[1]][["Bp"]] == 0)

      # if a basis function is paired -> the sibling basis function is set as unpaired
      if (new_basis[dropTerm][[1]][["status"]] == "paired") {

        # knot (t) and variable (v)
        check1 <- sapply(new_basis, "[[", "t") == new_basis[dropTerm][[1]][["t"]]
        check2 <- sapply(new_basis, "[[", "xi") == new_basis[dropTerm][[1]][["xi"]]

        # IDs for all basis functions in new_basis
        basis_ids <- sapply(new_basis, "[[", "id")

        # ID of the paired basis function
        paired_id <- basis_ids[check1 & check2]

        # ID of the sibling basis function
        sibling_id <- paired_id[!(id == paired_id)]

        # update new_basis
        new_basis[lapply(new_basis, "[[", "id") == sibling_id][[1]][["status"]] <- "unpaired"

      }

      # drop the selected basis function
      new_basis[dropTerm] <- NULL

      if (length(new_basis) == 1) {

        # bp_list
        new_Bp_list <- vector("list", nX)

        # B matrix
        new_B <- matrix(1, nrow = nrow(data))

        # vector of coefficients
        coefs <- unname(apply(data[, y, drop = FALSE], 2, max))

        # prediction
        y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

        for (out in 1:length(y)) {
          y_hat[, out] <- coefs[out]
        }

        LOF <- err_metric (
          y_obs = data[, y, drop = F],
          y_hat = y_hat,
          metric = metric,
          weight = 1 / dea_scores
        )

        GCV <- compute_gcv (
          y_obs = data[, y, drop = F],
          LOF = LOF,
          n_bf = ncol(new_B),
          kn_penalty = kn_penalty,
          kn_list = NULL,
          xi_degree = NULL
          )

      } else {

        # ============== #
        # Update Bp_list #
        # ============== #

        new_Bp_list <- update_Bp_list (
          basis = new_basis,
          nX = nX
          )

        # ============== #
        # Update It.list #
        # ============== #

        new_it_list <- set_intervals (
          Bp_list = new_Bp_list
          )

        # =============== #
        # Update B matrix #
        # =============== #

        new_B <- matrix(1, nrow = nrow(data))

        for (v in x) {

          for (side in c("paired", "right", "left")) {

            for (l in 1:length(new_Bp_list[[v]][[side]])) {

              bf_jp <- new_Bp_list[[v]][[side]][[l]][["Bp"]]
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

        try_estimation <- tryCatch (
          {
            coefs <- estimate_coefficients (
              model_type = model_type,
              B = new_B,
              Z = new_Z,
              y_obs = data[, y, drop = F],
              dea_scores = dea_scores,
              it_list = new_it_list,
              Bp_list = new_Bp_list,
              shape = shape
            )
          },

          error = function(e) {
            return("error")
          }
        )

        if (class(try_estimation[[1]]) == "character") next

        # ================ #
        # Update the knots #
        # ================ #

        knots <- update_knot_list (
          Bp_list = new_Bp_list
        )

        # =================== #
        # Compute error & GCV #
        # =================== #

        # prediction
        y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

        for (out in 1:length(y)) {
          y_hat[, out] <- new_B %*% coefs[, out, drop = F]
        }

        LOF <- err_metric (
          y_obs = data[, y, drop = F],
          y_hat = y_hat,
          metric = metric,
          weight = 1 / dea_scores
        )

        GCV <- compute_gcv (
          y_obs = data[, y, drop = F],
          LOF = LOF,
          n_bf = ncol(new_B),
          kn_penalty = kn_penalty,
          kn_list = knots,
          xi_degree = xi_degree
          )
      }

      if (LOF <= err) {
        best_basis <- new_basis
        best_Bp_list <- new_Bp_list
        best_B <- new_B
        err <- LOF
        best_gcv <- GCV
        best_knots <- knots
        best_coefs <- coefs
      }
    }

    # selected basis after removing
    basis <- best_basis

    # selected Bp_list after removing
    Bp_list <- best_Bp_list

    # Create p-th sub model
    p_model <- list (
      "id" = ncol(best_B),
      "B" = best_B,
      "LOF" = err,
      "GCV" = best_gcv,
      "t" = best_knots,
      "coefs" = best_coefs,
      "Bp_list" = best_Bp_list
    )

    # best model of p terms
    models <- append(models, list(p_model))

    # update number of terms
    terms <- terms - 1
  }

  return(models)
}

#' @title Compute Generalized Cross-Validation (GCV)
#'
#' @description
#'
#' This function computes the generalized cross-validation for the backward algorithm in Adaptive Constrained Enveloping Splines (ACES).
#'
#' @param y_obs
#' Vector of observed data.
#'
#' @param LOF
#' A \code{numeric} value representing the lack-of-fit of the model. This term quantifies how well the model fits the data.
#'
#' @param n_bf
#' A \code{integer} specifying the number of basis functions used in the model.
#'
#' @param kn_penalty
#' A \code{numeric} value representing the penalty applied for each knot in the model.
#'
#' @param kn_list
#'A \code{list} containing the set of knots used in the model.
#'
#' @param xi_degree
#' A \code{matrix} indicating the degree of each input variable.
#'
#' @return
#'
#'A \code{numeric} value representing the computed Generalized Cross-Validation (GCV) score.

compute_gcv <- function (
    y_obs,
    LOF,
    n_bf,
    kn_penalty,
    kn_list,
    xi_degree
    ) {

  # number of outputs
  nY <- ncol(y_obs)

  # number of observations
  N <- nrow(y_obs)

  # number of knots
  if (is.list(kn_list)) {

    # transform knot list to a data.frame
    knots_df <- do.call(rbind.data.frame, kn_list)

    # knot dimension
    knots_df$dim_knot <- xi_degree[2, match(knots_df$xi, xi_degree[1, ])]

    if (nrow(knots_df) == 0) {

      kn_num <- 0

    } else {

      # drop duplicated
      knots_df <- knots_df[!duplicated(knots_df[c("xi", "t")]), ]

      # number of knots
      kn_num <- sum(knots_df$dim_knot)
    }

  } else {

    if (is.null(kn_list)) kn_num <- 0

    kn_num <- sum(kn_list * xi_degree[2, ])

  }

  # cost-complexity measure
  ccm_val <- (1 - ((n_bf + kn_penalty * kn_num) / (nY * N))) ^ 2

  # generalized cross-validation
  GCV <- ifelse(ccm_val > nY * N, Inf, LOF / ccm_val)

  return(GCV)

}

#' @title Basis Functions to be Dropped
#'
#' @description
#' This function determines the basis function to be dropped from a set of intervals.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param conc
#' A \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @param nX
#' Number of inputs.
#'
#' @return
#' A \code{data.frame} indicating the candidate basis function to be dropped for each input.

basis_to_drop <- function (
    Bp_list,
    conc,
    nX
    ) {

  # set of intervals from Bp_list
  it_list <- set_intervals (
    Bp_list = Bp_list
    )

  # Initialize a data.frame for basis functions to drop
  bf_to_drop <- data.frame()

  for (v in 1:nX) {

    if (length(it_list[[v]]) != 0) {

      for (it in 1:length(it_list[[v]])) {

        interval <- it_list[[v]][[it]]

        if (conc) {

          # notice that concavity is lost in absence of a basis for a middle interval.
          # then, empty intervals are only available in the top right-side (last interval)
          # for example, if there is only one basis function in the first interval,
          # it can't be dropped.

          # check if it is not the last interval
          is_not_last_interval <- it != length(it_list[[v]])

          # check if there is only 1 basis function
          has_single_basis_function <- length(interval[["Bp"]][["Bp_xi"]]) == 1

          if (is_not_last_interval && has_single_basis_function) next

        }

        # add basis functions to drop
        for (j in it_list[[v]][[it]][["Bp"]][["Bp_xi"]]) {
          bf_to_drop <- rbind(bf_to_drop, c(v, j))
        }
      }
    }
  }

  if (nrow(bf_to_drop) != 0) {

    # drop duplicated values
    colnames(bf_to_drop) <- c("xi", "Bp_xi")
    bf_to_drop <- unique(bf_to_drop[order(bf_to_drop$xi, bf_to_drop$Bp_xi), ])

  }

  return(bf_to_drop)

}

#' @title Mapping between Basis Function ID and Matrix Column Index by Variable
#'
#' @description
#' This function creates a mapping between the ID of a basis function and its corresponding column index in the basis function matrix (B matrix) by variable.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param nX
#' Number of inputs.
#'
#' @return
#' A \code{matrix} that establishes the connection between the basis function ID and its column index in the basis function matrix by variable.

mapping_basis <- function (
    Bp_list,
    nX
    ) {

  # initialize matrix to store mapping
  id_Bp_xi <- data.frame (
    xi = NA,
    id = NA,
    Bp_xi = NA
    )

  k <- 1

  for (v in 1:nX) {
    for (side in c("paired", "right", "left")) {
      if (!is.null(Bp_list[[v]][[side]])) {
        for (l in 1:length(Bp_list[[v]][[side]])) {

          # select a Bp
          Bp <- Bp_list[[v]][[side]][[l]]

          # save the id for the selected Bp
          Bp_id <- Bp[["id"]]

          # save the matrix column index by the v-th variable
          Bp_xi <- Bp[["Bp_xi"]]

          if (side == "paired") {

            # id and matrix column connection
            id_Bp_xi[k:(k + 1), ] <- cbind(rep(v, 2), Bp_id, Bp_xi)

            k <- k + 2

          } else {

            id_Bp_xi[k, ] <- c(v, Bp_id, Bp_xi)

            k <- k + 1

          }
        }
      }
    }
  }

  return(id_Bp_xi)

}

#' @title Update knot_list
#'
#' @description
#' This function updates the knot_list during the backward algorithm.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @return
#' A \code{list} with the knots updated.

update_knot_list <- function (
    Bp_list
    ) {

  # number of inputs
  nX <- length(Bp_list)

  # initialize the knot list
  knots <- list()
  k <- 1

  for (v in 1:nX) {

    for (side in c("paired", "right", "left")) {

      if (!is.null(Bp_list[[v]][[side]])) {

        for (l in 1:length(Bp_list[[v]][[side]])) {

          # Bp_element
          Bp <- Bp_list[[v]][[side]][[l]]

          if (side == "paired") {

            knots[[k]] <- list (
              "xi" = v,
              "t"  = Bp[["t"]],
              "side" = "R",
              "status" = "paired"
            )

            knots[[k + 1]] <- list (
              "xi" = v,
              "t"  = Bp[["t"]],
              "side" = "L",
              "status" = "paired"
            )

            k <- k + 2

          } else if (side == "right") {

            knots[[k]] <- list (
              "xi" = v,
              "t"  = Bp[["t"]],
              "side" = "R",
              "status" = "unpaired"
            )

            k <- k + 1

          } else {

            knots[[k]] <- list (
              "xi" = v,
              "t"  = Bp[["t"]],
              "side" = "L",
              "status" = "unpaired"
            )

            k <- k + 1

          }
        }
      }
    }
  }

  return(knots)

}

#' @title Update Bp_list
#'
#' @description
#' This function updates the Bp_list during the backward algorithm.
#'
#' @param basis
#' A \code{list} containing the set of basis functions.
#'
#' @param nX
#' Number of inputs.
#'
#' @return
#' A \code{list} with the basis functions updated.

update_Bp_list <- function (
    basis, nX
    ) {

  # initialize the list of new basis functions
  new_Bp_list <- vector("list", nX)

  for (v in 1:nX) {

    new_Bp_list[[v]] <- list (
      "paired" = NULL,
      "right"  = NULL,
      "left"   = NULL
      )

  }

  # update the list of basis functions
  for (i in 2:length(basis)) {

    # id
    id <- basis[[i]][["id"]]

    # basis function
    Bp <- basis[[i]][["Bp"]]

    # variable
    xi <- basis[[i]][["xi"]]

    # knot
    kn <- basis[[i]][["t"]]

    # status
    status <- basis[[i]][["status"]]

    # side
    side <- basis[[i]][["side"]]

    if (status == "unpaired") {

      side_key <- ifelse(side == "R", "right", "left")

      new_Bp_list[[xi]][[side_key]] <- append (
        new_Bp_list[[xi]][[side_key]],
        list (
          list (
            "Bp" = Bp,
            "t"  = kn,
            "id" = id
            )
          )
        )

    } else {

      if (kn %in% sapply(new_Bp_list[[xi]][[status]], "[[", "t")) next

      new_Bp_list[[xi]][["paired"]] <- append (
        new_Bp_list[[xi]][["paired"]],
        list (
          list (
            "Bp" = cbind(Bp, basis[[i + 1]][["Bp"]]),
            "t"  = kn,
            "id" = c(id, id + 1)
          )
        )
      )
    }
  }

  # sort basis function by variable | side | knot
  for (v in 1:nX) {
    for (side in c("paired", "right", "left")) {
      if (!is.null(new_Bp_list[[v]][[side]])) {
        arrangement <- order(sapply(new_Bp_list[[v]][[side]], "[[", "t"))
        new_Bp_list[[v]][[side]] <- new_Bp_list[[v]][[side]][arrangement]
      }
    }
  }


  # number of column in matrix for ensuring concavity and monotonicity
  for (v in 1:nX) {

    Bp_xi <- 1

    for (side in c("paired", "right", "left")) {

      if (!is.null(new_Bp_list[[v]][[side]])) {

        for (l in 1:length(new_Bp_list[[v]][[side]])) {

          if (side == "paired") {

            new_Bp_list[[v]][[side]][[l]][["Bp_xi"]] <- c(Bp_xi, Bp_xi + 1)
            Bp_xi <- Bp_xi + 2

          } else {

            new_Bp_list[[v]][[side]][[l]][["Bp_xi"]] <- Bp_xi
            Bp_xi <- Bp_xi + 1

          }
        }
      }
    }
  }

  return(new_Bp_list)

}
