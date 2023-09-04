#' @title Backward Algorithm for Adaptive Constrained Enveloping Splines
#'
#' @description This function implements the Backward Algorithm in Adaptive Constrained Enveloping Splines to create a portfolio of sub-models by iteratively removing basis functions one by one.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param xi_degree Matrix indicating the degree of each input variable.
#' @param y Column indexes of output variables in \code{data}.
#' @param metric Lack-of-fit criterion to evaluate the model performance.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#' @param forward_model \code{list} containing the Forward ACES model.
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#' @param d An \code{integer} specifying the Generalized Cross Validation (GCV) penalty per knot.
#'
#' @return A \code{list} containing a portfolio of Adaptive Constrained Enveloping Splines sub-models.

aces_pruning <- function(
    data, xi_degree, y, metric, monotonicity, concavity, x0_y0, forward_model, Bp_list, d
    ) {

  # number of inputs
  nX <- length(Bp_list)

  # number of terms to drop
  terms <- length(forward_model[["basis"]]) - 1

  # knots sorted by variable and value
  knots <- update_knot_list (
    Bp_list = Bp_list
  )

  # B matrix
  B_mat <- forward_model[["Bmatx"]]

  # vector of coefficients from the forward ACES
  coefs <- forward_model[["basis"]][[terms + 1]][["coefs"]]

  # measures of performance
  LOF <- forward_model[["basis"]][[terms + 1]][["err"]]
  GCV <- compute_gcv (
    y_obs = data[, y, drop = F],
    y_hat = B_mat %*% coefs,
    metric = metric,
    B = B_mat,
    d = d,
    knots = knots,
    xi_degree = xi_degree
  )

  # set of ACES sub-models with p terms
  models <- list (
    list (
      "id"  = length(forward_model[["basis"]]),
      "B"   = B_mat,
      "LOF" = LOF,
      "GCV" = GCV,
      "t"   = knots,
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

    # Basis functions to be removed excluding:
    # - B1(X) = 1
    # - Basis function when appearing alone in an interval: excepting right-hand extreme interval
    bf_to_drop <- basis_to_drop (
      Bp_list = Bp_list,
      concavity = concavity,
      nX = nX
      )

    # Mapping between Basis Function ID and Matrix Column Index by Variable
    id_Bpxi <- mapping_basis (
      Bp_list = Bp_list,
      nX = nX
      )

    # Add id column to basis function to drop
    bf_to_drop <- merge(bf_to_drop, id_Bpxi, by = c("xi", "Bp_xi"))

    # Shuffle rows
    bf_to_drop <- bf_to_drop[sample(1:nrow(bf_to_drop)), ]

    if (nrow(bf_to_drop) == 0) break

    for (bf in 1:nrow(bf_to_drop)) {

      # ======================= #
      # Update the set of basis #
      # ======================= #

      new_basis <- basis

      # id to drop the basis function from new_basis
      id <- bf_to_drop[bf, "id"]

      # select the basis function to be dropped
      dropTerm <- sapply(new_basis, function(x) x[["id"]] == id)

      # Number of 0's in the basis function to be dropped (in case of same GCV)
      basis_0s <- sum(new_basis[dropTerm][[1]][["Bp"]] == 0)

      # If a basis function is paired -> the sibling basis function is set as unpaired
      if (new_basis[dropTerm][[1]][["status"]] == "paired") {

        # Knot (t) and variable (v)
        check1 <- sapply(new_basis, "[[", "t") == new_basis[dropTerm][[1]][["t"]]
        check2 <- sapply(new_basis, "[[", "xi") == new_basis[dropTerm][[1]][["xi"]]

        # IDs for all basis functions in new_basis
        basis_ids <- sapply(new_basis, "[[", "id")

        # ID of the paired basis function
        paired_id <- basis_ids[check1 & check2]

        # ID of the sibling basis function
        sibling_id <- paired_id[!(id == paired_id)]

        # Update new_basis
        new_basis[lapply(new_basis, "[[", "id") == sibling_id][[1]][["status"]] <- "unpaired"
      }

      # Drop the selected basis function
      new_basis[dropTerm] <- NULL

      if (length(new_basis) == 1) {

        new_Bp_list <- vector("list", nX)
        new_B <- matrix(1, nrow = nrow(data))
        coefs <- unname(apply(data[, y, drop = FALSE], 2, max))

        knots <- NULL

        # prediction
        y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

        for (out in 1:length(y)) {
          y_hat[, out] <- coefs[out]
        }

        LOF <- err_metric (
          y_obs = data[, y, drop = F],
          y_hat = y_hat,
          metric = metric
        )

        GCV <- compute_gcv (
          y_obs = data[, y, drop = F],
          y_hat = y_hat,
          metric = metric,
          B = new_B,
          d = d,
          knots = 0,
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
        for (v in 1:nX) {
          for (side in c("paired", "right", "left")) {
            for (l in 1:length(new_Bp_list[[v]][[side]])) {
              bf_jp <- new_Bp_list[[v]][[side]][[l]][["Bp"]]
              new_B <- cbind(new_B, bf_jp)
            }
          }
        }

        # ===================== #
        # Estimate coefficients #
        # ===================== #

        coefs <- estimate_coefficients (
          B = new_B,
          y = data[, y, drop = F],
          it_list = new_it_list,
          Bp_list = new_Bp_list,
          monotonicity = monotonicity,
          concavity = concavity,
          x0_y0 = x0_y0
          )

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
          metric = metric
        )

        GCV <- compute_gcv (
          y_obs = data[, y, drop = F],
          y_hat = y_hat,
          metric = metric,
          B = new_B,
          d = d,
          knots = knots,
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

    basis <- best_basis
    Bp_list <- best_Bp_list

    # Create P-model
    p_model <- list (
      "id"  = ncol(best_B),
      "B"   = best_B,
      "LOF" = err,
      "GCV" = best_gcv,
      "t"   = best_knots,
      "coefs" = best_coefs,
      "Bp_list" = best_Bp_list
    )

    # Best model of p terms
    models <- append(models, list(p_model))

    # Update number of terms
    terms <- terms - 1
  }

  return(models)
}

#' @title Compute Generalized Cross-Validation (GCV)
#'
#' @description This function computes the generalized cross-validation for the Backward Algorithm in Adaptive Constrained Enveloping Splines.
#'
#' @param y_obs Vector of observed data.
#' @param y_hat Vector of predicted values.
#' @param metric Lack-of-fit criterion to evaluate the model performance.
#' @param B Matrix of basis functions.
#' @param d Generalized Cross Validation (GCV) penalty per knot.
#' @param knots \code{list} with the set of knots.
#' @param xi_degree Matrix indicating the degree of each input variable.
#'
#' @return The Generalized Cross-Validation value.

compute_gcv <- function (
    y_obs, y_hat, metric, B, d, knots, xi_degree
    ) {

  # number of outputs
  nY <- ncol(y_obs)
  # number of observations
  N <- nrow(y_obs)

  # mean error
  lof_num <- err_metric (
    y_obs = y_obs,
    y_hat = y_hat,
    metric = metric
    )

  # number of knots
  if (is.list(knots)) {
    # knot list to data.frame
    knots_df <- do.call(rbind.data.frame, knots)
    # knot dimension
    knots_df$dimknot <- xi_degree[2, match(knots_df$xi, xi_degree[1, ])]

    if (nrow(knots_df) == 0) {
      nknots <- 0
    } else {
      # drop duplicated
      knots_df_drop_duplicated <- knots_df[!duplicated(knots_df[c("xi", "t")]), ]
      # number of knots
      nknots <- nrow(knots_df_drop_duplicated)
    }

  } else {
    nknots <- knots
  }

  # Cost-complexity measure + penalization of 0 coefficients
  ccm_val <- ncol(B) + d * nknots
  lof_den <- (1 - (ccm_val / (nY * N))) ^ 2

  # GCV
  if (ccm_val > nY * N) {
    GCV <- Inf
  } else {
    GCV <- lof_num / lof_den
  }

  return(GCV)
}

#' @title Basis Functions to be Dropped
#'
#' @description This function determines the basis function to be dropped from a set of intervals.
#'
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param nX Number of inputs.
#'
#' @return A \code{data.frame} indicating the candidate basis function to be dropped for each input.

basis_to_drop <- function (
    Bp_list, concavity, nX
    ) {

  # Set of intervals from Bp_list
  it_list <- set_intervals (
    Bp_list = Bp_list
    )

  # Initialize a data.frame
  bf_to_drop <- data.frame()
  for (v in 1:nX) {
    if (length(it_list[[v]]) != 0) {
      for (it in 1:length(it_list[[v]])) {

        if (concavity) {
          # notice that concavity is lost in absence of a basis for a middle interval.
          # then, empty intervals are only available in the top right-side (last interval)
          # for example, if in the first interval there is only one basis function,
          # it can't be dropped.

          # it != length(It.list[[v]]):
          # check if it is not the last interval
          check1 <- it != length(it_list[[v]])

          # length(It.list[[v]][[it]][["Bp"]][["Bp_xi"]]) == 1:
          # check if there is only 1 basis function
          check2 <- length(it_list[[v]][[it]][["Bp"]][["Bp_xi"]]) == 1

        } else {

          check1 <- check2 <- FALSE
        }

        if (check1 && check2) {
          next
        } else {
          for (j in it_list[[v]][[it]][["Bp"]][["Bp_xi"]]) {
            bf_to_drop <- rbind(bf_to_drop, c(v, j))
          }
        }
      }
    }
  }

  names(bf_to_drop) <- c("xi", "Bp_xi")

  # Drop duplicated values
  bf_to_drop <- unique(bf_to_drop[order(bf_to_drop$xi, bf_to_drop$Bp_xi), ])

  return(bf_to_drop)
}

#' @title Mapping between Basis Function ID and Matrix Column Index by Variable
#'
#' @description This function creates a mapping between the ID of a basis function and its corresponding column index in the basis function matrix (B matrix) by variable.
#'
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#' @param nX Number of inputs.
#'
#' @return A \code{matrix} that establishes the connection between the basis function ID and its column index in the basis function matrix by variable.

mapping_basis <- function (
    Bp_list, nX
    ) {

  id_Bp_xi <- data.frame(xi = NA, id = NA, Bp_xi = NA)
  k <- 1

  for (v in 1:nX) {
    for (side in c("paired", "right", "left")) {
      if (!is.null(Bp_list[[v]][[side]])) {
        for (l in 1:length(Bp_list[[v]][[side]])) {

          # Select a Bp
          Bp <- Bp_list[[v]][[side]][[l]]

          # Save the id for the selected Bp
          Bp_id <- Bp[["id"]]

          # Save the matrix column index by the "v" variable
          Bp_xi <- Bp[["Bp_xi"]]

          if (side == "paired") {
            # id and matrix column connection
            id_Bp_xi[k, ] <- c(v, Bp_id[[1]], Bp_xi[[1]])
            id_Bp_xi[k + 1, ] <- c(v, Bp_id[[2]], Bp_xi[[2]])

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
#' @description This function updates the knot_list during the backward algorithm.
#'
#' @param Bp_list A \code{list} containing the set of basis functions by input.
#'
#' @return A \code{list} with the knots updated.

update_knot_list <- function (
    Bp_list
    ) {

  # Number of inputs
  nX <- length(Bp_list)

  knots <- list(); k <- 1
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
#' @description This function updates the Bp_list during the backward algorithm.
#'
#' @param basis A \code{list} containing the set of basis functions.
#' @param nX Number of inputs.
#'
#' @return A \code{list} with the basis functions updated.

update_Bp_list <- function (
    basis, nX
    ) {

  new_Bp_list <- vector("list", nX)

  for (v in 1:nX) {
    new_Bp_list[[v]] <- list (
      "paired" = NULL,
      "right"  = NULL,
      "left"   = NULL
      )
    }

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
    side <- basis[[i]][["status"]]

    if (side == "unpaired" && basis[[i]][["side"]] == "R") {

      new_Bp_list[[xi]][["right"]] <- append (
        new_Bp_list[[xi]][["right"]],
        list (
          list (
            "Bp" = Bp,
            "t"  = kn,
            "id" = id
            )
          )
        )

    } else if (side == "unpaired" && basis[[i]][["side"]] == "L") {

      new_Bp_list[[xi]][["left"]] <- append (
        new_Bp_list[[xi]][["left"]],
        list (
          list (
            "Bp" = Bp,
            "t"  = kn,
            "id" = id
            )
          )
        )

    } else {

      if (kn %in% sapply(new_Bp_list[[xi]][[side]], "[[", "t")) next

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

  # Sort basis function by variable | side | knot
  for (v in 1:nX) {
    for (side in c("paired", "right", "left")) {
      if (!is.null(new_Bp_list[[v]][[side]])) {
        arrangement <- order(sapply(new_Bp_list[[v]][[side]], "[[", "t"))
        new_Bp_list[[v]][[side]] <- new_Bp_list[[v]][[side]][arrangement]
      }
    }
  }


  # Number of column in matrix for ensuring concavity and monotonicity
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

