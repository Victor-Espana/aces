#' @title Generate Side Knot Locations
#'
#' @description
#' This function creates the input space for locating side knots in the smoothing procedure.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param nX
#' Number of inputs in \code{data}.
#'
#' @param knots
#' A \code{data.frame} containing the knots from the backward step.
#'
#' @return
#' A \code{list} with the locations of the central and the side knots in each dimension.

side_knot_location <- function (
    data,
    nX,
    knots
    ) {

  # initialize a list to save side knot locations
  kn_side_loc <- vector("list", nX)

  for (v in 1:nX) {

    # get knots for the v-th variable
    kt_v <- sort(unique(knots[knots[, 1] == v, 2]))

    if (length(kt_v) == 0) next

    # add the initial and the end observation. They cannot be used as knots.
    kn_side_loc[[v]] <- c(min(data[, v]), kt_v, max(data[, v]))

    # compute the midpoints between central knots
    kn_side_loc[[v]] <- sort(c(kn_side_loc[[v]], kn_side_loc[[v]][- length(kn_side_loc[[v]])] + diff(kn_side_loc[[v]]) / 2))

  }

  return(kn_side_loc)

}

#' @title Generate a Suitable Triplet of Knots
#'
#' @description
#' This function generates a suitable triplet of knots for the smoothing procedure depending on the shape-constraints.
#'
#' @param knots
#' A \code{vector} with the initial locations for central and side knots.
#'
#' @param w
#' Hyperparameter for side knot distances in the cubic and quintic smoothing procedures.
#'
#' @param smoothing
#' Type of smoothing:
#' \itemize{
#' \item{\code{"cubic"}}: cubic smoothing procedure.
#' \item{\code{"quintic"}}: quintic smoothing procedure.
#' }
#'
#' @param monotonicity
#' A \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#'
#' @param concavity
#' A \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @return
#' A \code{vector} with a suitable triplet of knots.

set_triplet_knots <- function (
    knots,
    w,
    smoothing,
    mono,
    conc
    ) {

  # unconstrained knots
  t0 <- knots[1] # t-
  t1 <- knots[2] # t
  t2 <- knots[3] # t+

  # delta & epsilon
  d <- t1 - t0 # t  - t-
  e <- t2 - t1 # t+ - t

  # constrained knots
  if (mono | conc) {

    if (smoothing == "cubic") {

      ratio <- d / e

      if (ratio < 1 | ratio > 2) {
        # update t+
        e <- w * d
        t2 <- t1 + e
      }

    } else if (smoothing == "quintic") {

      ratio <- e / d

      if (ratio < 8/7 | ratio > 1.5) {
        # update t-
        d <- w * e
        t0 <- t1 - d

      }
    }
  }

  return(c(t0, t1, t2))

}

#' @title Fit a Smooth Adaptive Constrained Enveloping Splines using Cubic Functions
#'
#' @description
#'
#' This function performs smoothing on the Adaptive Constrained Enveloping Splines predictor while enforcing continuous first derivatives.
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
#' @param dea_eff
#' An indicator vector with 1s for efficient DMUs and 0s for inefficient DMUs.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function will estimate.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param kn_grid
#' A \code{data.frame} containing the set of knots from the backward step for smoothing.
#'
#' @param kn_side_loc
#' A \code{list} with the side knot locations for each input dimension.
#'
#' @param d
#' Generalized Cross Validation (GCV) penalty per knot.
#'
#' @param wc
#' Let `p` be the distance between the central knot and the right-side knot, and `v` be the distance between the central knot and the left-side knot during the smoothing procedure.A \code{numeric} value used for cubic smoothing \insertCite{friedman1991}{aces}. This parameter is defined as `v / p` and must be set between 1 and 2. If a \code{vector} is entered, the \code{wc} value that most reduced the lack-of-fit criterion is selected.
#'
#' @references
#'
#' \insertRef{friedman1991}{aces}
#'
#' @return
#'
#' A \code{list} containing information relative to the Smooth Adaptive Constrained Enveloping Splines through cubic functions.

cubic_aces <- function (
    data,
    x,
    y,
    dea_eff,
    model_type,
    metric,
    shape,
    kn_grid,
    kn_side_loc,
    d,
    wc
    ) {

  if (shape[["ptto"]] != FALSE) {
    data <- rbind(data, rep(0, ncol(data)))
  }

  # number of inputs
  nX <- length(x)

  # sample size
  N <- nrow(data)

  # select "wc" best distance between central and side knots
  w_list <- vector("list", length(wc))

  for (l in 1:length(wc)) {

    # new matrix of basis functions
    B <- matrix(data = rep(1, N), ncol = 1)

    # cubic knots
    cubic_knots <- vector("list", nX)

    # paired basis functions
    paired <- c()

    # not paired basis functions
    not_paired <- c()

    for (v in 1:nX) {

      if (is.null(kn_side_loc[[v]])) next

      # from first midpoint: position 2
      # to penultimate midpoint: position (-3)
      # step 2 to select midpoints

      for (i in seq(2, length(kn_side_loc[[v]]) - 3, 2)) {

        # select a central knot: position i + 1
        t <- kn_side_loc[[v]][i + 1]

        # side of that knot
        side <- kn_grid[kn_grid[, "t"] == t, "side"]

        # triplet of knots
        triplet <- set_triplet_knots (
          knots = kn_side_loc[[v]][i:(i + 2)],
          w = wc[l],
          smoothing = "cubic",
          mono = shape[["mono"]],
          conc = shape[["conc"]]
        )

        # update B with two new truncated cubic basis functions
        B <- create_cubic_basis (
          data = data,
          xi = v,
          knots = triplet,
          B = B,
          side = side
          )

        if (length(side) == 1) {

          not_paired <- c(not_paired, ncol(B))
          status <- "unpaired"
          side <- side

        } else {

          paired <- c(paired, (ncol(B) - 1):ncol(B))
          status <- "paired"
          side <- side

        }

        # update cubic knots
        cubic_knots[[v]] <- append (
          cubic_knots[[v]],
          list (
            list (
              t = triplet,
              status = status,
              side = side
            )
          )
        )
      }
    }

    if (shape[["mono"]] || shape[["conc"]]) {

      # number of paired basis functions
      n_pair <- sum(duplicated(kn_grid[, c("xi", "t")])) * 2

      # number of unpaired right side basis functions
      n_lsub <- nrow(kn_grid) - n_pair

      # reorder B to fit the positions of estimate_coefficients_smoothed
      B <- B[, c(1, paired, not_paired)]

      # estimate coefficients
      coefs <- estimate_coefficients_smoothed (
        model_type = model_type,
        B = B,
        y_obs = data[, y, drop = F],
        dea_eff = dea_eff,
        n_pair = n_pair,
        n_lsub = n_lsub,
        shape = shape
        )

    } else {

      coefs <- estimate_coefficients (
        model_type = model_type,
        B = B[, c(1, paired, not_paired)],
        y = data[, y, drop = F],
        dea_eff = dea_eff,
        it_list = NULL,
        Bp_list = NULL,
        shape = shape
      )

    }

    # prediction
    y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

    for (out in 1:length(y)) {
      y_hat[, out] <- B %*% coefs[, out, drop = F]
    }

    # compute lack-of-fit
    GCV <- compute_gcv (
      y_obs = data[, y, drop = F],
      y_hat = y_hat,
      metric = metric,
      n_bf = ncol(B),
      d = d,
      knots = sum(lengths(cubic_knots))
    )

    if (shape[["ptto"]] != FALSE) {
      B <- B[1:(nrow(B) - 1), ]
    }

    # save results of the smooth model
    w_list[[l]][["Bmatx"]] <- B
    w_list[[l]][["cubic_knots"]] <- cubic_knots
    w_list[[l]][["coefs"]] <- coefs
    w_list[[l]][["GCV"]] <- GCV
    w_list[[l]][["w"]] <- wc[l]
  }

  # lack-of-fit for each model
  GCVs <- sapply(w_list, function(x) x[["GCV"]])

  # model with minimum error
  aces_cubic <- w_list[[which.min(GCVs)]]

  return(aces_cubic)

}

#' @title Generate a New Pair of Cubic Basis Functions
#'
#' @description
#'
#' This function generates two new cubic basis functions from a variable and a triplet of knots.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param xi
#' Index of the variable that creates the new basis function(s).
#'
#' @param knots
#' Triplet of knots for creating the new basis function(s).
#'
#' @param B
#' A \code{matrix} of basis functions.
#'
#' @param side
#' Side of the new basis function inherited from the linear basis function.
#'
#' @return
#'
#' A \code{matrix} of basis functions updated with the new cubic basis functions.

create_cubic_basis <- function (
    data,
    xi,
    knots,
    B,
    side
    ) {

  # triplet of knots
  t0 <- knots[1]
  t1 <- knots[2]
  t2 <- knots[3]

  # delta & epsilon
  d <- t1 - t0 # t  - t-
  e <- t2 - t1 # t+ - t

  p1 <- (2 * e - d) / (e + d) ^ 2
  r1 <- (d - e) / (e + d) ^ 3

  p2 <- (2 * d - e) / (- e - d) ^ 2
  r2 <- (e - d) / (- e - d) ^ 3

  # side of the basis function: both or right
  if (length(side) == 2 || side == "R") {

    term1 <- p1 * (data[, xi] - t0) ^ 2
    term2 <- r1 * (data[, xi] - t0) ^ 3

    # C1
    C1 <- ifelse (
      data[, xi] <= t0,
      0,
      ifelse (
        data[, xi] > t0 & data[, xi] < t2,
        term1 + term2,
        data[, xi] - t1
        )
      )

    # add C1 to B matrix
    B <- cbind(B, C1)

  }

  # side of the basis function: both or right
  if (length(side) == 2 || side == "L") {

    term1 <- p2 * (data[, xi] - t2) ^ 2
    term2 <- r2 * (data[, xi] - t2) ^ 3

    # C2
    C2 <- ifelse (
      data[, xi] <= t0,
      t1 - data[, xi],
      ifelse (
        data[, xi] > t0 & data[, xi] < t2,
        term1 + term2,
        0
        )
      )

    # add C2 to B matrix
    B <- cbind(B, C2)

  }

  return(B)

}

#' @title Fit a Smooth Adaptive Constrained Enveloping Splines using Quintic Functions
#'
#' @description
#'
#' This function performs smoothing on the Adaptive Constrained Enveloping Splines predictor while enforcing continuous second derivatives.
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
#' @param dea_eff
#' Indicator vector with 1s for efficient DMUs and 0s for inefficient DMUs.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function will estimate.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param kn_grid
#' A \code{data.frame} containing the set of knots from the backward step for smoothing.
#'
#' @param kn_side_loc
#' A \code{list} with the side knot locations for each input dimension.
#'
#' @param d
#' Generalized Cross Validation (GCV) penalty per knot.
#'
#' @param wq
#' Let `p` be the distance between the central knot and the right-side knot, and `v` be the distance between the central knot and the left-side knot during the smoothing procedure. A \code{numeric} value used for quintic smoothing \insertCite{chen1999}{aces}. This parameter is defined as `p / v` and must be set between 8/7 and 1.5. If a \code{vector} is entered, the \code{wc} value that most reduced the lack-of-fit criterion is selected.
#'
#' @references
#'
#' \insertRef{chen1999}{aces}
#'
#' @importFrom Rdpack reprompt
#'
#' @return
#'
#' A \code{list} containing information relative to the Smooth Adaptive Constrained Enveloping Splines through quintic functions.

quintic_aces <- function (
    data,
    x,
    y,
    dea_eff,
    model_type,
    metric,
    shape,
    kn_grid,
    kn_side_loc,
    d,
    wq
    ) {

  if (shape[["ptto"]] != FALSE) {
    data <- rbind(data, rep(0, ncol(data)))
  }

  # sample size
  N <- nrow(data)

  # number of inputs
  nX <- length(x)

  # select "wq" best distance between central and side knots
  w_list <- vector("list", length(wq))

  for (l in 1:length(wq)) {

    # new matrix of basis functions
    B <- matrix(data = rep(1, N), ncol = 1)

    # quintic knots
    quintic_knots <- vector("list", nX)

    # paired basis functions
    paired <- c()

    # not paired basis functions
    not_paired <- c()

    for (v in 1:nX){
      if (is.null(kn_side_loc[[v]])) next

      # from first midpoint: position 2
      # to penultimate midpoint: position (-3)
      # step 2 to select midpoints

      for (i in seq(2, length(kn_side_loc[[v]]) - 3, 2)) {

        # select a central knot: position i + 1
        t <- kn_side_loc[[v]][i + 1]

        # side of that knot
        side <- kn_grid[kn_grid[, "t"] == t, "side"]

        # triplet of knots
        triplet <- set_triplet_knots (
          knots = kn_side_loc[[v]][i:(i + 2)],
          w = wq[l],
          smoothing = "quintic",
          mono = shape[["mono"]],
          conc = shape[["conc"]]
        )

        # update B with two new truncated quintic basis functions
        B <- create_quintic_basis (
          data = data,
          xi = v,
          knots = triplet,
          B = B,
          side = side
          )

        if (length(side) == 1) {

          not_paired <- c(not_paired, ncol(B))
          status <- "unpaired"
          side <- side

        } else {

          paired <- c(paired, (ncol(B) - 1):ncol(B))
          status <- "paired"
          side <- side

        }

        # update cubic knots
        quintic_knots[[v]] <- append (
          quintic_knots[[v]],
          list (
            list (
              t = triplet,
              status = status,
              side = side
            )
          )
        )

      }
    }

    if (shape[["mono"]] || shape[["conc"]]) {

      # number of paired basis functions
      n_pair <- sum(duplicated(kn_grid[, c("xi", "t")])) * 2

      # number of unpaired right side basis functions
      n_lsub <- nrow(kn_grid) - n_pair

      # reorder B to fit the positions of estimate_coefficients_smoothed
      B <- B[, c(1, paired, not_paired)]

      # estimate coefficients
      coefs <- estimate_coefficients_smoothed (
        model_type = model_type,
        B = B,
        y_obs = data[, y, drop = F],
        dea_eff = dea_eff,
        n_pair = n_pair,
        n_lsub = n_lsub,
        shape = shape
      )

    } else {

      coefs <- estimate_coefficients (
        model_type = model_type,
        B = B[, c(1, paired, not_paired)],
        y = data[, y, drop = F],
        dea_eff = dea_eff,
        it_list = NULL,
        Bp_list = NULL,
        shape = shape
      )

    }

    # prediction
    y_hat <- matrix(NA, nrow = N, ncol = length(y))

    for (out in 1:length(y)) {
      y_hat[, out] <- B %*% coefs[, out, drop = F]
    }

    # compute lack-of-fit
    GCV <- compute_gcv (
      y_obs = data[, y, drop = F],
      y_hat = y_hat,
      metric = metric,
      n_bf = ncol(B),
      d = d,
      knots = sum(lengths(quintic_knots))
    )

    if (shape[["ptto"]] != FALSE) {
      B <- B[1:(nrow(B) - 1), ]
    }

    # save results of the smooth model
    w_list[[l]][["Bmatx"]] <- B
    w_list[[l]][["quintic_knots"]] <- quintic_knots
    w_list[[l]][["coefs"]] <- coefs
    w_list[[l]][["GCV"]] <- GCV
    w_list[[l]][["w"]] <- wq[l]

  }

  # lack-of-fit for each model
  GCVs <- sapply(w_list, function(x) x[["GCV"]])

  # model with minimum error
  aces_quintic <- w_list[[which.min(GCVs)]]

  return(aces_quintic)

}

#' @title Generate a New Pair of Quintic Basis Functions
#'
#' @description
#'
#' This function generates two new quintic basis functions from a variable and a triplet of knots.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param xi
#' Index of the variable that creates the new basis function(s).
#'
#' @param knots
#' Triplet of knots for creating the new basis function(s).
#'
#' @param B
#' A \code{matrix} of basis functions.
#'
#' @param side
#' Side of the new basis function inherited from the linar basis function.
#'
#' @return
#'
#' A \code{matrix} of basis functions updated with the new quintic basis functions.

create_quintic_basis <- function (
    data,
    xi,
    knots,
    B,
    side
    ) {

  # triplet of knots
  t0 <- knots[1]
  t1 <- knots[2]
  t2 <- knots[3]

  d1 <- t2 - t1 # t+ - t
  d2 <- t1 - t0 # t  - t-
  d  <- t2 - t0 # t+ - t-

  alpha1 <- (6 * d1 - 4 * d2) / d ^ 3
  alpha2 <- (4 * d1 - 6 * d2) / d ^ 3

  beta1 <- (7 * d2 - 8 * d1) / d ^ 4
  beta2 <- (7 * d1 - 8 * d2) / d ^ 4

  gamma1 <- gamma2 <- (3 * d1 - 3 * d2) / d ^ 5

  # side of the basis function: both or right
  if (length(side) == 2 || side == "R") {

    term1 <- alpha1 * (data[, xi] - t0) ^ 3
    term2 <- beta1  * (data[, xi] - t0) ^ 4
    term3 <- gamma1 * (data[, xi] - t0) ^ 5

    # Q1
    Q1 <- ifelse (
      data[, xi] <= t0,
      0,
      ifelse (
        data[, xi] > t0 & data[, xi] < t2,
        term1 + term2 + term3,
        data[, xi] - t1
        )
      )

    # add Q1 to matrix B
    B  <- cbind(B, Q1)
  }

  # side of the basis function: both or right
  if (length(side) == 2 || side == "L") {

    term1 <- alpha2 * (data[, xi] - t2) ^ 3
    term2 <- beta2  * (data[, xi] - t2) ^ 4
    term3 <- gamma2 * (data[, xi] - t2) ^ 5

    # Q2
    Q2 <- ifelse (
      data[, xi] <= t0,
      t1 - data[, xi],
      ifelse (
        data[, xi] > t0 & data[, xi] < t2,
        term1 + term2 + term3,
        0
        )
      )

    # add Q2 to matrix B
    B  <- cbind(B, Q2)
  }

  return(B)

}
