#' @title Find Side-Knot Locations
#'
#' @description
#' Finds the neighboring data values available for side knots around each central
#' knot.
#'
#' @param data
#' A matrix containing the prepared model variables.
#'
#' @param nX
#' Number of inputs in \code{data}.
#'
#' @param knots
#' Data frame of knots retained after pruning.
#'
#' @return
#' A list of central and neighboring knot locations for each input.

side_knot_location <- function(
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
    kn_side_loc[[v]] <- sort(c(kn_side_loc[[v]], kn_side_loc[[v]][-length(kn_side_loc[[v]])] + diff(kn_side_loc[[v]]) / 2))
  }

  return(kn_side_loc)
}

#' @title Create a Smoothing Knot Triplet
#'
#' @description
#' Moves the side knots around a central knot according to the smoothing and shape
#' settings.
#'
#' @param knots
#' Numeric vector containing the left, central, and right knot locations.
#'
#' @param w
#' Side-knot distance parameter.
#'
#' @param smoothing
#' Smoothing type: \code{"cubic"} or \code{"quintic"}.
#'
#' @param mono
#' If \code{TRUE}, enforce non-decreasing monotonicity.
#'
#' @param conc
#' If \code{TRUE}, enforce concavity.
#'
#' @return
#' A numeric vector with the adjusted left, central, and right knots.

set_triplet_knots <- function(
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

      if (ratio < 8 / 7 | ratio > 1.5) {
        # update t-
        d <- w * e
        t0 <- t1 - d
      }
    }
  }

  return(c(t0, t1, t2))
}

#' @title Fit a Cubic-Smoothed ACES Model
#'
#' @description
#'
#' Replaces the linear hinges in a pruned ACES model with cubic basis functions
#' that join smoothly at the side knots and have continuous first derivatives.
#' Coefficients are re-estimated under the same envelopment and shape constraints,
#' and the candidate side-knot ratio with the lowest lack of fit is retained.
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
#' @param dea_scores
#' A matrix of output-oriented DEA-VRS scores, with one column per output.
#'
#' @param fdh_scores
#' A matrix of output-oriented FDH scores, with one column per output.
#'
#' @param metric
#' Character string specifying the lack-of-fit measure.
#'
#' @param shape
#' A list with logical elements \code{mono} and \code{conc}.
#'
#' @param kn_grid
#' Data frame of knots retained after pruning.
#'
#' @param kn_side_loc
#' Available side-knot locations for each input.
#'
#' @param kn_penalty
#' Penalty per knot used to compute GCV.
#'
#' @param xi_degree
#' A matrix that records the degree of each input in every expanded variable.
#'
#' @param wc
#' Cubic side-knot distance ratio \eqn{v / p}, where \eqn{v} and \eqn{p} are the
#' distances from the central knot to its left and right side knots. Values must
#' be between 1 and 2. If several values are supplied, the function selects the
#' one with the lowest lack of fit.
#'
#' @references
#'
#' \insertRef{friedman1991}{aces}
#'
#' @return
#'
#' A list containing the cubic-smoothed model, its knots, error, and GCV score.

cubic_aces <- function(
  data,
  x,
  y,
  dea_scores,
  fdh_scores,
  metric,
  shape,
  kn_grid,
  kn_side_loc,
  kn_penalty,
  xi_degree,
  wc
) {
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
        triplet <- set_triplet_knots(
          knots = kn_side_loc[[v]][i:(i + 2)],
          w = wc[l],
          smoothing = "cubic",
          mono = shape[["mono"]],
          conc = shape[["conc"]]
        )

        # update B with two new truncated cubic basis functions
        B <- create_cubic_basis(
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
        cubic_knots[[v]] <- append(
          cubic_knots[[v]],
          list(
            list(
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
      coefs <- estimate_coefficients_smoothed(
        B = B,
        y_obs = data[, y, drop = F],
        dea_scores = dea_scores,
        n_pair = n_pair,
        n_lsub = n_lsub,
        shape = shape
      )
    } else {
      coefs <- estimate_coefficients(
        B = B[, c(1, paired, not_paired)],
        y_obs = data[, y, drop = F],
        dea_scores = dea_scores,
        fdh_scores = fdh_scores,
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

    LOF <- err_metric(
      y_obs = data[, y, drop = F],
      y_hat = y_hat,
      metric = metric,
      weight = 1 / dea_scores
    )

    GCV <- compute_gcv(
      y_obs = data[, y, drop = F],
      LOF = LOF,
      n_bf = ncol(B),
      kn_penalty = kn_penalty,
      kn_list = lengths(cubic_knots),
      xi_degree = xi_degree
    )

    # save results of the smooth model
    w_list[[l]][["Bmatx"]] <- B
    w_list[[l]][["knots"]] <- cubic_knots
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

#' @title Add Cubic Basis Functions
#'
#' @description
#'
#' Evaluates cubic basis functions for one input and adds them to a basis matrix.
#'
#' @param data
#' A matrix containing the prepared model variables.
#'
#' @param xi
#' Index of the input used to create the basis functions.
#'
#' @param knots
#' Numeric vector containing the left, central, and right knots.
#'
#' @param B
#' Current basis matrix.
#'
#' @param side
#' Side inherited from the corresponding linear basis function.
#'
#' @return
#'
#' The basis matrix with the new cubic columns.

create_cubic_basis <- function(
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

  p1 <- (2 * e - d) / (e + d)^2
  r1 <- (d - e) / (e + d)^3

  p2 <- (2 * d - e) / (-e - d)^2
  r2 <- (e - d) / (-e - d)^3

  # side of the basis function: both or right
  if (length(side) == 2 || side == "R") {
    term1 <- p1 * (data[, xi] - t0)^2
    term2 <- r1 * (data[, xi] - t0)^3

    # C1
    C1 <- ifelse(
      data[, xi] <= t0,
      0,
      ifelse(
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
    term1 <- p2 * (data[, xi] - t2)^2
    term2 <- r2 * (data[, xi] - t2)^3

    # C2
    C2 <- ifelse(
      data[, xi] <= t0,
      t1 - data[, xi],
      ifelse(
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

#' @title Fit a Quintic-Smoothed ACES Model
#'
#' @description
#'
#' Replaces the linear hinges in a pruned ACES model with quintic basis functions
#' that join smoothly at the side knots and have continuous second derivatives.
#' Coefficients are re-estimated under the same envelopment and shape constraints,
#' and the candidate side-knot ratio with the lowest lack of fit is retained.
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
#' @param dea_scores
#' A matrix of output-oriented DEA-VRS scores, with one column per output.
#'
#' @param fdh_scores
#' A matrix of output-oriented FDH scores, with one column per output.
#'
#' @param metric
#' Character string specifying the lack-of-fit measure.
#'
#' @param shape
#' A list with logical elements \code{mono} and \code{conc}.
#'
#' @param kn_grid
#' Data frame of knots retained after pruning.
#'
#' @param kn_side_loc
#' Available side-knot locations for each input.
#'
#' @param kn_penalty
#' Penalty per knot used to compute GCV.
#'
#' @param xi_degree
#' A matrix that records the degree of each input in every expanded variable.
#'
#' @param wq
#' Quintic side-knot distance ratio used to place the two side knots around each
#' central knot. Values must be between 8/7 and 1.5. If several values are
#' supplied, the function selects the one with the lowest lack of fit.
#'
#' @references
#'
#' \insertRef{chen1999}{aces}
#'
#' @importFrom Rdpack reprompt
#'
#' @return
#'
#' A list containing the quintic-smoothed model, its knots, error, and GCV score.

quintic_aces <- function(
  data,
  x,
  y,
  dea_scores,
  fdh_scores,
  metric,
  shape,
  kn_grid,
  kn_side_loc,
  kn_penalty,
  xi_degree,
  wq
) {
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
        triplet <- set_triplet_knots(
          knots = kn_side_loc[[v]][i:(i + 2)],
          w = wq[l],
          smoothing = "quintic",
          mono = shape[["mono"]],
          conc = shape[["conc"]]
        )

        # update B with two new truncated quintic basis functions
        B <- create_quintic_basis(
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
        quintic_knots[[v]] <- append(
          quintic_knots[[v]],
          list(
            list(
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
      coefs <- estimate_coefficients_smoothed(
        B = B,
        y_obs = data[, y, drop = F],
        dea_scores = dea_scores,
        n_pair = n_pair,
        n_lsub = n_lsub,
        shape = shape
      )
    } else {
      coefs <- estimate_coefficients(
        B = B[, c(1, paired, not_paired)],
        y_obs = data[, y, drop = F],
        dea_scores = dea_scores,
        fdh_scores = fdh_scores,
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

    LOF <- err_metric(
      y_obs = data[, y, drop = F],
      y_hat = y_hat,
      metric = metric,
      weight = 1 / dea_scores
    )

    GCV <- compute_gcv(
      y_obs = data[, y, drop = F],
      LOF = LOF,
      n_bf = ncol(B),
      kn_penalty = kn_penalty,
      kn_list = lengths(quintic_knots),
      xi_degree = xi_degree
    )

    # save results of the smooth model
    w_list[[l]][["Bmatx"]] <- B
    w_list[[l]][["knots"]] <- quintic_knots
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

#' @title Add Quintic Basis Functions
#'
#' @description
#'
#' Evaluates quintic basis functions for one input and adds them to a basis matrix.
#'
#' @param data
#' A matrix containing the prepared model variables.
#'
#' @param xi
#' Index of the input used to create the basis functions.
#'
#' @param knots
#' Numeric vector containing the left, central, and right knots.
#'
#' @param B
#' Current basis matrix.
#'
#' @param side
#' Side inherited from the corresponding linear basis function.
#'
#' @return
#'
#' The basis matrix with the new quintic columns.

create_quintic_basis <- function(
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
  d <- t2 - t0 # t+ - t-

  alpha1 <- (6 * d1 - 4 * d2) / d^3
  alpha2 <- (4 * d1 - 6 * d2) / d^3

  beta1 <- (7 * d2 - 8 * d1) / d^4
  beta2 <- (7 * d1 - 8 * d2) / d^4

  gamma1 <- gamma2 <- (3 * d1 - 3 * d2) / d^5

  # side of the basis function: both or right
  if (length(side) == 2 || side == "R") {
    term1 <- alpha1 * (data[, xi] - t0)^3
    term2 <- beta1 * (data[, xi] - t0)^4
    term3 <- gamma1 * (data[, xi] - t0)^5

    # Q1
    Q1 <- ifelse(
      data[, xi] <= t0,
      0,
      ifelse(
        data[, xi] > t0 & data[, xi] < t2,
        term1 + term2 + term3,
        data[, xi] - t1
      )
    )

    # add Q1 to matrix B
    B <- cbind(B, Q1)
  }

  # side of the basis function: both or right
  if (length(side) == 2 || side == "L") {
    term1 <- alpha2 * (data[, xi] - t2)^3
    term2 <- beta2 * (data[, xi] - t2)^4
    term3 <- gamma2 * (data[, xi] - t2)^5

    # Q2
    Q2 <- ifelse(
      data[, xi] <= t0,
      t1 - data[, xi],
      ifelse(
        data[, xi] > t0 & data[, xi] < t2,
        term1 + term2 + term3,
        0
      )
    )

    # add Q2 to matrix B
    B <- cbind(B, Q2)
  }

  return(B)
}
