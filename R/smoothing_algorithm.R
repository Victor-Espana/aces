#' @title Generate Side Knot Locations
#'
#' @description This function creates the input space for locating side knots in the smoothing procedure.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param nX Number of inputs in \code{data}.
#' @param knots A \code{data.frame} containing the knots from the backward step.
#'
#' @return A \code{list} with the locations of the central and the side knots in each dimension.

side_knots_location <- function (
    data, nX, knots
    ) {

  # List for x_space
  x_space <- vector("list", nX)

  for (v in 1:nX) {

    # knots for the "v" variable
    kt_v <- sort(unique(knots[knots[, 1] == v, 2]))

    if (length(kt_v) == 0) next

    # add the initial and the end observation. They cannot be used as knots.
    x_space[[v]] <- c(min(data[, v]), kt_v, max(data[, v]))

    # calculate the midpoints between central knots
    x_space[[v]] <- sort(c(x_space[[v]], x_space[[v]][- length(x_space[[v]])] + diff(x_space[[v]]) / 2))
  }

  return(x_space)
}

#' @title Generate a Suitable Triplet of Knots
#'
#' @description This function generates a suitable triplet of knots for the smoothing procedure depending on the shape-constraints.
#'
#' @param knots A \code{vector} with the initial locations for central and side knots.
#' @param w Hyperparameter for the side knot distances in the cubic and quintic smoothing procedures.
#' @param smoothing Type of smoothing:
#' \itemize{
#' \item{\code{"cubic"}}: cubic smoothing procedure.
#' \item{\code{"quintic"}}: quintic smoothing procedure.
#' }
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @return A \code{vector} with a suitable triplet of knots.

set_triplet_knots <- function (
    knots, w, smoothing, monotonicity, concavity
    ) {

  # unconstrained knots
  t0 <- knots[1] # t-
  t1 <- knots[2] # t
  t2 <- knots[3] # t+

  # delta & epsilon
  d <- t1 - t0 # t  - t-
  e <- t2 - t1 # t+ - t

  # constrained knots
  if (monotonicity | concavity) {
    if (smoothing == "cubic") {
      ratio <- d / e
    } else {
      ratio <- e / d
    }

  if (smoothing == "cubic" && (ratio < 1 | ratio > 2)) {
    # update t+
    e <- w * d
    t2 <- t1 + e

  } else if (smoothing == "quintic" && (ratio < 8/7 | ratio > 1.5))
    # update t-
    d <- w * e
    t0 <- t1 - d
  }

  return(c(t0, t1, t2))
}

#' @title Fit a Smoothing Adaptive Constrained Enveloping Splines using Cubic Functions
#'
#' @description This function performs smoothing on the Adaptive Constrained Enveloping Splines predictor while enforcing continuous first derivatives.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param y Column output indexes in \code{data}.
#' @param nX Number of inputs in \code{data}.
#' @param metric The lack-of-fit criterion to evaluate the model performance.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#' @param knots A \code{data.frame} containing the set of knots from backward step for smoothing.
#' @param x_space A \code{list} with the central side knot locations in each dimension.
#' @param wc For the cubic smoothing procedure \insertCite{friedman1991}{aces}. Let `d` be the distance between the central knot and the left side knot, and let be `e` be the distance between the central knot and the right side knot. If the condition `1 < d / e < 2` is not satisfied, then `e = wc * d`. It must be set between 1 and 2. If a \code{vector} is entered, the \code{wc} value that most reduce the lack-of-fit criterion is selected.
#' @param d Generalized Cross Validation (GCV) penalty per knot.
#'
#' @references
#' \insertRef{friedman1991}{aces}
#'
#' @return A \code{list} containing the information of the Smoothing Adaptive Constrained Enveloping Splines through cubic functions.

cubic_aces <- function (
    data, y, nX, metric, monotonicity, concavity, x0_y0, knots, x_space, wc, d
    ) {

  # sample size
  N <- nrow(data)

  # select "wc" best distance between central and side knots
  w_list <- vector("list", length(wc))

  for (l in 1:length(wc)) {

    # new matrix of basis functions
    B <- matrix(data = rep(1, N), ncol = 1)

    # cubic knots
    cubic_knots <- vector("list", nX)

    paired <- c(); not_paired <- c()
    for (v in 1:nX) {
      if (is.null(x_space[[v]])) next

      # from first midpoint: position 2
      # to penultimate midpoint: position (-3)
      # step 2 to select midpoints
      for (i in seq(2, length(x_space[[v]]) - 3, 2)) {

        # select a central knot: position i + 1
        t <- x_space[[v]][i + 1]

        # side of that knot
        side <- knots[knots[, "t"] == t, "side"]

        # triplet of knots
        triplet <- set_triplet_knots (
          knots = x_space[[v]][i:(i + 2)],
          w = wc[l],
          smoothing = "cubic",
          monotonicity = monotonicity,
          concavity = concavity
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

    if (monotonicity || concavity) {

      # number of paired basis functions
      h <- sum(duplicated(knots[, c("xi", "t")])) * 2
      # number of unpaired right side basis functions
      r <- nrow(knots) - h

      # Reorder B to fit the positions of estimate_coefficients_smoothed
      B <- B[, c(1, paired, not_paired)]

      # estimate coefficients
      coefs <- estimate_coefficients_smoothed (
        B = B,
        y = data[, y, drop = F],
        h = h,
        r = r,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0
        )

    } else {

      coefs <- estimate_coefficients (
        B = B[, c(1, paired, not_paired)],
        y = data[, y, drop = F],
        it_list = NULL,
        Bp_list = NULL,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0
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
      B = B,
      d = d,
      knots = sum(lengths(cubic_knots))
    )

    # save results
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

#' @title Generate a new pair of Cubic Basis Functions
#'
#' @description This function generates two new cubic basis functions from a variable and a triplet of knots.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param xi Variable index of the variable that creates the new basis function(s).
#' @param knots Triplet of knots for creating the new basis function(s).
#' @param B A \code{matrix} of basis functions.
#' @param side Side of the new basis function.
#'
#' @return A \code{matrix} of basis functions updated with the new cubic basis functions.

create_cubic_basis <- function (
    data, xi, knots, B, side
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

  # both or right
  if (length(side) == 2 || side == "R") {

    term1 <- p1 * (data[, xi] - t0) ^ 2
    term2 <- r1 * (data[, xi] - t0) ^ 3

    # C1
    C1 <- ifelse(data[, xi] <= t0,
                 0,
                 (ifelse((data[, xi] > t0) & (data[, xi] < t2),
                         term1 + term2,
                         data[, xi] - t1)))


    B <- cbind(B, C1)
  }

  # both or left
  if (length(side) == 2 || side == "L") {

    term1 <- p2 * (data[, xi] - t2) ^ 2
    term2 <- r2 * (data[, xi] - t2) ^ 3

    # C2
    C2 <- ifelse(data[, xi] <= t0,
                 t1 - data[, xi],
                 (ifelse((data[, xi] > t0) & (data[, xi] < t2),
                         term1 + term2,
                         0)))

    B <- cbind(B, C2)
  }

  return(B)
}


#' @title Fit a Smoothing Adaptive Constrained Enveloping Splines using Quintic Functions
#'
#' @description This function performs smoothing on the Adaptive Constrained Enveloping Splines predictor while enforcing continuous second derivatives.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param y Column output indexes in \code{data}.
#' @param nX Number of inputs in \code{data}.
#' @param metric The lack-of-fit criterion to evaluate the model performance.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#' @param knots A \code{data.frame} containing the set of knots from backward step for smoothing.
#' @param x_space A \code{list} with the central side knot locations in each dimension.
#' @param wq For the quintic smoothing procedure \insertCite{chen1999}{aces}. Let `d` be the distance between the central knot and the left side knot, and let be `e` be the distance between the central knot and the right side knot. If the condition `8/7 < e / d < 1.5` is not satisfied, then `d = wq * e`. It must be set between 8/7 and 1.5. If a \code{vector} is entered, the \code{wc} value that most reduce the lack-of-fit criterion is selected.
#' @param d Generalized Cross Validation (GCV) penalty per knot.
#'
#' @references
#' \insertRef{chen1999}{aces}
#'
#' @importFrom Rdpack reprompt
#'
#' @return List containing the set of knots from backward (\code{knots}), the new quintic knots (\code{quintic_knots}) and the set of coefficients (\code{coefs}).

quintinc_aces <- function (
    data, y, nX, metric, monotonicity, concavity, x0_y0, knots, x_space, wq, d
    ) {

  # sample size
  N <- nrow(data)

  # select "wq" best distance between central and side knots
  w_list <- vector("list", length(wq))

  for (l in 1:length(wq)) {

    # new matrix of basis functions
    B <- matrix(data = rep(1, N), ncol = 1)

    # quintic knots
    quintic_knots <- vector("list", nX)

    paired <- c(); not_paired <- c()
    for (v in 1:nX){
      if (is.null(x_space[[v]])) next

      # from first midpoint: position 2
      # to penultimate midpoint: position (-3)
      # step 2 to select midpoints

      for (i in seq(2, length(x_space[[v]]) - 3, 2)) {

        # select a central knot: position i + 1
        t <- x_space[[v]][i + 1]

        # side of that knot
        side <- knots[knots[, "t"] == t, "side"]

        # triplet of knots
        triplet <- set_triplet_knots (
          knots = x_space[[v]][i:(i + 2)],
          w = wq[l],
          smoothing = "quintic",
          monotonicity = monotonicity,
          concavity = concavity
        )

        # update B with two new truncated quintci basis functions
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

    if (monotonicity || concavity) {

      # number of paired basis functions
      h <- sum(duplicated(knots[, c("xi", "t")])) * 2
      # number of unpaired right side basis functions
      r <- nrow(knots) - h

      # Reorder B to fit the positions of estimate_coefficients_smoothed
      B <- B[, c(1, paired, not_paired)]

      # estimate coefficients
      coefs <- estimate_coefficients_smoothed (
        B = B,
        y = data[, y, drop = F],
        h = h,
        r = r,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0
      )

    } else {

      coefs <- estimate_coefficients (
        B = B[, c(1, paired, not_paired)],
        y = data[, y, drop = F],
        it_list = NULL,
        Bp_list = NULL,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0
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
      B = B,
      d = d,
      knots = sum(lengths(quintic_knots))
    )

    # save results
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

#' @title Generate a new pair of Quintic Basis Functions
#'
#' @description This function generates two new quintic basis functions from a variable and a triplet of knots.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param xi Variable index of the variable that creates the new basis function(s).
#' @param knots Triplet of knots for creating the new basis function(s).
#' @param B A \code{matrix} of basis functions.
#' @param side Side of the new basis function.
#'
#' @return A \code{matrix} of basis functions updated with the new quintic basis functions.

create_quintic_basis <- function (
    data, xi, knots, B, side
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

  # both or right
  if (length(side) == 2 || side == "R") {

    term1 <- alpha1 * (data[, xi] - t0) ^ 3
    term2 <- beta1  * (data[, xi] - t0) ^ 4
    term3 <- gamma1 * (data[, xi] - t0) ^ 5

    # Q1
    Q1 <- ifelse(data[, xi] <= t0,
                     0,
                     (ifelse((data[, xi] > t0) & (data[, xi] < t2),
                             term1 + term2 + term3,
                             data[, xi] - t1)))
    B  <- cbind(B, Q1)
  }

  # both or left
  if (length(side) == 2 || side == "L") {

    term1 <- alpha2 * (data[, xi] - t2) ^ 3
    term2 <- beta2  * (data[, xi] - t2) ^ 4
    term3 <- gamma2 * (data[, xi] - t2) ^ 5

    # Q2
    Q2 <- ifelse(data[, xi] <= t0,
                 t1 - data[, xi],
                 (ifelse((data[, xi] > t0) & (data[, xi] < t2),
                         term1 + term2 + term3,
                         0)))
    B  <- cbind(B, Q2)
  }

  return(B)
}
