#' @title Smoothing Additive Adaptive Frontier Splines through Quintic Functions
#'
#' @description This function smoothes the AAFS predictor and imposes continuous second derivatives.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX number of inputs in \code{data}.
#' @param knots \code{data.frame} containing knots from Backward AAFS.
#' @param x.space \code{list} with the input space of each variable for smoothing.
#' @param y output indexes in \code{data}.
#' @param w The distance between the central knot and a side knot is at most \code{w} times the distance between that central knot and the other side knot for each truncated quintic basis function.
#'
#' @return List containing the set of knots from backward (\code{knots}), the new quintic knots (\code{quintic_knots}) and the set of coefficients (\code{coefs}).
QSAAFS <- function(data, nX, knots, x.space, y, w) {
  # Sample size
  N <- nrow(data)

  # Select best distance between central and side knots (w)
  w.list <- vector("list", length(w))

  for (l in 1:length(w)) {
    # New matrix of basis functions
    B <- matrix(data = rep(1, N), ncol = 1)

    # Cubic Knots
    Qknots <- vector("list", nX)

    paired <- c()
    not_paired <- c()
    for (xi in 1:nX){
      if (is.null(x.space[[xi]])) next
      # From first midpoint: position 2
      # To penultimate midpoint: position (-3)
      # Step 2 to select midpoints
      for (i in seq(2, length(x.space[[xi]]) - 3, 2)){

        # Select knot: position i + 1
        t <- x.space[[xi]][i + 1]

        # Sides of that knot
        side <- knots[knots[, "t"] == t, "side"]

        # Create two new truncated quintic functions
        # New constrained knots
        B_Qknots <- CreateQuinticBF(data, xi, x.space[[xi]][i:(i + 2)], B, side, w[l])

        B <- B_Qknots[["B"]]
        SQknots <- B_Qknots[["Qknots"]]

        if (length(side) == 1) {
          not_paired <- c(not_paired, ncol(B))
          status <- "unpaired"

        } else {
          paired <- c(paired, (ncol(B) - 1):ncol(B))
          status <- "paired"
        }

        # Update cubic knots
        Qknots[[xi]] <- append(Qknots[[xi]], list(list(t = SQknots, status = status)))
      }
    }

    # Estimate coefficients
    h <- sum(duplicated(knots[, "t"])) * 2
    r <- nrow(knots) - h

    # Reorder B to fit the positions of EstimCoeffsBackward [coefs paired, coefs unpaired]
    B <- B[, c(1, paired, not_paired)]

    coefs <- EstimCoeffsBackward(B, data[, y, drop = F], h, r)

    # Predictions
    y_hat <- matrix(NA, nrow = N, ncol = length(y))

    for (out in 1:length(y)) {
      y_hat[, out] <- B %*% coefs[, out]
    }

    err <- mean(unlist(me(data[, y, drop = F], y_hat[drop = F])))

    # Save results
    w.list[[l]][["B"]] <- B
    w.list[[l]][["Qknots"]] <- Qknots
    w.list[[l]][["w"]] <- w[l]
    w.list[[l]][["coefs"]] <- coefs
    w.list[[l]][["err"]] <- err
  }

  # Lack-of-fit at each model
  errors <- sapply(w.list, function(x) x[["err"]])

  # Model with minimum error
  SQAAFS <- w.list[[which.min(errors)]]

  return(SQAAFS)
}

#' @title Generate a new pair of Quintic Basis Functions
#'
#' @description This function generates two new quintic basis functions from a variable and a knot previously created during AAFS algorithm.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param xi Variable index of the new basis function(s).
#' @param knt Knots for creating the new basis function(s).
#' @param B Matrix of basis functions.
#' @param side Side of the basis function.
#' @param wl The distance between the central knot and a side knot is at most \code{wl} times the distance between that central knot and the other side knot. It must be between 8/7 and 1.5.
#'
#' @return Matrix of basis functions updated with the new basis functions.
CreateQuinticBF <- function(data, xi, knt, B, side, wl){

  # Unrestricted knots
  t0 <- knt[1] # t-
  t1 <- knt[2] # t
  t2 <- knt[3] # t+

  # Restricted deltas
  d2 <- t1 - t0 # t  - t-
  d1 <- wl * d2 # t+ - t
  t2 <- t1 + d1
  d  <- t2 - t0 # t+ - t-

  alpha1 <- (6 * d1 - 4 * d2) / d ^ 3
  alpha2 <- (4 * d1 - 6 * d2) / d ^ 3

  beta1 <- (7 * d2 - 8 * d1) / d ^ 4
  beta2 <- (7 * d1 - 8 * d2) / d ^ 4

  gamma1 <- gamma2 <- (3 * d1 - 3 * d2) / d ^ 5

  # Both or right
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


    B <- cbind(B, Q1)
  }

  # Both or left
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

    B <- cbind(B, Q2)
  }

  B_Qknots <- list(B = B, Qknots = c(t0, t1, t2))

  return(B_Qknots)
}

#' @title Smoothing Additive Adaptive Frontier Splines through Cubic Functions
#'
#' @description This function smoothes the AAFS predictor and imposes continuous first derivatives.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX number of inputs in \code{data}.
#' @param knots \code{data.frame} containing the knots from Backward AAFS.
#' @param x.space \code{list} with the input space of each variable for smoothing.
#' @param y output indexes in \code{data}.
#' @param w The distance between the central knot and a side knot is at most \code{w} times the distance between that central knot and the other side knot for each truncated cubic basis function.
#'
#' @return List containing the set of knots from backward (\code{knots}), the new cubic knots (\code{cubic_knots}) and the set of coefficients (\code{coefs}).
CSAAFS <- function(data, nX, knots, x.space, y, w) {
  # Sample size
  N <- nrow(data)

  # Select best distance between central and side knots (w)
  w.list <- vector("list", length(w))

  for (l in 1:length(w)) {
    # New matrix of basis functions
    B <- matrix(data = rep(1, N), ncol = 1)

    # Cubic Knots
    Cknots <- vector("list", nX)

    paired <- c()
    not_paired <- c()
    for (xi in 1:nX){
      if (is.null(x.space[[xi]])) next
      # From first midpoint: position 2
      # To penultimate midpoint: position (-3)
      # Step 2 to select midpoints
      for (i in seq(2, length(x.space[[xi]]) - 3, 2)){

        # Select knot: position i + 1
        t <- x.space[[xi]][i + 1]

        # Sides of that knot
        side <- knots[knots[, "t"] == t, "side"]

        # Create two new truncated cubic functions
        # New constrained knots
        B_Cknots <- CreateCubicBF(data, xi, x.space[[xi]][i:(i + 2)], B, side, w[l])

        B <- B_Cknots[["B"]]
        SCknots <- B_Cknots[["Cknots"]]

        if (length(side) == 1) {
          not_paired <- c(not_paired, ncol(B))
          status <- "unpaired"

        } else {
          paired <- c(paired, (ncol(B) - 1):ncol(B))
          status <- "paired"
        }

        # Update cubic knots
        Cknots[[xi]] <- append(Cknots[[xi]], list(list(t = SCknots, status = status)))
      }
    }

    # Estimate coefficients
    h <- sum(duplicated(knots[, "t"])) * 2
    r <- nrow(knots) - h

    # Reorder B to fit the positions of EstimCoeffsBackward [coefs paired, coefs unpaired]
    B <- B[, c(1, paired, not_paired)]

    coefs <- EstimCoeffsBackward(B, data[, y, drop = F], h, r)

    # Predictions
    y_hat <- matrix(NA, nrow = N, ncol = length(y))

    for (out in 1:length(y)) {
      y_hat[, out] <- B %*% coefs[, out]
    }

    err <- mean(unlist(me(data[, y, drop = F], y_hat[drop = F])))

    # Save results
    w.list[[l]][["B"]] <- B
    w.list[[l]][["Cknots"]] <- Cknots
    w.list[[l]][["w"]] <- w[l]
    w.list[[l]][["coefs"]] <- coefs
    w.list[[l]][["err"]] <- err
  }

  # Lack-of-fit at each model
  errors <- sapply(w.list, function(x) x[["err"]])

  # Model with minimum error
  CSAAFS <- w.list[[which.min(errors)]]

  return(CSAAFS)
}

#' @title Generate a new pair of Cubic Basis Functions
#'
#' @description This function generates two new cubic basis functions from a variable and a knot previously created during AAFS algorithm.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param xi Variable index of the new basis function(s).
#' @param knt Knots for creating the new basis function(s).
#' @param B Matrix of basis functions.
#' @param side Side of the basis function.
#' @param wl The distance between the central knot and a side knot is at most \code{wl} times the distance between that central knot and the other side knot. It must be between 1 and 2.
#'
#' @return Matrix of basis functions updated with the new basis functions.
CreateCubicBF <- function(data, xi, knt, B, side, wl){

  # Unrestricted knots
  t0 <- knt[1] # t-
  t1 <- knt[2] # t
  t2 <- knt[3] # t+

  # epsilon & delta
  d <- t1 - t0 # t  - t-
  e <- t2 - t1 # t+ - t

  if (d < e) {
    e  <- wl * d
    t2 <- t1 + e

  } else if (e < d) {
    d  <- wl * e
    t0 <- t1 - d
  }

  p1 <- (2 * e - d) / (e + d) ^ 2
  r1 <- (d - e) / (e + d) ^ 3

  p2 <- (2 * d - e) / (- e - d) ^ 2
  r2 <- (e - d) / (- e - d) ^ 3

  # Both or right
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

  # Both or left
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

  B_Cknots <- list(B = B, Cknots = c(t0, t1, t2))

  return(B_Cknots)
}

