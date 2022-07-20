#' @title Smoothing Multivariate Adaptive Frontier Splines through Quintic Functions
#'
#' @description This function smoothes the MAFS predictor and imposes continuous second derivatives.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX number of inputs in \code{data}.
#' @param knots \code{data.frame} containing knots from Backward MAFS.
#' @param y output indexes in \code{data}.
#'
#' @return List containing the set of knots from backward (\code{knots}), the new cubic knots (\code{cubic_knots}) and the set of coefficients (\code{alpha}).
QSmoothMAFS <- function(data, nX, knots, y) {
  N <- nrow(data)

  # New matrix of basis functions
  B <- matrix(data = rep(1, N), ncol = 1)

  # Quintic Knots
  Qknots <- vector("list", nX)

  paired <- c()
  not_paired <- c()
  for (xi in 1:nX){

    # Knots for the xi variable
    kt_xi <- sort(unique(knots[knots[, 1] == xi, 2]))

    if (length(kt_xi) == 0) next

    # Add the initial and the end observation. They cannot be used as knots.
    anova <- c(min(data[, xi]), kt_xi, max(data[, xi]))

    # Calculate Midpoints
    anova <- sort(c(anova, anova[- length(anova)] + diff(anova) / 2))

    # From first midpoint: position 2
    # To penultimate midpoint: position (-3)
    # Step 2 to select midpoints
    for (i in seq(2, length(anova) - 3, 2)){

      # Select knot: position i + 1
      t <- anova[i + 1]

      # Sides of that knot
      side <- knots[knots[, "t"] == t, "side"]

      # Create two new truncated quintic functions
      # New constrained knots
      B_Qknots <- CreateQuinticBF(data, xi, anova[i:(i + 2)], B, side)

      B       <- B_Qknots[["B"]]
      SQknots <- B_Qknots[["Qknots"]]

      if (length(side) == 1) {
        not_paired <- c(not_paired, ncol(B))
        status     <- "unpaired"

      } else {
        paired <- c(paired, (ncol(B) - 1):ncol(B))
        status <- "paired"

      }

      # Update quintic knots
      Qknots[[xi]] <- append(Qknots[[xi]], list(list(t = SQknots,
                                                     status = status)))
    }
  }

  # Estimate coefficients
  h <- sum(duplicated(knots[, "t"])) * 2
  r <- nrow(knots) - h

  # Reorder B to fit the positions of EstimCoeffsBackward [coefs paired, coefs unpaired]
  B <- B[, c(1, paired, not_paired)]

  alpha <- EstimCoeffsBackward(B, data[, y, drop = F], h, r)

  QSmoothMAFS <- list(knots  = knots,
                      Qknots = Qknots,
                      alpha  = alpha)

  return(QSmoothMAFS)
}

#' @title Generate a new pair of Quintic Basis Functions
#'
#' @description This function generates two new quintic basis functions from a variable and a knot previously created during MAFS algorithm.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param xi Variable index of the new basis function(s).
#' @param knt Knots for creating the new basis function(s).
#' @param B Matrix of basis functions.
#' @param side Side of the basis function.
#'
#' @return Matrix of basis functions updated with the new basis functions.
CreateQuinticBF <- function(data, xi, knt, B, side){

  # Unrestricted knots
  t0 <- knt[1] # t-
  t1 <- knt[2] # t
  t2 <- knt[3] # t+

  # Unrestricted deltas
  d  <- t2 - t0 # t+ - t-
  d1 <- t2 - t1 # t+ - t
  d2 <- t1 - t0 # t  - t-

  # New constrained knots and deltas
  if (d1 <= (8 / 7) * d2) {
    d2 <- (7 / 8) * d1
    t0 <- t1 - d2 # t-
  } else if (d1 >= (3 / 2) * d2) {
    d1 <- (3 / 2) * d2
    t2 <- d1 + t1 # t+
  }

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

#' @title Smoothing Multivariate Adaptive Frontier Splines through Cubic Functions
#'
#' @description This function smoothes the MAFS predictor and imposes continuous first derivatives.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX number of inputs in \code{data}.
#' @param knots \code{data.frame} containing knots from Backward MAFS.
#' @param y output indexes in \code{data}.
#'
#' @return List containing the set of knots from backward (\code{knots}), the new cubic knots (\code{cubic_knots}) and the set of coefficients (\code{alpha}).
CSmoothMAFS <- function(data, nX, knots, y) {
  N <- nrow(data)

  # New matrix of basis functions
  B <- matrix(data = rep(1, N), ncol = 1)

  # Quintic Knots
  Cknots <- vector("list", nX)

  paired <- c()
  not_paired <- c()
  for (xi in 1:nX){

    # Knots for the xi variable
    kt_xi <- sort(unique(knots[knots[, 1] == xi, 2]))

    if (length(kt_xi) == 0) next

    # Add the initial and the end observation. They cannot be used as knots.
    anova <- c(min(data[, xi]), kt_xi, max(data[, xi]))

    # Calculate Midpoints
    anova <- sort(c(anova, anova[- length(anova)] + diff(anova) / 2))

    # From first midpoint: position 2
    # To penultimate midpoint: position (-3)
    # Step 2 to select midpoints
    for (i in seq(2, length(anova) - 3, 2)){

      # Select knot: position i + 1
      t <- anova[i + 1]

      # Sides of that knot
      side <- knots[knots[, "t"] == t, "side"]

      # Create two new truncated quintic functions
      # New constrained knots
      B_Cknots <- CreateCubicBF(data, xi, anova[i:(i + 2)], B, side)

      B       <- B_Cknots[["B"]]
      SCknots <- B_Cknots[["Cknots"]]

      if (length(side) == 1) {
        not_paired <- c(not_paired, ncol(B))
        status     <- "unpaired"

      } else {
        paired <- c(paired, (ncol(B) - 1):ncol(B))
        status <- "paired"

      }

      # Update cubic Cknots
      Cknots[[xi]] <- append(Cknots[[xi]], list(list(t = SCknots,
                                                     status = status)))
    }
  }

  # Estimate coefficients
  h <- sum(duplicated(knots[, "t"])) * 2
  r <- nrow(knots) - h

  # Reorder B to fit the positions of EstimCoeffsBackward [coefs paired, coefs unpaired]
  B <- B[, c(1, paired, not_paired)]

  alpha <- EstimCoeffsBackward(B, data[, y, drop = F], h, r)

  CSmoothMAFS <- list(knots  = knots,
                      Cknots = Cknots,
                      alpha  = alpha)

  return(CSmoothMAFS)
}

#' @title Generate a new pair of Cubic Basis Functions
#'
#' @description This function generates two new cubic basis functions from a variable and a knot previously created during MAFS algorithm.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param xi Variable index of the new basis function(s).
#' @param knt Knots for creating the new basis function(s).
#' @param B Matrix of basis functions.
#' @param side Side of the basis function.
#'
#' @return Matrix of basis functions updated with the new basis functions.
CreateCubicBF <- function(data, xi, knt, B, side){

  # Unrestricted knots
  t0 <- knt[1] # t-
  t1 <- knt[2] # t
  t2 <- knt[3] # t+

  # epsilon & delta
  d <- t1 - t0 # t  - t-
  e <- 2 * d   # t+ - t

  # New constrained knots
  # To avoid t- <= 0
  t2 <- t1 + e # t+

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

