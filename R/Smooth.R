#' @title Smoothing Multivariate Adaptive Frontier Splines
#'
#' @description This function smoothes the MAFS predictor.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX number of inputs in \code{data}.
#' @param knots \code{data.frame} containing knots from Backward MAFS.
#' @param y output indexes in \code{data}.
#'
#' @return List containing the set of knots from backward (\code{knots}), the new cubic knots (\code{cubic_knots}) and the set of coefficients (\code{alpha}).
SmoothMAFS <- function(data, nX, knots, y) {
  N <- nrow(data)

  # New matrix of basis functions
  B <- matrix(data = rep(1, N), ncol = 1)

  # Cubic Knots
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

      # Create two new truncated cubic functions
      B <- CreateCubicBF(data, xi, anova[i:(i + 2)], B, side)

      if (length(side) == 1) {
        not_paired <- c(not_paired, ncol(B))
        status <- "not paired"

      } else {
        paired <- c(paired, (ncol(B) - 1):ncol(B))
        status <- "paired"

      }

      # Update cubic knots
      Cknots[[xi]] <- append(Cknots[[xi]], list(list(t = anova[i:(i + 2)],
                                                     status = status)))
    }
  }

  # Estimate coefficients
  h <- sum(duplicated(knots[, "t"])) * 2
  r <- nrow(knots) - h

  # Reorder B to fit the positions of EstimCoeffsBackward [coefs paired, coefs not paired]
  B <- B[, c(1, paired, not_paired)]

  alpha <- EstimCoeffsBackward(B, data[, y], h, r)

  Smooth_Backward_MAFS <- list(knots = knots,
                                cubic_knots = Cknots,
                                alpha = alpha)

  return(Smooth_Backward_MAFS)

}

#' @title Smoothing (Forward) Multivariate Adaptive Frontier Splines
#'
#' @description This function smoothes the Forward MAFS predictor.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX number of inputs in \code{data}.
#' @param knots \code{data.frame} containing knots from Forward MAFS.
#' @param y output indexes in \code{data}.
#'
#' @return List containing the set of knots from backward (\code{knots}), the new cubic knots (\code{cubic_knots}) and the set of coefficients (\code{alpha}).
SmoothForwardMAFS <- function(data, nX, knots, y) {

  N <- nrow(data)

  # New matrix of basis functions
  B <- matrix(data = rep(1, N), ncol = 1)

  # Cubic Knots
  Cknots <- vector("list", nX)

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

      # Sides of that knot. Always paired
      side <- c("R", "L")

      # Create two new truncated cubic functions
      B <- CreateCubicBF(data, xi, anova[i:(i + 2)], B, side)

      # Update cubic knots
      Cknots[[xi]] <- append(Cknots[[xi]], list(list(t = anova[i:(i + 2)],
                                                     status = "paired")))
    }
  }

  alpha <- EstimCoeffsForward(B, y)

  Smooth_Forward_MAFS <- list(knots = knots,
                               cubic_knots = Cknots,
                               alpha = alpha)

  return(Smooth_Forward_MAFS)
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

  t0 <- knt[1] # t-
  t1 <- knt[2] # t
  t2 <- knt[3] # t+

  # Both or right
  if (length(side) == 2 || side == "R") {
    p2 <- (2 * t2 + t0 - 3 * t1) / (t2 - t0) ^ 2 # p+
    r2 <- (2 * t1 - t2 - t0) / (t2 - t0) ^ 3 # r+

    # (x-t)+
    hinge1 <- ifelse(data[, xi] <= t0,
                     0,
                     (ifelse((data[, xi] > t0) & (data[, xi] < t2),
                             p2 * (data[, xi] - t0) ^ 2 + r2 * (data[, xi] - t0) ^ 3,
                             data[, xi] - t1)))


    B <- cbind(B, hinge1)

  }

  # Both or left
  if (length(side) == 2 || side == "L") {
    p0 <- (3 * t1 - 2 * t0 - t2) / (t0 - t2) ^ 2 # p-
    r0 <- (t0 + t2 - 2 * t1) / (t0 - t2) ^ 3 # r-

    # (t-x)+
    hinge2 <- ifelse(data[, xi] <= t0,
                     -(data[, xi] - t1),
                     (ifelse((data[, xi] > t0) & (data[, xi] < t2),
                             p0 * (data[, xi] - t2) ^ 2 + r0 * (data[, xi] - t2) ^ 3,
                             0)))

    B <- cbind(B, hinge2)
  }

  return(B)
}
