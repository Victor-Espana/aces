#' @title Mean Error
#'
#' @description This function computes the mean error between two numeric vectors.
#'
#' @param y Vector of actual data.
#' @param y_pred Vector of predicted values.
#'
#' @return Mean Error.
me <- function(y, y_pred){

  error <- vector("list", ncol(y))

  for (out in 1:ncol(y)) {
    error[out] <- round(sum(y_pred[, out] - y[, out]) / nrow(y), 4)
  }

  # error <- sum(y_pred - y) / (nrow(y) * ncol(y))

  return(error)
}

#' @title Add a New Pair of Basis Functions
#'
#' @description This function adds the pair of basis functions that produce the largest reduction in the lack-of-fit.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param ForwardModel \code{list} containing the set of basis functions (\code{BF}) and the B matrix (\code{B}).
#' @param knots.list \code{list} containing the set of selected knots.
#' @param Kp Maximum degree of interaction allowed.
#' @param L Minimum number of observations between two adjacent knots.
#' @param Le Minimum number of observations before the first and after the final knot.
#' @param knotsGrid Grid of virtual knots to perform AAFS.
#' @param err.min Forward algorithm minimum error in the given iteration.
#' @param err.idx Output with the maximum error.
#'
#' @return A \code{list} containing the matrix of basis functions (\code{B}), a \code{list} of basis functions (\code{BF}), a \code{list} of selected knots (\code{knots.list}) and the minimum given error (\code{err.min}).
AddBF <- function(data, x, y, ForwardModel,
                  knots.list, Kp, L , Le,
                  knotsGrid, err.min, err.idx) {

  # Samples in data
  N <- nrow(data)

  # Number of inputs
  nX <- length(x)

  # Set of basis functions
  BF.list <- ForwardModel[["BF"]]

  # Number of basis functions
  nBF <- length(BF.list)

  # Matrix of basis functions
  B <- ForwardModel[["B"]]

  signal <- 0 # improvement

  for (p in 1:nBF){

    # Basis function for the expansion
    bf <- BF.list[[p]]

    # =================== #
    # Condition I: degree #
    # =================== #

    # If the number of variables exceeds the maximum allowed degree,
    # the expansion cannot be carried out by this basis function.
    if(length(bf[['xi']]) >= Kp && !all(bf[['xi']] == -1)) next

    for (xi in 1:nX) {

      # The same xi cannot appear twice in a multivariate basis function.
      if (any(bf[["xi"]] == xi)) next

      # Create grid of knots
      knots <- setKnots(data, nX, xi, L, Le, knots.list, bf, knotsGrid)

      if (is.null(knots)) next

      knots <- sample(knots)

      for (i in 1:length(knots)) {
        # Update B
        New.B <- CreateBF(data, xi, knots[i], B, p)

        # Estimate coefficients
        coefs <- EstimCoeffsForward(New.B, data[, y, drop = F])

        # Predictions
        y_hat <- matrix(NA, nrow = N, ncol = length(y))

        for (out in 1:length(y)) {
          y_hat[, out] <- New.B %*% coefs[, out]
        }

        # mean error
        err <- me(data[, y, drop = F], y_hat[drop = F])

        if (err[[err.idx]] < err.min[[err.idx]]) {

          # Model has improved
          signal   <- 1
          err.min  <- err
          Best.B   <- New.B

          # index
          tindex <- which(knotsGrid[[xi]] == knots[i])

          # New pair of basis functions
          bf1 <- bf2 <-  bf

          # id
          bf1[["id"]] <- nBF + 1
          bf2[["id"]] <- nBF + 2

          # status
          bf1[["status"]] <- bf2[["status"]] <- "paired"

          # side
          bf1[["side"]] <- "R"
          bf2[["side"]] <- "L"

          # Bp
          bf1[['Bp']] <- New.B[, ncol(New.B) - 1]
          bf2[['Bp']] <- New.B[, ncol(New.B)]

          # xi
          if (all(bf[['xi']] == -1)){
            bf1[['xi']] <- bf2[['xi']] <- c(xi)

          } else {
            bf1[['xi']] <- bf2[['xi']] <- c(bf[['xi']], xi)
          }

          # t
          bf1[["t"]] <- bf2[["t"]] <- knots[i]

          # R
          bf1[['R']] <- bf2[['R']] <- err.min

          # coefficients
          bf1[['coefs']] <- bf2[['coefs']] <- coefs
        }
      }
    }
  }

  if (signal) {

    # Append new basis functions
    BF.list <- append(BF.list, list(bf1))
    BF.list <- append(BF.list, list(bf2))

    # Knots for each variable
    len.knt <- length(bf1[["xi"]])
    var.pos <- bf1[["xi"]][len.knt]

    knots.list[[var.pos]] <- append(knots.list[[bf1[['xi']][len.knt]]],
                                    list(list(t = bf1[["t"]], index = tindex)))

    return(list(Best.B, BF.list, knots.list, err.min))

  } else {

    return(list(B, BF.list, knots.list, err.min))
  }
}

#' @title Create the set of eligible knots
#'
#' @description This function generates a vector of knots to create a new pair of basis functions.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX \code{integer}. Number of inputs.
#' @param xi \code{integer}. Index of the variable that creates the new pair of basis functions.
#' @param L Minimum number of observations between two adjacent knots.
#' @param Le Minimum number of observations before the first and after the final knot.
#' @param knots.list \code{list} containing the set of selected knots.
#' @param bf \code{list}. Basis function for the expansion of the model.
#' @param knotsGrid Grid of virtual knots to perform AAFS.
#'
#' @importFrom dplyr %>%
#'
#' @return Numeric vector with the knots values.
setKnots <- function(data, nX, xi, L, Le, knots.list, bf, knotsGrid) {

  # Minimum number of observations between two adjacent knots
  if (length(L) > 1) {
    sp1 <- L[xi]
  } else {
    sp1 <- L
  }

  # Minimum number of observations before the first and after the final knot.
  if (length(Le) > 1) {
    sp2 <- Le[xi]
  } else {
    sp2 <- Le
  }

  # Observations in the space of the basis function.
  nBF <- bf[['Bp']] != 0

  # Minimum span for Friedman's approach
  if (sp1 == - 1) {
    sp1 <- ceiling(- log2(- (1 / (nX * sum(nBF))) * log(1 - 0.05)) / 2.5)
  }

  # ===================== #
  # Condition III: L & Le #
  # ===================== #

  # Boolean: possible knots
  index <- rep(TRUE, length(knotsGrid[[xi]]))

  # Set to FALSE first and last "Le" observations
  if (sp2 != 0) {
    index[c(1:sp2, (length(index) - sp2 + 1):length(index))] <- FALSE
  }

  # Set indexes to FALSE based on minspan
  if (!is.null(knots.list[[xi]])) {
    # Used indexes
    used.idx <- sapply(knots.list[[xi]], "[[", "index")
    setL.idx <- c()

    for (idx in 1:length(used.idx)) {
      # L indexes of distance between knots
      L.index  <- c((used.idx[idx] - sp1):(used.idx[idx] + sp1))
      setL.idx <- append(setL.idx, L.index)
    }

    # Indexes must be greater than 0 and lower than N
    setL.idx <- unique(sort(setL.idx[setL.idx > 0 & setL.idx <= sum(nBF)]))
    index[setL.idx] <- FALSE
  }

  # Knots
  if (all(index == FALSE)) {
    knots <- NULL
  } else {
    # Minimum value of observations in BF
    minim <- min(data[nBF, xi])
    # Maximum value of observations in BF
    maxim <- max(data[nBF, xi])

    # Exclude observations according to minspan and endspan
    varKnots <- knotsGrid[[xi]][order(knotsGrid[[xi]])][index]

    # Knots between minim and maxim
    knots <- varKnots[varKnots >= minim & varKnots <= maxim]
  }

  if (length(knots) <= 1) knots <- NULL

  return(knots)
}

#' @title Generate a new pair of Basis Functions
#'
#' @description This function generates two new basis functions from a variable and a knot.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param xi \code{integer}. Variable index of the new basis function(s).
#' @param knt Knot for creating the new basis function(s).
#' @param B \code{matrix} of basis functions on which the new pair of functions is added.
#' @param p \code{integer}. Parent basis function index.
#'
#' @return Matrix of basis functions (\code{B}) updated with the new pair of basis functions.
CreateBF <- function(data, xi, knt, B, p) {

  # Create (xi-t)+ and (t-xi)+
  hinge1 <- ifelse(data[, xi] > knt, data[, xi] - knt, 0)
  hinge2 <- ifelse(data[, xi] < knt, knt - data[, xi], 0)

  # two new basis functions
  bf1 <- B[, p] * hinge1
  bf2 <- B[, p] * hinge2

  # update B
  B <- cbind(B, bf1, bf2)

  return(B)
}

#' @title Estimate Coefficients in Additive Adaptive Frontier Splines during the Forward Procedure.
#'
#' @description This function solves a Linear Programming Problem to obtain a set of coefficients that impose concavity and monotonocity in the estimator.
#'
#' @param B \code{matrix} of basis functions.
#' @param y Output \code{matrix} in data.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return \code{vector} with the coefficients estimated.
EstimCoeffsForward <- function(B, y){

  n <- nrow(B)
  p <- ncol(B)

  coefs <- matrix(NA, nrow = p, ncol = ncol(y))

  for (out in 1:ncol(y)) {

    # Select variable
    y.ind <- y[, out]

    # vars: c(coef_0, coef_1, ..., coef_P, e_1, ... , e_n)
    objVal <- c(rep(0, p), rep(1, n))

    # = #
    # A #
    # = #
    # Equality constraints: y_hat - e = y
    Amat1 <- cbind(B, diag(rep(- 1, n), n))

    # Concavity
    Amat2 <- matrix(0, nrow = (p - 1) / 2, ncol = p + n)

    pairs <- 1:p
    pairs <- pairs[lapply(pairs, "%%", 2) == 0]

    for (i in pairs){
      Amat2[i / 2, c(i, i + 1)] <- - 1
    }

    # Increasing monotony
      a0 <- rep(0, p - 1)
    alps <- diag(x = c(1, - 1), p - 1)
    mat0 <- matrix(0, p - 1, n)

    Amat3 <- cbind(a0, alps, mat0)

    Amat  <- rbind(Amat1, Amat2, Amat3)

    # = #
    # b #
    # = #
    bvec <- c(y.ind, rep(0, (p - 1) / 2), rep(0, p - 1))

    # Direction of inequality
    dirs <- c(rep("==", n), rep(">=", nrow(Amat) - n))

    # Bounds
    bnds <- list(lower = list(ind = 1L:p, val = rep(- Inf, p)))

    # Solve
    sols <- Rglpk_solve_LP(objVal, Amat, dirs, bvec, bnds)

    coefs[, out] <- sols$solution[1:p]
  }

  return(coefs)
}
