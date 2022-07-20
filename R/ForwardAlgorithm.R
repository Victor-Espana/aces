#' @title Mean Squared Error
#'
#' @description This function computes the mean squared error between two numeric vectors.
#'
#' @param y Vector of actual data.
#' @param y_pred Vector of predicted values.
#'
#' @return Mean Squared Error.
mse <- function(y, y_pred){

  error <- sum(y_pred - y) / (nrow(y) * ncol(y))
  return(round(error, 4))
}

#' @title Add a New Pair of Basis Functions
#'
#' @description This function adds the best pair of basis functions to the model.
#'
#' @param data data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param ForwardModel \code{list} containing the set of basis functions and the B matrix.
#' @param knots.list \code{list} containing the set of selected knots.
#' @param Kp Maximum degree of interaction allowed.
#' @param L Minimum number of observations between two adjacent knots.
#' @param Le Minimum number of observations before the first and after the final knot.
#' @param knotsGrid Grid of virtual knots to perform MAFS. Can be set to \code{NULL}.
#' @param linpreds \code{logical}. If \code{TRUE}, predictors can enter linearly.
#' @param err.min Minimun error in the split.
#'
#' @return A \code{list} containing the matrix of basis functions (\code{B}), a \code{list} of basis functions (\code{BF}), a \code{list} of selected knots (\code{knots.list}) and the minimun error (\code{err.min}).
AddBF <- function(data, x, y, ForwardModel, knots.list,
                  Kp, L, Le, knotsGrid, linpreds, err.min) {

  N  <- nrow(data)
  nX <- length(x)

  # Set of basis functions
  BF.list <- ForwardModel[["BF"]]
  nBF     <- length(BF.list)

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
    # the division cannot be carried out by this basis function.
    if(length(bf[['xi']]) >= 1 && !all(bf[['xi']] == -1)) next # 1 = Kp

    for (xi in 1:nX) {

      knots <- setKnots(data, nX, xi, L, Le, knots.list, bf, knotsGrid, linpreds)

      if (is.null(knots)) next

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

        # mse
        err <- mse(data[, y, drop = F], y_hat[drop = F])

        if (err < err.min) {

          # Model has improved
          signal  <- 1
          err.min <- err
          Best.B  <- New.B

          # index
          if (is.null(knotsGrid)) {
            tindex <- which(data[, xi] == knots[i])

          } else {
            tindex <- which(knotsGrid[[xi]] == knots[i])
          }

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

          # alpha
          bf1[['alpha']] <- bf2[['alpha']] <- coefs
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

#' @title Create the set of knots
#'
#' @description This function generates the vector of knots to create a new pair of basis functions.
#'
#' @param data data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param nX \code{integer}. Number of inputs.
#' @param xi \code{integer}. Index of the variable that creates the basis function.
#' @param L Minimum number of observations between two adjacent knots.
#' @param Le Minimum number of observations before the first and after the final knot.
#' @param knots.list \code{list} containing the set of selected knots.
#' @param bf \code{list}. Basis function for the expansion of the model.
#' @param knotsGrid Grid of virtual knots to perform MAFS. Can be set to \code{NULL}.
#' @param linpreds \code{logical}. If \code{TRUE}, predictors can enter linearly.
#'
#' @importFrom dplyr %>%
#'
#' @return Numeric vector with the knots values.
setKnots <- function(data, nX, xi, L, Le, knots.list, bf, knotsGrid, linpreds) {

  # Minimum number of observations between two adjacents knots
  if (length(L) > 1) {
    sp1 <- L[xi]
  } else {
    sp1 <- L
  }

  if (length(Le) > 1) {
    sp2 <- Le[xi]
  } else {
    sp2 <- Le
  }

  # ===================== #
  # Condition II: product #
  # ===================== #

  # The same xi cannot appear twice in a product.
  if (any(bf[["xi"]] == xi)) knots <- NULL

  # Observations in the space.
  if (is.null(knotsGrid)) {
    index <- bf[['Bp']] != 0
    N     <- length(index)

  } else {
    index      <- knotsGrid[[xi]] != 0
    N          <- length(index)
    data.index <- bf[['Bp']] != 0
  }

  # Minimum span for Friedman approach
  if (sp1 == - 1) {
    sp1 <- ceiling(- log2(- (1 / (nX * sum(index))) * log(1 - 0.05)) / 2.5)
  }

  # ===================== #
  # Condition III: L & Le #
  # ===================== #

  # Set to FALSE first and last "Le" observations
  index[c(1:sp2, (length(index) - sp2 + 1):length(index))] <- FALSE

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
    setL.idx        <- unique(sort(setL.idx[setL.idx > 0 & setL.idx <= N]))
    index[setL.idx] <- FALSE
  }

  if (all(index == FALSE)) knots <- NULL

  # Knots
  if(is.null(knotsGrid)) {
    knots <- data[order(data[, xi]), xi][index] %>% unique()

  } else {
    minim    <- min(data[order(data[, xi]), xi][data.index])
    maxim    <- max(data[order(data[, xi]), xi][data.index])

    varKnots <- knotsGrid[[xi]][index]

    knots    <- varKnots[varKnots >= minim & varKnots <= maxim]
  }

  if (linpreds == TRUE) {
    knots <- c(0, knots)
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
#' @return Matrix of basis function (\code{B}) updated with the new basis functions.
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

#' @title Estimate Coefficients in Multivariate Adaptive Frontier Splines during Forward Procedure.
#'
#' @description This function solves a Linear Programming Problem to obtain a set of coefficients.
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

  alpha <- matrix(NA, nrow = p, ncol = ncol(y))

  for (out in 1:ncol(y)) {

    # Select variable
    y.ind <- y[, out]

    # vars: c(alpha_0, alpha_1, ..., alpha_P, e_1, ... , e_n)
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
    a0   <- rep(0, p - 1)
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

    alpha[, out] <- sols$solution[1:p]
  }

  return(alpha)
}
