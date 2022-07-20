#' @title Backward Algorithm for Multivariate Adaptive Frontier Splines.
#'
#' @description This function creates a portfolio of sub-models by removing basis functions one to one.
#'
#' @param data \code{Data.frame} or \code{matrix} containing the variables in the model.
#' @param y Column output index in \code{data}.
#' @param ForwardModel \code{list} containing the Forward MAFS model.
#' @param d Integer. Generalized Cross Validation (GCV) penalty per knot. If set to \code{-1}, \code{GCV = RSS / n}.
#'
#' @return \code{list} containing a portfolio of MAFS sub-models.
PruningMAFS <- function(data, y, ForwardModel, d) {

  # terms
  terms <- length(ForwardModel[["BF"]]) - 1

  # set of MAFS models with p terms
  models <- list(list(
    id    = NA,
    B     = NA,
    LOF   = NA,
    t     = NA,
    alpha = NA
    ))

  # ==
  # id
  # ==
  models[[1]][["id"]] <- length(ForwardModel[["BF"]])

  # ==
  # B
  # ==
  models[[1]][["B"]] <- ForwardModel[["B"]]

  # ==
  # knots
  # ==
  knots <- lapply(ForwardModel[["BF"]][- 1], function(x) list(
    xi   = x[["xi"]],
    t    = x[["t"]],
    side = x[["side"]]
  ))

  models[[1]][["t"]] <- knots

  # ==
  # LOF
  # ==
  # Predictions
  y_hat <- ForwardModel[["B"]] %*% ForwardModel[["alpha"]]

  # mse
  num <- mse(data[, y, drop = F], y_hat[drop = F])

  if (d != - 1) {
    len.knots <- length(unique(lapply(knots, function(x) list(
      xi = x[["xi"]],
      t = x[["t"]])
    )))

    CM  <- ncol(ForwardModel[["B"]]) + d * len.knots
    den <- (1 - (CM / nrow(data))) ^ 2

  } else {
    den <- 1
  }

  models[[1]][["LOF"]] <- num / den

  # ==
  # alpha
  # ==
  models[[1]][["alpha"]] <- ForwardModel[["alpha"]]

  # ======= #
  # Pruning #
  # ======= #

      B <- ForwardModel[["B"]]
  SetBF <- ForwardModel[["BF"]]

  while (terms > 0){

    # error
    err <- Inf

    # Basis functions to be removed excluding:
    # - B1(X) = 1
    # - Right-side hinge functions when appearing unpaired
    # - Left-side hinge function from reflected pairs
    NBFs <- unlist(lapply(SetBF, function(x) if((x["status"] == "paired" && x[["side"]] == "R") ||
                                                (x["status"] == "unpaired" && x[["side"]] == "L"))
                                                 x[['id']]))

    for (bf in NBFs) {
      # Update SetBF
      New.SetBF <- SetBF

      # Drop element
      dropTerm <- sapply(New.SetBF, function(x) x[["id"]] == bf)

      # Update B
      # Ensure matrix format when B is just a vector
      New.B <- B[, !dropTerm]

      # After Forward Algorithm: paired | R --> Always even number
      if (!bf %% 2){
        New.SetBF[lapply(New.SetBF, "[[", "id") == bf + 1][[1]][["status"]] <- "unpaired"
      }

      # Drop Basis Function
      New.SetBF[dropTerm] <- NULL

      # Number of paired basis functions
      h <- length(unlist(lapply(New.SetBF, function(x) if((x["status"] == "paired"))
                                                           x[['id']])))

      # Number of left-side unpaired basis functions
      r <- length(unlist(lapply(New.SetBF, function(x) if((x["status"] == "unpaired"))
                                                           x[['id']])))

      # Estimate coefficients
      # Predict y_hat
      if (terms > 1){

        not.paired <- (1:terms)[sapply(New.SetBF, function(x) x[["status"]] == "unpaired")]
        paired     <- (1:terms)[sapply(New.SetBF, function(x) x[["status"]] != "unpaired")]
        colsOrder  <- c(paired, not.paired)

        coefs <- EstimCoeffsBackward(New.B[, colsOrder], data[, y, drop = F], h, r)

        y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

        for (out in 1:length(y)) {
          y_hat[, out] <- New.B[, colsOrder] %*% coefs[, out]
        }

      } else {
        colsOrder <- 1

        # Last iteration. Need to transform New.B from vector to matrix
        coefs <- matrix(apply(data[, y, drop = FALSE], 2, max), ncol = length(y))
        New.B <- matrix(New.B)

        y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

        for (out in 1:length(y)) {
          y_hat[, out] <- New.B %*% coefs[, out, drop = F]
        }
      }

      # mse
      num <- mse(data[, y, drop = F], y_hat[drop = F])

      # Knots
      knots <- lapply(New.SetBF[colsOrder][- 1],
                      function(x) if(all(x[["xi"]] != - 1)) list(
                        xi     = x[['xi']],
                        t      = x[['t']],
                        side   = x[["side"]],
                        status = x[["status"]]
                      ))
      # LOF
      if (d != - 1) {

        len.knots <- length(unique(lapply(knots, function(x) list(
          xi = x[["xi"]],
           t = x[["t"]])
          )))

        CM  <- ncol(New.B) + d * len.knots
        den <- (1 - (CM / nrow(data))) ^ 2

      } else {
        den <- 1
      }

      LOF <- num / den

      if (LOF < err) {
        Best.SetBF <- New.SetBF[colsOrder]
            Best.B <- New.B[, colsOrder]
               err <- LOF
        Best.knots <- knots
        Best.coefs <- coefs
      }
    }

    # Update set of basis functions and matrix of basis functions
    SetBF <- Best.SetBF
        B <- Best.B

    # Create P-model
    P_Model <- list(
      id    = ncol(Best.B),
      B     = Best.B,
      LOF   = err,
      t     = Best.knots,
      alpha = Best.coefs
      )

    # Best model of p terms
    models <- append(models, list(P_Model))

    # Update number of terms
    terms <- terms - 1
  }

  return(models)
}

#' @title Estimate Coefficients in Multivariate Adaptive Frontier Splines during Backward Procedure.
#'
#' @description This function solves a Linear Programming Problem to obtain a set of coefficients.
#'
#' @param B \code{matrix} of basis functions.
#' @param y Output \code{vector} in data.
#' @param h Number of coefficients associated with paired basis functions.
#' @param r Number of coefficients associated with unpaired left-side basis functions.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return \code{vector} with the coefficients estimated.
EstimCoeffsBackward <- function(B, y, h, r){

  n <- nrow(B)
  p <- ncol(B) # (p = 1 + h + r)

  alpha <- matrix(NA, nrow = p, ncol = ncol(y))

  for (out in 1:ncol(y)) {

    # Select the variable
    y.ind <- y[, out]

    # vars: c(tau_0, alpha_1, beta_2, ... , alpha_(H-1), beta_H, w_1, ..., w_R, e_1, ... , e_n)
    objVal <- c(rep(0, p), rep(1, n))

    # = #
    # A #
    # = #
    # Equality constraints: y_hat - e = y
    Amat1 <- cbind(B, diag(rep(-1, n), n))

    # Concavity [paired basis functions]
    Amat2 <- matrix(0, nrow = h / 2, ncol = p + n)

    if (h > 0) {
      for (i in seq(2, h, 2)){
        Amat2[i / 2, c(i, i + 1)] <- - 1
      }
    }

    # Increasing monotony [paired basis functions]
    Amat3 <- cbind(rep(0, h), diag(x = c(1, -1), h), matrix(0, h, r), matrix(0, h, n))

    # Concavity & increasing monotony [left-side basis functions]
    Amat4 <- cbind(matrix(0, r, 1 + h), diag(- 1, r), matrix(0, r, n))

    Amat  <- rbind(Amat1, Amat2, Amat3, Amat4)

    # = #
    # b #
    # = #
    bvec <- c(y.ind, rep(0, h / 2), rep(0, r), rep(0, h))

    # Direction of inequality
    dir <- c(rep("==", n), rep(">=", nrow(Amat) - n))

    # Bounds
    bounds <- list(lower = list(ind = 1L:p, val = rep(-Inf, p)))

    # Solve
    sol <- Rglpk_solve_LP(objVal, Amat, dir, bvec, bounds)

    alpha[, out] <- sol$solution[1:p]
  }

  return(alpha)
}

#' @title Add interaction of variables
#'
#' @description This function adds interaction of variables in the model.
#'
#' @param data \code{Data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param model MAFS model from backward algorithm.
#' @param knots Knots from (backward) MAFS algorithm.
#' @param err.red Minimum reduced error rate for the addition of two new basis functions. Default is \code{0.01}.
#' @param Kp Maximum degree of interaction allowed. Default is \code{1}.
#' @param d Integer. Generalized Cross Validation (GCV) penalty per knot. If set to \code{-1}, \code{GCV = RSS / n}.
#'
#' @return \code{vector} with the coefficients estimated.
InteractionModel <- function(data, x, y, model, knots, err.red, Kp, d){

  # Current B matrix
  B <- model[["B"]]

  # Current MSE
  CM  <- ncol(B) + d * length(unique(knots$t))
  den <- (1 - (CM / nrow(data))) ^ 2

  err <- model[["LOF"]] * den

  # Nº of terms
  terms <- ncol(B)

  # Paired knots
  h <- sum(knots[, "status"] == "paired")
  r <- sum(knots[, "status"] == "unpaired")

  # Unpaired cols
  ucols <- (terms - r + 1):terms

  # Unpaired knots
  uknots <- knots[knots[, "status"] == "unpaired", ]

  # NewB to update B
  NewB <- B

  # New model
  for (col in ucols) {
    for (kn in 1:nrow(uknots)) {
      # knots (based on observed data)
      knotvar   <- uknots[kn, "xi"]
      knotval   <- uknots[kn, "t"]

      knotsGrid <- data[data[, knotvar] <= knotval, ]

      for (var in x[x != knotvar]) {
        for (i in 1:nrow(knotsGrid)) {

          knt <- knotsGrid[i, var]

          NewB[, col] <- B[, col] * ifelse(data[, var] < knt, knt - data[, var], 0)

          # New model
          coefs <- EstimCoeffsBackward(NewB, data[, y, drop = F], h, r)

          # New predictions
          y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

          for (out in 1:length(y)) {
            y_hat[, out] <- NewB %*% coefs[, out, drop = F]
          }

          # mse
          new.err <- mse(data[, y, drop = F], y_hat[drop = F])

          if (new.err < err * (1 - err.red)) {

            # Update model
            model[["B"]] <- NewB

            # Update knot
            model[["t"]][[col - 1]][["xi"]] <- c(knotvar, var)
            model[["t"]][[col - 1]][["t"]]  <- c(knotval, knt)
          }
        }
      }
    }
  }
}

