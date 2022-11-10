#' @title Backward Algorithm for Additive Adaptive Frontier Splines.
#'
#' @description This function creates a portfolio of sub-models by removing basis functions one to one.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param y Column output index in \code{data}.
#' @param AAFS.Forward \code{list} containing the Forward AAFS model.
#' @param d Integer. Generalized Cross Validation (GCV) penalty per knot.
#'
#' @return \code{list} containing a portfolio of Additive Adaptive Frontier Splines sub-models.
PruningAAFS <- function(data, y, AAFS.Forward, d) {

  # terms
  terms <- length(AAFS.Forward[["BF"]]) - 1

  # set of AAFS models with p terms
  models <- list(list(
       id = NA,
        B = NA,
      GCV = NA,
        t = NA,
    coefs = NA
    ))

  # ==
  # id
  # ==
  models[[1]][["id"]] <- length(AAFS.Forward[["BF"]])

  # ==
  # B
  # ==
  models[[1]][["B"]] <- AAFS.Forward[["B"]]

  # ==
  # knots
  # ==
  knots <- lapply(AAFS.Forward[["BF"]][- 1], function(x) list(
    xi   = x[["xi"]],
    t    = x[["t"]],
    side = x[["side"]]
  ))

  models[[1]][["t"]] <- knots

  # ==
  # coefficients
  # ==
  models[[1]][["coefs"]] <- AAFS.Forward[["BF"]][[terms + 1]][["coefs"]]

  # ==
  # GCV
  # ==
  models[[1]][["GCV"]] <- computeGCV(data, y, AAFS.Forward[["B"]], models[[1]][["coefs"]], d, knots)

  # ======= #
  # Pruning #
  # ======= #
      B <- AAFS.Forward[["B"]]
  SetBF <- AAFS.Forward[["BF"]]

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
      if (terms > 1){

        not.paired <- (1:terms)[sapply(New.SetBF, function(x) x[["status"]] == "unpaired")]
            paired <- (1:terms)[sapply(New.SetBF, function(x) x[["status"]] != "unpaired")]
         colsOrder <- c(paired, not.paired)

        coefs <- EstimCoeffsBackward(New.B[, colsOrder], data[, y, drop = F], h, r)
        New.B <- New.B[, colsOrder]

      } else {
        colsOrder <- 1

        # Last iteration. Need to transform New.B from vector to matrix
        coefs <- matrix(apply(data[, y, drop = FALSE], 2, max), ncol = length(y))
        New.B <- matrix(New.B)
      }

      # GCV
      knots <- lapply(New.SetBF[colsOrder][- 1],
                      function(x) if(all(x[["xi"]] != - 1)) list(
                          xi = x[['xi']],
                           t = x[['t']],
                        side = x[["side"]],
                      status = x[["status"]]))

      GCV <- computeGCV(data, y, New.B, coefs, d, knots)

      if (GCV < err) {
        Best.SetBF <- New.SetBF[colsOrder]
            Best.B <- New.B
               err <- GCV
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
      GCV   = err,
      t     = Best.knots,
      coefs = Best.coefs
      )

    # Best model of p terms
    models <- append(models, list(P_Model))

    # Update number of terms
    terms <- terms - 1
  }

  return(models)
}

#' @title Estimate Coefficients in Additive Adaptive Frontier Splines during Backward Procedure.
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

  coefs <- matrix(NA, nrow = p, ncol = ncol(y))

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
    dirs <- c(rep("==", n), rep(">=", nrow(Amat) - n))

    # Bounds
    bnds <- list(lower = list(ind = 1L:p, val = rep(-Inf, p)))

    # Solve
    sols <- Rglpk_solve_LP(objVal, Amat, dirs, bvec, bnds)

    coefs[, out] <- sols$solution[1:p]
  }

  return(coefs)
}

#' @title Compute Generalized Cross-Validation
#'
#' @description This function computes the generalized cross-validation for the backward procedure.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param y Column output index in \code{data}.
#' @param B \code{matrix} of basis functions.
#' @param coefs \code{vector} of coefficients in the model.
#' @param d Generalized Cross Validation (GCV) penalty per knot.
#' @param knots \code{list} of selected knots.
#'
#' @return Generalized Cross-Validation value
computeGCV <- function(data, y, B, coefs, d, knots) {

  # Predictions
  y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

  for (out in 1:length(y)) {
    y_hat[, out] <- B %*% coefs[, out, drop = F]
  }

  # Mean error
  LOF.numer <- mean(unlist(me(data[, y, drop = F], y_hat[drop = F])))

  # Cost-complexity cost
  len.knots <- length(unique(lapply(knots, function(x) list(xi = x[["xi"]], t = x[["t"]]))))
  CM <- ncol(B) + d * len.knots
  LOF.denom <- (1 - (CM / nrow(data))) ^ 2

  GCV <- LOF.numer / LOF.denom

  return(GCV)
}

