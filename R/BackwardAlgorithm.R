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
       id = NA,
        B = NA,
      LOF = NA,
        t = NA,
    alpha = NA))

  # ==
  # id
  # ==
  models[[1]][["id"]] <- length(ForwardModel[["BF"]])

  # ==
  # B matrix
  # ==
  models[[1]][["B"]] <- ForwardModel[["B"]]

  # ==
  # knots
  # ==
  knots <- lapply(ForwardModel[["BF"]][- 1], function(x) list(xi = x[["xi"]],
                                                               t = x[["t"]],
                                                            side = x[["side"]]))
  models[[1]][["t"]] <- knots

  # ==
  # LOF
  # ==
  # Predictions
  y_hat <- ForwardModel[["B"]] %*% ForwardModel[["alpha"]]

  # mse
  num <- mse(data[, y], y_hat)

  if (d != - 1) {

    lgt_knots <- length(unique(lapply(knots, function(x) list(xi = x[["xi"]],
                                                              t = x[["t"]]))))

    CM <- ncol(ForwardModel[["B"]]) + d * lgt_knots
    den <- (1 - (CM / nrow(data))) ^ 2

  } else {
    den <- 1
  }

  LOF <- num / den

  models[[1]][["LOF"]] <- LOF

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
                                                (x["status"] == "not paired" && x[["side"]] == "L"))
                                                 x[['id']]))

    for (bf in NBFs) {
      # Update SetBF
      New_SetBF <- SetBF

      # Drop element
      dropTerm <- sapply(New_SetBF, function(x) x[["id"]] == bf)

      # Update B
      # Ensure matrix format when B is just a vector
      New_B <- B[, !dropTerm]

      # After Forward Algorithm: paired | R --> Always even number
      if (!bf %% 2){
        New_SetBF[lapply(New_SetBF, "[[", "id") == bf + 1][[1]][["status"]] <- "not paired"
      }

      # Drop Basis Function
      New_SetBF[dropTerm] <- NULL

      # Count pair of basis functions
      h <- length(unlist(lapply(New_SetBF, function(x) if((x["status"] == "paired"))
                                                           x[['id']])))

      # Count left-side basis functions
      r <- length(unlist(lapply(New_SetBF, function(x) if((x["status"] == "not paired"))
                                                           x[['id']])))
      # Estimate coefficients
      # Predict y_hat
      if (terms > 1){

        not_paired <- (1:terms)[sapply(New_SetBF, function(x) x[["status"]] == "not paired")]
            paired <- (1:terms)[sapply(New_SetBF, function(x) x[["status"]] != "not paired")]
         colsOrder <- c(paired, not_paired)

             coefs <- EstimCoeffsBackward(New_B[, colsOrder], data[, y], h, r)
             y_hat <- New_B[, colsOrder] %*% coefs

      } else {

        colsOrder <- 1

        # Last iteration. Need to transform New_B from vector to matrix
        coefs <- max(data[, y])
        New_B <- matrix(New_B)
        y_hat <- New_B %*% coefs
      }

      # mse
      num <- mse(data[, y], y_hat)

      # Knots
      knots <- lapply(New_SetBF[colsOrder][- 1],
                      function(x) if(all(x[["xi"]] != - 1)) list(xi = x[['xi']],
                                                                  t = x[['t']],
                                                               side = x[["side"]]))
      # LOF
      if (d != - 1) {

        lgt_knots <- length(unique(lapply(knots, function(x) list(xi = x[["xi"]],
                                                                   t = x[["t"]]))))

               CM <- ncol(New_B) + d * lgt_knots
              den <- (1 - (CM / nrow(data))) ^ 2

      } else {
        den <- 1
      }

      LOF <- num / den

      if (LOF < err) {
        Best_SetBF <- New_SetBF[colsOrder]
            Best_B <- New_B[, colsOrder]
               err <- LOF
        Best_knots <- knots
        Best_coefs <- coefs
      }
    }

    # Update set of basis functions and matrix of basis functions
    SetBF <- Best_SetBF
        B <- Best_B

    # Create P-model
    P_Model <- list(
              id = ncol(Best_B),
               B = Best_B,
             LOF = err,
               t = Best_knots,
           alpha = Best_coefs
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
#' @description This function solves a Quadratic Programming Problem to obtain a set of coefficients.
#'
#' @param B \code{matrix} of basis functions.
#' @param y Output \code{vector} in data.
#' @param h Number of coefficients associated with pairs of basis functions.
#' @param r Number of coefficients associated with left-side basis functions.
#'
#' @importFrom quadprog solve.QP
#'
#' @return \code{vector} with the coefficients estimated.
EstimCoeffsBackward <- function(B, y, h, r){

  n <- nrow(B)
  p <- ncol(B) # (p = 1 + h + r)

  # vars: c(tau_0, alpha_1, beta_2, ... , alpha_(H-1), beta_H, w_1, ..., w_R, e_1, ... , e_n)

  # ================= #
  # D: Quadratic part #
  # ================= #
  # Identity matrix
  Dmat <- diag(rep(1, p + n))
  # Near 0 for alpha
  Dmat[1:p, 1:p] <- diag(rep(1e-12, p), p)

  # ============== #
  # d: Linear part #
  # ============== #
  dvec = rep(0, p + n)

  # = #
  # A #
  # = #
  # Equality constraints: y_hat - e = y
  Amat1 <- cbind(B, diag(rep(-1, n), n))

  # e_n >= 0
  Amat2 <- cbind(matrix(0, n, p), diag(rep(1, n)))

  # Concavity [paired basis functions]
  Amat3 <- matrix(0, nrow = h / 2, ncol = p + n)

  if (h > 0) {
    for (i in seq(2, h, 2)){
      Amat3[i / 2, c(i, i + 1)] <- - 1
    }
  }

  # Increasing monotony [paired basis functions]
  Amat4 <- cbind(rep(0, h), diag(x = c(1, -1), h), matrix(0, h, r), matrix(0, h, n))

  # Concavity & increasing monotony [left-side basis functions]
  Amat5 <- cbind(matrix(0, r, 1 + h), diag(- 1, r), matrix(0, r, n))

  Amat <- rbind(Amat1, Amat2, Amat3, Amat4, Amat5)

  # = #
  # b #
  # = #
  bvec <- c(y,  rep(0, n), rep(0, h / 2), rep(0, r), rep(0, h))

  s <- solve.QP(Dmat = Dmat, dvec = dvec,
                Amat = t(Amat), bvec = bvec,
                meq = n)

  alpha <- s$solution[1:p]

  return(alpha)
}

