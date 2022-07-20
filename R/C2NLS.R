#' @title Corrected Concave Nonparametric Least Squares
#'
#' @description This function performs Corrected Concave Nonparametric Least Squares.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param X input indexes in \code{data}.
#' @param y output indexes in \code{data}.
#'
#' @return A \code{C2NLS} object.
#'
#' @export
C2NLS <- function(data,x,y) {

   n <- nrow(data)
  nX <- length(x)

  # Corrected Nonparametric Least Squares
  df_CNLS <- CNLS(data,x,y)

  epsil <- df_CNLS$epsilon - max(df_CNLS$epsilon)
  alpha <- df_CNLS$alpha + max(df_CNLS$epsilon)

  beta  <- matrix(nrow = n, ncol = nX)

  for (i in 1:n) {
    for (j in 1:nX) {
      beta[i, j] = df_CNLS[i, j + 2] # f, alpha, beta, nepsil
    }
  }

  f <- alpha + rowSums(beta * data[, x])

  C2NLS <- C2NLS_object(data, x, y, alpha, beta, epsil, f)

  return(C2NLS)
}

#' @title Corrected Nonparametric Least Squares
#'
#' @description This function performs Corrected Nonparametric Least Squares.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param X input indexes in \code{data}.
#' @param y output indexes in \code{data}.
#'
#' @importFrom dplyr mutate_if
#' @importFrom quadprog solve.QP
#'
#' @return \code{data.frame} with the solutions of the optimization problem solved.
CNLS <- function(data, x, y) {

   n <- dim(data)[1]
  nX <- length(x)

  # num_vars = n + n * nX + n
  # vars: c(alpha_0, ..., alpha_N, beta_1, ..., beta_n, e_1, ... , e_n)
  colnames <- c(paste("alpha_", 1:n, sep = ""),
                paste(paste("beta_", sort(rep(1:n, nX)), sep = ""), rep(1:nX, n), sep = "_"),
                paste("ep_", 1:n, sep = ""))

  # ================= #
  # D: Quadratic part #
  # ================= #

  # Identity matrix (1 for epsilons)
  Dmat <- diag(rep(1e-12, 2 * n + n * nX))

  # Near 0 for alpha and beta
  Dmat[(n + n * nX + 1):(2 * n + n * nX), (n + n * nX + 1):(2 * n + n * nX)] <- diag(rep(1, n), n)
  colnames(Dmat) <- colnames

  # ============== #
  # d: Linear part #
  # ============== #
  dvec <- rep(0, 2 * n + n * nX)

  # ================================== #
  # A: matrix defining the constraints #
  # ================================== #

  # R1 --> alpha_i + beta_i * x_i + ep_i = y_i
  # n restricciones de tipo R1
  Amat1 <- cbind(diag(rep(1, n), n), matrix(rep(diag(0, n), nX), n), diag(rep(1, n), n))
  colnames(Amat1) <- colnames

  cont <- n + 1
  for (i in 1:n) {
    for (j in 1:nX) {
      Amat1[i,cont] <- data[i, j]
      cont <- cont + 1
    }
  }

  # R2 --> alpha_i - alpha_h + beta_i * x_i - beta_h * x_i + ep_i - ep_h = 0
  # n * n restricciones tipo R2
  Amat2 <- c()
  for (i in 1:n) {
    Amat2_N <- cbind(diag(rep(1, n), n), matrix(rep(diag(0, n), nX), n), diag(rep(0, n), n))
    Amat2_N[, i] <- Amat2_N[, i] - 1
    cont <- n + 1
    for (k in 1:n) {
      for (j in 1:nX) {
        Amat2_N[k,cont] <- data[i, j]
        cont <- cont + 1
      }
    }
    for (k in 1:n) {
      for (j in 1:nX) {
        Amat2_N[k, n + j + (i - 1) * nX] <- Amat2_N[k, n + j + (i - 1) * nX] - data[i, j]
      }
    }
    Amat2 <- rbind(Amat2, Amat2_N)
  }

  # R3 --> b_i >= 0
  # n restricciones de tipo R3
  Amat3 <- cbind(matrix(data = 0, ncol = n, nrow = n * nX), diag(rep(1, n * nX), n * nX),
                 matrix(data = 0, ncol = n, nrow = n * nX))

  # n + n * n + n * nX
  Amat <- rbind(Amat1, Amat2, Amat3)
  colnames(Amat) <- colnames

  # ================================ #
  # b: right size of the constraints #
  # ================================ #
  bvec <- c(data[, y], rep(0, n * n), rep(0, n * nX))

  # ============================ #
  # Solve the quadratic problems #
  # ============================ #
  s <- solve.QP(Dmat = Dmat, dvec = dvec,
                Amat = t(Amat), bvec = bvec,
                meq = n)

  # Solutions
  alpha <- s$solution[1:n]
  beta  <- matrix(nrow = n, ncol = nX)
  cont  <- n + 1

  for (i in 1:n) {
    for (j in 1:nX) {
      beta[i, j] <- s$solution[cont]
      cont <- cont + 1
    }
  }

  ep <- tail(s$solution, n)

  df <- data.frame(alpha = alpha, beta = beta, epsilon = ep) %>%
    mutate_if(is.numeric, round, digits = 10)

  f  <- alpha + rowSums(beta * data[, x])
  df <- cbind(f, df)

  return(df)
}

#' @title Create a C2NLS object
#'
#' @description This function saves information about the Corrected Concave Nonparametric Least Squares model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param alpha alpha variable.
#' @param beta beta variable.
#' @param epsilon epsilon variable.
#' @param f alpha + beta * x + e.
#'
#' @return A \code{C2NLS} object.
C2NLS_object <- function(data, x, y, alpha, beta, epsilon, f) {

  C2NLS_object <- list("data" = list(data = data,
                                     x = x,
                                     y = y,
                                     input_names = names(data)[x],
                                     output_names = names(data)[y]),
                       "alpha" = alpha,
                       "beta"  = beta,
                       "epsilon" = epsilon,
                       "f" = f)

  class(C2NLS_object) <- "C2NLS"

  return(C2NLS_object)
}

#' @title Model Prediction for Corrected Concave Nonparametric Least Squares
#'
#' @description This function predicts the expected output by a \code{C2NLS} object.
#'
#' @param object A \code{C2NLS} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#'
#' @return \code{data.frame} with the predicted values.
#'
#' @export
predict.C2NLS <- function(object, newdata, x) {

  # Check that both data samples (train and predict) have the same number of inputs
  train_x <- object[["data"]][["x"]]
  if (length(x) != length(train_x)) {
    stop("Different variable names in training and test sets.")
  }

  n_pred <- nrow(newdata)

  # Get alpha and betas from training
  n_train <- nrow(object[["data"]][["data"]])
  alpha   <- object[["alpha"]]
  beta    <- object[["beta"]]
  epsilon <- object[["epsilon"]]

  # Predict: f(x) = min {alpha_i + beta_i * x}
  f <- matrix(nrow = n_pred, ncol = 1)

  for (k in 1:n_pred) {
    min <- Inf
    for (i in 1:n_train) {
      if (length(train_x) == 1) { # mono-input
        comb_actual = alpha[k] + beta[k] * newdata[k, x]
      } else { # multi-input
        comb_actual = alpha[k] + rowSums(beta[k, ] * newdata[k, x])
      }
      if (comb_actual < min) {
        min <- comb_actual
      }
    }
    f[k] = min
  }

  return(data.frame(x = newdata[, x], f = f))
}
