#' @title Model Prediction for Adaptive Constrained Enveloping Splines (ACES).
#'
#' @description This function predicts the expected output by an \code{aces} object.
#'
#' @param object An \code{aces} object.
#' @param newdata A \code{data.frame} containing the input and netput variables to predict on.
#' @param x Input indexes in \code{newdata}.
#' @param method Model for prediction.
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_cubic"}}: Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_quintic"}}: Quintic Smoothed Adaptive Constrained Enveloping Splines model.
#' }
#'
#' @return A \code{data.frame} with the predicted values through the Adaptive Constrained Enveloping Splines model.
#'
#' @export

predict.aces <- function (
    object, newdata, x, method = "aces"
    ) {

  # strategy to predict y: "all" or "individual"
  if (length(names(object)) == 1 && names(object) == "y_all") {
    y_type <- "all"
  } else {
    y_type <- "individual"
  }

  if (y_type == "all") {
    nY <- length(object[["y_all"]][["data"]][["y"]])
  } else {
    nY <- length(object)
  }

  # number of models;
  # all: 1 model
  # individual: 1 model for each output
  models <- length(object)

  # output predictions
  y_hat <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

  for (m in 1:models) {
    # select a model
    model <- object[[m]]

    # look for netputs if y_type = "individual"
    z <- model[["data"]][["z"]]

    # check if training and test names are equal
    tr_names <- model[["data"]][["xnames"]]
    ts_names <- colnames(newdata)[c(x, z)]

    if (!identical(sort(tr_names), sort(ts_names))) {
      stop("Different variable names in training data and newdata.")
    }

    # data in [x, z, y] format with interaction of variables included
    data <- set_interactions (
      data = newdata,
      x = x,
      y = NULL,
      z = z,
      degree = model[["control"]][["degree"]]
    )

    # sample size
    N <- nrow(data)

    if (method == "aces_forward") {
      aces_model <- model[["methods"]][["aces_forward"]]
      knots <- aces_model[["knots"]]

    } else if (method == "aces") {
      aces_model <- model[["methods"]][["aces"]]
      knots <- aces_model[["knots"]]

    } else if (method == "aces_cubic") {
      aces_model <- model[["methods"]][["aces_cubic"]]
      knots <- aces_model[["cubic_knots"]]

    } else if (method == "aces_quintic") {
      aces_model <- model[["methods"]][["aces_quintic"]]
      knots <- aces_model[["quintic_knots"]]

    } else {
      stop("Not available method Please, check help(\"predict\")")
    }

    # matrix of basis function
    B <- set_B (
      newdata = data,
      model = aces_model,
      knots = knots,
      method = method
      )

    # prediction
    if (y_type == "all") {
      for (out in 1:nY) {
        y_hat[, out] <- pmax(0, B %*% aces_model[["coefs"]][, out, drop = F])
      }
      names(y_hat) <- paste(model[["data"]][["ynames"]], "_pred", sep = "")
    } else {
      y_hat[, m] <- pmax(0, B %*% aces_model[["coefs"]])
      names(y_hat)[m] <- paste(model[["data"]][["ynames"]], "_pred", sep = "")
    }
  }

  return(y_hat)
}

#' @title Build (B) Matrix of Basis Functions
#'
#' @description This function builds the (B) matrix of basis functions for prediction given a model.
#'
#' @param newdata A \code{data.frame} containing the input and netput variables to predict on.
#' @param model A \code{list} with data of a model.
#' @param knots Set of knots of the \code{model}.
#' @param method Type of model for prediction.
#'
#' @return The (B) Matrix of Basis Functions

set_B <- function (
    newdata, model, knots, method
    ) {

  # Sample size
  N <- nrow(newdata)
  # Initialize B matrix
  B <- matrix(rep(1, N), nrow = N)

  if (method == "aces_forward") {

    for (i in 1:nrow(knots)) {
      # Variable and knot
      v <- knots[i, 1]
      t <- knots[i, 2]

      # Basis functions for the new data
      hinge1 <- pmax(0, newdata[, v] - t)
      hinge2 <- pmax(0, t - newdata[, v])

      # Update B matrix
      B <- cbind(B, hinge1, hinge2)
    }

  } else if (method == "aces") {

    for (i in 1:nrow(knots)) {
      # Variable and knot
      v <- knots[i, 1]
      t <- knots[i, 2]

      # Basis functions for the new data
      if (knots[i, 3] == "R") {
        hinge1 <- pmax(0, newdata[, v] - t)

        # Update B
        B <- cbind(B, hinge1)

      } else {
        hinge2 <- pmax(0, t - newdata[, v])

        # Update B
        B <- cbind(B, hinge2)
      }
    }

  } else if (method == "aces_cubic") {

    for (status in c("paired", "unpaired")) {
      for (v in 1:length(knots)) { # variable
        if (is.null(knots[[v]])) next # variable is not used
        for (t in 1:length(knots[[v]])) {
          if (knots[[v]][[t]][["status"]] != status) next # keep the order

          # knot side
          side <- knots[[v]][[t]][["side"]]

          # cubic knot
          Ct <- knots[[v]][[t]][["t"]]

          # t-; t; t+
          t0 <- Ct[1]; t1 <- Ct[2]; t2 <- Ct[3]

          d <- t1 - t0 # t  - t-
          e <- t2 - t1 # t+ - t

          p1 <- (2 * e - d) / (e + d) ^ 2
          r1 <- (d - e) / (e + d) ^ 3

          p2 <- (2 * d - e) / (- e - d) ^ 2
          r2 <- (e - d) / (- e - d) ^ 3

          if (status == "paired" || (status == "unpaired" && side == "R")) {

            term1 <- p1 * (newdata[, v] - t0) ^ 2
            term2 <- r1 * (newdata[, v] - t0) ^ 3

            # C1
            C1 <- ifelse(newdata[, v] <= t0,
                         0,
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2,
                                 newdata[, v] - t1)))
            B  <- cbind(B, C1)
          }

          if (status == "paired" || (status == "unpaired" && side == "L")) {

            term1 <- p2 * (newdata[, v] - t2) ^ 2
            term2 <- r2 * (newdata[, v] - t2) ^ 3

            # C2
            C2 <- ifelse(newdata[, v] <= t0,
                         t1 - newdata[, v],
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2,
                                 0)))
            B  <- cbind(B, C2)
          }
        }
      }
    }

  } else if (method == "aces_quintic") {

    for (status in c("paired", "unpaired")) {
      for (v in 1:length(knots)) { # variable
        if (is.null(knots[[v]])) next # variable is not used
        for (t in 1:length(knots[[v]])) {
          if (knots[[v]][[t]][["status"]] != status) next # keep the order

          # knot side
          side <- knots[[v]][[t]][["side"]]

          # quintic knot
          Qt <- knots[[v]][[t]][["t"]]

          # t-; t; t+
          t0 <- Qt[1]; t1 <- Qt[2]; t2 <- Qt[3]

          d  <- t2 - t0 # t+ - t-
          d1 <- t2 - t1 # t+ - t
          d2 <- t1 - t0 # t  - t-

          alpha1 <- (6 * d1 - 4 * d2) / d ^ 3
          alpha2 <- (4 * d1 - 6 * d2) / d ^ 3

          beta1 <- (7 * d2 - 8 * d1) / d ^ 4
          beta2 <- (7 * d1 - 8 * d2) / d ^ 4

          gamma1 <- gamma2 <- (3 * d1 - 3 * d2) / d ^ 5

          if (status == "paired" || (status == "unpaired" && side == "R")) {

            term1 <- alpha1 * (newdata[, v] - t0) ^ 3
            term2 <- beta1  * (newdata[, v] - t0) ^ 4
            term3 <- gamma1 * (newdata[, v] - t0) ^ 5

            # Q1
            Q1 <- ifelse(newdata[, v] <= t0,
                         0,
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2 + term3,
                                 newdata[, v] - t1)))

            B  <- cbind(B, Q1)
          }

          if (status == "paired" || (status == "unpaired" && side == "L")) {

            term1 <- alpha2 * (newdata[, v] - t2) ^ 3
            term2 <- beta2  * (newdata[, v] - t2) ^ 4
            term3 <- gamma2 * (newdata[, v] - t2) ^ 5

            # Q2
            Q2 <- ifelse(newdata[, v] <= t0,
                         t1 - newdata[, v],
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2 + term3,
                                 0)))
            B  <- cbind(B, Q2)
          }
        }
      }
    }

  } else {
    stop("Not available method Please, check help(\"predict\")")
  }

  return(B)
}
