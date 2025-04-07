#' @title Model Prediction for Adaptive Constrained Enveloping Splines (ACES).
#'
#' @description
#'
#' This function predicts the expected output using an \code{aces} object.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param newdata
#' A \code{data.frame} containing the input and netput variables to predict on.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param method
#' Model for prediction:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_cubic"}}: Cubic Smoothed Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_quintic"}}: Quintic Smoothed Adaptive Constrained Enveloping Splines.
#' }
#'
#' @return
#'
#' A \code{data.frame} with the predicted values through the Adaptive Constrained Enveloping Splines model.
#'
#' @export

predict.aces <- function (
    object,
    newdata,
    x,
    method = "aces"
    ) {

  # number of outputs
  nY <- length(object[["data"]][["y"]])

  # determine error type
  error_type <- object[["control"]][["error_type"]]

  # output predictions
  y_hat_aux <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

  # check if training and test names are equal
  tr_names <- object[["data"]][["xnames"]]
  ts_names <- colnames(newdata)[x]

  if (!identical(sort(tr_names), sort(ts_names))) {
    stop("Different variable names in training data and newdata.")
  }

  # data in [x, y] format with interaction and / or transformation of variables included
  data <- prepare_data (
    data = newdata,
    x = x,
    y = NULL,
    max_degree = object[["control"]][["max_degree"]],
    error_type = error_type
  )

  # technology
  tecno <- object[["technology"]][[method]]

  if (!method %in% c("aces_forward", "aces", "aces_cubic", "aces_quintic")) {
    stop("Not available method. Please, check help(\"predict\")")
  }

  # model
  aces_model <- object[["methods"]][[method]]

  # set of knots
  knots <- aces_model[["knots"]]

  # matrix of basis function
  B <- set_Bmat (
    newdata = data,
    model = aces_model,
    knots = knots,
    method = method
    )

  if (error_type == "add") {

    for (out in 1:nY) {
      y_hat_aux[, out] <- pmax(0, B %*% aces_model[["coefs"]][, out, drop = F])
    }

  } else {

    for (out in 1:nY) {
      y_hat_aux[, out] <- B %*% aces_model[["coefs"]][, out, drop = F]
      y_hat_aux[, out] <- exp(y_hat_aux[, out])
    }

  }

  # compute DEA scores
  scores <- rad_out (
    tech_xmat = tecno[["xmat"]],
    tech_ymat = tecno[["ymat"]],
    eval_xmat = as.matrix(newdata[, x]),
    eval_ymat = as.matrix(y_hat_aux),
    convexity = TRUE,
    returns = "variable"
  )

  # remove unfeasibilities
  scores[scores < -10000 | scores > 10000] <- NA

  # predictions
  y_hat <- as.data.frame(y_hat_aux * scores)

  names(y_hat) <- paste(object[["data"]][["ynames"]], "_pred", sep = "")

  return(y_hat)

}

#' @title Model Prediction for Random Forest Adaptive Constrained Enveloping Splines (RF-ACES)
#'
#' @description
#'
#' This function predicts the expected output by an \code{rf_aces} object.
#'
#' @param object
#' A \code{rf_aces} object.
#'
#' @param newdata
#' A \code{data.frame} containing the input and netput variables to predict on.
#'
#' @param x
#' Input indexes in \code{newdata}.
#'
#' @param method
#' Model for prediction:
#' \itemize{
#' \item{\code{"rf_aces"}}: Random Forest Adaptive Constrained Enveloping Splines.
#' \item{\code{"rf_aces_cubic"}}: Random Forest Cubic Smoothed Adaptive Constrained Enveloping Splines.
#' \item{\code{"rf_aces_quintic"}}: Random Forest Quintic Smoothed Adaptive Constrained Enveloping Splines.
#' }
#'
#' @return
#'
#' A \code{data.frame} with the predicted values through the Random Forest Adaptive Constrained Enveloping Splines model.
#'
#' @export

predict.rf_aces <- function (
    object,
    newdata,
    x,
    method = "rf_aces"
    ) {

  # number of outputs
  nY <- length(object[["data"]][["y"]])

  # determine error type
  error_type <- object[["control"]][["error_type"]]

  # check if training and test names are equal
  tr_names <- object[["data"]][["xnames"]]
  ts_names <- colnames(newdata)[c(x)]

  if (!identical(sort(tr_names), sort(ts_names))) {
    stop("Different variable names in training data and newdata.")
  }

  # data in [x, y] format with interaction of variables included
  data <- prepare_data (
    data = newdata,
    x = x,
    y = NULL,
    max_degree = object[["control"]][["max_degree"]],
    error_type = error_type
  )

  # technology
  tecno <- object[["technology"]][[method]]

  # number of models in Random Forest
  RF_models <- length(object[["forest"]])

  # list of predictions for each model
  y_hat_RF <- vector("list", RF_models)

  for (t in 1:RF_models) {

    # output predictions
    y_hat <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

    # select an element from the forest
    model <- object[["forest"]][[t]]

    # method
    aces_model <- model[["methods"]][[method]]
    knots <- aces_model[["knots"]]

    # matrix of basis function
    B <- set_Bmat (
      newdata = data,
      model = aces_model,
      knots = knots,
      method = method
    )

    if (error_type == "add") {

      for (out in 1:nY) {
        y_hat[, out] <- pmax(0, B %*% aces_model[["coefs"]][, out, drop = F])
      }

    } else {

      for (out in 1:nY) {
        y_hat[, out] <- B %*% aces_model[["coefs"]][, out, drop = F]
        y_hat[, out] <- exp(y_hat[, out])
      }

    }

    y_hat_RF[[t]] <- as.data.frame(y_hat)

  }

  # point estimation
  y_hat_aux <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

  # mean prediction
  for (var in 1:nY) {

    # select "var" variable for each data.frame
    rf_estimation <- lapply(y_hat_RF, function(df) df[, var])

    # transform to matrix
    matrix_var <- do.call(cbind, rf_estimation)

    # mean predictions
    y_hat_aux[, var] <- rowMeans(matrix_var, na.rm = TRUE)

  }

  # compute DEA scores
  scores <- rad_out (
    tech_xmat = tecno[["xmat"]],
    tech_ymat = tecno[["ymat"]],
    eval_xmat = as.matrix(newdata[, x]),
    eval_ymat = as.matrix(y_hat_aux),
    convexity = TRUE,
    returns = "variable"
  )

  # remove unfeasibilities
  scores[scores < -10000 | scores > 10000] <- NA

  # predictions
  y_hat <- as.data.frame(y_hat_aux * scores)

  names(y_hat) <- paste(object[["data"]][["ynames"]], "_pred", sep = "")

  return(y_hat)

}

#' @title Build (B) Matrix of Basis Functions
#'
#' @description
#' This function builds the (B) matrix of basis functions for prediction given a model.
#'
#' @param newdata
#' A \code{data.frame} containing the input and netput variables to predict on.
#'
#' @param model
#' A \code{list} with data of a model.
#'
#' @param knots
#' Set of knots of the \code{model}.
#'
#' @param method
#' Type of model for prediction.
#'
#' @return
#' The (B) Matrix of Basis Functions

set_Bmat <- function (
    newdata,
    model,
    knots,
    method
    ) {

  # sample size
  N <- nrow(newdata)

  # initialize B matrix
  B <- matrix(rep(1, N), nrow = N)

  if (method %in% c("aces_forward", "rf_aces")) {

    for (i in 1:nrow(knots)) {

      # variable
      v <- knots[i, 1]

      # knot
      t <- knots[i, 2]

      # basis functions for the new data
      hinge1 <- pmax(0, newdata[, v] - t)
      hinge2 <- pmax(0, t - newdata[, v])

      # update B matrix
      B <- cbind(B, hinge1, hinge2)
    }

  } else if (method == "aces") {

    for (i in 1:nrow(knots)) {

      # variable
      v <- knots[i, 1]

      # knot
      t <- knots[i, 2]

      # basis functions for the new data
      if (knots[i, 3] == "R") {

        # right-hand basis function
        hinge1 <- pmax(0, newdata[, v] - t)

        # update B matrix
        B <- cbind(B, hinge1)

      } else {

        # left-hand basis function
        hinge2 <- pmax(0, t - newdata[, v])

        # update B matrix
        B <- cbind(B, hinge2)

      }
    }

  } else if (method %in% c("aces_cubic", "rf_aces_cubic")) {

    for (status in c("paired", "unpaired")) {

      # select a variable
      for (v in 1:length(knots)) {

        # the variable is not used
        if (is.null(knots[[v]])) next

        # select a knot
        for (t in 1:length(knots[[v]])) {

          # keep the order: first paired and then unpaired
          if (knots[[v]][[t]][["status"]] != status) next

          # knot side
          side <- knots[[v]][[t]][["side"]]

          # cubic knot
          Ct <- knots[[v]][[t]][["t"]]

          # t-
          t0 <- Ct[1]

          # t
          t1 <- Ct[2]

          # t+
          t2 <- Ct[3]

          # t  - t-
          d <- t1 - t0

          # t+ - t
          e <- t2 - t1

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

  } else if (method %in% c("aces_quintic", "rf_aces_quintic")) {

    for (status in c("paired", "unpaired")) {

      # select a variable
      for (v in 1:length(knots)) {

        # the variable is not used
        if (is.null(knots[[v]])) next

        # select a knot
        for (t in 1:length(knots[[v]])) {

          # keep the order: first paired and then unpaired
          if (knots[[v]][[t]][["status"]] != status) next

          # knot side
          side <- knots[[v]][[t]][["side"]]

          # quintic knot
          Qt <- knots[[v]][[t]][["t"]]

          # t-
          t0 <- Qt[1]

          # t
          t1 <- Qt[2]

          # t+
          t2 <- Qt[3]

          # t+ - t-
          d  <- t2 - t0

          # t+ - t
          d1 <- t2 - t1

          # t  - t-
          d2 <- t1 - t0

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

    stop("Not available method. Please, check help(\"predict\")")

  }

  return(B)

}
