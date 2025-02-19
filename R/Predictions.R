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
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_cubic"}}: Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_quintic"}}: Quintic Smoothed Adaptive Constrained Enveloping Splines model.
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
  y_hat <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

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

  # sample size
  N <- nrow(data)

  if (!method %in% c("aces_forward", "aces", "aces_cubic", "aces_quintic")) {
    stop("Not available method. Please, check help(\"predict\")")
  }

  # model
  aces_model <- object[["methods"]][[method]]

  # set of knots
  knots <- aces_model[["knots"]]

  # technology
  tecno <- object[["technology"]][[method]]

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

  # compute DEA scores
  scores <- rad_out (
    tech_xmat = tecno[["xmat"]],
    tech_ymat = tecno[["ymat"]],
    eval_xmat = as.matrix(newdata[, x]),
    eval_ymat = as.matrix(y_hat),
    convexity = TRUE,
    returns = "variable"
  )

  # predictions
  y_hat <- as.data.frame(y_hat * scores)

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
#' \item{\code{"rf_aces"}}: Random Forest Adaptive Constrained Enveloping Splines model.
#' \item{\code{"rf_aces_cubic"}}: Random Forest Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"rf_aces_quintic"}}: Random Forest Quintic Smoothed Adaptive Constrained Enveloping Splines model.
#' }
#'
#' @param stochastic_pred
#' This parameter should be kept as \code{FALSE} (default) if an additive model is estimated. If a stochastic model is estimated, one of the following elements should be specified:
#' \itemize{
#'   \item{"avg"}: average-practice production function.
#'   \item{"mom"}: moments-estimation production function.
#'   \item{"pse"}: pseudolikelihood-estimation production function.
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
    method = "rf_aces",
    stochastic_pred = FALSE
    ) {

  # number of outputs
  nY <- length(object[["models"]][[1]][["data"]][["y"]])

  # determine model type
  model_type <- object[["models"]][[1]][["control"]][["model_type"]]

  if (model_type == "stochastic" && !stochastic_pred %in% c("avg", "mom", "pse")) {
    stop("A stochastic prediction must be selected if a stochastic model is estimated.")
  }

  # determine error type
  error_type <- object[["models"]][[1]][["control"]][["error_type"]]

  # number of models in Random Forest
  RF_models <- length(object[["models"]])

  # list of predictions for each model
  y_hat_RF <- vector("list", RF_models)

  for (b in 1:RF_models) {

    # output predictions
    y_hat <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

    # select a model
    model <- object[["models"]][[b]]

    # check if training and test names are equal
    tr_names <- model[["data"]][["xnames"]]
    ts_names <- colnames(newdata)[c(x)]

    if (!identical(sort(tr_names), sort(ts_names))) {
      stop("Different variable names in training data and newdata.")
    }

    # data in [x, y] format with interaction of variables included
    data <- prepare_data (
      data = newdata,
      x = x,
      y = NULL,
      degree = model[["control"]][["degree"]],
      error_type = error_type
    )

    # sample size
    N <- nrow(data)

    if (method == "rf_aces") {
      aces_model <- model[["methods"]][["rf_aces"]]
      knots <- aces_model[["knots"]]

    } else if (method == "rf_aces_cubic") {
      aces_model <- model[["methods"]][["rf_aces_cubic"]]
      knots <- aces_model[["cubic_knots"]]

    } else if (method == "rf_aces_quintic") {
      aces_model <- model[["methods"]][["rf_aces_quintic"]]
      knots <- aces_model[["quintic_knots"]]

    } else {
      stop("Not available method Please, check help(\"predict\")")

    }

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

      # compute DEA scores
      scores <- rad_out (
        tech_xmat = aces_model[["technology"]][["xmat"]],
        tech_ymat = aces_model[["technology"]][["ymat"]],
        eval_xmat = as.matrix(newdata[, x]),
        eval_ymat = as.matrix(y_hat),
        convexity = TRUE,
        returns = "variable"
      )

      # predictions
      y_hat <- as.data.frame(y_hat * scores)

    } else {

      for (out in 1:nY) {
        y_hat[, out] <- B %*% aces_model[["coefs"]][, out, drop = F]
        y_hat[, out] <- exp(y_hat[, out])
      }

    }

    if (model_type == "stochastic") {

      # standard deviation for inefficiency term
      std_u <- switch (
        stochastic_pred,
        "avg" = 0,
        "mom" = aces_model[["stochastic"]][["mom"]][["std_u"]],
        "pse" = aces_model[["stochastic"]][["pse"]][["std_u"]]
      )

      adjustment_factor <- std_u * sqrt(2 / pi)

      # production function
      if (error_type == "add") {

        y_hat <- y_hat + adjustment_factor

      } else {

        y_hat <- y_hat * exp(adjustment_factor)

      }
    }
  }


  #     # prediction
  #     if (y_type == "all") {
  #
  #       for (out in 1:nY) {
  #         y_hat[, out] <- pmax(0, B %*% aces_model[["coefs"]][, out, drop = F])
  #       }
  #
  #       names(y_hat) <- paste(model[["data"]][["ynames"]], "_pred", sep = "")
  #
  #     } else {
  #
  #       y_hat[, m] <- pmax(0, B %*% aces_model[["coefs"]])
  #
  #       names(y_hat)[m] <- paste(model[["data"]][["ynames"]], "_pred", sep = "")
  #
  #     }
  #   }
  #
  #   if (error_type == "mul") {
  #     # change to original scale
  #     y_hat <- exp(y_hat)
  #   }
  #
  #   y_hat <- rad_out (
  #     tech_xmat = as.matrix(newdata[, x]),
  #     tech_ymat = as.matrix(y_hat),
  #     eval_xmat = as.matrix(newdata[, x]),
  #     eval_ymat = as.matrix(y_hat),
  #     convexity = TRUE,
  #     returns = "variable"
  #   ) * y_hat
  #
  #   y_hat_RF[[b]] <- y_hat
  #
  # }

  # point estimation
  y_hat_aux <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

  # interval estimation
  inter_estim <- vector("list", nY)

  # mean prediction
  for (var in 1:nY) {

    # select "var" variable for each data.frame
    rf_estimation <- lapply(y_hat_RF, function(df) df[, var])

    # add to interval estimation
    inter_estim[[nY]] <- matrix(do.call(cbind, rf_estimation), ncol = 1)

    # transform to matrix
    matrix_var <- do.call(cbind, rf_estimation)

    # mean predictions
    y_hat_aux[, var] <- rowMeans(matrix_var, na.rm = TRUE)

  }

  names(y_hat_aux) <- names(y_hat)
  names(inter_estim) <- names(y_hat)

  y_hat <- y_hat_aux

  return(y_hat)

}

#' @title Model Prediction for Adaptive Constrained Enveloping Splines (ACES).
#'
#' @description
#'
#' This function predicts the expected output by an \code{aces} object.
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
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_cubic"}}: Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_quintic"}}: Quintic Smoothed Adaptive Constrained Enveloping Splines model.
#' }
#'
#' @param stochastic_pred
#' This parameter should be kept as \code{FALSE} (default) if an additive model is estimated. If a stochastic model is estimated, one of the following elements should be specified:
#' \itemize{
#'   \item{"avg"}: average-practice production function.
#'   \item{"mom"}: moments-estimation production function.
#'   \item{"pse"}: pseudolikelihood-estimation production function.
#' }
#'
#' @return
#'
#' A \code{data.frame} with the predicted values through the Adaptive Constrained Enveloping Splines model.
#'
#' @export

predict.s_aces <- function (
    object,
    newdata,
    x,
    method = "aces",
    stochastic_pred = FALSE
) {

  # number of outputs
  nY <- length(object[["data"]][["y"]])

  # determine model type
  model_type <- object[["control"]][["model_type"]]

  if (model_type == "stochastic" && !stochastic_pred %in% c("avg", "mom", "pse")) {
    stop("A stochastic prediction must be selected if a stochastic model is estimated.")
  }

  # determine error type
  error_type <- object[["control"]][["error_type"]]

  # output predictions
  y_hat <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

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

  # sample size
  N <- nrow(data)

  if (method == "aces_forward") {
    aces_model <- object[["methods"]][["aces_forward"]]
    knots <- aces_model[["knots"]]
    tecno <- object[["technology"]][["aces_forward"]]

  } else if (method == "aces") {
    aces_model <- object[["methods"]][["aces"]]
    knots <- aces_model[["knots"]]
    tecno <- object[["technology"]][["aces"]]

  } else if (method == "aces_cubic") {
    aces_model <- object[["methods"]][["aces_cubic"]]
    knots <- aces_model[["cubic_knots"]]
    tecno <- object[["technology"]][["aces_cubic"]]

  } else if (method == "aces_quintic") {
    aces_model <- object[["methods"]][["aces_quintic"]]
    knots <- aces_model[["quintic_knots"]]
    tecno <- object[["technology"]][["aces_quintic"]]

  } else {
    stop("Not available method. Please, check help(\"predict\")")

  }

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

  # compute DEA scores
  scores <- rad_out (
    tech_xmat = tecno[["xmat"]],
    tech_ymat = tecno[["ymat"]],
    eval_xmat = as.matrix(newdata[, x]),
    eval_ymat = as.matrix(y_hat),
    convexity = TRUE,
    returns = "variable"
  )

  # predictions
  y_hat <- as.data.frame(y_hat * scores)

  if (model_type == "stochastic") {

    # standard deviation for inefficiency term
    std_u <- switch (
      stochastic_pred,
      "avg" = 0,
      "mom" = aces_model[["stochastic"]][["mom"]][["std_u"]],
      "pse" = aces_model[["stochastic"]][["pse"]][["std_u"]]
    )

    adjustment_factor <- std_u * sqrt(2 / pi)

    # production function
    if (error_type == "add") {

      y_hat <- y_hat + adjustment_factor

    } else {

      y_hat <- y_hat * exp(adjustment_factor)

    }
  }

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
