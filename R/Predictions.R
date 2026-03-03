#' @title Compute Efficient Targets using Adaptive Constrained Enveloping Splines (ACES)
#'
#' @description
#' Projects a set of Decision Making Units (DMUs) onto the estimated production frontier. This function calculates the efficient input/output targets required for each DMU to achieve full efficiency, based on the specified ACES model and distance measure.
#'
#' @param eval_data
#' A \code{data.frame} or a \code{matrix} containing the DMUs to be evaluated.
#'
#' @param x
#' Column indexes of input variables in \code{eval_data}.
#'
#' @param y
#' Column indexes of output variables in \code{eval_data}.
#'
#' @param relevant
#' A \code{logical} indicating if only relevant variables should be included in the technology definition.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param method
#' Model prediction method used to compute predictions of inputs and obtain a new vector of outputs for \code{eval_data}:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_cubic"}}: Cubic Smooth Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_quintic"}}: Quintic Smooth Adaptive Constrained Enveloping Splines.
#' }
#'
#' @param measure
#' Mathematical programming model to calculate scores:
#' \itemize{
#' \item{\code{rad_out}} Output-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{rad_inp}} Input-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{ddf}}     Directional distance function proposed by \insertCite{chambers1998;textual}{aces}.
#' \item{\code{rsl_out}} Output-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{rsl_inp}} Input-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{wam_mip}} Measure of Inefficiency Proportions proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{wam_nor}} Normalized Weighted Additive Model proposed by \insertCite{lovell1995;textual}{aces}.
#' \item{\code{wam_ram}} Range Adjusted Measure proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{wam_bam}} Bounded Adjusted Measure proposed by \insertCite{cooper2011;textual}{aces}.
#' }
#'
#' @param returns
#' Type of returns to scale:
#' \itemize{
#' \item{\code{"constant"}} Constant returns to scale.
#' \item{\code{"variable"}} Variable returns to scale (default).
#' }
#'
#' @param direction
#' Direction of the vector to project on the frontier. Only applied if \code{measure = "ddf"}. A \code{matrix} or \code{data.frame} with \code{n} rows (number of DMUs to be evaluated) and \code{nX + nY} columns, containing the direction of the input variables followed by the direction of the output variables in the same order as they appear in the data.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{banker1984}{aces} \cr \cr
#' \insertRef{chambers1998}{aces} \cr \cr
#' \insertRef{fare1978}{aces} \cr \cr
#' \insertRef{cooper1999}{aces} \cr \cr
#' \insertRef{lovell1995}{aces} \cr \cr
#' \insertRef{cooper2011}{aces}
#'
#' @importFrom dplyr summarise %>% mutate_if
#' @importFrom stats median quantile sd
#'
#' @return
#' A \code{data.frame} with input and output targets computed through an Adaptive Constrained Enveloping Splines model.
#'
#' @export
#'
get_targets <- function (
    eval_data,
    x,
    y,
    relevant = FALSE,
    object,
    method = "aces",
    measure = "rad_out",
    returns = "variable",
    direction = NULL
    ) {

  # handle errors:
  display_errors_scores (
    data = eval_data,
    x = x,
    y = y,
    object = object,
    method = method,
    measure = measure,
    returns = returns,
    direction = direction
  )

  # number of inputs
  nX <- length(x)

  # number of outputs
  nY <- length(y)

  # =================== #
  # Data for technology #
  # =================== #

  # matrix of inputs
  tech_xmat <- as.matrix(object[["technology"]][[method]][["xmat"]])

  # auxiliar variable for using only relevant variables
  rel_x <- NULL

  if (relevant) {

    # set of knots
    knots <- object[["methods"]][[method]][["knots"]]

    # variable degree
    xi_degree <- object[["control"]][["xi_degree"]]
    colnames(xi_degree) <- names(object[["control"]][["kn_grid"]])

    # check participating variables
    participating_vars <- intersect(xi_degree[1, ], unique(knots$xi))

    # names of participant variables
    participating_vars_names <- colnames(xi_degree)[participating_vars]

    # extract individual variables assuming interaction variable names use "_"
    split_vars <- unique(unlist(strsplit(participating_vars_names, "_")))

    # get the column indices for the unique variables
    rel_x <- sort(match(split_vars, colnames(xi_degree)))

    # update technology
    tech_xmat <- as.matrix(tech_xmat[, rel_x, drop = FALSE])

  }

  # matrix of outputs
  tech_ymat <- as.matrix(object[["technology"]][[method]][["ymat"]])

  # ======================= #
  # Data for evaluated DMUs #
  # ======================= #

  # matrix of inputs
  if (relevant && !is.null(rel_x)) {
    eval_xmat <- as.matrix(eval_data[, x[rel_x]])
  } else {
    eval_xmat <- as.matrix(eval_data[, x])
  }

  # matrix of outputs
  eval_ymat <- as.matrix(eval_data[, y])

  # scaling setup
  scaling <- object[["control"]][["scale"]]

  if (!is.null(scaling) && scaling$is_scaled) {

    if (relevant) {
      sx <- scaling$mean_x[rel_x]
    } else {
      sx <- scaling$mean_x
    }

    sy <- scaling$mean_y

    # apply scaling
    eval_xmat <- sweep(eval_xmat, 2, sx, "/")
    eval_ymat <- sweep(eval_ymat, 2, sy, "/")

    if (!is.null(direction) && (is.matrix(direction) || is.data.frame(direction))) {

      dir_x <- as.matrix(direction[, 1:nX])
      dir_y <- as.matrix(direction[, (nX + 1):(nX + nY)])

      if (relevant) {
        if (ncol(dir_x) == length(scaling$mean_x)) {
          dir_x <- dir_x[, rel_x, drop = FALSE]
        }
      }

      # apply scaling
      dir_x_scaled <- sweep(dir_x, 2, sx, "/")
      dir_y_scaled <- sweep(dir_y, 2, sy, "/")

      # rebuild direction vectors
      direction <- cbind(dir_x_scaled, dir_y_scaled)

    }

  }

  # ========== #
  # Get Scores #
  # ========== #

  if (measure == "rad_out") {

    scores <- rad_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1]

  } else if (measure == "rad_inp") {

    scores <- rad_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1]

  } else if (measure == "ddf") {

    scores <- ddf (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      direction = direction,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1]

  } else if (measure == "rsl_out") {

    scores <- rsl_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1:nY]

  } else if (measure == "rsl_inp") {

    scores <- rsl_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1:nX]

  } else if (grepl("wam", measure)) {

    scores <- wam (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      weights = measure,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1:(nX + nY)]

  }

  # =========== #
  # Get Targets #
  # =========== #

  if (measure == "rad_out") {
    x_hat <- eval_xmat
    y_hat <- eval_ymat * scores

  } else if (measure == "rad_inp") {
    x_hat <- eval_xmat * scores
    y_hat <- eval_ymat

  } else if (measure == "ddf") {
    x_hat <- eval_xmat - scores * direction[, 1:nX]
    y_hat <- eval_ymat + scores * direction[, (nX + 1):(nX + nY)]

  } else if (measure == "rsl_out") {
    x_hat <- eval_xmat
    y_hat <- eval_ymat * scores

  } else if (measure == "rsl_inp") {
    x_hat <- eval_xmat * scores
    y_hat <- eval_ymat

  } else if (grepl("wam", measure)) {

    # extract scaled slacks directly from the linear programming solution
    slacks_x <- scores[, 1:nX, drop = FALSE]
    slacks_y <- scores[, (nX + 1):(nX + nY), drop = FALSE]

    x_hat <- eval_xmat - slacks_x
    y_hat <- eval_ymat + slacks_y

  }

  if (!is.null(scaling) && scaling$is_scaled) {

    # restore the targets to their original scale
    x_hat <- sweep(x_hat, 2, sx, "*")
    y_hat <- sweep(y_hat, 2, sy, "*")

  }

  # save names
  colnames(x_hat) <- paste(object[["data"]][["xnames"]][rel_x], "_hat", sep = "")
  colnames(y_hat) <- paste(object[["data"]][["ynames"]], "_hat", sep = "")

  # return predictions
  targets <- data.frame (
    x_hat,
    y_hat
  )

  return(targets)

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

  # check if training and test names are equal
  tr_names <- object[["data"]][["xnames"]]
  ts_names <- colnames(newdata)[c(x)]

  if (!identical(sort(tr_names), sort(ts_names))) {
    stop("Different variable names in training data and newdata.")
  }

  # =========================== #
  # 1. SCALING INPUTS           #
  # =========================== #

  # scaling parameters
  scaling <- object[["control"]][["scale"]]

  # data to work with
  data_work <- newdata

  if (!is.null(scaling) && scaling$is_scaled) {

    # original inputs
    raw_x_pred <- as.matrix(data_work[, x])

    # scaled inputs
    scaled_x_pred <- sweep(raw_x_pred, 2, scaling$mean_x, "/")

    # update data
    data_work[, x] <- scaled_x_pred

  }

  # =========================== #
  # 2. PREDICTION               #
  # =========================== #

  # data in [x, y] format with interaction of variables included
  data <- set_data (
    data = data_work,
    x = x,
    y = NULL,
    max_degree = object[["control"]][["max_degree"]]
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

    # prediction
    for (out in 1:nY) {
      y_hat[, out] <- pmax(0, B %*% aces_model[["coefs"]][, out, drop = F])
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

  # =========================== #
  # 3. RESTORE ORIGINAL SCALE   #
  # =========================== #

  if (!is.null(scaling) && scaling$is_scaled) {
    y_hat <- sweep(y_hat, 2, scaling$mean_y, "*")
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
