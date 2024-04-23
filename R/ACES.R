#' @title Fit an Adaptive Constrained Enveloping Splines (ACES) model
#'
#' @description
#'
#' This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity. Both stochastic (as StoNED) and envelopment (as DEA) versions are available. These estimations are based upon the adaptation of the machine learning technique known as Multivariate Adaptive Regression Splines (MARS) developed by \insertCite{friedman1991;textual}{aces}. An adaptation of Random Forest \insertCite{breiman2001}{aces} is also included. For details, see \insertCite{espana2024;textual}{aces}
#'
#' @name aces
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the variables in the model.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param y_type
#' A \code{character} string that determines the prediction approach for \code{y}. It can be either:
#' \itemize{
#' \item{\code{"ind"}} to predict one output at a time and use the rest of outputs as additional inputs (netputs).
#' \item{\code{"all"}} to predict all the outputs at the same time with the original set of inputs.
#' }
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function will estimate. It can be either:
#' \itemize{
#' \item{\code{"env"}}: The model fits an enveloping production frontier.
#' \item{\code{"sto"}}: The model fits a stochastic production frontier.
#' }
#'
#' @param error_type
#' A \code{character} string specifying the error structure that the function will use when fitting the model. It can be either:
#' \itemize{
#' \item{\code{"add"}}: The model assumes an additive error structure.
#' \item{\code{"mul"}}: The model assumes a multiplicative error structure.
#' }
#'
#' @param RF
#' A \code{list} indicating if a bootstrap aggregation methodology as in Random Forest must be used with the following items:
#' \itemize{
#' \item{\code{apply}}: A \code{logical} indicating if a bagging methodology as in Random Forest must be used.
#' \item{\code{sample}}: A \code{numeric} indicating the sample size for bagging.
#' \item{\code{models}}: A \code{numeric} indicating the number of models for bagging.
#' \item{\code{nvars}}: An \code{integer} indicating the number of variables randomly chosen at each split.
#' \item{\code{oob_red}}: A \code{numeric} value specifying the minimum improvement ratio in the moving average of the out-of-bag error over a period of 10 compared to the previous period for the addition of a new model to the ensamble. Default is \code{0.001}.
#' }
#'
#' @param mul_BF
#' A \code{list} specifying the maximum degree of the BFs and the cost of introducing a multivariate BF:
#' \itemize{
#' \item{\code{degree}}: Either a \code{list} with the input indexes for interaction of variables or a \code{numeric} value that determines the maximum degree of interaction between variables. Basis functions products are constrained to contain factors involving distinct variables to ensure interpretability and avoid multicollinearity.
#' \item{\code{hd_cost}}: A \code{numeric} value specifying the minimum percentage of improvement over the best 1-degree basis function to incorporate a higher degree basis function. Default is \code{0.20}.
#' }
#'
#' @param metric
#' A \code{list} specifying the lack-of-fit criterion to evaluate the model performance.
#' \itemize{
#' \item{\code{"mae"}}: Mean Absolute Error.
#' \item{\code{"mape"}}: Mean Absolute Percentage Error.
#' \item{\code{"mse"}}: Mean Squared Error.
#' \item{\code{"msle"}}: Mean Squared Logarithmic Error.
#' \item{\code{"rmse"}}: Root Mean Squared Error.
#' \item{\code{"nrmse1"}}: Normalized Root Mean Squared Error (using mean).
#' \item{\code{"nrmse2"}}: Normalized Root Mean Squared Error (using range).
#' }
#'
#' @param shape
#' A \code{list} with the following items:
#' \itemize{
#' \item{\code{mon}}: A \code{logical} value indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator. If \code{TRUE}, the estimator is constrained to be non-decreasing.
#' \item{\code{con}}: A \code{logical} value indicating whether to enforce the constraint of concavity in the estimator. If \code{TRUE}, the estimator is constrained to be concave.
#' \item{\code{ori}}: A \code{logical} value indicating whether the estimator should satisfy f(0) = 0. If \code{TRUE}, the estimator is constrained to pass through the origin.
#' }
#'
#' @param nterms
#' A positive \code{integer} specifying the maximum number of terms created before pruning. Default is \code{50}.
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree basis functions. Default is \code{0.01}.
#'
#' @param kn_grid
#' Either a \code{-1} (default) the use the original approach of \insertCite{friedman1991;textual}{aces} based on the observed data, or a a \code{list} with the grid of knots to perform ACES. Each element of the \code{list} contains a vector with the knots of one variable (e.g., the second element of the list contains the knot values of the second variable and so on).
#'
#' @param minspan
#' A \code{numeric} value specifying the minimum number of observations between two adjacent knots. It can be one of the following:
#' \itemize{
#' \item{\code{minspan = -2}}: Computed as in \insertCite{zhang1994;textual}{aces}.
#' \item{\code{minspan = -1}}: Computed as in \insertCite{friedman1991;textual}{aces}.
#' \item{\code{minspan = +m}}: A user-specified positive integer.
#' }
#'
#' @param endspan
#' A \code{numeric} value specifying the minimum number of observations before the first and after the final knot. It can be one of the following:
#' \itemize{
#' \item{\code{endspan = -2}}: Computed as in \insertCite{zhang1994;textual}{aces}.
#' \item{\code{endspan = -1}}: Computed as in \insertCite{friedman1991;textual}{aces}.
#' \item{\code{endspan = +m}}: A user-specified positive integer.
#' }
#'
#' @param kn_penalty
#' A positive \code{numeric} value specifying the Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}.
#'
#' @param smoothing
#' A \code{list} specifying conditions for the smoothing procedure:
#' \itemize{
#' \item{\code{wc}}:  A \code{numeric} value used for the cubic smoothing procedure \insertCite{friedman1991}{aces}. This parameter adjusts the distance `e` between the central knot and the right side knot based on the distance `d` between the central knot and the left side knot. If the condition `1 < d / e <  2` is nos satisfied, then `e` is adjusted to be `wc * d`. This parameter must be set between 1 and 2. If a \code{vector} is entered, the \code{wc} value that most reduced the lack-of-fit criterion is selected..
#' \item{\code{wq}}: A \code{numeric} value used for the quintic smoothing procedure \insertCite{chen1999}{aces}. This parameter adjusts the distance `d` between the central knot and the left side knot based on the distance `e` between the central knot and the right side knot. If the condition `8/7 < e / d <  1.5` is nos satisfied, then `d` is adjusted to be `wq * e`. This parameter must be set between 8/7 and 1.5. If a \code{vector} is entered, the \code{wc} value that most reduced the lack-of-fit criterion is selected.
#' }
#'
#' @details
#'
#' This function generates a production frontier satisfying that adheres to certain classical production theory axioms, such as monotonicity and concavity. Users can be choose between enveloping or stochastic versions of the production functions. The algorithm comprises two procedures:
#'
#' 1. A forward selection algorithm that creates a ser of linear basis functions, which may initially overfit the training data.
#'
#' 2. A backward elimination algorithm that remove basis that do not significantly contribute to the model's performance.
#'
#' 3. Two smoothing procedures are available: one using cubic functions and another using quintic functions.
#'
#' If the stochastic version is selected, the frontier's shape is estimated without imposing enveloping constraints on the observations. In the second stage, the expected value of inefficiency is estimated using the residuals obtained from the first stage. This procedure is inspired by \insertCite{kuosmanen2012;textual}{aces}.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{breiman2001}{aces} \cr \cr
#' \insertRef{zhang1994}{aces} \cr \cr
#' \insertRef{chen1999}{aces} \cr \cr
#' \insertRef{kuosmanen2012}{aces}
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc
#' @importFrom extraDistr rsign
#'
#' @return
#'
#' An \code{aces} object.
#'
#' @export

aces <- function (
    data,
    x,
    y,
    y_type = "ind",
    model_type = "env",
    error_type = "add",
    RF = list (
      "apply" = FALSE,
      "sample" = nrow(data),
      "models" = 100,
      "nvars" = ceiling(length(x) / 3),
      "oob_red" = 0.001
      ),
    mul_BF = list (
      "degree" = 1,
      "hd_cost" = 0.05
    ),
    metric = "mse",
    shape = list (
      "mon" = TRUE,
      "con" = TRUE,
      "ori" = FALSE
      ),
    nterms = 50,
    err_red = 0.01,
    kn_grid = - 1,
    minspan = - 1,
    endspan = - 1,
    kn_penalty = 2,
    smoothing = list (
      "wc" = seq(1, 2, length.out = 5),
      "wq" = seq(8 / 7, 1.5, length.out = 5)
    )
  ) {

  # possible error messages:
  display_errors (
    caller = "aces",
    data = data,
    x = x,
    y = y,
    y_type = y_type,
    model_type = model_type,
    error_type = error_type,
    nvars = RF[["nvars"]],
    degree = mul_BF[["degree"]],
    metric = metric,
    nterms = nterms,
    err_red = err_red,
    hd_cost = mul_BF[["hd_cost"]],
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    d = kn_penalty,
    wc = smoothing[["wc"]],
    wq = smoothing[["wq"]],
    object = NULL,
    measure = NULL,
    returns = NULL,
    direction = NULL,
    digits = NULL
  )

  # list with individual ACES models
  ACES <- list()

  if (!RF[["apply"]]) {
    if (y_type == "all") {

      ACES[[1]] <- aces_algorithm (
        data = data,
        x = x,
        y = y,
        z = NULL,
        y_type = y_type,
        model_type = model_type,
        error_type = error_type,
        degree = mul_BF[["degree"]],
        metric = metric,
        shape = list (
          "mon" = shape[["mon"]],
          "con" = shape[["con"]],
          "ori" = shape[["ori"]]
        ),
        nterms = nterms,
        err_red = err_red,
        hd_cost = mul_BF[["hd_cost"]],
        minspan = minspan,
        endspan = endspan,
        kn_grid = kn_grid,
        d = kn_penalty,
        wc = smoothing[["wc"]],
        wq = smoothing[["wq"]]
      )

      # name
      names(ACES)[1] <- "y_all"

    } else {

      for (out in 1:length(y)) {

        # output indicator
        indx <- rep(FALSE, length(y))

        # select the output
        indx[out] <- TRUE

        # outputs as inputs
        nets <- y[!indx]
        outs <- y[indx]

        ACES[[out]] <- aces_algorithm (
          data = data,
          x = x,
          y = outs,
          z = nets,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          degree = mul_BF[["degree"]],
          metric = metric,
          shape = list (
            "mon" = shape[["mon"]],
            "con" = shape[["con"]],
            "ori" = shape[["ori"]]
          ),
          nterms = nterms,
          err_red = err_red,
          hd_cost = mul_BF[["hd_cost"]],
          minspan = minspan,
          endspan = endspan,
          kn_grid = kn_grid,
          d = kn_penalty,
          wc = smoothing[["wc"]],
          wq = smoothing[["wq"]]
        )

        # name
        names(ACES)[out] <- names(data)[y[out]]
      }
    }

  } else {

    # number of models to train
    number_models <- RF[["models"]]

    # sample size for each model
    sample_size <- RF[["sample"]]

    # number of outputs
    nY <- length(y)

    # RF-ACES models
    RF_ACES <- vector("list", number_models)

    # ProgressBar
    pb <- txtProgressBar(min = 0, max = number_models, style = 3)

    # out-of-bag predictions
    oob_pred <- vector("list", number_models)

    # moving-average for the out-of-bag
    oob_vec <- c()
    mov_avg_oob_05 <- c()
    mov_avg_oob_10 <- c()
    mov_avg_oob_25 <- c()

    stopping_condition_counter <- 0

    for (m in 1:length(RF_ACES)) {

      # update ProgressBar
      setTxtProgressBar(pb, m)

      # index of samples in the model
      sample_bag <- sample (
        1:nrow(data),
        size = sample_size,
        replace = TRUE
      )

      # data for the m-model
      data_bag <- data[sample_bag, ]

      # out of bag data for the m-model
      inb_idxs <- 1:nrow(data) %in% sample_bag
      oob_idxs <- which(!(1:nrow(data) %in% sample_bag))
      data_oob <- data[!inb_idxs, ]

      if (y_type == "all") {

        RF_ACES[[m]][[1]] <- rf_aces_algorithm (
          data = data_bag,
          x = x,
          y = y,
          z = NULL,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          degree = mul_BF[["degree"]],
          metric = metric,
          shape = list (
            "mon" = shape[["mon"]],
            "con" = shape[["con"]],
            "ori" = shape[["ori"]]
          ),
          nterms = nterms,
          nvars = RF[["nvars"]],
          err_red = err_red,
          hd_cost = mul_BF[["hd_cost"]],
          minspan = minspan,
          endspan = endspan,
          kn_grid = kn_grid,
          wc = smoothing[["wc"]],
          wq = smoothing[["wq"]]
        )

        # name
        names(RF_ACES[[m]])[1] <- "y_all"

        # predictions
        y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

        for (out in 1:length(y)) {
          Bmatx <- RF_ACES[[m]][[1]][["methods"]][["aces_forward"]][["Bmatx"]]
          coefs <- RF_ACES[[m]][[1]][["methods"]][["aces_forward"]][["coefs"]]
          y_hat[, out] <- Bmatx %*% coefs[, out]
        }

        y_hat <- y_hat[oob_idxs, ]

        oob_pred[[m]] <- matrix (
          c (
            "y_hat" = y_hat,
            "idx" = oob_idxs
            ),
          ncol = nY + 1
        )

        # OOB performance
        oob_mse <- compute_oob (
          data = data,
          y = y,
          oob_idxs = oob_idxs,
          oob_pred = oob_pred,
          model = m
        )

        # out-of-bag error
        RF_ACES[[m]][[1]][["OOB"]] <- oob_mse

        # moving average for OOB (from the forward algorithm)
        oob_vec <- c(oob_vec, RF_ACES[[m]][["y_all"]][["OOB"]])
        mov_avg_oob_05 <- c(mov_avg_oob_05, mean(tail(oob_vec, 05)))
        mov_avg_oob_10 <- c(mov_avg_oob_10, mean(tail(oob_vec, 10)))
        mov_avg_oob_25 <- c(mov_avg_oob_25, mean(tail(oob_vec, 25)))

        if (m > min(nrow(data), 200)) {
          if (oob_vec[m] >= oob_vec[m - 1]) {

            if (oob_vec[m] >= mov_avg_oob_05[m] * (1 - RF[["oob_red"]])) {
              c05 <- TRUE
            } else {
              c05 <- FALSE
            }

            if (oob_vec[m] > mov_avg_oob_10[m] * (1 - RF[["oob_red"]])) {
              c10 <- TRUE
            } else {
              c10 <- FALSE
            }

            if (oob_vec[m] > mov_avg_oob_25[m] * (1 - RF[["oob_red"]])) {
              c25 <- TRUE
            } else {
              c25 <- FALSE
            }

            if (c05 + c10 + c25 == 3) {
              stopping_condition_counter <- stopping_condition_counter + 1
            }

            if (stopping_condition_counter == 5) {
              cat("\n")
              print(paste("Out-of-bag error stabilized after", m, "models"))
              break
            }
          } else {
            stopping_condition_counter == 0
          }
        }

    } else {

      for (out in 1:length(y)) {

        # output indicator
        indx <- rep(FALSE, length(y))

        # select the output
        indx[out] <- TRUE

        # outputs as inputs
        nets <- y[!indx]
        outs <- y[indx]

        RF_ACES[[m]][[out]] <- rf_aces_algorithm (
          data = data_bag,
          x = x,
          y = outs,
          z = nets,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          degree = mul_BF[["degree"]],
          metric = metric,
          shape = list (
            "mon" = shape[["mon"]],
            "con" = shape[["con"]],
            "ori" = shape[["ori"]]
          ),
          nterms = nterms,
          nvars = RF[["nvars"]],
          err_red = err_red,
          hd_cost = mul_BF[["hd_cost"]],
          minspan = minspan,
          endspan = endspan,
          kn_grid = kn_grid,
          wc = smoothing[["wc"]],
          wq = smoothing[["wq"]]
        )

        # name
        names(RF_ACES[[m]])[out] <- names(data)[y[out]]

      }

      # predictions
      y_hat <- matrix(NA, nrow = nrow(data), ncol = length(y))

      for (out in 1:length(y)) {
        Bmatx <- RF_ACES[[m]][[out]][["methods"]][["aces_forward"]][["Bmatx"]]
        coefs <- RF_ACES[[m]][[out]][["methods"]][["aces_forward"]][["coefs"]]
        y_hat[, out] <- Bmatx %*% coefs
      }

      y_hat <- y_hat[oob_idxs, ]

      oob_pred[[m]] <- matrix (
        c (
          "y_hat" = y_hat,
          "idx" = oob_idxs
        ),
        ncol = nY + 1
      )

      # OOB performance
      oob_mse <- compute_oob (
        data = data,
        y = y,
        oob_idxs = oob_idxs,
        oob_pred = oob_pred,
        model = m
      )

      # out-of-bag error
      for (out in 1:length(y)) {
        RF_ACES[[m]][[out]][["OOB"]] <- oob_mse
      }

      # moving average for OOB (from the forward algorithm)
      oob_vec <- c(oob_vec, RF_ACES[[m]][["y_all"]][["OOB"]])
      mov_avg_oob_05 <- c(mov_avg_oob_05, mean(tail(oob_vec, 05)))
      mov_avg_oob_10 <- c(mov_avg_oob_10, mean(tail(oob_vec, 10)))
      mov_avg_oob_25 <- c(mov_avg_oob_25, mean(tail(oob_vec, 25)))

      if (m > min(nrow(data), 200)) {
        if (oob_vec[m] >= oob_vec[m - 1]) {

          if (oob_vec[m] >= mov_avg_oob_05[m] * (1 - RF[["oob_red"]])) {
            c05 <- TRUE
          } else {
            c05 <- FALSE
          }

          if (oob_vec[m] > mov_avg_oob_10[m] * (1 - RF[["oob_red"]])) {
            c10 <- TRUE
          } else {
            c10 <- FALSE
          }

          if (oob_vec[m] > mov_avg_oob_25[m] * (1 - RF[["oob_red"]])) {
            c25 <- TRUE
          } else {
            c25 <- FALSE
          }

          if (c05 + c10 + c25 == 3) {
            stopping_condition_counter <- stopping_condition_counter + 1
          }

          if (stopping_condition_counter == 5) {
            cat("\n")
            print(paste("Out-of-bag error stabilized after", m, "models"))
            break
          }
        } else {
          stopping_condition_counter == 0
        }
      }
    }
  }

  # random-forest hyperparameters
  rf_list = list (
    "apply" = RF[["apply"]],
    "sample" = RF[["sample"]],
    "models" = m,
    "nvars" = RF[["nvars"]],
    "oob_red" = RF[["oob_red"]]
  )

  ACES <- RF_ACES[1:m]

  for (j in seq_along(ACES)) {
    for (out in 1:length(ACES[[1]])) {
      ACES[[j]][[out]][["control"]] <- append (
        ACES[[j]][[out]][["control"]],
        list("RF" = rf_list)
      )
    }
  }

  }

  # type of object
  if (!RF[["apply"]]) {
    class(ACES) <- "aces"

  } else {
    class(ACES) <- "rf_aces"
  }

  return(ACES)
}

#' @title Algorithm of Adaptive Constrained Enveloping Splines (ACES).
#'
#' @description
#'
#' This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity. Both stochastic (as StoNED) and envelopment (as DEA) versions are available. These estimations are based upon the adaptation of the machine learning technique known as Multivariate Adaptive Regression Splines (MARS) developed by \insertCite{friedman1991;textual}{aces}.
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the variables in the model.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param z
#' Column indexes of netput variables in \code{data} (outputs evaluated as inputs). These variables must be considered as output when computing predictions or efficiency scores.
#'
#' @param y_type
#' A \code{character} string that determines the prediction approach for \code{y}.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function will estimate.
#'
#' @param error_type
#' A \code{character} string specifying the error structure that the function will use when fitting the model.
#'
#' @param degree
#' Maximum degree of interaction between variables.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param nterms
#' Maximum number of terms created before pruning.
#'
#' @param err_red
#' Minimum reduced error rate for the addition of a new pair of 1-degree basis functions.
#'
#' @param hd_cost
#' Minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.
#'
#' @param minspan
#' Minimum number of observations between two adjacent knots.
#'
#' @param endspan
#' Minimum number of observations before the first and after the final knot.
#'
#' @param kn_grid
#' Design of the grid of knots to perform ACES
#'
#' @param d
#' Generalized Cross Validation (GCV) penalty per knot.
#'
#' @param wc
#' Hyperparameter for the side knot distances in the cubic smoothing procedure \insertCite{friedman1991}{aces}.
#'
#' @param wq
#' Hyperparameter for the side knot distances in the quintic smoothing procedure \insertCite{chen1999}{aces}.
#'
#' @references
#'
#' \insertRef{zhang1994}{aces} \cr \cr
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{chen1999}{aces}
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc
#'
#' @return
#'
#' An \code{aces} object.
#'
#' @export

aces_algorithm <- function (
    data,
    x,
    y,
    z,
    y_type,
    model_type,
    error_type,
    degree,
    metric,
    shape,
    nterms,
    err_red,
    hd_cost,
    minspan,
    endspan,
    kn_grid,
    d,
    wc,
    wq
    ) {

  # shape constraints
  monotonicity <- shape[["mon"]]
  concavity    <- shape[["con"]]
  origin       <- shape[["ori"]]

  # accelerate the algorithm by selecting the efficient DMUs in DEA to impose only
  # the envelope on these points.
  # only valid under monotonicity and concavity constraints

  if (monotonicity && concavity) {
    dea_scores <- rad_out (
      tech_xmat = as.matrix(data[, x]),
      tech_ymat = as.matrix(data[, y]),
      eval_xmat = as.matrix(data[, x]),
      eval_ymat = as.matrix(data[, y]),
      convexity = TRUE,
      returns = "variable"
    )

    dea_eff <- c(1:nrow(data))[dea_scores <= 1.001]

  } else {

    dea_eff <- c(1:nrow(data))
  }

  # original input indexes
  inps <- x

  # original output indexes
  outs <- y

  # original netput indexes
  nets <- z

  # original set of DMUs
  dmus <- data

  # change to logarithmic scale if multiplicative error type
  if (error_type == "mul") {
    data[, c(y)] <- log(data[, c(y)])
  }

  # origin
  if (origin && error_type == "add") {
    origin <- "0"

  } else if (origin && error_type == "mul") {
    origin <- "1"

  } else {
    origin <- "FALSE"

  }

  # data in [x, z, y] format with interaction of variables included
  data <- set_interactions (
    data = data,
    x = x,
    y = y,
    z = z,
    degree = degree
    )

  # samples in data
  N <- nrow(data)

  # reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  # number of inputs / outputs as inputs
  nX <- length(x)

  # mumber of outputs
  nY <- length(y)

  # matrix with:
  # row 1: the index of the variable
  # row 2: the degree of the variable
  xi_degree <- matrix (
    c(x, rep(0, length(x))),
    byrow = TRUE,
    nrow = 2,
    ncol = length(x)
    )

  if (!is.list(degree)) {

    v <- 0

    for (i in 1:degree) {
      combs <- combn(1:(length(inps) + length(nets)), i)

      for (k in 1:ncol(combs)) {
        v <- v + 1
        xi_degree[2, v] <- i
      }
    }

  } else {
    xi_degree[2, inps] <- 1

    for (k in 1:length(degree)) {
      xi_degree[2, length(inps) + k] <- length(degree[[k]])
    }
  }

  # ===================== #
  #   FORWARD ALGORITHM   #
  # ===================== #

  if (model_type == "env") {
    y_hat <- matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, max)

  } else {
    y_hat <- matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, mean)

  }

  # basis function
  #     id: index
  # status: intercept / paired / unpaired
  #   side: E (entire) / R (right) / L (left)
  #     Bp: basis function
  #     xi: variable for splitting
  #      t: knot for splitting
  #      R: mean error between true data and predicted data (B %*% coefs)
  #    GCV: generalized cross validation
  #   GRSq: generalized residual of squares
  #  coefs: regression coefficients

  bf <- list (
    "id" = 1,
    "status" = "intercept",
    "side" = "E",
    "Bp" = rep(1, N),
    "xi" = c(- 1),
    "t" = c(- 1),
    "R" = err_metric(data[, y, drop = F], y_hat, metric),
    "GCV" = compute_gcv(data[, y, drop = F], y_hat, metric, 1, 1, 0, xi_degree),
    "GRSq" = 0,
    "coefs" = unname(apply(data[, y, drop = FALSE], 2, max))
  )

  # set of knots. It saves indexes of data used as knots.
  kn_list <- vector("list", nX)

  # set of basis functions by variable.
  Bp_list <- vector("list", nX)

  for (xi in 1:nX) {
    Bp_list[[xi]] <- list (
      "paired" = NULL,
      "right" = NULL,
      "left" = NULL
      )
  }

  # set of basis functions (bf_set) and the matrix of basis functions (B)
  aces_forward <- list (
    "bf_set" = list(bf),
    "B" = matrix(rep(1, N))
    )

  # error of the first basis function
  err <- bf[["R"]]

  # minimum span (minspan) and end span (endspan)
  L_Le <- compute_span (
    data = data,
    minspan = minspan,
    endspan = endspan,
    nX = nX
    )

  # minimum span
  L  <- L_Le[[1]]

  # end span
  Le <- L_Le[[2]]

  # set the grid of knots
  kn_grid <- set_kn_grid (
    data = data,
    nX = nX,
    inps = length(inps),
    kn_grid = kn_grid
    )

  # initial error
  err_min <- err

  while(length(aces_forward[["bf_set"]]) + 2 < nterms) {

    # negative GRSq
    last_bf <- aces_forward[["bf_set"]][[length(aces_forward[["bf_set"]])]]
    if (last_bf[["GRSq"]] < 0) break

    # add 2 new basis functions to the model:
    B_bf_knt_err <- add_basis_function (
      data = data,
      x = x,
      y = y,
      xi_degree = xi_degree,
      dea_eff = dea_eff,
      model_type = model_type,
      metric = metric,
      forward_model = aces_forward,
      Bp_list = Bp_list,
      monotonicity = monotonicity,
      concavity = concavity,
      origin = origin,
      kn_list = kn_list,
      kn_grid = kn_grid,
      L = L,
      Le = Le,
      err_min = err,
      hd_cost = hd_cost
      )

    if (!is.list(B_bf_knt_err)) {
      break
    } else {
      new_err <- B_bf_knt_err[[5]]
    }

    # t_old <- length(aces_forward[["bf_set"]])
    # t_new <- length(aces_forward[["bf_set"]]) + 2
    # prop_err_red <- round(100 - round(new_err[1] * 100 / err[1]), 2)
    # print(paste("From:", t_old, "To:", t_new, "-->", prop_err_red, "% of reduction"))

    if (new_err[1] < err[1] * (1 - err_red[1])) {

      # update B
      aces_forward[["B"]] <- B_bf_knt_err[[1]]

      # update basis functions
      aces_forward[["bf_set"]] <- B_bf_knt_err[[2]]

      # update the knots list
      kn_list <- B_bf_knt_err[[3]]

      # update the Bp list
      Bp_list <- B_bf_knt_err[[4]]

      # update error
      err <- new_err

    } else {
      break

    }
  }

  # set of knots from forward algorithm
  var <- c()
  knt <- c()
  sts <- c()

  for (v in 1:nX) {
    for (side in c("paired", "right", "left")) {
      if (!is.null(Bp_list[[v]][[side]])) {
        # knots in variable "xi"
        knt_xi <- sapply(Bp_list[[v]][[side]], "[[", "t")

        # knots
        knt <- c(knt, knt_xi)
        # variable
        var <- c(var, rep(v, length(knt_xi)))
        # status
        sts <- c(sts, rep(side, length(knt_xi)))
      }
    }
  }

  kn_forward <- data.frame (
    xi = var,
    t = knt,
    status = sts
    )

  # ==
  # forward aces
  # ==

  aces_forward = list (
    "basis" = aces_forward[["bf_set"]],
    "Bmatx" = aces_forward[["B"]],
    "knots" = kn_forward,
    "coefs" = rev(aces_forward[["bf_set"]])[[1]][["coefs"]]
  )

  # ====================== #
  #   BACKWARD ALGORITHM   #
  # ====================== #

  aces_submodels <- aces_pruning (
      data = data,
      y = y,
      xi_degree = xi_degree,
      dea_eff = dea_eff,
      model_type = model_type,
      metric = metric,
      forward_model = aces_forward,
      Bp_list = Bp_list,
      monotonicity = monotonicity,
      concavity = concavity,
      origin = origin,
      d = d
    )

  # generalized cross-validation for each model
  GCVs <- sapply(aces_submodels, function(x) x[["GCV"]])

  # model with minimum error (excluding the model without knots)
  aces_backward <- aces_submodels[[which.min(GCVs[1:(length(aces_submodels) - 1)])]]

  # set of surviving knots
  kn_backward <- do.call(rbind.data.frame, aces_backward[["t"]])

  # sort the knots by "xi", "status", "side", "t"
  knots_backward_order <- order (
    kn_backward$xi,
    kn_backward$status,
    ifelse(kn_backward$status == "paired", kn_backward$t, desc(kn_backward$side)),
    ifelse(kn_backward$status == "paired", desc(kn_backward$side), kn_backward$t)
    )

  # update the set of knots
  kn_backward <- kn_backward[knots_backward_order, ]

  # ==
  # aces
  # ==

  aces <- list (
    "aces_submodels" = aces_submodels,
    "Bmatx" = aces_backward[["B"]],
    "knots" = kn_backward,
    "coefs" = aces_backward[["coefs"]],
    "GCV" = aces_backward[["GCV"]]
  )

  # ======================= #
  #   SMOOTHING PROCEDURE   #
  # ======================= #

  # sub-models sorted by gcv value
  aces_submodels_gcv <- aces_submodels[order(sapply(aces_submodels, "[[", "GCV"))]

  # initialize a list with smoothed sub-models
  aces_smoothed_submodels <- vector("list", length(aces_submodels_gcv))

  for (s in 1:length(aces_submodels_gcv)) {

    # select a model to be smoothed
    aces_smoothed <- aces_submodels_gcv[[s]]

    # skip if there are not knots
    if (is.null(aces_smoothed[["t"]])) next

    # data.frame of knots
    kn_smoothed <- do.call(rbind.data.frame, aces_smoothed[["t"]])

    # if monotonicity is required:
    # 1- wc in (1, 2) & wq in (8/7, 1.5)

    # If concavity is required:
    # 1- wc in (1, 2) & wq in (8/7, 1.5)
    # 2- unpaired right basis functions are not allowed

    # check for right-side unpaired basis functions
    check1 <- kn_smoothed$side == "R"
    check2 <- kn_smoothed$status == "unpaired"

    if (concavity && max(check1 + check2) == 2) {
      next

    } else {
      aces_smoothed_submodels[[s]][["Model"]] <- aces_smoothed
      aces_smoothed_submodels[[s]][["Knots"]] <- kn_smoothed

      }
    }

    # initialize list of cubic aces models
    cubic_aces_models <- vector("list", length(aces_smoothed_submodels))

    # initialize list of quintic aces models
    quintic_aces_models <- vector("list", length(aces_smoothed_submodels))

    for (m in 1:length(aces_smoothed_submodels)) {
      if (is.null(aces_smoothed_submodels[[m]])) {
        next

      } else {

        # select a model
        aces_smoothed <- aces_smoothed_submodels[[m]][["Model"]]

        # select a set of knots
        kn_smoothed <- aces_smoothed_submodels[[m]][["Knots"]]

        # generate the input space for side knots location
        kn_side_loc <- side_knots_location (
          data = data,
          nX = nX,
          knots = kn_smoothed
        )

        # ==
        # smoothing cubic aces
        # ==

        cubic_aces_models[[m]] <- cubic_aces (
          data = data,
          x = x,
          y = y,
          dea_eff = dea_eff,
          model_type = model_type,
          metric = metric,
          monotonicity = monotonicity,
          concavity = concavity,
          origin = origin,
          kn_grid = kn_smoothed,
          kn_side_loc = kn_side_loc,
          d = 0.25,
          wc = wc
        )

        # ==
        # smoothing quintic aces
        # ==

        quintic_aces_models[[m]] <- quintic_aces (
          data = data,
          x = x,
          y = y,
          dea_eff = dea_eff,
          model_type = model_type,
          metric = metric,
          monotonicity = monotonicity,
          concavity = concavity,
          origin = origin,
          kn_grid = kn_smoothed,
          kn_side_loc = kn_side_loc,
          d = 0.25,
          wq = wq
        )
      }
    }

    # GCVs of cubic models
    aces_cubic_gcvs <- sapply (
      cubic_aces_models,
      function (x) {
        ifelse (
          is.null(x[["GCV"]]),
          Inf,
          x[["GCV"]]
          )
        }
      )

    min_gcv <- which.min(aces_cubic_gcvs)

    # cubic aces
    aces_cubic <- cubic_aces_models[[min_gcv]]

    # GCVs of quintic models
    aces_quintic_gcvs <- sapply (
      quintic_aces_models,
      function (x) {
        ifelse (
          is.null(x[["GCV"]]),
          Inf,
          x[["GCV"]]
          )
        }
      )

    min_gcv <- which.min(aces_quintic_gcvs)

    # quintic aces
    aces_quintic <- quintic_aces_models[[min_gcv]]

    # ==================================================== #
    #  STOCHASTIC ADAPTIVE CONSTRAINED ENVELOPING SPLINES  #
    # ==================================================== #

    if (model_type == "sto") {
      for (model in c("aces_forward", "aces", "aces_cubic", "aces_quintic")) {

        # =========
        # Residuals
        # =========

        if (model == "aces_forward") {
          y_hat <- aces_forward[["Bmatx"]] %*% aces_forward[["coefs"]]

        } else if (model == "aces") {
          y_hat <- aces[["Bmatx"]] %*% aces[["coefs"]]

        } else if (model == "aces_cubic") {
          y_hat <- aces_cubic[["Bmatx"]] %*% aces_cubic[["coefs"]]

        } else {
          y_hat <- aces_quintic[["Bmatx"]] %*% aces_quintic[["coefs"]]

        }

        if (error_type == "mul") {
          mean_pred <- exp(y_hat)

        } else {
          mean_pred <- y_hat
        }

        mean_pred <- rad_out (
          tech_xmat = as.matrix(data[, inps]),
          tech_ymat = as.matrix(mean_pred),
          eval_xmat = as.matrix(data[, inps]),
          eval_ymat = as.matrix(mean_pred),
          convexity = TRUE,
          returns = "variable"
        ) * mean_pred

        if (error_type == "mul") {
          mean_pred <- log(mean_pred)

        } else {
          mean_pred <- mean_pred
        }

        residuals <- data[, y] - mean_pred

        # H0: residuals are symmetrically distributed
        # H1: residuals are negatively skewed

        sqrt_b1_btp <- c()
        for (m in 1:100) {
          # number of residulas
          nresid <- length(residuals)

          # stage 1: re-centered residuals
          recentered_residuals <- residuals - mean(residuals)

          # stage 2: symmetrized residuals
          random_signs <- rsign(nresid)
          symmetrized_residuals <- random_signs * recentered_residuals

          # stage 3: re-sample
          resample_idx <- sample(1:nresid, nresid, replace = TRUE)
          bootstrapped <- symmetrized_residuals[resample_idx]

          # stage 4: compute the sqrt_b1 statistic for the bootstrapped sample
          M2_btp <- mean((bootstrapped - mean(bootstrapped)) ^ 2)
          M3_btp <- mean((bootstrapped - mean(bootstrapped)) ^ 3)

          sqrt_b1_btp[m] <- M3_btp / M2_btp ^ (3 / 2)
        }

        # 1 - alpha percentile
        critical_value <- quantile(sqrt_b1_btp, 0.05, na.rm = TRUE)

        # the sqrt_b1 test from a symmetric data sample takes a value lower than
        # critical_value the 95% of the times.

        # compute the sqrt_b1 statistic
        M2 <- mean((residuals - mean(residuals)) ^ 2)
        M3 <- mean((residuals - mean(residuals)) ^ 3)

        if (M3 > 0) {
          M3 <- - 0.0001
        }

        sqrt_b1 <- M3 / M2 ^ (3 / 2)

        if (sqrt_b1 > critical_value) {
          residual_shape <- "symmetric"
          # warning (
          #   "The statistical evidence is not sufficient to reject the hypothesis
          # that the residuals are symmetrically distributed."
          # )
        } else {
          residual_shape <- "left-skewed"

        }

        # ==
        # Efficiency estimation by the method of moments
        # ==

        # standard deviation for inefficiency term
        std_u_mm <- (M3 / (sqrt(2 / pi) * (1 - 4 / pi))) ^ (1 / 3)

        # standard deviation for error term
        std_v_mm <- sqrt(M2 - ((pi - 2) / pi) * std_u_mm ^ 2)

        if (is.nan(std_v_mm)) {
          std_v_mm <- 0.0001
        }

        # conditional efficiency
        eff_score_mm <- cond_eff (
          data = data,
          y = y,
          mean_pred = mean_pred,
          std_ineff = std_u_mm,
          std_error = std_v_mm,
          error_type = error_type
        )

        # ==
        # Efficiency estimation by pseudolikelihood
        # ==

        sigma_hat <- function (l) {
          sqrt(mean(residuals ^ 2) / (1 - ((2 * l ^ 2) / (pi * (1 + l ^ 2)))))
        }

        epsilon_hat <- function (l) {
          residuals - (sqrt(2) * l * sigma_hat(l)) / (pi * (1 + l ^ 2)) ^ (1 / 2)
        }

        log_like <- function(l) {
          # log-likelihood function
          logmv <- - N * log(sigma_hat(l)) +
            sum(pnorm((- epsilon_hat(l) * l) / sigma_hat(l), log.p = TRUE)) -
            0.5 * (sum(epsilon_hat(l) ^ 2) / sigma_hat(l) ^ 2)

          return(- sum(logmv))
        }

        lambda_hat <- optim (
          c(1), log_like,
          method = c("L-BFGS-B"),
          lower = 0
        )[["par"]]

        # standard deviation for error term
        std_v_pl <- (sigma_hat(lambda_hat) ^ 2 / (1 + lambda_hat ^ 2)) ^ (1 / 2)

        # standard deviation for inefficiency term
        std_u_pl <- std_v_pl * lambda_hat

        if (std_v_pl < 0) {
          std_v_pl <- 0.0001
        }

        # conditional efficiency
        eff_score_pl <- cond_eff (
          data = data,
          y = y,
          mean_pred = mean_pred,
          std_ineff = std_u_pl,
          std_error = std_v_pl,
          error_type = error_type
        )

        # results
        stochastic_analysis <- list (
          "mean_pred" = mean_pred,

          "residuals" = residuals,
          "residual_shape" = residual_shape,

          "mom" = list (
            "std_u" = std_u_mm,
            "std_v" = std_v_mm,
            "eff_score" = eff_score_mm
          ),

          "pse" = list (
            "std_u" = std_u_pl,
            "std_v" = std_v_pl,
            "eff_score" = eff_score_pl
          )
        )

        if (model == "aces_forward") {
          aces_forward[["sto"]] <- stochastic_analysis

        } else if (model == "aces") {
          aces[["sto"]] <- stochastic_analysis

        } else if (model == "aces_cubic") {
          aces_cubic[["sto"]] <- stochastic_analysis

        } else {
          aces_quintic[["sto"]] <- stochastic_analysis

        }
      }
    }

    # =========== #
    # ACES OBJECT #
    # =========== #

    ACES <- aces_object (
      data = dmus,
      x = inps,
      y = outs,
      z = nets,
      y_type = y_type,
      model_type = model_type,
      error_type = error_type,
      degree = degree,
      metric = metric,
      monotonicity = monotonicity,
      concavity = concavity,
      origin = origin,
      nterms = ncol(aces_forward[["Bmatx"]]),
      err_red = err_red,
      hd_cost = hd_cost,
      minspan = minspan,
      endspan = endspan,
      kn_grid = kn_grid,
      d = d,
      wc = wc,
      wq = wq,
      aces_forward = aces_forward,
      aces = aces,
      aces_cubic = aces_cubic,
      aces_quintic = aces_quintic
    )

    return(ACES)
}

#' @title Create an aces object
#'
#' @description
#'
#' This function saves information about the Adaptive Constrained Enveloping Splines (ACES) model.
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the variables in the model.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param z
#' Column indexes of netput variables in \code{data} (outputs evaluated as inputs).
#'
#' @param y_type
#' A \code{character} string that determines the prediction approach for \code{y}.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function will estimate.
#'
#' @param error_type
#' A \code{character} string specifying the error structure that the function will use when fitting the model.
#'
#' @param degree
#' Maximum degree of interaction between variables.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.
#'
#' @param monotonicity
#' A \code{logical} value indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#'
#' @param concavity
#' A \code{logical} value indicating whether to enforce the constraint of concavity in the estimator.
#'
#' @param origin
#' A \code{logical} value indicating whether the estimator should satisfy f(0) = 0.
#'
#' @param nterms
#' Maximum number of terms created before pruning.
#'
#' @param err_red
#' Minimum reduced error rate for the addition of a new pair of 1-degree basis functions.
#'
#' @param hd_cost
#' Minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.
#'
#' @param minspan
#' Minimum number of observations between two adjacent knots.
#'
#' @param endspan
#' Minimum number of observations before the first and after the final knot.
#'
#' @param kn_grid
#' Design of the grid of knots to perform ACES
#'
#' @param d
#' Generalized Cross Validation (GCV) penalty per knot.
#'
#' @param wc
#' Hyperparameter for the side knot distances in the cubic smoothing procedure \insertCite{friedman1991}{aces}.
#'
#' @param wq
#' Hyperparameter for the side knot distances in the quintic smoothing procedure \insertCite{chen1999}{aces}.
#'
#' @param aces_forward
#' A \code{list} containing the forward step of the Adaptive Constrained Enveloping Splines model.
#'
#' @param aces
#' A \code{list} containing the Adaptive Constrained Enveloping Splines model.
#'
#' @param aces_cubic
#' A \code{list} containing the smoothed version of Adaptive Constrained Enveloping Splines through cubic basis functions.
#'
#' @param aces_quintic
#' A \code{list} containing the smoothed version of Adaptive Constrained Enveloping Splines through quintic basis functions.
#'
#' @references
#'
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{chen1999}{aces}
#'
#' @return
#'
#' An \code{aces} object.

aces_object <- function (
    data,
    x,
    y,
    z,
    y_type,
    model_type,
    error_type,
    degree,
    metric,
    monotonicity,
    concavity,
    origin,
    nterms,
    err_red,
    hd_cost,
    minspan,
    endspan,
    kn_grid,
    d,
    wc,
    wq,

    aces_forward,
    aces,
    aces_cubic,
    aces_quintic
    ) {

  object <- list()

  object[["data"]] <- list (
    "df" = data,
    "x" = x,
    "z" = z,
    "y" = y,
    "xnames" = colnames(data)[c(x, z)],
    "ynames" = colnames(data)[y],
    "y_type" = y_type,
    "rownames" = rownames(data)
  )

  object[["control"]] <- list (
    "model_type" = model_type,
    "error_type" = error_type,
    "degree" = degree,
    "metric" = metric,
    "monotonicity" = monotonicity,
    "concavity" = concavity,
    "origin" = as.logical(origin),
    "nterms" = nterms,
    "err_red" = err_red,
    "hd_cost" = hd_cost,
    "minspan" = minspan,
    "endspan" = endspan,
    "kn_grid" = kn_grid,
    "d" = d,
    "wc" = wc,
    "wq" = wq
  )

  object[["methods"]] <- list (
    "aces_forward" = aces_forward,
    "aces" = aces,
    "aces_cubic" = aces_cubic,
    "aces_quintic" = aces_quintic
  )

  return(object)
}

#' @title Error Messaging in aces functions.
#'
#' @description
#'
#' This function displays error messages if hyperparameters are bad introduced.
#'
#' @param caller
#' A \code{character} string specifying the function that calls \code{display_errors}.
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the variables in the model.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param y_type
#' A \code{character} string that determines the prediction approach for \code{y}. It can be either: \code{"ind"} or \code{"all"}.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function will estimate. It can be either: \code{"env"} or \code{"sto"}.
#'
#' @param error_type
#'  A \code{character} string specifying the error structure that the function will use when fitting the model. It can be either: \code{"add"} or \code{"mul"}.
#'
#' @param nvars
#' An \code{integer} indicating the number of variables randomly chosen at each split in RF-ACES.
#'
#' @param degree
#' Either a \code{list} with the input indexes for interaction of variables or a \code{numeric} value that determines the maximum degree of interaction between variables.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance. It can be: \code{"mae"}, \code{"mape"}, \code{"mse"}, \code{"rmse"}, \code{"nrmse1"} or \code{"nrmse2"}.
#'
#' @param nterms
#' A positive \code{integer} specifying the maximum number of terms created before pruning
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree basis functions.
#'
#' @param hd_cost
#' A \code{numeric} value specifying the minimum percentage of improvement over the best 1-degree basis function to incorporate a higher degree basis function.
#'
#' @param minspan
#' A \code{numeric} value specifying the minimum number of observations between two adjacent knots. It can be: \code{"-2"}, \code{"-1"} or \code{"m"}.
#'
#' @param endspan
#' A \code{numeric} value specifying the minimum number of observations before the first and after the final knot. It can be: \code{"-2"}, \code{"-1"} or \code{"m"}.
#'
#' @param kn_grid
#' Grid of knots to perform ACES. It can be: \code{-1} or a \code{list}.
#'
#' @param d
#' A positive \code{numeric} value specifying the Generalized Cross Validation (GCV) penalty per knot.
#'
#' @param wc
#' A numeric value used for the cubic smoothing procedure.
#'
#' @param wq
#' A numeric value used for the quintic smoothing procedure.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param measure
#' Mathematical programming model to compute the efficiency scores.
#'
#' @param returns
#' Type of returns to scale for computing the efficiency scores.
#'
#' @param direction
#' Vector direction for DMU projection in Directional Distance Function when computing the efficiency scores.
#'
#' @param digits
#' Number of digits to round efficiency scores.
#'
#' @importFrom stats na.omit
#'
#' @return
#' This function return error messages if hyperparameters are incorrectly specified.

display_errors <- function (
    caller,
    data,
    x,
    y,
    y_type,
    model_type,
    error_type,
    nvars,
    degree,
    metric,
    nterms,
    err_red,
    hd_cost,
    minspan,
    endspan,
    kn_grid,
    d,
    wc,
    wq,
    object,
    measure,
    returns,
    direction,
    digits
    ) {

  if (caller == "aces") {

    # data is a matrix or a data.frame
    if (is.list(data) && !is.data.frame(data)) {
      stop("data must be a data.frame or a matrix")
    }

    # x in data
    tryCatch(data[, x], error = function(e) {
      message("Index values from x are not in data.")
    })

    # y in data
    tryCatch(data[, y], error = function(e) {
      message("Index values from y are not in data.")
    })

    # variables classes are valid
    if (!all(sapply(data[, c(x, y)], is.numeric))) {
      stop("data variables must be numeric. Please, check sapply(data, is.numeric))")
    }

    # NA values
    if (any(is.na(data[, c(x, y)]))) {
      data <- na.omit(data[, c(x, y)])
      warning("Rows with NA values have been omitted .\n")
    }

    # y_type must be "ind" or "all"
    if (!is.null(y_type) && !y_type %in% c("ind", "all")) {
      stop("Not available y_type. Please, check help(\"aces\")")
    }

    # model_type must be "env" or "sto"
    if (!is.null(model_type) && !model_type %in% c("env", "sto")) {
      stop("Not available model_type. Please, check help(\"aces\")")
    }

    # error_type must be "add" or "mul"
    if (!is.null(error_type) && !error_type %in% c("add", "mul")) {
      stop("Not available error_type. Please, check help(\"aces\")")
    }

    # nvars lower than number of inputs
    if (y_type == "ind") {
      nets <- length(y) - 1
    } else {
      nets <- 0
    }

    if (!is.null(nvars) && nvars > (length(x) + nets)) {
      stop("nvars must be lower than the number of inputs.")
    }

    # degree must be a valid number
    if (!is.list(degree)) {
      case1 <- y_type == "ind" && degree > length(x) + length(y) - 1
      case2 <- y_type == "all" && degree > length(x)

      if (case1 || case2) {
        stop("degree must be lower than the number of inputs.")
      }
    }

    # the lack-of-fit criterion must be a valid measure
    if (!metric %in% c("mae", "mape", "mse", "msle", "rmse", "nrmse1", "nrmse2")) {
      stop(paste(metric, "is not available. Please, check help(\"aces\")"))
    }

    # nterms must be a positive integer
    if (floor(nterms) != nterms) {
      stop("nterms must be a positive integer.")
    }

    # err_red must be between 0 and 1
    if (!is.null(err_red) && !(err_red >= 0 && err_red <= 1)) {
      stop("err_red must be between 0 and 1.")
    }

    # hd_cost must be between 0 and 1
    if (!(hd_cost >= 0 && hd_cost <= 1)) {
      stop("hd_cost must be between 0 and 1.")
    }

    # minspan must -2, -1 or a positive integer
    if (minspan < -2 || floor(minspan) != minspan) {
      stop("minspan must be - 2, - 1 or a positive integer.")
    }

    # endspan must -2, -1 or a positive integer
    if (endspan < -2 || floor(endspan) != endspan) {
      stop("endspan must be - 2, - 1 or a positive integer.")
    }

    # the kn_grid must have the same length that x
    if (is.list(kn_grid) && length(kn_grid) != length(x)) {
      stop ("If kn_grid is entered, it must have the same length that x.")
    }

    # d must be a semi-positive number
    if (!is.null(d) && d < 0) {
      stop("d must be greater than 0.")
    }

    # wc must be a valid number to ensure shape-constraints:
    if (!all(wc >= 1 & wc <= 2)) {
      stop("wc must be between 1 and 2")
    }

    # wq must be a valid number to ensure shape-constraints:
    if (!all(wq >= 8 / 7 & wq <= 1.5)) {
      stop("wq must be between 8/7 and 1.5")
    }

  } else {

    if (!class(object) %in% c("aces", "rf_aces")) {
      stop(paste(deparse(substitute(object)), "must be an aces or and rf_aces object."))
    }

    if (!is.null(measure) && !measure %in% c("rad_out", "rad_inp", "ddf", "rsl_out", "rsl_inp", "wam")) {
      stop(paste(measure, "is not available. Please, check help(\"aces_scores\")"))
    }

    if (!is.null(returns) && !returns %in% c("constant", "variable")) {
      stop(paste(returns, "is not available. Please, check help(\"aces_scores\")"))
    }

    if (!is.null(direction) && !direction %in% c("mean", "briec")) {
      stop(paste(direction, "is not available. Please, check help(\"aces_scores\")"))
    }

    if (digits < 0) {
      stop("digits must be greater than 0.")
    }

  }
}

#' @title Create Additional Inputs through Interaction Between Variables.
#'
#' @description This function creates additional inputs through the interaction of the original variables and returns a matrix with the data in a [x,y] format.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param x Column indexes of input variables in \code{data}.
#' @param y Column indexes of output variables in \code{data}.
#' @param z Column indexes of netput variables in \code{data} (outputs evaluated as inputs).
#' @param degree A \code{list} with the input indexes for interaction of variables or a \code{numeric} value that determines the maximum degree of interaction between variables. Basis functions products are constrained to contain factors involving distinct variables to ensure interpretability and avoid multicollinearity.
#'
#' @return A \code{matrix} in a [x,y] format with the interaction of variables included.

set_interactions <- function (
    data,
    x,
    y,
    z,
    degree
    ) {

  # 1. Transform the output(s) as input(s) through the inverse: netput(s)
  data[, z] <- 1 / data[, z]

  # 2. Interaction effects
  if (is.list(degree) || degree > 1) {
    if (!is.list(degree)) {
      # Create a list with all the possible combinations between 1 and as much l(x) + l(y) elements
      max_degree <- degree
      degree <- list()
      for (i in 2:max_degree) {
        combs <- combn(1:(length(x) + length(z)), i)
        for (col in 1:ncol(combs)) {
          degree <- append(degree, list(combs[, col]))
        }
      }
    }

    # Number of additional variables
    IVars <- length(degree)

    # New x indexes
    x <- c(x, z, (ncol(data) + 1):(ncol(data) + IVars))

    # Create the new variables
    for (p in 1:IVars) {
      vars <- x[degree[[p]]]
      name_vars <- colnames(data)[vars]
      name <- paste(name_vars, collapse = "_")
      data[, name] <- apply(data[, vars], 1, prod)
    }
  } else {
    x <- c(x, z)
  }

  # 3. Data in the correct order
  data <- data[, c(x, y)]

  return(as.matrix(data))
}

#' @title Error metric for model evaluation.
#'
#' @description
#'
#' This function computes an error metric for model evaluation.
#'
#' @param y_obs
#' Vector of observed data.
#'
#' @param y_hat
#' Vector of predicted values.
#'
#' @param metric
#' Lack-of-fit criterion to evaluate the model performance:
#' \itemize{
#' \item{\code{mae}}: Mean Absolute Error
#' \item{\code{mape}}: Mean Absolute Percentage Error
#' \item{\code{mse}}: Mean Squared Error
#' \item{\code{rmse}}: Root Mean Squared Error
#' \item{\code{nrmse1}}: Normalized Root Mean Squared Error (using mean)
#' \item{\code{nrmse2}}: Normalized Root Mean Squared Error (using range)
#' }
#'
#' @return
#' The calculated error metric.

err_metric <- function (
    y_obs,
    y_hat,
    metric
    ) {

  # samples in data
  N <- nrow(y_obs)

  # number of outputs
  nY <- ncol(y_obs)

  if (all(y_obs > 0) && any(y_hat < 0) ) {
    # do not "allow" negative predictions (negative outputs)
    error <- Inf

  } else if (metric == "mae") {
    # mean absolute error
    devtn <- abs(y_hat - y_obs)
    error <- sum(devtn) / (N * nY)

  } else if (metric == "mape") {
    # mean absolute percentage error
    devtn <- abs(y_hat - y_obs) / y_obs
    error <- sum(devtn) / (N * nY) * 100

  } else if (metric == "mse") {
    # mean squared error
    devtn <- (y_hat - y_obs) ^ 2
    error <- sum(devtn) / (N * nY)

  } else if (metric == "msle") {
    # mean squared logarithmic error
    devtn <- (log(y_hat + 1) - log(y_obs + 1)) ^ 2
    error <- sum(devtn) / (N * nY)

  } else if (metric == "rmse") {
    # root mean squared error
    devtn <- (y_hat - y_obs) ^ 2
    error <- sqrt(sum(devtn) / (N * nY))

  } else if (metric == "nrmse1") {
    # normalized root mean squared error by the mean
    devtn <- (y_hat - y_obs) ^ 2
    error <- sqrt(sum(devtn) / (N * nY)) / mean(y_obs)

  } else {
    # calculate the mean of column-wise maximums and minimums in y
    ymax <- mean(apply(y_obs, 2, max))
    ymin <- mean(apply(y_obs, 2, min))

    # normalized root mean squared error by the range
    devtn <- (y_hat - y_obs) ^ 2
    error <- sqrt(sum(devtn) / (N * nY)) / (ymax - ymin)
  }

  return(error)
}

#' @title Compute Minimum and End Spans
#'
#' @description This function computes the minimum span (i.e., the minimum number of observations between two adjacent knots) and the end span (i.e., the minimum number of observations before the first and after the final knot).
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param minspan Minimum number of observations between two adjacent knots.
#' \itemize{
#' \item{\code{minspan = -2}} Computed as in \insertCite{zhang1994;textual}{aces}.
#' \item{\code{minspan = -1}} Computed as in \insertCite{friedman1991;textual}{aces}.
#' \item{\code{minspan = +m}}
#' }
#' @param endspan Minimum number of observations before the first and after the final knot.
#' \itemize{
#' \item{\code{endspan = -2}} Computed as in \insertCite{zhang1994;textual}{aces}.
#' \item{\code{endspan = -1}} Computed as in \insertCite{friedman1991;textual}{aces}.
#' \item{\code{endspan = +m}}
#' }
#' @param nX Number of inputs.
#'
#' @return A \code{list} with the minimum span and the end span.

compute_span <- function(data, minspan, endspan, nX) {

  # Sample size
  N <- nrow(data)

  # Minimum span (L)
  if (minspan == - 2) {
    # Zhang approach
    L <- c()

    for (var in 1:nX) {
      # Top 3 values
      max3 <- data[order(data[, var], decreasing = TRUE), var][1:3]
      # Bottom 3 values
      min3 <- data[order(data[, var]), var][1:3]

      m1 <- - (max3[1] - min3[1]) / (2.5 * (N - 1)) * log2(- (1 / N) * log(0.95))
      m2 <- (1 / N) * sum(max3 - min3)

      # Lvar limited for the 10% of the DMUs
      Lvar <- floor(min(N * 0.10, max(m1, m2)))

      L <- c(L, Lvar)
    }

  } else if (minspan == - 1) {
    # Friedman approach (this value must be computed later)
    L <- - 1

  } else {
    L <- min(N * 0.10, minspan)
  }

  # End span (Le)
  if (endspan == - 2) {
    # Zhang approach
    Le <- c()

    for (var in 1:nX) {
      max3 <- data[order(data[, var], decreasing = TRUE), var][1:3]
      min3 <- data[order(data[, var]), var][1:3]

      m1 <- - (max3[1] - min3[1]) / (2.5 * (N - 1)) * log2(- (1 / N) * log(0.95))
      m2 <- (1 / N) * sum(max3 - min3)

      Levar <- floor(min(N * 0.10, max(m1, m2)))

      Le <- c(Le, Levar)
    }

  } else if (endspan == - 1) {
    # Friedman approach
    Le <- floor(min(N * 0.1, 3 - log2(0.05 / nX)))

  } else {
    Le <- min(N * 0.1, endspan)

  }

  return(list(L, Le))
}

#' @title Set the Grid of Knots
#'
#' @description This function sets the grid of knots to perform ACES.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param nX Number of inputs (with interactions and netputs).
#' @param inps Number of original inputs (without interactions).
#' @param kn_grid Grid of knots to perform ACES.
#'
#' @return A \code{list} with the available vector of knots for each variable.

set_kn_grid <- function(data, nX, inps, kn_grid) {

  # Case 1: kn_grid is provided (list) and new variables are created (nX > inputs):
    # expand the kn_grid list.
  # Case 2: kn_grid is provided (list) and new variables are not created (nX = inputs):
    # keep the same kn_grid list.
  # Case 3: kn_grid is not provided:
    # create the kn_grid list.

  if (is.list(kn_grid)) {
    if (nX > inps) {
      # New variables (through interactions and netputs)
      NewVars <- nX - inps

      for (iv in 1:NewVars) {
        # variable index
        varIndx <- nX - NewVars + iv
        # variable name
        varName <- colnames(data)[varIndx]
        # variable data
        varData <- data[, varIndx]

        # length of the maximum grid
        max_len_grid <- max(sapply(kn_grid, length))

        # grid of knots for the VarName
        kn_grid[[varName]] <- seq (
          from = min(varData),
          to = max(varData),
          length.out = max_len_grid
          )
      }

    } else {
      kn_grid <- kn_grid
    }

  } else {
    kn_grid <- lapply(1:nX, function(i) data[, i])
  }

  return(kn_grid)
}
