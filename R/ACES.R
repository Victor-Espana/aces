#' @title Fit an Adaptive Constrained Enveloping Splines (ACES) Model
#'
#' @description
#'
#' This function estimates a production frontier that satisfies classical production theory axioms, such as monotonicity and concavity. Both stochastic and deterministic versions are available. The estimations are based on the adaptation of the Multivariate Adaptive Regression Splines (MARS) technique developed by \insertCite{friedman1991;textual}{aces}. An adaptation of Random Forest \insertCite{breiman2001}{aces} is also included. For details, see \insertCite{espana2024;textual}{aces}
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
#' A \code{character} string specifying the nature of the production frontier to estimate. Options are:
#' \itemize{
#'   \item{\code{"ind"}}: Predict one output at a time, using the rest of the outputs as additional inputs (netputs). Netputs are treated as inputs for prediction and as outputs for efficiency estimation.
#'   \item{\code{"all"}}: Predict all the outputs simultaneously with the original set of inputs.
#' }
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier to estimate. Options are:
#' \itemize{
#'   \item{\code{"env"}}: Fit an enveloping production frontier.
#'   \item{\code{"sto"}}: Fit a stochastic production frontier.
#' }
#'
#' @param error_type
#' A \code{character} string specifying the error structure to use. Options are:
#' \itemize{
#'   \item{\code{"add"}}: Additive error structure.
#'   \item{\code{"mul"}}: Multiplicative error structure.
#' }
#'
#' @param RF
#' A \code{list} indicating if a bootstrap aggregation methodology like Random Forest should be used. Items include:
#' \itemize{
#'   \item{\code{apply}}: A \code{logical} indicating if bagging should be applied.
#'   \item{\code{sample}}: A \code{numeric} indicating the sample size for bagging.
#'   \item{\code{models}}: A \code{numeric} indicating the number of models for bagging.
#'   \item{\code{nvars}}: An \code{integer} indicating the number of variables randomly chosen at each split.
#'   \item{\code{oob_red}}: A \code{numeric} specifying the minimum improvement ratio in the moving average of the out-of-bag error over a period of 10 compared to the previous period for the addition of a new model to the ensamble. Default is \code{0.001}.
#' }
#'
#' @param mul_BF
#' A \code{list} specifying the maximum degree of basis functions (BFs) and the cost of introducing a multivariate BF. Items include:
#' \itemize{
#' \item{\code{degree}}: A \code{list} with input indexes for interaction of variables, or a \code{numeric} specifying the maximum degree of interaction. Basis functions products are constrained to contain factors involving distinct variables to ensure interpretability and avoid multicollinearity.
#' \item{\code{hd_cost}}:  A \code{numeric} specifying the minimum percentage improvement over the best 1-degree BF to incorporate a higher degree BF. Default is \code{0.05}.
#' }
#'
#' @param metric
#' A \code{list} specifying the lack-of-fit criterion to evaluate the model performance. Options are:
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
#' A \code{list} specifying shape constraints for the estimator. Items include:
#' \itemize{
#'   \item{\code{mono}}: A \code{logical} indicating if non-decreasing monotonicity should be enforced.
#'   \item{\code{conc}}: A \code{logical} indicating if concavity should be enforced.
#'   \item{\code{ptto}}: A \code{logical} indicating if the estimator should satisfy \code{f(0) = 0}.
#' }
#'
#' @param nterms
#' A positive \code{integer} specifying the maximum number of terms created before pruning. Default is \code{50}.
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree basis functions. Default is \code{0.01}.
#'
#' @param kn_grid
#' Either a \code{-1} (default) to use the original approach of \insertCite{friedman1991;textual}{aces} based on the observed data, or a a \code{list} with the grid of knots for performing ACES. Each element of the \code{list} contains a vector with the knots of one variable (e.g., the first element of the list contains the knot values of the first variable and so on).
#'
#' @param minspan
#' A \code{numeric} specifying the minimum number of observations between two adjacent knots. Options are:
#' \itemize{
#' \item{\code{minspan = -2}}: Computed as in \insertCite{zhang1994;textual}{aces}.
#' \item{\code{minspan = -1}}: Computed as in \insertCite{friedman1991;textual}{aces}.
#' \item{\code{minspan = +m}}: A user-specified positive integer.
#' }
#'
#' @param endspan
#' A \code{numeric} specifying the minimum number of observations before the first and after the final knot. Options are:
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
#' Let `p` be the distance between the central knot and the right-side knot, and `v` be the distance between the central knot and the left-side knot during the smoothing procedure. A \code{list} specifying conditions for the smoothing procedure. Items include:
#' \itemize{
#' \item{\code{wc}}:  A \code{numeric} value used for cubic smoothing \insertCite{friedman1991}{aces}. This parameter is defined as `v / p` and must be set between 1 and 2. If a \code{vector} is entered, the \code{wc} value that most reduced the lack-of-fit criterion is selected.
#' \item{\code{wq}}: A \code{numeric} value used for quintic smoothing \insertCite{chen1999}{aces}. This parameter is defined as `p / v` and must be set between 8/7 and 1.5. If a \code{vector} is entered, the \code{wc} value that most reduced the lack-of-fit criterion is selected.
#' }
#'
#' @details
#'
#' This function generates a production frontier adhering to classical production theory axioms, such as monotonicity and concavity. Users can choose between enveloping or stochastic versions of the production functions. The algorithm comprises three main procedures:
#'
#' 1. A forward selection algorithm that creates a ser of linear basis functions, which may initially overfit the training data.
#'
#' 2. A backward elimination algorithm that remove basis that do not significantly contribute to the model's performance.
#'
#' 3. Two smoothing procedures are available: one using cubic functions and another using quintic functions.
#'
#' If the stochastic version is selected, the frontier's shape is estimated without imposing enveloping constraints on the observations. In the second stage, the expected value of inefficiency is estimated using the residuals obtained from the first stage.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{breiman2001}{aces} \cr \cr
#' \insertRef{zhang1994}{aces} \cr \cr
#' \insertRef{chen1999}{aces}
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
      "mono" = TRUE,
      "conc" = TRUE,
      "ptto" = FALSE
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
    convexity = NULL,
    returns = NULL,
    direction = NULL,
    digits = NULL
  )

  # adapt "pass through origin" hyperparameter based on error_type
  if (shape[["ptto"]] && error_type == "add") {
    shape[["ptto"]] <- "0"

  } else if (shape[["ptto"]] && error_type == "mul") {
    shape[["ptto"]] <- "1"

  }

  # list with individual ACES models
  ACES <- list()

  if (!RF[["apply"]]) {

    if (y_type == "all") {

      ACES[[1]] <- aces_algorithm (
        data = data,
        inps = x,
        outs = y,
        nets = NULL,
        y_type = y_type,
        model_type = model_type,
        error_type = error_type,
        degree = mul_BF[["degree"]],
        metric = metric,
        shape = list (
          "mono" = shape[["mono"]],
          "conc" = shape[["conc"]],
          "ptto" = shape[["ptto"]]
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
          inps = x,
          outs = outs,
          nets = nets,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          degree = mul_BF[["degree"]],
          metric = metric,
          shape = list (
            "mono" = shape[["mono"]],
            "conc" = shape[["conc"]],
            "ptto" = shape[["ptto"]]
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

    # Progress Bar
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

      # update Progress Bar
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
          inps = x,
          outs = y,
          nets = NULL,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          degree = mul_BF[["degree"]],
          metric = metric,
          shape = list (
            "mono" = shape[["mono"]],
            "conc" = shape[["conc"]],
            "ptto" = shape[["ptto"]]
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

        # early stopping RF-ACES based on moving average
        oob_vec <- c(oob_vec, RF_ACES[[m]][["y_all"]][["OOB"]])
        mov_avg_oob_05 <- c(mov_avg_oob_05, mean(tail(oob_vec, 05)))
        mov_avg_oob_10 <- c(mov_avg_oob_10, mean(tail(oob_vec, 10)))
        mov_avg_oob_25 <- c(mov_avg_oob_25, mean(tail(oob_vec, 25)))

        if (m > min(nrow(data), 200) && oob_vec[m] >= oob_vec[m - 1]) {

          stopping_condition_counter <- ifelse (
            oob_vec[m] >= mov_avg_oob_05[m] * (1 - RF[["oob_red"]]) &
            oob_vec[m] >= mov_avg_oob_10[m] * (1 - RF[["oob_red"]]) &
            oob_vec[m] >= mov_avg_oob_25[m] * (1 - RF[["oob_red"]]),
            stopping_condition_counter + 1,
            0
          )

        if (stopping_condition_counter == 5) {

          cat("\n")
          print(paste("Out-of-bag error stabilized after", m, "models"))
          break

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
          inps = x,
          outs = outs,
          nets = nets,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          degree = mul_BF[["degree"]],
          metric = metric,
          shape = list (
            "mono" = shape[["mono"]],
            "conc" = shape[["conc"]],
            "ptto" = shape[["ptto"]]
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

      # early stopping RF-ACES based on moving average
      oob_vec <- c(oob_vec, RF_ACES[[m]][["y1"]][["OOB"]])
      mov_avg_oob_05 <- c(mov_avg_oob_05, mean(tail(oob_vec, 05)))
      mov_avg_oob_10 <- c(mov_avg_oob_10, mean(tail(oob_vec, 10)))
      mov_avg_oob_25 <- c(mov_avg_oob_25, mean(tail(oob_vec, 25)))

      if (m > min(nrow(data), 200) && oob_vec[m] >= oob_vec[m - 1]) {

        stopping_condition_counter <- ifelse (
          oob_vec[m] >= mov_avg_oob_05[m] * (1 - RF[["oob_red"]]) &
            oob_vec[m] >= mov_avg_oob_10[m] * (1 - RF[["oob_red"]]) &
            oob_vec[m] >= mov_avg_oob_25[m] * (1 - RF[["oob_red"]]),
          stopping_condition_counter + 1,
          0
        )

        if (stopping_condition_counter == 5) {

          cat("\n")
          print(paste("Out-of-bag error stabilized after", m, "models"))
          break

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
  class(ACES) <- ifelse(RF[["apply"]], "rf_aces", "aces")

  return(ACES)

}

#' @title Algorithm of Adaptive Constrained Enveloping Splines (ACES).
#'
#' @description
#'This function implements the Adaptive Constrained Enveloping Splines (ACES) algorithm, which estimates a production frontier satisfying classical production theory axioms like monotonicity and concavity. It offers both stochastic and envelopment versions. These estimations are based on the adaptation of the Multivariate Adaptive Regression Splines (MARS) technique developed by \insertCite{friedman1991;textual}{aces}. For details, see \insertCite{espana2024;textual}{aces}
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the variables in the model.
#'
#' @param inps
#' Column indexes of input variables in \code{data}.
#'
#' @param outs
#' Column indexes of output variables in \code{data}.
#'
#' @param nets
#' Column indexes of netput (outputs evaluated as inputs) variables in \code{data}. These variables are treated as inputs during prediction computation and as outputs when computing efficiency scores.
#'
#' @param y_type
#' A \code{character} string that determines the prediction approach for \code{y}.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the production frontier that the function estimates.
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
#' Grid design for knots placement in ACES.
#'
#' @param d
#' Penalty per knot for Generalized Cross Validation.
#'
#' @param wc
#' Hyperparameter for side knot distances in the cubic smoothing procedure.
#'
#' @param wq
#' Hyperparameter for the side knot distances in the quintic smoothing procedure.
#'
#' @references
#'
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{espana2024}{aces} \cr \cr
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
    inps,
    outs,
    nets,
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

  # accelerate the algorithm by selecting the efficient DMUs in DEA to impose only
  # the envelope on these points.
  # only valid under monotonicity and concavity constraints

  if (shape[["mono"]] && shape[["conc"]]) {

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

  # save a copy of the original data
  dmus <- data

  # data in [x, z, y] format with interaction and / or transformation of
  # variables included
  data <- prepare_data (
    data = data,
    x = inps,
    y = outs,
    z = nets,
    degree = degree,
    error_type = error_type
    )

  # samples in data
  N <- nrow(data)

  # reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - length(outs))
  y <- (length(x) + 1):ncol(data)

  # number of inputs / outputs as inputs
  nX <- length(x)

  # number of outputs
  nY <- length(y)

  # matrix with:
  # row 1: the index of the variable
  # row 2: the degree of the variable (netput = 1)
  xi_degree <- matrix (
    c(x, rep(1, length(x))),
    byrow = TRUE,
    nrow = 2,
    ncol = length(x)
    )

  if (!is.list(degree)) {

    v <- 0

    for (i in 1:degree) {

      combs <- combn(1:length(inps), i)

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

  } else if (model_type == "sto") {
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
  kn_grid <- set_knots_grid (
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
      shape = shape,
      kn_list = kn_list,
      kn_grid = kn_grid,
      L = L,
      Le = Le,
      err_min = err,
      hd_cost = hd_cost
      )

    if (!is.list(B_bf_knt_err)) break

    # new best error
    new_err <- B_bf_knt_err[[5]]

    # update model
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
      shape = shape,
      d = d
    )

  # generalized cross-validation for each model
  GCVs <- sapply(aces_submodels, function(x) x[["GCV"]])

  # model with minimum error (excluding the model without knots)
  aces_backward <- aces_submodels[[which.min(GCVs[1:(length(aces_submodels) - 1)])]]

  # set of surviving knots
  kn_backward <- do.call(rbind.data.frame, aces_backward[["t"]])

  # sort the knots by "xi", "status", "side", "t"
  knots_backward_order <- with (
    kn_backward, {
    order(xi, status, ifelse(status == "paired", t, desc(side)), ifelse(status == "paired", desc(side), t))
    }
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

    # transform to data.frame
    kn_smoothed <- do.call(rbind.data.frame, aces_smoothed[["t"]])

    # if monotonicity is required:
    # 1- wc in (1, 2) and wq in (8/7, 1.5)

    # If concavity is required:
    # 1- wc in (1, 2) and wq in (8/7, 1.5)
    # 2- unpaired right basis functions are not allowed

    # check for right-side unpaired basis functions
    check1 <- kn_smoothed$side == "R"
    check2 <- kn_smoothed$status == "unpaired"

    if (shape[["conc"]] && max(check1 + check2) == 2) {

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
        kn_side_loc <- side_knot_location (
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
          shape = shape,
          kn_grid = kn_smoothed,
          kn_side_loc = kn_side_loc,
          d = d,
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
          shape = shape,
          kn_grid = kn_smoothed,
          kn_side_loc = kn_side_loc,
          d = d,
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
      shape = shape,
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

#' @title Create an aces Object
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
#' A \code{character} string specifying the nature of the production frontier that the function estimates.
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
    shape,
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
    "shape" = shape,
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

#' @title Prepare Data for Fitting
#'
#' @description
#' This function prepares the data for model fitting by generating additional input variables through interactions between variables. It also performs any necessary transformations, such as changing to a logarithmic scale if the error type is multiplicative. It returns a matrix in [x, z, y] format, where x represents input variables, z represents netput variables, and y represents output variables.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param z
#' Column indexes of netput variables in \code{data}. These variables are not considered for interaction with other variables.
#'
#' @param degree
#'  Maximum degree of interaction between variables. It can be a \code{list} of input indexes for interactions or a \code{numeric} value determining the maximum degree of interaction.
#'
#' @param error_type
#' A \code{character} string specifying the error structure when fitting the model.
#'
#' @return
#' A \code{matrix} in a [x, z, y] format with variable interactions and / or transformations included.

prepare_data <- function (
    data,
    x,
    y,
    z,
    degree,
    error_type
    ) {

  # 1. change to logarithmic scale if the error_type is multiplicative
  if (error_type == "mul") {
    data[, c(y)] <- log(data[, c(y)])
  }

  # 2. transform the output(s) as input(s) through the opposite: netput(s)
  data[, z] <- - data[, z]

  # 3. generate interaction effects
  if (is.list(degree) || degree > 1) {

    if (!is.list(degree)) {

      # create a list with all the possible combinations between 1 and as much
      # len(x) elements
      max_degree <- degree
      degree <- list()

      for (i in 2:max_degree) {

        combs <- combn(1:length(x), i)

        for (col in 1:ncol(combs)) {
          degree <- append(degree, list(combs[, col]))
        }

      }
    }

    # number of additional variables
    IVars <- length(degree)

    # new x indexes
    new_x <- c(x, (ncol(data) + 1):(ncol(data) + IVars))

    # create the new variables
    for (p in 1:IVars) {

      # select the variables
      vars <- new_x[degree[[p]]]

      # name the new variable
      name_vars <- colnames(data)[vars]
      name <- paste(name_vars, collapse = "_")

      # create the new variable
      data[, name] <- apply(data[, vars], 1, prod)

    }

  } else {

    new_x <- x

  }

  # 4. data correctly sorted
  data <- data[, c(new_x, z, y)]

  return(as.matrix(data))

}

#' @title Error Metric for Model Evaluation.
#'
#' @description
#' Computes an error metric for model evaluation based on observed and predicted values.
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

    # compute the mean of column-wise maximums and minimums in y
    ymax <- mean(apply(y_obs, 2, max))
    ymin <- mean(apply(y_obs, 2, min))

    # normalized root mean squared error by the range
    devtn <- (y_hat - y_obs) ^ 2
    error <- sqrt(sum(devtn) / (N * nY)) / (ymax - ymin)
  }

  return(error)
}

#' @title Compute Minimum and End Span
#'
#' @description
#' This function computes the minimum span, which is the minimum number of observations between two adjacent knots, and the end span, which is the minimum number of observations before the first knot and after the final knot.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param minspan
#' A \code{numeric} value or vector specifying the minimum number of observations between two adjacent knots. The following options are available:
#' \itemize{
#'   \item{\code{minspan = -2}}: Computed according to the method proposed by \insertCite{zhang1994;textual}{aces}.
#'   \item{\code{minspan = -1}}: Computed according to the method proposed by \insertCite{friedman1991;textual}{aces}.
#'   \item{\code{minspan = +m}}: A positive integer specifying the exact number of observations.
#' }
#'
#' @param endspan
#' A \code{numeric} value or vector specifying the minimum number of observations before the first knot and after the final knot. The following options are available:
#' \itemize{
#'   \item{\code{endspan = -2}}: Computed according to the method proposed by \insertCite{zhang1994;textual}{aces}.
#'   \item{\code{endspan = -1}}: Computed according to the method proposed by \insertCite{friedman1991;textual}{aces}.
#'   \item{\code{endspan = +m}}: A positive integer specifying the exact number of observations.
#' }
#'
#' @param nX
#' Number of input variables.
#'
#' @references
#'
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{zhang1994}{aces}
#'
#' @return
#' A \code{list} with two components:
#' \itemize{
#'   \item{\code{minspan}}: The computed minimum span.
#'   \item{\code{endspan}}: The computed end span.
#' }

compute_span <- function (
    data,
    minspan,
    endspan,
    nX
    ) {

  # sample size
  N <- nrow(data)

  # minimum span (L)
  if (minspan == - 2) { # Zhang approach

    L <- numeric(nX)

    # fixed log_factor
    log_factor <- log2(- (1 / N) * log(0.95))

    for (var in 1:nX) {

      # sorted variable
      sorted_var <- sort(data[, var])

      # 3 highest values
      max3 <- tail(sorted_var, 3)

      # 3 lowest values
      min3 <- head(sorted_var, 3)

      m1 <- - (max3[1] - min3[1]) / (2.5 * (N - 1)) * log_factor
      m2 <- (1 / N) * sum(max3 - min3)

      # Lvar limited for the 10% of the DMUs
      L[var] <- floor(min(N * 0.10, max(m1, m2)))

    }

  } else if (minspan == - 1) { # Friedman approach (this value must be computed later)

    L <- - 1

  } else {

    L <- min(N * 0.10, minspan)

  }

  # end span (Le)
  if (endspan == - 2) { # Zhang approach

    Le <- numeric(nX)

    # fixed log_factor
    log_factor <- log2(- (1 / N) * log(0.95))

    for (var in 1:nX) {

      # sorted variable
      sorted_var <- sort(data[, var])

      # 3 highest values
      max3 <- tail(sorted_var, 3)

      # 3 lowest values
      min3 <- head(sorted_var, 3)

      m1 <- - (max3[1] - min3[1]) / (2.5 * (N - 1)) * log_factor
      m2 <- (1 / N) * sum(max3 - min3)

      Le[var] <- floor(min(N * 0.10, max(m1, m2)))

    }

  } else if (endspan == - 1) { # Friedman approach

    Le <- floor(min(N * 0.1, 3 - log2(0.05 / nX)))

  } else {

    Le <- min(N * 0.1, endspan)

  }

  return(list(L, Le))

}

#' @title Set the Grid of Knots
#'
#' @description
#' This function sets the grid of knots to perform Adaptive Concave Estimation for Stochastic Frontier (ACES).
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param nX
#' Number of inputs (including interactions and netputs).
#'
#' @param inps
#' Number of original inputs (excluding interactions).
#'
#' @param kn_grid
#' Grid of knots to perform ACES. If not provided, the function creates a grid of knots for each variable.
#'
#' @return
#' A \code{list} with the available vector of knots for each variable.

set_knots_grid <- function (
    data,
    nX,
    inps,
    kn_grid
    ) {

  # Case 1: kn_grid is provided (list) and new variables are created (nX > inputs):
    # expand the kn_grid list.

  # Case 2: kn_grid is provided (list) and new variables are not created (nX = inputs):
    # keep the same kn_grid list.

  # Case 3: kn_grid is not provided:
    # create the kn_grid list.

  if (is.list(kn_grid)) { # if kn_grid is provided

    if (nX > inps) {

      # Number of new variables (through interactions and netputs)
      new_vars <- nX - inps

      for (v in seq_len(new_vars)) {

        # variable index
        var_idx <- nX - new_vars + v

        # variable name
        var_name <- colnames(data)[var_idx]

        # variable data
        var_data <- data[, var_idx]

        # length of the maximum grid
        max_len_grid <- max(sapply(kn_grid, length))

        # grid of knots for the new variable
        kn_grid[[varName]] <- seq (
          from = min(var_data),
          to = max(var_data),
          length.out = max_len_grid
          )
      }

    }

  } else { # if kn_grid is not provided, create it

    kn_grid <- lapply(1:nX, function(i) data[, i])

    # names
    names(kn_grid) <- colnames(data)[1:nX]

  }

  return(kn_grid)

}
