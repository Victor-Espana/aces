#' @title Algorithm of Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).
#'
#' @description
#' This function implements the Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) algorithm, which estimates a production frontier satisfying classical production theory axioms like monotonicity and concavity. It offers both stochastic and envelopment versions. These estimations are based on the adaptation of the Multivariate Adaptive Regression Splines (MARS) technique developed by \insertCite{friedman1991;textual}{aces} and the Random Forest methodology \insertCite{breiman2001}{aces}. For details, see \insertCite{espana2024;textual}{aces}.
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
#' @param nvars
#' An \code{integer} indicating the number of variables randomly chosen at each split.
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
#' @param wc
#' Hyperparameter for side knot distances in the cubic smoothing procedure.
#'
#' @param wq
#' Hyperparameter for the side knot distances in the quintic smoothing procedure.
#'
#' @references
#'
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{breiman2001}{aces}
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc
#'
#' @return
#'
#' An \code{aces} object.
#'
#' @export

rf_aces_algorithm <- function (
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
    nvars,
    err_red,
    hd_cost,
    minspan,
    endspan,
    kn_grid,
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

  # original set of DMUs
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

  # samples size
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
  rf_aces <- list (
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

  while(length(rf_aces[["bf_set"]]) + 2 < nterms) {

    # negative GRSq
    last_bf <- rf_aces[["bf_set"]][[length(rf_aces[["bf_set"]])]]
    if (last_bf[["GRSq"]] < 0) break

    while (last_bf[["id"]] == 1) {

      # random input for bagging
      rnd_x <- sample(x, nvars)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function (
        data = data,
        x = rnd_x,
        y = y,
        xi_degree = xi_degree,
        dea_eff = dea_eff,
        model_type = model_type,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        L = L,
        Le = Le,
        err_min = err,
        hd_cost = hd_cost
      )

      if (is.list(B_bf_knt_err)) break

    }

    if (last_bf[["id"]] != 1) {

      # random input for bagging
      rnd_x <- sample(x, nvars)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function (
        data = data,
        x = rnd_x,
        y = y,
        xi_degree = xi_degree,
        dea_eff = dea_eff,
        model_type = model_type,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        L = L,
        Le = Le,
        err_min = err,
        hd_cost = hd_cost
      )

    }

    if (!is.list(B_bf_knt_err)) break

    # new best error
    new_err <- B_bf_knt_err[[5]]

    # update model
    if (last_bf[["id"]] == 1 || new_err[1] < err[1] * (1 - err_red[1])) {

        # update B
        rf_aces[["B"]] <- B_bf_knt_err[[1]]

        # update basis functions
        rf_aces[["bf_set"]] <- B_bf_knt_err[[2]]

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
  # rf-aces
  # ==

  rf_aces = list (
    "basis" = rf_aces[["bf_set"]],
    "Bmatx" = rf_aces[["B"]],
    "knots" = kn_forward,
    "coefs" = rev(rf_aces[["bf_set"]])[[1]][["coefs"]]
  )

  # ======================= #
  #   SMOOTHING PROCEDURE   #
  # ======================= #

  # select a model to be smoothed
  rf_aces_smoothed <- rf_aces

  # data.frame of knots
  kn_smoothed <- rf_aces_smoothed[["knots"]]
  kn_smoothed$side <- "R"

  kn_smoothed_left <- kn_smoothed
  kn_smoothed_left$side <- "L"

  kn_smoothed <- rbind (
    kn_smoothed,
    kn_smoothed_left
    )

  # generate the input space for side knots location
  kn_side_loc <- side_knot_location (
    data = data,
    nX = nX,
    knots = kn_smoothed
  )

  # ==
  # smoothing cubic aces
  # ==

  aces_cubic <- cubic_aces (
    data = data,
    x = x,
    y = y,
    dea_eff = dea_eff,
    model_type = model_type,
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    d = 0,
    wc = wc
  )

  aces_quintic <- quintic_aces (
    data = data,
    x = x,
    y = y,
    dea_eff = dea_eff,
    model_type = model_type,
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    d = 0,
    wq = wq
  )

  # =========== #
  # ACES OBJECT #
  # =========== #

  RF_ACES <- aces_object (
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
    nterms = ncol(rf_aces[["Bmatx"]]),
    err_red = err_red,
    hd_cost = hd_cost,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    d = NULL,
    wc = wc,
    wq = wq,
    aces_forward = rf_aces,
    aces = NULL,
    aces_cubic = aces_cubic,
    aces_quintic = aces_quintic
  )

  return(RF_ACES)

}

#' @title Compute Out-of-Bag Error in Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).
#'
#' @description
#'
#' This function calculates the Out-Of-Bag (OOB) error for a given Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) with a certain number of base learners.
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
#' @param object
#' An \code{aces} that has been trained using bagging as in Random Forest.
#'
#' @param oob_data
#' A \code{data.frame} containing the out-of-bag samples. These are the sample that were not used during the training of the ACES model and will be used to evaluate its performance.
#'
#' @param inb_idxs
#' A \code{vector} of integers representing the indices of the samples that are in the bag.
#'
#' @param oob_pred
#' A \code{list} with the Out-Of-Bag predictions for each model.
#'
#' @param model
#' Index indicating the model of the ensamble.
#'
#' @param method
#' Model for prediction:
#' \itemize{
#' \item{\code{"aces_forward"}}: Random Forest Forward Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_cubic"}}: Random Forest Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_quintic"}}: Random Forest Quintic Smoothed Adaptive Constrained Enveloping Splines model.
#' }
#'
#' @return
#'
#' This function returns the Out-of-Bag error metric, which provides an estimate of the model's prediction error on unseen data.

compute_oob <- function (
    data,
    y,
    oob_idxs,
    oob_pred,
    model
    ) {

  # number of outputs
  nY <- length(y)

  # sample size
  N <- nrow(data)

  # combine predictions from all models up to the current model
  model_pred <- as.data.frame(do.call(rbind, oob_pred[1:model]))

  # convert the predictions to a data frame and split by indices
  model_pred <- lapply (
    split (
      model_pred[, 1:nY, drop = F],
      model_pred[, ncol(model_pred)]
      ),
    function(t) apply(t, 2, mean)
    )

  # calculate the mean prediction for each index
  model_pred <- do.call(rbind, model_pred)

  # calculate the OOB mean squared error
  oob_idxs <- as.numeric(rownames(model_pred))
  oob_mse <- sum((data[oob_idxs, y] - model_pred[, 1:nY, drop = F]) ^ 2) / (length(oob_idxs) * nY)

  return(oob_mse)

}

#' @title Compute Variable Importance in Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).
#'
#' @description
#'
#' This function computes a metric of variable importance for a given Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) with a certain number of base learners.
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
#' @param object
#' A \code{rf_aces} object.
#'
#' @param repeats
#' Number of times the variable importance procedure is repeated.
#'
#' @return
#'
#' This function returns a metric of variable importance for each input variable

compute_variable_importance <- function (
    data,
    x,
    y,
    object,
    repeats = 1
    ) {

  # ACES meta-data
  control_features <- object[[1]][[1]][["control"]]

  model_type <- control_features[["model_type"]]
  error_type <- control_features[["error_type"]]
  y_type <- object[[1]][[1]][["data"]][["y_type"]]
  degree <- control_features[["degree"]]
  metric <- control_features[["metric"]]
  mon <- control_features[["monotonicity"]]
  con <- control_features[["concavity"]]
  ori <- as.logical(control_features[["origin"]])
  err_red <- control_features[["err_red"]]
  hd_cost <- control_features[["hd_cost"]]
  minspan <- control_features[["minspan"]]
  endspan <- control_features[["endspan"]]
  wc <- control_features[["wc"]]
  wq <- control_features[["wq"]]

  # RF-ACES meta-data
  RF_features <- object[[1]][[1]][["control"]][["RF"]]

  RF_sample <- RF_features[["sample"]]
  RF_models <- RF_features[["models"]]
  RF_nvars <- RF_features[["nvars"]]
  RF_oob_red <- RF_features[["oob_red"]]

  # number of estimated technologies
  nrankings <- length(object[[1]])

  # ranking
  rankings <- vector("list", nrankings)

  # sample size
  N <- nrow(data)

  for (k in 1:nrankings) {

    # select a model
    data_info <- object[[1]][[k]][["data"]]

    # look for netputs if y_type = "individual"
    z <- data_info[["z"]]

    # output index
    y <- data_info[["y"]]

    # data in [x, z, y] format with interaction of variables included
    data_model <- set_interactions (
      data = data,
      x = x,
      y = NULL,
      z = z,
      degree = degree
    )

    # add output data to input data
    data_model <- cbind(data_model, data[, y, drop = FALSE])

    # reorder index 'x' and 'y' in data
    x <- 1:(ncol(data_model) - length(y))
    y <- (length(x) + 1):ncol(data_model)

    # number of inputs
    nX <- length(x)

    # initialize matrix of variable importance
    mat_varimp <- matrix (
      0,
      nrow = nX
    )
    rownames(mat_varimp) <- colnames(data_model)[x]
    colnames(mat_varimp) <- "importance"

    # out-of-bag error of the complete model
    oob <- object[[length(object)]][[1]][["OOB"]]

    # out-of-bag excluding the j-th input variable
    oob_j <- c()

    # compute oob shuffling each variable "j"
    for (j in 1:nX) {
      for (n in 1:repeats) {
        data_shuffled <- data_model
        xvar_shuffled <- data_model[sample(1:N), x[j]]
        data_shuffled[, x[j]] <- xvar_shuffled

        model_j <- aces (
          data = data_shuffled,
          x = x,
          y = y,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          RF = list (
            "apply" = TRUE,
            "sample" = RF_sample,
            "models" = RF_models,
            "nvars" = RF_nvars,
            "oob_red" = RF_oob_red
          ),
          mul_BF = list (
            "degree" = degree,
            "hd_cost" = 0
          ),
          metric = "mse",
          shape = list (
            "mon" = mon,
            "con" = con,
            "ori" = ori
          ),
          nterms = nrow(data),
          err_red = err_red,
          kn_grid = - 1,
          minspan = minspan,
          endspan = endspan,
          kn_penalty = NULL,
          smoothing = list (
            "wc" = wc,
            "wq" = wq
          )
        )

        # out-of-bag error for the model with the jth variable shuffled
        oob_j <- c(oob_j, model_j[[length(model_j)]][[1]][["OOB"]])
      }

      # add oob_j to the matrix of variable importance
      mat_varimp[j, "importance"] <- round(100 * ((mean(oob_j) - oob) / mean(oob_j)), 2)
    }

    rankings[[k]] <- mat_varimp

    return(rankings)
}








}
