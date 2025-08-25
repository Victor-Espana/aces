#' @title Fit a Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) model.
#'
#' @description
#'
#' This function estimates a deterministic production frontier that adheres to classical production theory axioms, such as monotonicity and concavity. The estimation is based on adaptations of the Multivariate Adaptive Regression Splines (MARS) technique, developed by \insertCite{friedman1991;textual}{aces} and Random Forest, introduced by \insertCite{breiman2001;textual}{aces}. For comprehensive details on the methodology and implementation, please refer to \insertCite{espana2024;textual}{aces} and \insertCite{espana2024rf;textual}{aces}.
#'
#' @name rf_aces
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
#' @param quick_aces
#' A \code{logical} indicating if the fast version of ACES should be employed.
#'
#' @param error_type
#' A \code{character} string specifying the error structure to use. Options are:
#' \itemize{
#'   \item{\code{"add"}}: Additive error structure.
#'   \item{\code{"mul"}}: Multiplicative error structure.
#' }
#'
#' @param mul_BF
#' A \code{list} specifying the maximum degree of basis functions (BFs) and the cost of introducing a higher-degree BF. Items include:
#' \itemize{
#' \item{\code{max_degree}}: A \code{list} with input indexes for interaction of variables, or a \code{numeric} specifying the maximum degree of interaction. BFs products are constrained to contain factors involving distinct variables to ensure interpretability and avoid multicollinearity.
#' \item{\code{inter_cost}}: A \code{numeric} specifying the minimum percentage improvement over the best 1-degree BF to incorporate a higher degree BF. Default is \code{0.05}.
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
#' @param learners
#' An \code{integer} indicating the number of models for bagging.
#'
#' @param bag_size
#' An \code{integer} indicating the number of samples to draw from \code{data} to train each base estimator.
#'
#' @param max_feats
#' An \code{integer} indicating the number of variables randomly chosen at each split in RF-ACES.
#'
#' @param early_stopping
#' A \code{list} specifying the early stopping criteria based on the out-of-bag (OOB) error to unnecessary training by stopping when the OOB error no longer improves. The list should contain the following elements:
#' \itemize{
#'   \item{\code{ma_window}}: An \code{integer} specifying the size of the moving average window used to smooth the OOB error. A larger value makes the stopping condition more stable but less reactive to short-term fluctuations. Default is 10.
#'   \item{\code{tolerance}}: An \code{integer} indicating the number of consecutive iterations without significant improvement in the OOB error before stopping the training. Higher values make the model more tolerant to small fluctuations. Default is 5.
#' }
#'
#' @param max_terms
#' A positive \code{integer} specifying the maximum number of terms created during the forward step. Default is \code{50}.
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree basis functions. Default is \code{0.01}.
#'
#' @param kn_grid
#' Specifies the grid of knots to be used in the ACES algorithm. Accepts two options:
#' \itemize{
#'    \item{\code{-1}} (default): Uses the original approach by \insertCite{friedman1991;textual}{aces}, where the knots are automatically selected based on the observed data.
#'    \item{\code{list}}: A user-defined grid of knots for each variable. This must be a \code{list} where each element corresponds to a vector of knots for a specific variable (e.g., the first element of the \code{list} contains the knot values for the first variable, the second element contains the knots for the second variable, and so on). This option allows for greater control over the knot placement when customizing the estimation process.
#' }
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
#' @references
#'
#' \insertRef{espana2024rf}{aces} \cr \cr
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{breiman2001}{aces} \cr \cr
#' \insertRef{zhang1994}{aces} \cr
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc group_by across all_of summarise
#' @importFrom extraDistr rsign
#'
#' @return
#'
#' A \code{rf_aces} object.
#'
#' @export

rf_aces <- function (
    data,
    x,
    y,
    quick_aces = TRUE,
    error_type = "add",
    mul_BF = list (
      "max_degree" = 1,
      "inter_cost" = 0.05
    ),
    metric = "mse",
    shape = list (
      "mono" = TRUE,
      "conc" = TRUE,
      "ptto" = FALSE
      ),
    learners = 100,
    bag_size = nrow(data),
    max_feats = length(x) / 3,
    early_stopping = list (
      "ma_window" = 10,
      "tolerance" = 5
    ),
    max_terms = 50,
    err_red = 0.01,
    kn_grid = - 1,
    minspan = - 1,
    endspan = - 1
    ) {

  # possible error messages:
  display_errors_rf_aces (
    data = data,
    x = x,
    y = y,
    quick_aces = quick_aces,
    error_type = error_type,
    max_degree = mul_BF[["max_degree"]],
    inter_cost = mul_BF[["inter_cost"]],
    metric = metric,
    learners = learners,
    bag_size = bag_size,
    max_feats = max_feats,
    max_terms = max_terms,
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid
  )

  if (shape[["ptto"]]) {
    data <- rbind(data, rep(1e-10, ncol(data)))
  }

  # RF-ACES models
  RF_ACES <- vector("list", learners)

  # Progress Bar
  pb <- txtProgressBar(min = 0, max = learners, style = 3)

  # out-of-bag predictions
  oob_pred <- vector("list", learners)

  # moving-average for the out-of-bag
  oob_trend <- c()
  ma_window <- early_stopping[["ma_window"]]
  tolerance <- early_stopping[["tolerance"]]

  stopping_condition_counter <- 0

  for (m in 1:length(RF_ACES)) {

    # update Progress Bar
    setTxtProgressBar(pb, m)

    # index of samples in the model
    sample_bag <- sample (
      1:nrow(data),
      size = bag_size,
      replace = TRUE
    )

    # data for the m-model
    data_bag <- data[sample_bag, ]

    # out of bag data for the m-model
    inb_idxs <- 1:nrow(data) %in% sample_bag
    oob_idxs <- which(!(1:nrow(data) %in% sample_bag))
    data_oob <- data[oob_idxs, ]

    RF_ACES[[m]] <- rf_aces_algorithm (
      data = data_bag,
      x_vars = x,
      y_vars = y,
      quick_aces = quick_aces,
      error_type = error_type,
      max_degree = mul_BF[["max_degree"]],
      inter_cost = mul_BF[["inter_cost"]],
      metric = metric,
      shape = list (
        "mono" = shape[["mono"]],
        "conc" = shape[["conc"]],
        "ptto" = shape[["ptto"]]
      ),
      max_feats = max_feats,
      max_terms = max_terms,
      err_red = err_red,
      minspan = minspan,
      endspan = endspan,
      kn_grid = kn_grid
    )

    # predictions
    Bmatx <- RF_ACES[[m]][["methods"]][["rf_aces"]][["Bmatx"]]
    coefs <- RF_ACES[[m]][["methods"]][["rf_aces"]][["coefs"]]
    y_hat <- Bmatx %*% coefs

    # select out-of-bag indices
    y_hat <- y_hat[oob_idxs, ]

    oob_pred[[m]] <- matrix (
      c (
        "y_hat" = y_hat,
        "idx" = oob_idxs
        ),
      ncol = length(y) + 1
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
    RF_ACES[[m]][["OOB"]] <- oob_mse

    # early stopping RF-ACES based on moving average
    oob_trend <- c(oob_trend, oob_mse)

    if (m > ma_window) {

      # compute moving average over the last `ma_window` iterations
      oob_ma <- mean(tail(oob_trend, ma_window))

      # current oob greater than moving average over last oob
      if (oob_mse >= oob_ma) {
        stopping_condition_counter <- stopping_condition_counter + 1

      } else {
        stopping_condition_counter <- 0

      }
    }

    if (stopping_condition_counter == tolerance) {

      cat("\n")
      print(paste("Out-of-bag error stabilized after", m, "models"))
      break
    }

  }

  # random-forest hyperparameters
  rf_list = list (
    learners = learners,
    bag_size = bag_size,
    max_feats = max_feats,
    early_stopping = early_stopping
  )

  RF_ACES <- list("forest" = RF_ACES[1:m])

  # ========== #
  # Technology #
  # ========== #

  xmat_rf_aces <- xmat_rf_aces_cubic <- xmat_rf_aces_quintic <- NULL
  ymat_rf_aces <- ymat_rf_aces_cubic <- ymat_rf_aces_quintic <- NULL

  # save technology for Random Forest
  for (i in seq_along(RF_ACES[["forest"]])) {

    xmat_rf_aces <- rbind (
      xmat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["xmat"]]
      )

    ymat_rf_aces <- rbind (
      ymat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["ymat"]]
    )

    xmat_rf_aces_cubic <- rbind (
      xmat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["xmat"]]
    )

    ymat_rf_aces_cubic <- rbind (
      ymat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["ymat"]]
    )

    xmat_rf_aces_quintic <- rbind (
      xmat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["xmat"]]
    )

    ymat_rf_aces_quintic <- rbind (
      ymat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["ymat"]]
    )

  }

  data_sorted <- as.data.frame (
    set_data (
      data = data,
      x = x,
      y = y,
      max_degree = 1,
      error_type = "add"
    )
  )

  inp_cols <- names(data_sorted)[1:length(x)]
  out_cols <- names(data_sorted)[(length(x) + 1):(length(x) + length(y))]

  colnames(xmat_rf_aces) <- colnames(xmat_rf_aces_cubic) <- colnames(xmat_rf_aces_quintic) <- inp_cols
  colnames(ymat_rf_aces) <- colnames(ymat_rf_aces_cubic) <- colnames(ymat_rf_aces_quintic) <- out_cols

  rf_aces_technology <- data.frame(xmat_rf_aces, ymat_rf_aces) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), mean, .names = "{col}"), .groups = 'drop')

  rf_aces_cubic_technology <- data.frame(xmat_rf_aces_cubic, ymat_rf_aces_cubic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), mean, .names = "{col}"), .groups = 'drop')

  rf_aces_quintic_technology <- data.frame(xmat_rf_aces_quintic, ymat_rf_aces_quintic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), mean, .names = "{col}"), .groups = 'drop')

  RF_ACES[["technology"]] <- list (
    "rf_aces" = list (
      "xmat" = as.matrix(rf_aces_technology[, x]),
      "ymat" = as.matrix(rf_aces_technology[, y])
    ),
    "rf_aces_cubic" = list (
      "xmat" = as.matrix(rf_aces_cubic_technology[, x]),
      "ymat" = as.matrix(rf_aces_cubic_technology[, y])
    ),
    "rf_aces_quintic" = list (
      "xmat" = as.matrix(rf_aces_quintic_technology[, x]),
      "ymat" = as.matrix(rf_aces_quintic_technology[, y])
    )
  )

  # ====== #
  # Object #
  # ====== #

  RF_ACES[["data"]] <- list (
    "df" = data,
    "x" = x,
    "y" = y,
    "xnames" = colnames(data)[x],
    "ynames" = colnames(data)[y],
    "rownames" = rownames(data)
  )

  RF_ACES[["control"]] <- list (
    "quick_aces" = quick_aces,
    "error_type" = error_type,
    "max_degree" = mul_BF[["max_degree"]],
    "inter_cost" = mul_BF[["inter_cost"]],
    "xi_degree" = RF_ACES[["models"]][[1]][["control"]][["xi_degree"]],
    "metric" = metric,
    "shape" = shape,
    "learners" = rf_list[["learners"]],
    "bag_size" = rf_list[["bag_size"]],
    "max_feats" = rf_list[["max_feats"]],
    "early_stopping" = rf_list[["early_stopping"]],
    "max_terms" = max_terms,
    "err_red" = err_red,
    "minspan" = minspan,
    "endspan" = endspan,
    "kn_grid" = kn_grid,
    "kn_penalty" = NULL,
    "wc" = NULL,
    "wq" = NULL
  )

  # change the order
  RF_ACES <- RF_ACES[c("data", "control", "forest", "technology")]

  # type of object
  class(RF_ACES) <- "rf_aces"

  return(RF_ACES)

}


#' @title Algorithm of Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).
#'
#' @description
#' This function implements the Random Forest Adaptive Constrained Enveloping Splines (ACES) algorithm.
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the variables in the model.
#'
#' @param x_vars
#' Column indexes of input variables in \code{data}.
#'
#' @param y_vars
#' Column indexes of output variables in \code{data}.
#'
#' @param quick_aces
#' A \code{logical} indicating if the fast version of ACES should be employed.
#'
#' @param error_type
#' A \code{character} string specifying the error structure when fitting the model.
#'
#' @param max_degree
#' Maximum degree of interaction between variables.
#'
#' @param inter_cost
#' Minimum percentage of improvement over the best 1 degree BF to incorporate a higher degree BF.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion employed to evaluate the model performance.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param max_feats
#' An \code{integer} indicating the number of variables randomly chosen at each split.
#'
#' @param max_terms
#' Maximum number of terms created before pruning.
#'
#' @param err_red
#' Minimum reduced error rate for the addition of a new pair of 1-degree basis functions.
#'
#' @param minspan
#' Minimum number of observations between two adjacent knots.
#'
#' @param endspan
#' Minimum number of observations before the first and after the final knot.
#'
#' @param kn_grid
#' Grid design for knots placement in RF-ACES.
#
#' @importFrom dplyr desc
#'
#' @return
#'
#' A \code{rf-aces} object.
#'
#' @export

rf_aces_algorithm <- function (
    data,
    x_vars,
    y_vars,
    quick_aces,
    error_type,
    max_degree,
    inter_cost,
    metric,
    shape,
    max_feats,
    max_terms,
    err_red,
    minspan,
    endspan,
    kn_grid
    ) {

  # save a copy of the original data
  DMUs <- data

  # data in [x, z, y] format with interaction of variables included
  data <- set_data (
    data = data,
    x = x_vars,
    y = y_vars,
    max_degree = max_degree,
    error_type = error_type
  )

  # samples size
  N <- nrow(data)

  # set "x", "z" and "y" indexes in data
  x <- 1:(ncol(data) - length(y_vars))
  y <- (length(x) + 1):ncol(data)

  # maximum features before and after interaction variables
  max_feats <- max_feats / length(x_vars)
  max_feats <- max_feats * length(x)

  # set number of inputs and outputs
  nX <- length(x)
  nY <- length(y)

  # variable importance
  var_imp <- matrix (
    rep(0, nX),
    nrow = 1
  )
  colnames(var_imp) <- colnames(data)[1:nX]

  x_drop <- c()

  # remove variables with low correlation
  if (quick_aces) {

    # Spearman’s Rank Correlation
    spearman_corr <- cor(data, method = "spearman")
    spearman_corr <- spearman_corr[1:length(x), (length(x) + 1):ncol(data)]

    # Kendall’s Tau
    kendall_corr <- cor(data, method = "kendall")
    kendall_corr <- kendall_corr[1:length(x), (length(x) + 1):ncol(data)]

    # compute threshold for removing variables
    t1_quantile <- quantile(spearman_corr[spearman_corr > 0], probs = 0.2)
    threshold_1 <- min(0.1, t1_quantile)

    t2_quantile <- quantile(kendall_corr[kendall_corr > 0], probs = 0.2)
    threshold_2 <- min (0.1, t2_quantile)

    # iterate over the variables
    for (j in 1:nX) {

      # check for both Spearman and Kendall correlations
      if (spearman_corr[j] < threshold_1 & kendall_corr[j] < threshold_2) {

        var_imp[1, j] <- - 1

      }
    }

    relevant_variables <- c()

    for (col in colnames(var_imp)) {

      if (var_imp[1, col] == 0) {

        variables <- unlist(strsplit(col, "_"))
        relevant_variables <- unique(c(relevant_variables, variables))

      }
    }

    relevant_variables <- intersect(colnames(data), relevant_variables)

    # inputs removed
    x_drop <- x_vars[!which(colnames(data) %in% relevant_variables)]

  }

  # table of scores
  table_scores <- matrix (
    ncol = nY + 1,
    nrow = nrow(data),
    dimnames = list(NULL, c("y_all", paste("y", 1:nY, sep = "")))
  ) %>% as.data.frame()

  x_filtered <- if (length(x_drop) == 0) x_vars else x_vars[- x_drop]

  table_scores[, 1] <- rad_out (
    tech_xmat = as.matrix(DMUs[, x_filtered]),
    tech_ymat = as.matrix(DMUs[, y_vars]),
    eval_xmat = as.matrix(DMUs[, x_filtered]),
    eval_ymat = as.matrix(DMUs[, y_vars]),
    convexity = TRUE,
    returns = "variable"
  )[, 1]

  for (out in 1:nY) {

    table_scores[, 1 + out] <- rad_out (
      tech_xmat = as.matrix(DMUs[, x_filtered]),
      tech_ymat = as.matrix(DMUs[, y_vars[out]]),
      eval_xmat = as.matrix(DMUs[, x_filtered]),
      eval_ymat = as.matrix(DMUs[, y_vars[out]]),
      convexity = TRUE,
      returns = "variable"
    )[, 1]

  }

  # weights for error metrics based on DEA
  dea_scores <-  table_scores[, 2:ncol(table_scores)]

  # matrix with:
  # row 1: the index of the variable
  # row 2: the degree of the variable
  xi_degree <- matrix (
    c(x, rep(1, length(x))),
    byrow = TRUE,
    nrow = 2,
    ncol = length(x)
  )

  if (!is.list(max_degree)) {

    v <- 0

    for (i in 1:max_degree) {

      combs <- combn(1:length(x_vars), i)

      for (k in 1:ncol(combs)) {
        v <- v + 1
        xi_degree[2, v] <- i
      }

    }

  } else {

    xi_degree[2, x_vars] <- 1

    for (k in 1:length(max_degree)) {
      xi_degree[2, length(x_vars) + k] <- length(max_degree[[k]])
    }

  }

  # ===================== #
  #   FORWARD ALGORITHM   #
  # ===================== #

  y_hat <- matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, max)

  # lack-of-fit
  LOF <- err_metric (
    y_obs = data[, y, drop = F],
    y_hat = y_hat,
    metric = metric,
    weight = 1 / dea_scores
  )

  # basis function
  #     id: index
  # status: intercept / paired / unpaired
  #   side: E (entire) / R (right) / L (left)
  #     Bp: basis function
  #     xi: variable for splitting
  #      t: knot for splitting
  #      R: mean error between true data and predicted data (B %*% coefs)
  #  coefs: regression coefficients

  bf <- list (
    "id" = 1,
    "status" = "intercept",
    "side" = "E",
    "Bp" = rep(1, N),
    "xi" = c(- 1),
    "t" = c(- 1),
    "LOF" = LOF,
    "coefs" = unname(apply(data[, y, drop = FALSE], 2, max))
  )

  # set of knots to save indexes of data used as knots
  kn_list <- vector("list", nX)

  # set of basis functions by variable
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
  err <- bf[["LOF"]]

  # set the grid of knots
  kn_grid <- set_knots_grid (
    data = data,
    n_input_1 = length(x_vars),
    n_input_2 = nX,
    kn_grid = kn_grid,
    quick_aces = quick_aces,
    dea_scores = table_scores[, 1]
  )

  # minimum span (minspan) and end span (endspan)
  L_Le <- compute_span (
    kn_grid = kn_grid,
    minspan = minspan,
    endspan = endspan,
    n_input = nX
  )

  # list to save technologies created through RF-ACES
  technology <- list()

  # initial error
  err_min <- err

  while(length(rf_aces[["bf_set"]]) + 2 < max_terms) {

    last_bf <- rf_aces[["bf_set"]][[length(rf_aces[["bf_set"]])]]

    while (last_bf[["id"]] == 1) {

      # random input for bagging
      random_x <- sample(x, max_feats)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function (
        data = data,
        x = random_x,
        y = y,
        xi_degree = xi_degree,
        inter_cost = inter_cost,
        model_type = "envelopment",
        dea_scores = dea_scores,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        span = c(L_Le[[1]], L_Le[[2]]),
        err_min = err,
        var_imp = var_imp,
        quick_aces = quick_aces
      )

      if (is.list(B_bf_knt_err)) break

    }

    if (last_bf[["id"]] != 1) {

      # random input for bagging
      random_x <- sample(x, max_feats)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function (
        data = data,
        x = random_x,
        y = y,
        xi_degree = xi_degree,
        inter_cost = inter_cost,
        model_type = "envelopment",
        dea_scores = dea_scores,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        span = c(L_Le[[1]], L_Le[[2]]),
        err_min = err,
        var_imp = var_imp,
        quick_aces = quick_aces
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

  # generate technology
  technology[["rf_aces"]] <- generate_technology (
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = rf_aces[["Bmatx"]] %*% rf_aces[["coefs"]],
    error_type = error_type,
    table_scores = table_scores,
    ptto = shape[["ptto"]]
  )

  # ======================= #
  #   SMOOTHING PROCEDURE   #
  # ======================= #

  # select a model to be smoothed
  rf_aces_smoothed <- rf_aces

  # distance between knots in smooth models
  wc <- seq(1, 2, length.out = 5)
  wq <- seq(8 / 7, 1.5, length.out = 5)

  # data.frame of knots
  kn_smoothed_right <- rf_aces_smoothed[["knots"]]
  kn_smoothed_right$side <- "R"

  kn_smoothed_left <- rf_aces_smoothed[["knots"]]
  kn_smoothed_left$side <- "L"

  kn_smoothed <- rbind (
    kn_smoothed_right,
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

  rf_aces_cubic <- cubic_aces (
    data = data,
    x = x,
    y = y,
    dea_scores = dea_scores,
    model_type = "envelopment",
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    kn_penalty = 1,
    xi_degree = xi_degree,
    wc = wc
  )

  # generate technology
  technology[["rf_aces_cubic"]] <- generate_technology (
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = rf_aces_cubic[["Bmatx"]] %*% rf_aces_cubic[["coefs"]],
    error_type = error_type,
    table_scores = table_scores,
    ptto = shape[["ptto"]]
  )

  rf_aces_quintic <- quintic_aces (
    data = data,
    x = x,
    y = y,
    dea_scores = dea_scores,
    model_type = "envelopment",
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    kn_penalty = 1,
    xi_degree = xi_degree,
    wq = wq
  )

  # generate technology
  technology[["rf_aces_quintic"]] <- generate_technology (
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = rf_aces_quintic[["Bmatx"]] %*% rf_aces_quintic[["coefs"]],
    error_type = error_type,
    table_scores = table_scores,
    ptto = shape[["ptto"]]
  )

  # =========== #
  # ACES OBJECT #
  # =========== #

  RF_ACES <- aces_object (
    data = DMUs,
    x = x_vars,
    y = y_vars,
    quick_aces = quick_aces,
    error_type = error_type,
    max_degree = max_degree,
    inter_cost = inter_cost,
    xi_degree = xi_degree,
    metric = metric,
    shape = shape,
    max_terms = ncol(rf_aces[["Bmatx"]]),
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = NULL,
    wc = rf_aces_cubic[["w"]],
    wq = rf_aces_quintic[["w"]],
    aces_forward = rf_aces,
    aces = NULL,
    aces_cubic = rf_aces_cubic,
    aces_quintic = rf_aces_quintic,
    technology = technology
  )

  RF_ACES[["data"]] <- NULL

  RF_ACES[["control"]] <- list (
    "xi_degree" = RF_ACES[["control"]][["xi_degree"]],
    "wc" = RF_ACES[["control"]][["wc"]],
    "wq" = RF_ACES[["control"]][["wq"]]
  )

  RF_ACES[["methods"]][["aces"]] <- NULL
  names(RF_ACES[["methods"]]) <- c (
    "rf_aces",
    "rf_aces_cubic",
    "rf_aces_quintic"
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

  # calculate the OOB-MSE: index
  oob_idxs <- as.numeric(rownames(model_pred))

  # calculate the OOB-MSE: value
  oob_mse <- sum((data[oob_idxs, y] - model_pred[, 1:nY, drop = F]) ^ 2) /
    (length(oob_idxs) * nY)

  return(oob_mse)

}

#' @title Compute Variable Importance for Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).
#'
#' @description
#' Computes a robust measure of variable importance for a fitted Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) model. Importance scores are obtained by evaluating the increase in prediction error when permuting each input variable as in \insertCite{breiman2001;textual}{aces}.
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
#' @references
#' \insertRef{breiman2001}{aces} \cr
#'
#' @return
#' This function returns a metric of variable importance for each input variable.

rf_aces_varimp <- function (
    data,
    x,
    y,
    object,
    repeats = 1
    ) {

  # RF-ACES configuration
  control_features <- object[["control"]]

  quick_aces <- control_features[["quick_aces"]]
  error_type <- control_features[["error_type"]]
  max_degree <- control_features[["max_degree"]]
  inter_cost <- control_features[["inter_cost"]]
  metric <- control_features[["metric"]]
  shape <- control_features[["shape"]]
  learners <- control_features[["learners"]]
  bag_size <- control_features[["bag_size"]]
  max_feats <- control_features[["max_feats"]]
  early_stopping <- control_features[["early_stopping"]]
  max_terms <- control_features[["max_terms"]]
  err_red <- control_features[["err_red"]]
  kn_grid <- control_features[["kn_grid"]]
  minspan <- control_features[["minspan"]]
  endspan <- control_features[["endspan"]]

  # initialize matrix of variable importance
  mat_varimp <- matrix (
    0,
    nrow = length(x)
  )

  rownames(mat_varimp) <- colnames(data)[x]
  colnames(mat_varimp) <- "importance"

  # out-of-bag error of the complete model
  oob <- object[["forest"]][[length(object[["forest"]])]][["OOB"]]

  # out-of-bag excluding the j-th input variable
  oob_j <- c()

  # compute oob shuffling each variable "j"
  for (j in 1:length(x)) {
    for (n in 1:repeats) {

      data_shuffled <- data
      xvar_shuffled <- data[sample(1:nrow(data)), x[j]]
      data_shuffled[, x[j]] <- xvar_shuffled

      model_j <- rf_aces (
        data = data_shuffled,
        x = x,
        y = y,
        quick_aces = quick_aces,
        error_type = error_type,
        mul_BF = list (
          "max_degree" = max_degree,
          "inter_cost" = inter_cost
        ),
        metric = metric,
        shape = shape,
        learners = learners,
        bag_size = bag_size,
        max_feats = max_feats,
        early_stopping = early_stopping,
        max_terms = max_terms,
        err_red = err_red,
        kn_grid = kn_grid,
        minspan = minspan,
        endspan = endspan
        )

        # out-of-bag error for the model with the jth variable shuffled
        oob_j <- c(oob_j, model_j[["forest"]][[length(model_j[["forest"]])]][["OOB"]])

      }

      # add oob_j to the matrix of variable importance
      mat_varimp[j, "importance"] <- round(100 * ((mean(oob_j) - oob) / mean(oob_j)), 2)

  }

  ranking <- mat_varimp[order(mat_varimp[, 1], decreasing = TRUE), ]

  return(rankings)

}

#' @title Percentile-Based Prediction for Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).
#'
#' @description
#' Predicts the optimal output level under the RF-ACES model using a specified percentile of the conditional distribution, rather than the mean. As in DEA, inputs are held fixed, and the output is predicted via Random Forests by selecting a given percentile (e.g., the 5th percentile) of the conditional distribution. This approach enables modeling conservative or optimistic production frontiers by adjusting the percentile parameter.
#'
#' @param object
#' A \code{rf_aces} object.
#'
#' @param newdata
#' A \code{data.frame} or \code{matrix} containing the variables in the model.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param p
#' A numeric value between 0 and 1 specifying the conditional percentile to be used for prediction.
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
#' A numeric vector of predicted output values for each observation, based on the specified percentile.

rf_aces_p_predict <- function (
    object,
    newdata,
    x,
    y,
    p,
    method
    ) {

  RF_ACES <- object

  xmat_rf_aces <- xmat_rf_aces_cubic <- xmat_rf_aces_quintic <- NULL
  ymat_rf_aces <- ymat_rf_aces_cubic <- ymat_rf_aces_quintic <- NULL

  # save technology for Random Forest
  for (i in seq_along(RF_ACES[["forest"]])) {

    xmat_rf_aces <- rbind (
      xmat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["xmat"]]
    )

    ymat_rf_aces <- rbind (
      ymat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["ymat"]]
    )

    xmat_rf_aces_cubic <- rbind (
      xmat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["xmat"]]
    )

    ymat_rf_aces_cubic <- rbind (
      ymat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["ymat"]]
    )

    xmat_rf_aces_quintic <- rbind (
      xmat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["xmat"]]
    )

    ymat_rf_aces_quintic <- rbind (
      ymat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["ymat"]]
    )

  }

  data_sorted <- as.data.frame (
    set_data (
      data = newdata,
      x = x,
      y = y,
      max_degree = 1,
      error_type = "add"
    )
  )

  inp_cols <- names(data_sorted)[1:length(x)]
  out_cols <- names(data_sorted)[(length(x) + 1):(length(x) + length(y))]

  colnames(xmat_rf_aces) <- colnames(xmat_rf_aces_cubic) <- colnames(xmat_rf_aces_quintic) <- inp_cols
  colnames(ymat_rf_aces) <- colnames(ymat_rf_aces_cubic) <- colnames(ymat_rf_aces_quintic) <- out_cols

  rf_aces_technology <- data.frame(xmat_rf_aces, ymat_rf_aces) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), ~ quantile(.x, probs = p), .names = "{col}"), .groups = 'drop')

  rf_aces_cubic_technology <- data.frame(xmat_rf_aces_cubic, ymat_rf_aces_cubic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), ~ quantile(.x, probs = p), .names = "{col}"), .groups = 'drop')

  rf_aces_quintic_technology <- data.frame(xmat_rf_aces_quintic, ymat_rf_aces_quintic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), ~ quantile(.x, probs = p), .names = "{col}"), .groups = 'drop')

  RF_ACES[["technology"]] <- list (
    "rf_aces" = list (
      "xmat" = as.matrix(rf_aces_technology[, x]),
      "ymat" = as.matrix(rf_aces_technology[, y])
    ),
    "rf_aces_cubic" = list (
      "xmat" = as.matrix(rf_aces_cubic_technology[, x]),
      "ymat" = as.matrix(rf_aces_cubic_technology[, y])
    ),
    "rf_aces_quintic" = list (
      "xmat" = as.matrix(rf_aces_quintic_technology[, x]),
      "ymat" = as.matrix(rf_aces_quintic_technology[, y])
    )
  )

  # number of outputs
  nY <- length(RF_ACES[["data"]][["y"]])

  # determine error type
  error_type <- RF_ACES[["control"]][["error_type"]]

  # data in [x, y] format with interaction of variables included
  data <- set_data (
    data = newdata,
    x = x,
    y = NULL,
    max_degree = RF_ACES[["control"]][["max_degree"]],
    error_type = error_type
  )

  # technology
  tecno <- RF_ACES[["technology"]][[method]]

  # number of models in Random Forest
  RF_models <- length(RF_ACES[["forest"]])

  # list of predictions for each model
  y_hat_RF <- vector("list", RF_models)

  for (t in 1:RF_models) {

    # output predictions
    y_hat <- as.data.frame(matrix(NA, nrow = nrow(newdata), ncol = nY))

    # select an element from the forest
    model <- RF_ACES[["forest"]][[t]]

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

  names(y_hat) <- paste(RF_ACES[["data"]][["ynames"]], "_pred", sep = "")

  return(y_hat)

}

