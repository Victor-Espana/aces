#' @title Fit an RF-ACES Model
#'
#' @description
#'
#' Estimates maximum attainable outputs with an ensemble of randomized ACES
#' learners and uses their results to build an aggregate production technology.
#' Each learner is trained on a bootstrap sample and considers a random subset
#' of inputs, following the main ideas of Random Forests.
#'
#' @details
#' When \code{shape$mono} and \code{shape$conc} are both \code{TRUE}, each
#' learner estimates a non-decreasing and concave production function for every
#' output. If either restriction is relaxed, its spline models are intermediate
#' predictors of maximum attainable output. In both cases, RF-ACES constructs
#' learner-specific technologies and combines their reference points. Their
#' envelopment supplies the production-set structure used for efficiency
#' measurement. Aggregating the learners can reduce the variability of a single
#' ACES fit. See
#' \insertCite{espana2024rf;textual}{aces} for methodological details.
#'
#' @name rf_aces
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the model variables.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param scale_data
#' If \code{TRUE}, divide each input and output by its mean before fitting. This
#' can improve solver convergence.
#'
#' @param quick_aces
#' If \code{TRUE}, use Quick ACES. It removes an input only when its Spearman and
#' Kendall correlations fall below their thresholds for every output, reduces
#' the knot grid around efficient DMUs, and allocates later candidates according
#' to recent error reductions. Examining fewer basis functions can substantially
#' shorten fitting time for larger problems. The internal reduction history is
#' retained in each learner and can be summarized with \code{aces_varimp()}.
#'
#' @param mul_BF
#' A \code{list} with two elements:
#' \itemize{
#'   \item{\code{max_degree}: A positive integer giving the maximum interaction
#'   degree, or a list of input-index vectors defining the interactions to use.}
#'   \item{\code{inter_cost}: The minimum relative improvement over the best
#'   first-degree basis function required to add a higher-degree basis function.
#'   The default, \code{0.05}, means 5 percent.}
#' }
#'
#' @param metric
#' A character string specifying the lack-of-fit measure:
#' \itemize{
#'   \item{\code{"mae"}: Mean absolute error.}
#'   \item{\code{"mape"}: Mean absolute percentage error.}
#'   \item{\code{"mse"}: Mean squared error.}
#'   \item{\code{"msle"}: Mean squared logarithmic error.}
#'   \item{\code{"rmse"}: Root mean squared error.}
#'   \item{\code{"nrmse1"}: Root mean squared error divided by the mean.}
#'   \item{\code{"nrmse2"}: Root mean squared error divided by the range.}
#' }
#'
#' @param shape
#' A \code{list} controlling the restrictions imposed during spline estimation:
#' \itemize{
#'   \item{\code{mono}: Enforce non-decreasing monotonicity.}
#'   \item{\code{conc}: Enforce concavity.}
#' }
#' When both are \code{TRUE}, each learner estimates a production function for
#' every output. Otherwise, the spline models estimate maximum attainable outputs
#' before the learner-specific technology is constructed.
#'
#' @param learners
#' Positive integer giving the maximum number of ACES learners in the forest.
#' Fitting may stop earlier when the OOB stopping rule is satisfied.
#'
#' @param bag_size
#' Positive integer giving the number of observations sampled with replacement
#' for each learner. Values smaller than the sample size increase diversity
#' between learners.
#'
#' @param max_feats
#' Number of inputs randomly considered by each learner. Smaller values increase
#' randomization; larger values make learners more similar to standard ACES
#' fits.
#'
#' @param early_stopping
#' A \code{list} that controls early stopping based on out-of-bag (OOB) error:
#' \itemize{
#'   \item{\code{ma_window}: Size of the moving-average window used to smooth
#'   OOB error. Larger windows are more stable but react more slowly.}
#'   \item{\code{tolerance}: Number of consecutive iterations without
#'   improvement before fitting stops. Larger values allow more learners to be
#'   fitted before stopping.}
#' }
#'
#' @param max_terms
#' Positive integer giving the maximum number of terms created during the
#' forward step. The algorithm may stop earlier when the improvement falls
#' below \code{err_red}.
#'
#' @param err_red
#' Minimum relative error reduction required to add a new pair of first-degree
#' basis functions. Larger values produce simpler models and stop the forward
#' step sooner.
#'
#' @param kn_grid
#' Knot candidates. Use \code{-1} to select eligible knots from the observed
#' input values following the standard ACES procedure. Alternatively, supply a
#' list with one numeric vector of candidates for each original input.
#'
#' @param minspan
#' Minimum number of observations between adjacent knots. Use \code{-2} for
#' the rule in \insertCite{zhang1994;textual}{aces}, \code{-1} for the rule in
#' \insertCite{friedman1991;textual}{aces}, or a positive integer.
#'
#' @param endspan
#' Minimum number of observations before the first knot and after the last knot.
#' Use \code{-2} for the rule in \insertCite{zhang1994;textual}{aces}, \code{-1}
#' for the rule in \insertCite{friedman1991;textual}{aces}, or a positive integer.
#'
#' @references
#'
#' \insertRef{espana2024rf}{aces} \cr \cr
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{espana2025}{aces} \cr \cr
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
#' An object of class \code{rf_aces}. It contains the training data and controls,
#' all fitted learners, OOB information, and the aggregated technology reference
#' points used for prediction and efficiency measurement.
#'
#' @export

rf_aces <- function(
  data,
  x,
  y,
  scale_data = TRUE,
  quick_aces = FALSE,
  mul_BF = list(
    "max_degree" = 1,
    "inter_cost" = 0.05
  ),
  metric = "mse",
  shape = list(
    "mono" = TRUE,
    "conc" = TRUE
  ),
  learners = 100,
  bag_size = nrow(data),
  max_feats = length(x) / 3,
  early_stopping = list(
    "ma_window" = 10,
    "tolerance" = 5
  ),
  max_terms = 50,
  err_red = 0.01,
  kn_grid = -1,
  minspan = -1,
  endspan = -1
  ) {

  # ====================== #
  # DISPLAY ERRORS RF-ACES #
  # ====================== #

  display_errors_rf_aces(
    data = data,
    x = x,
    y = y,
    scale_data = scale_data,
    quick_aces = quick_aces,
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

  # Preserve the exact public fitting specification. It is required by LOCO
  # importance, which refits the complete forest after removing one input.
  fit_spec <- list(
    scale_data = scale_data,
    quick_aces = quick_aces,
    mul_BF = mul_BF,
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

  # =================== #
  #    SCALING SETUP    #
  # =================== #

  # Factores por defecto (Neutros)
  sx <- rep(1, length(x))
  sy <- rep(1, length(y))

  # Copia de trabajo
  data_algo <- data

  if (scale_data) {
    # 1. Extraer matrices crudas
    raw_x <- as.matrix(data[, x])
    raw_y <- as.matrix(data[, y])

    # 2. Calcular factores (Media)
    sx <- colMeans(raw_x, na.rm = TRUE)
    sx[sx == 0] <- 1

    sy <- colMeans(raw_y, na.rm = TRUE)
    sy[sy == 0] <- 1

    # 3. Aplicar escalado
    data_algo[, x] <- sweep(raw_x, 2, sx, "/")
    data_algo[, y] <- sweep(raw_y, 2, sy, "/")
  }

  # ========= #
  # ALGORITHM #
  # ========= #

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
    sample_bag <- sample(
      1:nrow(data),
      size = bag_size,
      replace = TRUE
    )

    # data for the m-model
    data_bag <- data_algo[sample_bag, ]

    # out of bag data for the m-model
    inb_idxs <- 1:nrow(data_algo) %in% sample_bag
    oob_idxs <- which(!(1:nrow(data_algo) %in% sample_bag))
    data_oob <- data_algo[oob_idxs, ]

    RF_ACES[[m]] <- rf_aces_algorithm(
      data = data_bag,
      x_vars = x,
      y_vars = y,
      quick_aces = quick_aces,
      max_degree = mul_BF[["max_degree"]],
      inter_cost = mul_BF[["inter_cost"]],
      metric = metric,
      shape = list(
        "mono" = shape[["mono"]],
        "conc" = shape[["conc"]]
      ),
      max_feats = max_feats,
      max_terms = max_terms,
      err_red = err_red,
      minspan = minspan,
      endspan = endspan,
      kn_grid = kn_grid
    )

    # predictions for OOB observations
    aces_model <- RF_ACES[[m]][["methods"]][["rf_aces"]]
    coefs <- aces_model[["coefs"]]
    knots <- aces_model[["knots"]]

    oob_expanded <- set_data(
      data = data_oob,
      x = x,
      y = y,
      max_degree = mul_BF[["max_degree"]]
    )

    x_algo <- 1:(ncol(oob_expanded) - length(y))

    B_oob <- set_Bmat(
      newdata = oob_expanded[, x_algo, drop = FALSE],
      model = aces_model,
      knots = knots,
      method = "rf_aces"
    )

    y_hat <- B_oob %*% coefs

    oob_pred[[m]] <- matrix(
      c(
        "y_hat" = y_hat,
        "idx" = oob_idxs
      ),
      ncol = length(y) + 1
    )

    # OOB performance
    oob_mse <- compute_oob(
      data = data_algo,
      y = y,
      oob_idxs = oob_idxs,
      oob_pred = oob_pred,
      model = m
    )

    # out-of-bag error
    RF_ACES[[m]][["OOB"]] <- oob_mse

    # store bag indices for variable importance
    RF_ACES[[m]][["sample_bag"]] <- sample_bag

    # early stopping RF-ACES based on moving average
    if (m > ma_window) {
      # compute moving average over the previous `ma_window` iterations
      oob_ma <- mean(tail(oob_trend, ma_window))

      # current oob greater than or equal to moving average
      if (oob_mse >= oob_ma) {
        stopping_condition_counter <- stopping_condition_counter + 1
      } else {
        stopping_condition_counter <- 0
      }
    }

    oob_trend <- c(oob_trend, oob_mse)

    if (stopping_condition_counter == tolerance) {
      cat("\n")
      print(paste("Out-of-bag error stabilized after", m, "models"))
      break
    }
  }

  # random-forest hyperparameters
  rf_list <- list(
    learners = m,
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
    xmat_rf_aces <- rbind(
      xmat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["xmat"]]
    )

    ymat_rf_aces <- rbind(
      ymat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["ymat"]]
    )

    xmat_rf_aces_cubic <- rbind(
      xmat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["xmat"]]
    )

    ymat_rf_aces_cubic <- rbind(
      ymat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["ymat"]]
    )

    xmat_rf_aces_quintic <- rbind(
      xmat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["xmat"]]
    )

    ymat_rf_aces_quintic <- rbind(
      ymat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["ymat"]]
    )
  }

  data_sorted <- as.data.frame(
    set_data(
      data = data,
      x = x,
      y = y,
      max_degree = 1
    )
  )

  inp_cols <- names(data_sorted)[1:length(x)]
  out_cols <- names(data_sorted)[(length(x) + 1):(length(x) + length(y))]

  colnames(xmat_rf_aces) <- colnames(xmat_rf_aces_cubic) <- colnames(xmat_rf_aces_quintic) <- inp_cols
  colnames(ymat_rf_aces) <- colnames(ymat_rf_aces_cubic) <- colnames(ymat_rf_aces_quintic) <- out_cols

  rf_aces_technology <- data.frame(xmat_rf_aces, ymat_rf_aces) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), mean, .names = "{col}"), .groups = "drop")

  rf_aces_cubic_technology <- data.frame(xmat_rf_aces_cubic, ymat_rf_aces_cubic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), mean, .names = "{col}"), .groups = "drop")

  rf_aces_quintic_technology <- data.frame(xmat_rf_aces_quintic, ymat_rf_aces_quintic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), mean, .names = "{col}"), .groups = "drop")

  RF_ACES[["technology"]] <- list(
    "rf_aces" = list(
      "xmat" = as.matrix(rf_aces_technology[, 1:length(x)]),
      "ymat" = as.matrix(rf_aces_technology[, (length(x) + 1):ncol(rf_aces_technology)])
    ),
    "rf_aces_cubic" = list(
      "xmat" = as.matrix(rf_aces_cubic_technology[, 1:length(x)]),
      "ymat" = as.matrix(rf_aces_cubic_technology[, (length(x) + 1):ncol(rf_aces_cubic_technology)])
    ),
    "rf_aces_quintic" = list(
      "xmat" = as.matrix(rf_aces_quintic_technology[, 1:length(x)]),
      "ymat" = as.matrix(rf_aces_quintic_technology[, (length(x) + 1):ncol(rf_aces_quintic_technology)])
    )
  )

  # ====== #
  # Object #
  # ====== #

  RF_ACES[["data"]] <- list(
    "df" = data,
    "x" = x,
    "y" = y,
    "xnames" = colnames(data)[x],
    "ynames" = colnames(data)[y],
    "rownames" = rownames(data)
  )

  RF_ACES[["control"]] <- list(
    "quick_aces" = quick_aces,
    "max_degree" = mul_BF[["max_degree"]],
    "inter_cost" = mul_BF[["inter_cost"]],
    "xi_degree" = RF_ACES[["forest"]][[1]][["control"]][["xi_degree"]],
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
    "kn_grid" = RF_ACES[["forest"]][[1]][["control"]][["kn_grid"]],
    "kn_penalty" = NULL,
    "psi" = 0.05,
    "wc" = NULL,
    "wq" = NULL,
    scale = list(
      is_scaled = scale_data,
      mean_x = sx,
      mean_y = sy
    ),
    fit_spec = fit_spec
  )

  # change the order
  RF_ACES <- RF_ACES[c("data", "control", "forest", "technology")]

  # type of object
  class(RF_ACES) <- "rf_aces"

  return(RF_ACES)
}


#' @title Run One RF-ACES Learner
#'
#' @description
#' Fits one randomized ACES learner for use in an RF-ACES forest. The function
#' expands interactions, samples the allowed inputs, estimates maximum attainable
#' outputs, refines the output vectors, and constructs the learner-specific
#' technology.
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the model variables.
#'
#' @param x_vars
#' Column indexes of input variables in \code{data}.
#'
#' @param y_vars
#' Column indexes of output variables in \code{data}.
#'
#' @param quick_aces
#' If \code{TRUE}, use the Quick ACES input filtering, knot reduction, and
#' adaptive candidate allocation rules.
#'
#' @param max_degree
#' Maximum interaction degree, or a list of input-index vectors defining the
#' interactions to use.
#'
#' @param inter_cost
#' Minimum relative improvement over the best first-degree basis function
#' required to add a higher-degree basis function.
#'
#' @param metric
#' Character string specifying the lack-of-fit measure.
#'
#' @param shape
#' A list with logical elements \code{mono} and \code{conc}. When both are
#' \code{TRUE}, estimate a non-decreasing and concave production function for
#' each output; otherwise, use the spline models as intermediate predictors of
#' maximum attainable output before constructing the learner-specific technology.
#'
#' @param max_feats
#' Number of inputs randomly considered by the learner.
#'
#' @param max_terms
#' Maximum number of terms created during the forward step.
#'
#' @param err_red
#' Minimum relative error reduction required to add a new pair of first-degree
#' basis functions.
#'
#' @param minspan
#' Minimum number of observations between two adjacent knots.
#'
#' @param endspan
#' Minimum number of observations before the first knot and after the last knot.
#'
#' @param kn_grid
#' Knot candidates. Use \code{-1} for automatic selection or supply a list.
#'
#' @importFrom dplyr desc
#'
#' @return
#'
#' An \code{aces} object representing one RF-ACES learner.
#'
#' @export

rf_aces_algorithm <- function(
  data,
  x_vars,
  y_vars,
  quick_aces,
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
  data <- set_data(
    data = data,
    x = x_vars,
    y = y_vars,
    max_degree = max_degree
  )

  # samples size
  N <- nrow(data)

  # set "x", "z" and "y" indexes in data
  x <- 1:(ncol(data) - length(y_vars))
  y <- (length(x) + 1):ncol(data)

  # maximum features before and after interaction variables
  max_feats <- max_feats / length(x_vars)
  max_feats <- max(round(max_feats * length(x), 0), 1)

  # set number of inputs and outputs
  nX <- length(x)
  nY <- length(y)

  # ================== #
  # VARIABLE FILTERING #
  # ================== #

  # Best candidate reduction for each prepared input in the current iteration.
  # NA means that the input was not evaluated in that iteration.
  var_imp <- matrix(
    rep(NA_real_, nX),
    nrow = 1
  )
  colnames(var_imp) <- colnames(data)[1:nX]

  quick_keep <- rep(TRUE, nX)

  x_filtered <- x_vars

  # remove variables with low correlation
  if (quick_aces) {
    correlation_filter <- quick_aces_correlation_filter(data, x, y)
    quick_keep <- unname(correlation_filter$keep)

    x_filtered <- quick_aces_retained_inputs(
      keep = quick_keep,
      x_vars = x_vars,
      max_degree = max_degree
    )
  }

  # table of scores
  table_scores <- matrix(
    ncol = nY + 1,
    nrow = nrow(data),
    dimnames = list(NULL, c("y_all", paste("y", 1:nY, sep = "")))
  ) %>% as.data.frame()

  # ========== #
  # DEA SCORES #
  # ========== #

  table_scores[, 1] <- rad_out(
    tech_xmat = as.matrix(DMUs[, x_filtered]),
    tech_ymat = as.matrix(DMUs[, y_vars]),
    eval_xmat = as.matrix(DMUs[, x_filtered]),
    eval_ymat = as.matrix(DMUs[, y_vars]),
    convexity = TRUE,
    returns = "variable"
  )[, 1]

  for (out in 1:nY) {
    table_scores[, 1 + out] <- rad_out(
      tech_xmat = as.matrix(DMUs[, x_filtered]),
      tech_ymat = as.matrix(DMUs[, y_vars[out]]),
      eval_xmat = as.matrix(DMUs[, x_filtered]),
      eval_ymat = as.matrix(DMUs[, y_vars[out]]),
      convexity = TRUE,
      returns = "variable"
    )[, 1]
  }

  # weights for error metrics based on DEA
  dea_scores <- table_scores[, 2:ncol(table_scores)]

  # ========== #
  # FDH SCORES #
  # ========== #

  table_scores[, 1] <- rad_out(
    tech_xmat = as.matrix(DMUs[, x_filtered]),
    tech_ymat = as.matrix(DMUs[, y_vars]),
    eval_xmat = as.matrix(DMUs[, x_filtered]),
    eval_ymat = as.matrix(DMUs[, y_vars]),
    convexity = FALSE,
    returns = "variable"
  )[, 1]

  for (out in 1:nY) {
    table_scores[, 1 + out] <- rad_out(
      tech_xmat = as.matrix(DMUs[, x_filtered]),
      tech_ymat = as.matrix(DMUs[, y_vars[out]]),
      eval_xmat = as.matrix(DMUs[, x_filtered]),
      eval_ymat = as.matrix(DMUs[, y_vars[out]]),
      convexity = FALSE,
      returns = "variable"
    )[, 1]
  }

  # weights for error metrics based on DEA
  fdh_scores <- table_scores[, 2:ncol(table_scores), drop = FALSE]

  # ==================== #
  # VARIABLE INTERACTION #
  # ==================== #

  # matrix with:
  # row 1: the index of the variable
  # row 2: the degree of the variable
  xi_degree <- matrix(
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
    xi_degree[2, seq_along(x_vars)] <- 1

    for (k in 1:length(max_degree)) {
      xi_degree[2, length(x_vars) + k] <- length(max_degree[[k]])
    }
  }

  # ===================== #
  #   FORWARD ALGORITHM   #
  # ===================== #

  y_hat <- matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, max)

  # lack-of-fit
  LOF <- err_metric(
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
  #    LOF: mean error between true data and predicted data (B %*% coefs)
  #  coefs: regression coefficients

  bf <- list(
    "id" = 1,
    "status" = "intercept",
    "side" = "E",
    "Bp" = rep(1, N),
    "xi" = c(-1),
    "t" = c(-1),
    "LOF" = LOF,
    "coefs" = unname(apply(data[, y, drop = FALSE], 2, max))
  )

  # set of knots to save indexes of data used as knots
  kn_list <- vector("list", nX)

  # set of basis functions by variable
  Bp_list <- vector("list", nX)

  for (xi in 1:nX) {
    Bp_list[[xi]] <- list(
      "paired" = NULL,
      "right" = NULL,
      "left" = NULL
    )
  }

  # set of basis functions (bf_set) and the matrix of basis functions (B)
  rf_aces <- list(
    "bf_set" = list(bf),
    "B" = matrix(rep(1, N))
  )

  # error of the first basis function
  err <- bf[["LOF"]]

  # set the grid of knots
  kn_grid <- set_knots_grid(
    data = data,
    n_input_1 = length(x_vars),
    n_input_2 = nX,
    kn_grid = kn_grid,
    quick_aces = quick_aces,
    dea_scores = table_scores[, 1]
  )

  # minimum span (minspan) and end span (endspan)
  L_Le <- compute_span(
    kn_grid = kn_grid,
    minspan = minspan,
    endspan = endspan,
    n_input = nX
  )

  # list to save technologies created through RF-ACES
  technology <- list()

  # initial error
  err_min <- err

  # number of accepted forward iterations
  iter <- 0L

  while (length(rf_aces[["bf_set"]]) + 2 < max_terms) {
    last_bf <- rf_aces[["bf_set"]][[length(rf_aces[["bf_set"]])]]

    while (last_bf[["id"]] == 1) {
      # random input for bagging
      random_x <- sample(x, max_feats)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function(
        data = data,
        x = random_x,
        y = y,
        xi_degree = xi_degree,
        inter_cost = inter_cost,
        dea_scores = dea_scores,
        fdh_scores = fdh_scores,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        span = c(L_Le[[1]], L_Le[[2]]),
        err_min = err,
        var_imp = var_imp,
        quick_keep = quick_keep,
        quick_aces = quick_aces
      )

      if (is.list(B_bf_knt_err)) break
    }

    if (last_bf[["id"]] != 1) {
      # random input for bagging
      random_x <- sample(x, max_feats)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function(
        data = data,
        x = random_x,
        y = y,
        xi_degree = xi_degree,
        inter_cost = inter_cost,
        dea_scores = dea_scores,
        fdh_scores = fdh_scores,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        span = c(L_Le[[1]], L_Le[[2]]),
        err_min = err,
        var_imp = var_imp,
        quick_keep = quick_keep,
        quick_aces = quick_aces
      )
    }

    if (!is.list(B_bf_knt_err)) break

    # new best error
    new_err <- B_bf_knt_err[[5]]

    # update model
    if (last_bf[["id"]] == 1 || new_err[1] < err[1] * (1 - err_red[1])) {
      iter <- iter + 1L

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

      # retain all candidate reductions, not only the variable selected for the
      # split. The next row belongs to the next accepted iteration.
      var_imp <- B_bf_knt_err[[6]]
      var_imp <- rbind(var_imp, rep(NA_real_, nX))
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

  kn_forward <- data.frame(
    xi = var,
    t = knt,
    status = sts
  )

  # ==
  # rf-aces
  # ==

  variable_reduction <- if (iter > 0L) {
    var_imp[seq_len(iter), , drop = FALSE]
  } else {
    var_imp[FALSE, , drop = FALSE]
  }
  rownames(variable_reduction) <- paste0(
    "iteration_",
    seq_len(nrow(variable_reduction))
  )

  rf_aces <- list(
    "basis" = rf_aces[["bf_set"]],
    "Bmatx" = rf_aces[["B"]],
    "knots" = kn_forward,
    "coefs" = rev(rf_aces[["bf_set"]])[[1]][["coefs"]],
    "variable_reduction" = variable_reduction
  )

  # generate technology
  var_names <- colnames(DMUs[, c(x_vars, y_vars)])
  technology[["rf_aces"]] <- generate_technology(
    var_names = var_names,
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = rf_aces[["Bmatx"]] %*% rf_aces[["coefs"]],
    table_scores = table_scores
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

  kn_smoothed <- rbind(
    kn_smoothed_right,
    kn_smoothed_left
  )

  # generate the input space for side knots location
  kn_side_loc <- side_knot_location(
    data = data,
    nX = nX,
    knots = kn_smoothed
  )

  # ==
  # smoothing cubic aces
  # ==

  rf_aces_cubic <- cubic_aces(
    data = data,
    x = x,
    y = y,
    dea_scores = dea_scores,
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    kn_penalty = 1,
    xi_degree = xi_degree,
    wc = wc
  )

  # generate technology
  technology[["rf_aces_cubic"]] <- generate_technology(
    var_names = var_names,
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = rf_aces_cubic[["Bmatx"]] %*% rf_aces_cubic[["coefs"]],
    table_scores = table_scores
  )

  rf_aces_quintic <- quintic_aces(
    data = data,
    x = x,
    y = y,
    dea_scores = dea_scores,
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    kn_penalty = 1,
    xi_degree = xi_degree,
    wq = wq
  )

  # generate technology
  technology[["rf_aces_quintic"]] <- generate_technology(
    var_names = var_names,
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = rf_aces_quintic[["Bmatx"]] %*% rf_aces_quintic[["coefs"]],
    table_scores = table_scores
  )

  # =========== #
  # ACES OBJECT #
  # =========== #

  RF_ACES <- aces_object(
    data = DMUs,
    x = x_vars,
    y = y_vars,
    quick_aces = quick_aces,
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
    psi = NULL,
    wc = rf_aces_cubic[["w"]],
    wq = rf_aces_quintic[["w"]],
    aces_forward = rf_aces,
    aces = NULL,
    aces_cubic = rf_aces_cubic,
    aces_quintic = rf_aces_quintic,
    technology = technology
  )

  RF_ACES[["data"]] <- NULL

  RF_ACES[["control"]] <- list(
    "xi_degree" = RF_ACES[["control"]][["xi_degree"]],
    "kn_grid" = RF_ACES[["control"]][["kn_grid"]],
    "wc" = RF_ACES[["control"]][["wc"]],
    "wq" = RF_ACES[["control"]][["wq"]]
  )

  RF_ACES[["methods"]][["aces"]] <- NULL
  names(RF_ACES[["methods"]]) <- c(
    "rf_aces",
    "rf_aces_cubic",
    "rf_aces_quintic"
  )

  return(RF_ACES)
}

#' @title Compute RF-ACES Out-of-Bag Error
#'
#' @description
#'
#' Combines the OOB predictions available up to a given learner and computes
#' their mean squared error against the observed outputs. Repeated predictions
#' for the same observation are averaged before the error is calculated.
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the model variables.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param oob_idxs
#' Row indexes of observations with at least one OOB prediction. The function
#' updates this value from \code{oob_pred} before computing the error.
#'
#' @param oob_pred
#' A list of OOB predictions. Each element contains predicted outputs and the
#' corresponding row index.
#'
#' @param model
#' Index of the last learner to include.
#'
#' @return
#'
#' A single numeric OOB mean squared error.

compute_oob <- function(
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
  model_pred <- lapply(
    split(
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
  oob_mse <- sum((data[oob_idxs, y] - model_pred[, 1:nY, drop = F])^2) /
    (length(oob_idxs) * nY)

  return(oob_mse)
}

#' @title Compute RF-ACES Variable Importance
#'
#' @description
#' Measures the predictive importance of each original input in an RF-ACES
#' model. Each input is permuted in the out-of-bag (OOB) observations, and the
#' resulting increase in prediction error is averaged across learners and
#' repetitions. Larger positive values indicate greater importance.
#' This legacy interface is retained for compatibility; use
#' \code{aces_varimp(object, importance = "permutation")} for the unified
#' RF-ACES implementation and access to forward and LOCO importance.
#'
#' @details
#' The error measure is the one used to fit \code{object}. For each learner, the
#' function compares its error before and after permuting one input. Interaction
#' variables are then rebuilt, so the result includes the input's main effect
#' and every interaction that contains it. The function always evaluates the
#' linear \code{"rf_aces"} learner, even when a smoothed method is used elsewhere.
#'
#' With \code{normalize = FALSE}, the score is 100 times the mean signed change
#' in the selected error measure; it is not a percentage increase.
#' Positive values mean that permutation worsened prediction, values near zero
#' indicate little change, and negative values mean that permutation improved
#' prediction. With \code{normalize = TRUE}, negative scores are first shifted
#' so that the smallest is zero, and the largest score is scaled to 100. The
#' normalized values therefore provide a relative ranking only.
#'
#' The scores measure predictive sensitivity, not statistical significance or a
#' causal effect. Strongly correlated inputs can substitute for one another and
#' may therefore receive smaller or less stable permutation scores. Raw scores
#' should be compared only when the same error measure and data scale are used.
#'
#' This function is diagnostic. Its ranking is not used by \code{get_scores()}
#' or \code{get_targets()} when \code{relevant = TRUE}; those functions use the
#' inputs represented in selected basis functions instead.
#'
#' @param data
#' Training data used to fit \code{object}, with the same rows, row order, and
#' model variables.
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
#' Positive integer giving the number of permutations for each input. Results
#' are averaged across repetitions.
#'
#' @param normalize
#' If \code{TRUE}, shift negative scores to start at zero and scale the largest
#' score to 100. If \code{FALSE}, return the signed error increase multiplied by
#' 100.
#'
#' @references
#' \insertRef{breiman2001}{aces} \cr
#'
#' @export
#'
#' @return
#' A data frame with columns \code{variable} and \code{importance}, sorted from
#' highest to lowest importance. There is one row for each original input.

rf_aces_varimp <- function(
  data,
  x,
  y,
  object,
  repeats = 1,
  normalize = TRUE
) {
  # number of outputs
  nY <- length(y)

  # number of trees
  n_trees <- length(object[["forest"]])

  # metric used during training
  metric <- object[["control"]][["metric"]]

  # max_degree for building the data matrix
  max_degree <- object[["control"]][["max_degree"]]

  # sample size
  N <- nrow(data)

  # scaling setup
  scaling <- object[["control"]][["scale"]]
  data_work <- data

  if (!is.null(scaling) && scaling$is_scaled) {
    data_work[, x] <- sweep(as.matrix(data[, x]), 2, scaling$mean_x, "/")
    data_work[, y] <- sweep(as.matrix(data[, y]), 2, scaling$mean_y, "/")
  }

  # data in [x, z, y] format with interaction of variables included
  data_algo <- set_data(
    data = data_work,
    x = x,
    y = y,
    max_degree = max_degree
  )

  # indexes in the expanded data
  x_algo <- 1:(ncol(data_algo) - nY)
  y_algo <- (length(x_algo) + 1):ncol(data_algo)

  # DEA scores for weighting (computed once on the full sample, unscaled data)
  dea_scores <- matrix(NA, nrow = N, ncol = nY)

  for (out in 1:nY) {
    dea_scores[, out] <- rad_out(
      tech_xmat = as.matrix(data[, x]),
      tech_ymat = as.matrix(data[, y[out]]),
      eval_xmat = as.matrix(data[, x]),
      eval_ymat = as.matrix(data[, y[out]]),
      convexity = TRUE,
      returns = "variable"
    )[, 1]
  }

  # initialize importance accumulator
  mat_varimp <- matrix(0, nrow = length(x))
  rownames(mat_varimp) <- colnames(data)[x]
  colnames(mat_varimp) <- "importance"

  for (n in 1:repeats) {
    # for each tree, compute baseline OOB error and permuted OOB error
    tree_baseline_err <- rep(NA, n_trees)
    tree_perm_err <- matrix(NA, nrow = n_trees, ncol = length(x))

    for (t in 1:n_trees) {
      tree <- object[["forest"]][[t]]
      aces_model <- tree[["methods"]][["rf_aces"]]
      knots <- aces_model[["knots"]]
      coefs <- aces_model[["coefs"]]

      # OOB indices
      sample_bag <- tree[["sample_bag"]]
      if (is.null(sample_bag)) {
        oob_idxs <- 1:N
      } else {
        oob_idxs <- which(!(1:N %in% sample_bag))
      }

      if (length(oob_idxs) == 0) next

      # OOB data in expanded format
      oob_data <- data_algo[oob_idxs, , drop = FALSE]
      y_oob <- as.matrix(oob_data[, y_algo, drop = FALSE])

      # DEA weights for OOB observations
      weight_oob <- 1 / dea_scores[oob_idxs, , drop = FALSE]

      # build B matrix for OOB observations
      B_oob <- set_Bmat(
        newdata = oob_data[, x_algo, drop = FALSE],
        model = aces_model,
        knots = knots,
        method = "rf_aces"
      )

      # baseline predictions
      y_hat_oob <- B_oob %*% coefs

      # baseline error for this tree using the configured metric
      tree_baseline_err[t] <- err_metric(
        y_obs = y_oob,
        y_hat = y_hat_oob,
        metric = metric,
        weight = weight_oob
      )

      # permuted predictions for each variable
      for (j in 1:length(x)) {
        # permute the j-th input in OOB data (scaled space)
        data_perm <- data_work[oob_idxs, , drop = FALSE]
        data_perm[, x[j]] <- data_perm[sample(nrow(data_perm)), x[j]]

        # rebuild the expanded data with the permuted variable
        data_perm_algo <- set_data(
          data = data_perm,
          x = x,
          y = y,
          max_degree = max_degree
        )

        # build B matrix with permuted OOB data
        B_perm <- set_Bmat(
          newdata = data_perm_algo[, x_algo, drop = FALSE],
          model = aces_model,
          knots = knots,
          method = "rf_aces"
        )

        # permuted predictions
        y_hat_perm <- B_perm %*% coefs

        # permuted error
        tree_perm_err[t, j] <- err_metric(
          y_obs = y_oob,
          y_hat = y_hat_perm,
          metric = metric,
          weight = weight_oob
        )
      }
    }

    # importance = mean increase in error across trees (per variable)
    for (j in 1:length(x)) {
      diffs <- tree_perm_err[, j] - tree_baseline_err
      diffs[!is.finite(diffs)] <- NA
      mat_varimp[j, 1] <- mat_varimp[j, 1] + mean(diffs, na.rm = TRUE)
    }
  }

  # average over repeats
  mat_varimp[, 1] <- 100 * mat_varimp[, 1] / repeats

  # ranking of variable importance
  var_ord <- order(mat_varimp[, 1], decreasing = TRUE, na.last = TRUE)

  ranking <- data.frame(
    variable = rownames(mat_varimp)[var_ord],
    importance = round(as.numeric(mat_varimp[var_ord, 1]), 2),
    row.names = NULL,
    check.names = FALSE
  )

  if (normalize) {
    imp <- ranking[, "importance"]

    mn <- suppressWarnings(min(imp[is.finite(imp)], na.rm = TRUE))
    if (is.finite(mn) && mn < 0) {
      imp <- imp - mn
    }

    mx <- suppressWarnings(max(imp[is.finite(imp)], na.rm = TRUE))
    if (is.finite(mx) && mx > 0) {
      ranking[, "importance"] <- round(100 * imp / mx, 2)
    } else {
      warning("Normalization skipped: non-positive or non-finite max after shift.")
      ranking[, "importance"] <- round(imp, 2)
    }
  }

  return(ranking)
}

#' @title Predict an RF-ACES Percentile
#'
#' @description
#' Predicts outputs from a chosen percentile of the learner-specific frontier
#' predictions instead of their mean. Lower percentiles provide more
#' conservative predictions, while higher percentiles move toward the upper
#' part of the estimated conditional distribution.
#'
#' @param object
#' A \code{rf_aces} object.
#'
#' @param newdata
#' A \code{data.frame} or \code{matrix} containing the new observations.
#'
#' @param x
#' Column indexes of input variables in \code{newdata}.
#'
#' @param y
#' Column indexes of output variables in \code{newdata}.
#'
#' @param p
#' Percentile to predict, as a number between 0 and 1.
#'
#' @param method
#' RF-ACES model to use:
#' \itemize{
#' \item{\code{"rf_aces"}}: Random Forest Adaptive Constrained Enveloping Splines.
#' \item{\code{"rf_aces_cubic"}}: Cubic-smoothed RF-ACES.
#' \item{\code{"rf_aces_quintic"}}: Quintic-smoothed RF-ACES.
#' }
#'
#' @return
#' A numeric vector or matrix of percentile predictions.

rf_aces_p_predict <- function(
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
    xmat_rf_aces <- rbind(
      xmat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["xmat"]]
    )

    ymat_rf_aces <- rbind(
      ymat_rf_aces,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces"]][["ymat"]]
    )

    xmat_rf_aces_cubic <- rbind(
      xmat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["xmat"]]
    )

    ymat_rf_aces_cubic <- rbind(
      ymat_rf_aces_cubic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_cubic"]][["ymat"]]
    )

    xmat_rf_aces_quintic <- rbind(
      xmat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["xmat"]]
    )

    ymat_rf_aces_quintic <- rbind(
      ymat_rf_aces_quintic,
      RF_ACES[["forest"]][[i]][["technology"]][["rf_aces_quintic"]][["ymat"]]
    )
  }

  data_sorted <- as.data.frame(
    set_data(
      data = newdata,
      x = x,
      y = y,
      max_degree = 1
    )
  )

  inp_cols <- names(data_sorted)[1:length(x)]
  out_cols <- names(data_sorted)[(length(x) + 1):(length(x) + length(y))]

  colnames(xmat_rf_aces) <- colnames(xmat_rf_aces_cubic) <- colnames(xmat_rf_aces_quintic) <- inp_cols
  colnames(ymat_rf_aces) <- colnames(ymat_rf_aces_cubic) <- colnames(ymat_rf_aces_quintic) <- out_cols

  rf_aces_technology <- data.frame(xmat_rf_aces, ymat_rf_aces) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), ~ quantile(.x, probs = p), .names = "{col}"), .groups = "drop")

  rf_aces_cubic_technology <- data.frame(xmat_rf_aces_cubic, ymat_rf_aces_cubic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), ~ quantile(.x, probs = p), .names = "{col}"), .groups = "drop")

  rf_aces_quintic_technology <- data.frame(xmat_rf_aces_quintic, ymat_rf_aces_quintic) %>%
    group_by(across(all_of(inp_cols))) %>%
    summarise(across(all_of(out_cols), ~ quantile(.x, probs = p), .names = "{col}"), .groups = "drop")

  RF_ACES[["technology"]] <- list(
    "rf_aces" = list(
      "xmat" = as.matrix(rf_aces_technology[, x]),
      "ymat" = as.matrix(rf_aces_technology[, y])
    ),
    "rf_aces_cubic" = list(
      "xmat" = as.matrix(rf_aces_cubic_technology[, x]),
      "ymat" = as.matrix(rf_aces_cubic_technology[, y])
    ),
    "rf_aces_quintic" = list(
      "xmat" = as.matrix(rf_aces_quintic_technology[, x]),
      "ymat" = as.matrix(rf_aces_quintic_technology[, y])
    )
  )

  # number of outputs
  nY <- length(RF_ACES[["data"]][["y"]])

  # data in [x, y] format with interaction of variables included
  data <- set_data(
    data = newdata,
    x = x,
    y = NULL,
    max_degree = RF_ACES[["control"]][["max_degree"]]
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
    B <- set_Bmat(
      newdata = data,
      model = aces_model,
      knots = knots,
      method = method
    )

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
  scores <- rad_out(
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
