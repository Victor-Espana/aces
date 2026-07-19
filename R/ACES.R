#' @title Fit an ACES Model
#'
#' @description
#'
#' Estimates maximum attainable outputs with Adaptive Constrained Enveloping
#' Splines (ACES) and uses them to build a production technology. ACES adapts
#' Multivariate Adaptive Regression Splines (MARS) to frontier estimation.
#'
#' @details
#' The forward step builds a flexible spline model for all outputs, and the
#' backward step selects a more parsimonious set of basis functions. When
#' \code{shape$mono} and \code{shape$conc} are both \code{TRUE}, the model for
#' each output is a non-decreasing and concave production function. If either
#' restriction is relaxed, the spline model is an intermediate predictor that
#' moves observed outputs toward their estimated maximum production capacity.
#'
#' ACES always constructs a production technology. For multiple outputs, the
#' predicted components are first refined using the agreement between joint and
#' output-specific radial scores. The original inputs and refined output vectors
#' then define the reference points used for efficiency measurement. Their
#' convex, free-disposal envelopment supplies the final production-set structure
#' under the selected returns-to-scale assumption. See
#' \insertCite{espana2024;textual}{aces} and
#' \insertCite{espana2025;textual}{aces} for methodological details.
#'
#' @name aces
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
#' shorten fitting time for larger problems. The best reduction found for every
#' prepared input in each accepted forward iteration is retained and can be
#' summarized with \code{aces_varimp()}.
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
#' When both are \code{TRUE}, ACES estimates a production function for each
#' output. Otherwise, the spline model is used to estimate maximum attainable
#' outputs before the final production technology is constructed.
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
#' @param kn_penalty
#' Non-negative penalty per knot used to compute generalized cross-validation
#' (GCV) during pruning. Larger values favor models with fewer knots.
#'

#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{espana2025}{aces} \cr \cr
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{zhang1994}{aces} \cr
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc group_by across all_of summarise
#' @importFrom extraDistr rsign
#'
#' @return
#'
#' An object of class \code{aces}. It contains the original data and controls,
#' the forward and pruned models, the cubic- and quintic-smoothed models, and the
#' reference points that define the estimated technology for each available
#' method.
#'
#' @export

aces <- function(
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
  max_terms = 50,
  err_red = 0.01,
  kn_grid = -1,
  minspan = -1,
  endspan = -1,
  kn_penalty = 2
) {

  # =================== #
  # DISPLAY ERRORS ACES #
  # =================== #

  display_errors_aces(
    data = data,
    x = x,
    y = y,
    scale_data = scale_data,
    quick_aces = quick_aces,
    max_degree = mul_BF[["max_degree"]],
    inter_cost = mul_BF[["inter_cost"]],
    metric = metric,
    max_terms = max_terms,
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = kn_penalty
  )

  # =================== #
  #    SCALING SETUP    #
  # =================== #

  # scale factor for inputs
  sx <- rep(1, length(x))

  # scale factor for outputs
  sy <- rep(1, length(y))

  # copy of data for fitting
  data_algo <- data

  if (scale_data) {
    raw_x <- as.matrix(data[, x])
    raw_y <- as.matrix(data[, y])

    # calculate factors (Mean)
    sx <- colMeans(raw_x, na.rm = TRUE)
    sx[sx == 0] <- 1

    sy <- colMeans(raw_y, na.rm = TRUE)
    sy[sy == 0] <- 1

    # apply scaling
    data_algo[, x] <- sweep(raw_x, 2, sx, "/")
    data_algo[, y] <- sweep(raw_y, 2, sy, "/")
  }

  # ========= #
  # ALGORITHM #
  # ========= #

  ACES <- aces_algorithm(
    data = data_algo,
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
    max_terms = max_terms,
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = kn_penalty
  )

  # add scaling options
  ACES[["control"]][["scale"]] <- list(
    is_scaled = scale_data,
    mean_x = sx,
    mean_y = sy
  )

  # Preserve the user-facing fitting specification so diagnostics that require
  # refitting (notably LOCO importance) can reproduce the original setup.
  ACES[["control"]][["fit_spec"]] <- list(
    scale_data = scale_data,
    quick_aces = quick_aces,
    mul_BF = mul_BF,
    metric = metric,
    shape = shape,
    max_terms = max_terms,
    err_red = err_red,
    kn_grid = kn_grid,
    minspan = minspan,
    endspan = endspan,
    kn_penalty = kn_penalty
  )

  # type of object
  class(ACES) <- "aces"

  return(ACES)
}

#' @title Screen Inputs by Rank Correlation
#'
#' @description
#' Applies the correlation filter used by Quick ACES. An input is retained when
#' its Spearman or Kendall correlation reaches the corresponding threshold for
#' at least one output.
#'
#' @details
#' Each threshold is the smaller of 0.1 and the 20th percentile of the positive
#' correlations for that measure. An input is removed only when both measures
#' fall below their thresholds for every output. Undefined correlations, such as
#' those produced by a constant input, do not satisfy a threshold.
#'
#' @param data
#' A matrix containing the prepared inputs followed by the outputs.
#'
#' @param x
#' Column indexes of prepared inputs in \code{data}.
#'
#' @param y
#' Column indexes of outputs in \code{data}.
#'
#' @return
#' A list containing the logical vector \code{keep}, the Spearman and Kendall
#' correlation matrices, and their two filtering thresholds.
#'
#' @keywords internal

quick_aces_correlation_filter <- function(data, x, y) {
  input_data <- data[, x, drop = FALSE]
  output_data <- data[, y, drop = FALSE]

  spearman_corr <- suppressWarnings(stats::cor(
    input_data,
    output_data,
    method = "spearman"
  ))

  kendall_corr <- suppressWarnings(stats::cor(
    input_data,
    output_data,
    method = "kendall"
  ))

  correlation_threshold <- function(correlations) {
    positive <- correlations[
      is.finite(correlations) & correlations > 0
    ]

    if (length(positive) == 0) {
      return(0.1)
    }

    min(
      0.1,
      as.numeric(stats::quantile(positive, probs = 0.2, names = FALSE))
    )
  }

  thresholds <- c(
    spearman = correlation_threshold(spearman_corr),
    kendall = correlation_threshold(kendall_corr)
  )

  keep <- vapply(seq_len(nrow(spearman_corr)), function(j) {
    reaches_threshold <-
      (is.finite(spearman_corr[j, ]) &
        spearman_corr[j, ] >= thresholds[["spearman"]]) |
      (is.finite(kendall_corr[j, ]) &
        kendall_corr[j, ] >= thresholds[["kendall"]])

    any(reaches_threshold)
  }, logical(1))

  names(keep) <- colnames(input_data)

  list(
    keep = keep,
    spearman = spearman_corr,
    kendall = kendall_corr,
    thresholds = thresholds
  )
}

#' @title Map Quick ACES Filtering to Original Inputs
#'
#' @description
#' Returns the original input columns retained by the Quick ACES correlation
#' filter. An original input is retained when either its own term or an
#' interaction containing it passes the filter.
#'
#' @param keep
#' A logical vector indicating which prepared inputs passed the correlation
#' filter. Original inputs must appear first, followed by interaction terms.
#'
#' @param x_vars
#' Column indexes of the original inputs in the unprepared data.
#'
#' @param max_degree
#' Maximum interaction degree, or a list of input-index vectors defining the
#' interactions used to prepare the data.
#'
#' @return
#' An integer vector containing the retained original input-column indexes.
#'
#' @keywords internal

quick_aces_retained_inputs <- function(keep, x_vars, max_degree) {
  n_inputs <- length(x_vars)

  if (length(keep) < n_inputs) {
    stop("'keep' must contain one value for every original input.", call. = FALSE)
  }

  retained <- as.logical(keep[seq_len(n_inputs)])
  n_interactions <- length(keep) - n_inputs

  if (n_interactions == 0) {
    return(x_vars[retained])
  }

  if (is.list(max_degree)) {
    interaction_groups <- lapply(max_degree, function(group) {
      match(group, x_vars)
    })
  } else {
    interaction_groups <- list()

    for (degree in 2:max_degree) {
      interaction_groups <- c(
        interaction_groups,
        utils::combn(seq_len(n_inputs), degree, simplify = FALSE)
      )
    }
  }

  if (length(interaction_groups) != n_interactions ||
      anyNA(unlist(interaction_groups))) {
    stop(
      "Prepared interactions do not match 'max_degree'.",
      call. = FALSE
    )
  }

  retained_interactions <- which(
    keep[n_inputs + seq_len(n_interactions)]
  )

  for (interaction in retained_interactions) {
    retained[interaction_groups[[interaction]]] <- TRUE
  }

  x_vars[retained]
}

.aces_varimp_training_data <- function(object) {
  data <- as.data.frame(object[["data"]][["df"]], check.names = FALSE)

  # RF-ACES keeps the original data in the public object; a plain ACES object
  # keeps the data used by the algorithm and therefore needs to be unscaled.
  if (inherits(object, "rf_aces")) return(data)

  x <- object[["data"]][["x"]]
  y <- object[["data"]][["y"]]
  scaling <- object[["control"]][["scale"]]

  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    data[, x] <- sweep(
      as.matrix(data[, x, drop = FALSE]),
      2,
      scaling$mean_x,
      "*"
    )
    data[, y] <- sweep(
      as.matrix(data[, y, drop = FALSE]),
      2,
      scaling$mean_y,
      "*"
    )
  }

  data
}

.aces_varimp_predict <- function(object, data, x, method) {
  if (inherits(object, "rf_aces")) {
    return(.rf_aces_varimp_predict(object, data, x, method))
  }

  fitted_method <- object[["methods"]][[method]]

  if (is.null(fitted_method)) {
    stop("The fitted object does not contain method '", method, "'.", call. = FALSE)
  }

  working_data <- as.data.frame(data, check.names = FALSE)
  scaling <- object[["control"]][["scale"]]

  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    working_data[, x] <- sweep(
      as.matrix(data[, x, drop = FALSE]),
      2,
      scaling$mean_x[match(x, object[["data"]][["x"]])],
      "/"
    )
  }

  expanded <- set_data(
    data = working_data,
    x = x,
    y = NULL,
    max_degree = object[["control"]][["max_degree"]]
  )

  basis <- set_Bmat(
    newdata = expanded,
    model = fitted_method,
    knots = fitted_method[["knots"]],
    method = method
  )

  estimate <- pmax(0, basis %*% as.matrix(fitted_method[["coefs"]]))
  estimate <- matrix(
    estimate,
    nrow = nrow(data),
    ncol = length(object[["data"]][["y"]])
  )

  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    estimate <- sweep(estimate, 2, scaling$mean_y, "*")
  }

  estimate
}

.rf_aces_varimp_expanded <- function(object, data, x) {
  working_data <- as.data.frame(data, check.names = FALSE)
  scaling <- object[["control"]][["scale"]]

  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    scale_index <- match(x, object[["data"]][["x"]])
    if (anyNA(scale_index)) {
      stop("The prediction inputs do not match the fitted RF-ACES inputs.", call. = FALSE)
    }
    working_data[, x] <- sweep(
      as.matrix(working_data[, x, drop = FALSE]),
      2,
      scaling$mean_x[scale_index],
      "/"
    )
  }

  set_data(
    data = working_data,
    x = x,
    y = NULL,
    max_degree = object[["control"]][["max_degree"]]
  )
}

.rf_aces_varimp_predict <- function(object, data, x, method) {
  expanded <- .rf_aces_varimp_expanded(object, data, x)
  n_outputs <- length(object[["data"]][["y"]])
  estimate <- matrix(0, nrow = nrow(data), ncol = n_outputs)
  n_learners <- length(object[["forest"]])

  if (n_learners == 0L) {
    stop("The RF-ACES object contains no learners.", call. = FALSE)
  }

  for (learner in object[["forest"]]) {
    fitted_method <- learner[["methods"]][[method]]
    if (is.null(fitted_method)) {
      stop("The fitted RF-ACES learners do not contain method '", method, "'.", call. = FALSE)
    }
    basis <- set_Bmat(
      newdata = expanded,
      model = fitted_method,
      knots = fitted_method[["knots"]],
      method = method
    )
    learner_estimate <- pmax(0, basis %*% as.matrix(fitted_method[["coefs"]]))
    estimate <- estimate + matrix(
      learner_estimate,
      nrow = nrow(data),
      ncol = n_outputs
    )
  }

  estimate <- estimate / n_learners
  scaling <- object[["control"]][["scale"]]
  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    estimate <- sweep(estimate, 2, scaling$mean_y, "*")
  }
  estimate
}

.rf_aces_varimp_oob_predict <- function(object, data, x, method) {
  expanded <- .rf_aces_varimp_expanded(object, data, x)
  n_observations <- nrow(data)
  n_outputs <- length(object[["data"]][["y"]])
  prediction_sum <- matrix(0, nrow = n_observations, ncol = n_outputs)
  prediction_count <- integer(n_observations)

  for (learner in object[["forest"]]) {
    sample_bag <- learner[["sample_bag"]]
    if (is.null(sample_bag)) {
      stop(
        "OOB importance is unavailable because bag membership is not stored; refit RF-ACES with the current version.",
        call. = FALSE
      )
    }
    oob <- which(!seq_len(n_observations) %in% unique(sample_bag))
    if (length(oob) == 0L) next

    fitted_method <- learner[["methods"]][[method]]
    if (is.null(fitted_method)) {
      stop("The fitted RF-ACES learners do not contain method '", method, "'.", call. = FALSE)
    }
    basis <- set_Bmat(
      newdata = expanded[oob, , drop = FALSE],
      model = fitted_method,
      knots = fitted_method[["knots"]],
      method = method
    )
    learner_estimate <- pmax(0, basis %*% as.matrix(fitted_method[["coefs"]]))
    prediction_sum[oob, ] <- prediction_sum[oob, , drop = FALSE] + matrix(
      learner_estimate,
      nrow = length(oob),
      ncol = n_outputs
    )
    prediction_count[oob] <- prediction_count[oob] + 1L
  }

  rows <- which(prediction_count > 0L)
  if (length(rows) == 0L) {
    stop("No observation has an OOB prediction.", call. = FALSE)
  }
  prediction <- matrix(NA_real_, nrow = n_observations, ncol = n_outputs)
  prediction[rows, ] <- prediction_sum[rows, , drop = FALSE] / prediction_count[rows]

  scaling <- object[["control"]][["scale"]]
  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    prediction[rows, ] <- sweep(
      prediction[rows, , drop = FALSE],
      2,
      scaling$mean_y,
      "*"
    )
  }

  list(prediction = prediction, rows = rows, counts = prediction_count)
}

.rf_aces_forward_varimp <- function(object, normalize, type) {
  histories <- lapply(object[["forest"]], function(learner) {
    learner[["methods"]][["rf_aces"]][["variable_reduction"]]
  })
  if (length(histories) == 0L || any(vapply(histories, is.null, logical(1)))) {
    stop(
      "Forward-search reductions are not stored in every learner; refit RF-ACES with the current version.",
      call. = FALSE
    )
  }

  variable_names <- unique(unlist(lapply(histories, colnames), use.names = FALSE))
  if (length(variable_names) == 0L) {
    stop("Stored RF-ACES reduction histories have no input names.", call. = FALSE)
  }

  align_history <- function(history) {
    history <- as.matrix(history)
    storage.mode(history) <- "double"
    aligned <- matrix(
      NA_real_,
      nrow = nrow(history),
      ncol = length(variable_names),
      dimnames = list(rownames(history), variable_names)
    )
    common <- intersect(colnames(history), variable_names)
    aligned[, common] <- history[, common, drop = FALSE]
    aligned
  }
  histories <- lapply(histories, align_history)

  if (type == "iterations") {
    long_history <- lapply(seq_along(histories), function(i) {
      history <- as.data.frame(histories[[i]], check.names = FALSE)
      data.frame(
        learner = i,
        iteration = seq_len(nrow(history)),
        history,
        check.names = FALSE,
        row.names = NULL
      )
    })
    return(do.call(rbind, long_history))
  }

  learner_scores <- vapply(
    histories,
    function(history) colSums(history, na.rm = TRUE),
    numeric(length(variable_names))
  )
  if (is.null(dim(learner_scores))) {
    learner_scores <- matrix(
      learner_scores,
      nrow = length(variable_names),
      ncol = length(histories)
    )
  }
  score <- rowMeans(learner_scores)
  iterations_evaluated <- Reduce(
    `+`,
    lapply(histories, function(history) colSums(!is.na(history)))
  )
  learners_evaluated <- Reduce(
    `+`,
    lapply(histories, function(history) as.integer(colSums(!is.na(history)) > 0L))
  )

  if (normalize) {
    maximum <- suppressWarnings(max(score, na.rm = TRUE))
    score <- if (is.finite(maximum) && maximum > 0) 100 * score / maximum else score * 0
  }

  ranking <- data.frame(
    variable = variable_names,
    importance = round(as.numeric(score), if (normalize) 2 else 6),
    learners_evaluated = as.integer(learners_evaluated),
    iterations_evaluated = as.integer(iterations_evaluated),
    check.names = FALSE,
    row.names = NULL
  )
  ranking <- ranking[order(ranking$importance, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  rownames(ranking) <- NULL
  ranking
}

.aces_varimp_training_weights <- function(data, x, y) {
  scores <- matrix(NA_real_, nrow = nrow(data), ncol = length(y))

  for (out in seq_along(y)) {
    scores[, out] <- rad_out(
      tech_xmat = as.matrix(data[, x, drop = FALSE]),
      tech_ymat = as.matrix(data[, y[out], drop = FALSE]),
      eval_xmat = as.matrix(data[, x, drop = FALSE]),
      eval_ymat = as.matrix(data[, y[out], drop = FALSE]),
      convexity = TRUE,
      returns = "variable",
      type = "objective"
    )[, 1]
  }

  weights <- 1 / scores
  if (any(!is.finite(weights))) {
    stop("Non-finite training weights prevent variable-importance calculation.", call. = FALSE)
  }

  weights
}

.aces_varimp_comparison_ranking <- function(
  variable_names,
  baseline_loss,
  comparison_loss,
  normalize
) {
  if (length(baseline_loss) == 1L) {
    baseline_loss <- rep(baseline_loss, length(variable_names))
  }
  if (length(baseline_loss) != length(variable_names)) {
    stop("Internal error: baseline losses do not match the inputs.", call. = FALSE)
  }
  loss_increase <- comparison_loss - baseline_loss
  importance <- loss_increase
  importance[!is.finite(importance)] <- NA_real_

  if (normalize) {
    importance <- pmax(importance, 0)
    maximum <- suppressWarnings(max(importance, na.rm = TRUE))

    if (is.finite(maximum) && maximum > 0) {
      importance <- 100 * importance / maximum
    } else {
      importance[is.finite(importance)] <- 0
    }
  }

  ranking <- data.frame(
    variable = variable_names,
    importance = round(as.numeric(importance), if (normalize) 2 else 6),
    loss_increase = round(as.numeric(loss_increase), 6),
    baseline_loss = round(as.numeric(baseline_loss), 6),
    comparison_loss = round(as.numeric(comparison_loss), 6),
    check.names = FALSE,
    row.names = NULL
  )

  ranking <- ranking[
    order(ranking$importance, decreasing = TRUE, na.last = TRUE),
    ,
    drop = FALSE
  ]
  rownames(ranking) <- NULL
  ranking
}

.aces_loco_spec <- function(object, dropped_input) {
  original_x <- object[["data"]][["x"]]
  retained_x <- original_x[original_x != dropped_input]

  if (length(retained_x) == 0) {
    stop("LOCO cannot remove the only input in the model.", call. = FALSE)
  }

  spec <- object[["control"]][["fit_spec"]]
  used_fallback <- is.null(spec)

  if (used_fallback) {
    scaling <- object[["control"]][["scale"]]
    common_spec <- list(
      scale_data = !is.null(scaling) && isTRUE(scaling$is_scaled),
      quick_aces = object[["control"]][["quick_aces"]],
      mul_BF = list(
        max_degree = object[["control"]][["max_degree"]],
        inter_cost = object[["control"]][["inter_cost"]]
      ),
      metric = object[["control"]][["metric"]],
      shape = object[["control"]][["shape"]],
      max_terms = object[["control"]][["max_terms"]],
      err_red = object[["control"]][["err_red"]],
      kn_grid = -1,
      minspan = object[["control"]][["minspan"]],
      endspan = object[["control"]][["endspan"]]
    )
    if (inherits(object, "rf_aces")) {
      spec <- c(
        common_spec,
        list(
          learners = object[["control"]][["learners"]],
          bag_size = object[["control"]][["bag_size"]],
          max_feats = object[["control"]][["max_feats"]],
          early_stopping = object[["control"]][["early_stopping"]]
        )
      )
    } else {
      spec <- c(
        common_spec,
        list(kn_penalty = object[["control"]][["kn_penalty"]])
      )
    }
  }

  max_degree <- spec$mul_BF$max_degree
  if (is.list(max_degree)) {
    max_degree <- Filter(
      function(group) !dropped_input %in% group && all(group %in% retained_x),
      max_degree
    )
    if (length(max_degree) == 0) max_degree <- 1L
  } else {
    max_degree <- max(1L, min(max_degree, length(retained_x)))
  }
  spec$mul_BF$max_degree <- max_degree

  if (is.list(spec$kn_grid)) {
    if (length(spec$kn_grid) == length(original_x)) {
      spec$kn_grid <- spec$kn_grid[original_x != dropped_input]
    } else {
      spec$kn_grid <- -1
      used_fallback <- TRUE
    }
  }

  if (inherits(object, "rf_aces")) {
    spec$max_feats <- max(1L, min(spec$max_feats, length(retained_x)))
  }

  list(x = retained_x, spec = spec, used_fallback = used_fallback)
}

.aces_varimp_control <- function(importance, control, model_family = "aces") {
  if (!is.list(control)) {
    stop("'control' must be a named list.", call. = FALSE)
  }

  default_method <- if (identical(model_family, "rf_aces")) "rf_aces" else "aces"
  default_evaluation <- if (identical(model_family, "rf_aces")) "oob" else NULL

  defaults <- switch(
    importance,
    forward = list(type = "overall"),
    permutation = list(
      method = default_method,
      repeats = 10L,
      eval_data = default_evaluation,
      seed = NULL
    ),
    loco = list(
      method = default_method,
      eval_data = default_evaluation,
      verbose = FALSE
    )
  )

  if (length(control) == 0) return(defaults)

  control_names <- names(control)
  if (is.null(control_names) || any(control_names == "") ||
      anyDuplicated(control_names)) {
    stop("Every element of 'control' must have a unique name.", call. = FALSE)
  }

  unknown <- setdiff(control_names, names(defaults))
  if (length(unknown) > 0) {
    stop(
      "Unknown control option for importance = '", importance, "': ",
      paste(unknown, collapse = ", "), ". Allowed options are: ",
      paste(names(defaults), collapse = ", "), ".",
      call. = FALSE
    )
  }

  utils::modifyList(defaults, control, keep.null = TRUE)
}

#' @title Compute ACES and RF-ACES Variable Importance
#'
#' @description
#' Computes ACES or RF-ACES variable importance from forward-search reductions,
#' input permutation, or leave-one-covariate-out (LOCO) refitting.
#'
#' @details
#' With \code{importance = "forward"}, at every accepted forward iteration each
#' ACES learner retains the largest relative training-loss reduction obtained by
#' any candidate pair associated with each prepared input. The raw ACES score is
#' the sum of those reductions. For RF-ACES, the sums are computed within each
#' learner and averaged across all learners; an input not evaluated by a learner
#' contributes zero to that learner. This calculation does not refit the model.
#'
#' With \code{importance = "permutation"}, each original input is permuted while
#' the fitted model is held fixed. Prepared interactions are rebuilt, so the
#' result captures the main effect and every interaction involving that input.
#' The raw score is the mean increase in the fitted loss across permutations.
#'
#' With \code{importance = "loco"}, ACES or the complete RF-ACES forest is
#' refitted after removing each original input and every interaction that
#' contains it. The raw score is the increase in loss relative to the complete
#' model. LOCO therefore performs one additional model fit per original input
#' and can be computationally expensive, especially for RF-ACES. RF-ACES LOCO
#' refits are randomized; use enough learners and call \code{set.seed()} before
#' \code{aces_varimp()} when reproducibility is required.
#'
#' For ACES, permutation and LOCO use the original training observations when
#' \code{control$eval_data = NULL}. For RF-ACES, \code{eval_data} is either an
#' external data set or \code{"oob"}, its default. An OOB prediction for an
#' observation is the
#' mean only of the learners that did not use that observation in their bootstrap
#' sample. In OOB LOCO, each full-versus-reduced comparison uses the intersection
#' of observations with predictions from both forests. Supplying a data frame or
#' matrix evaluates either model on external data; LOCO models are still fitted
#' only to the original training observations. Training and OOB losses use the
#' fitted loss metric and DEA-based training weights. External losses use equal
#' observation weights.
#'
#' Normalized permutation and LOCO scores set negative loss increases to zero
#' and scale the largest remaining value to 100. Raw scores retain their sign
#' and are meaningful only for a common loss metric and output scale. None of
#' the three approaches measures statistical significance or causality.
#'
#' @param object
#' A fitted \code{aces} or \code{rf_aces} object.
#'
#' @param normalize
#' If \code{TRUE}, scale the largest importance to 100. If \code{FALSE}, return
#' the raw cumulative reduction for forward importance or the signed loss
#' increase for permutation and LOCO.
#'
#' @param importance
#' Importance approach: \code{"forward"}, \code{"permutation"}, or
#' \code{"loco"}.
#'
#' @param control
#' A named list containing only the options used by the selected approach:
#' \itemize{
#'   \item Forward: \code{type}, either \code{"overall"} or
#'   \code{"iterations"}.
#'   \item Permutation: \code{method}, \code{repeats}, \code{eval_data}, and
#'   \code{seed}.
#'   \item LOCO: \code{method}, \code{eval_data}, and \code{verbose}.
#' }
#' For ACES, the available fitted methods are \code{"aces_forward"},
#' \code{"aces"}, \code{"aces_cubic"}, and \code{"aces_quintic"}; the default
#' is \code{"aces"}. For RF-ACES they are \code{"rf_aces"},
#' \code{"rf_aces_cubic"}, and \code{"rf_aces_quintic"}; the default is
#' \code{"rf_aces"}. For ACES, \code{eval_data} may be a data frame or matrix,
#' or \code{NULL} for training evaluation. For RF-ACES it may be a data frame or
#' matrix, or \code{"oob"}. An unknown option produces an error instead of being
#' silently ignored.
#'
#' @return
#' Forward importance returns the ranking or iteration history requested in
#' \code{control$type}. Permutation and LOCO return a data frame with the original input
#' name, importance, signed loss increase, baseline loss, and comparison loss.
#'
#' @export

aces_varimp <- function(
  object,
  normalize = TRUE,
  importance = c("forward", "permutation", "loco"),
  control = list()
) {
  if (!inherits(object, "aces") && !inherits(object, "rf_aces")) {
    stop("'object' must be an aces or rf_aces object.", call. = FALSE)
  }

  if (!is.logical(normalize) || length(normalize) != 1 || is.na(normalize)) {
    stop("'normalize' must be TRUE or FALSE.", call. = FALSE)
  }

  model_family <- if (inherits(object, "rf_aces")) "rf_aces" else "aces"
  importance <- match.arg(importance)
  control <- .aces_varimp_control(importance, control, model_family)

  if (importance == "forward") {
    type <- match.arg(control$type, c("overall", "iterations"))

    if (model_family == "rf_aces") {
      return(.rf_aces_forward_varimp(object, normalize, type))
    }

    reduction <- object[["methods"]][["aces_forward"]][["variable_reduction"]]

    if (is.null(reduction)) {
      stop(
        "Forward-search reductions are not stored in this object; refit it with the current version of aces.",
        call. = FALSE
      )
    }

    reduction <- as.matrix(reduction)
    storage.mode(reduction) <- "double"

    if (type == "iterations") {
      history <- as.data.frame(reduction, check.names = FALSE)
      return(data.frame(
        iteration = seq_len(nrow(reduction)),
        history,
        check.names = FALSE,
        row.names = NULL
      ))
    }

    score <- colSums(reduction, na.rm = TRUE)
    iterations_evaluated <- colSums(!is.na(reduction))

    if (normalize) {
      maximum <- suppressWarnings(max(score, na.rm = TRUE))
      if (is.finite(maximum) && maximum > 0) {
        score <- 100 * score / maximum
      } else {
        score[] <- 0
      }
    }

    variable_names <- colnames(reduction)
    if (is.null(variable_names)) {
      variable_names <- paste0("x", seq_along(score))
    }

    ranking <- data.frame(
      variable = variable_names,
      importance = round(as.numeric(score), if (normalize) 2 else 6),
      iterations_evaluated = as.integer(iterations_evaluated),
      check.names = FALSE,
      row.names = NULL
    )

    ranking <- ranking[
      order(ranking$importance, decreasing = TRUE, na.last = TRUE),
      ,
      drop = FALSE
    ]
    rownames(ranking) <- NULL
    return(ranking)
  }

  available_methods <- if (model_family == "rf_aces") {
    c("rf_aces", "rf_aces_cubic", "rf_aces_quintic")
  } else {
    c("aces_forward", "aces", "aces_cubic", "aces_quintic")
  }
  method <- match.arg(control$method, available_methods)
  eval_data <- control$eval_data

  evaluation_is_oob <- is.character(eval_data) && length(eval_data) == 1L &&
    !is.na(eval_data) && identical(tolower(eval_data), "oob")
  if (is.character(eval_data) && !evaluation_is_oob) {
    stop("Character 'control$eval_data' must be exactly 'oob'.", call. = FALSE)
  }
  if (evaluation_is_oob && model_family != "rf_aces") {
    stop("'control$eval_data = \"oob\"' is available only for RF-ACES.", call. = FALSE)
  }
  if (model_family == "rf_aces" && is.null(eval_data)) {
    stop(
      "For RF-ACES, 'control$eval_data' must be 'oob' or a data frame or matrix.",
      call. = FALSE
    )
  }

  if (importance == "permutation") {
    repeats <- control$repeats
    seed <- control$seed

    if (!is.numeric(repeats) || length(repeats) != 1 || is.na(repeats) ||
        repeats < 1 || repeats != floor(repeats)) {
      stop("'control$repeats' must be a positive integer.", call. = FALSE)
    }

    if (!is.null(seed)) {
      if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) ||
          seed != floor(seed)) {
        stop("'control$seed' must be NULL or a single integer.", call. = FALSE)
      }
      had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (had_random_seed) {
        previous_random_seed <- get(".Random.seed", envir = .GlobalEnv)
      }

      on.exit({
        if (had_random_seed) {
          assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)

      set.seed(seed)
    }
  } else {
    verbose <- control$verbose
    if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
      stop("'control$verbose' must be TRUE or FALSE.", call. = FALSE)
    }
  }

  training_data <- .aces_varimp_training_data(object)
  x <- object[["data"]][["x"]]
  y <- object[["data"]][["y"]]
  variable_names <- object[["data"]][["xnames"]]
  metric <- object[["control"]][["metric"]]

  if (evaluation_is_oob) {
    truth <- as.matrix(training_data[, y, drop = FALSE])
    weights <- .aces_varimp_training_weights(training_data, x, y)
    baseline_oob <- .rf_aces_varimp_oob_predict(
      object = object,
      data = training_data,
      x = x,
      method = method
    )

    if (importance == "permutation") {
      baseline_rows <- baseline_oob$rows
      baseline_loss <- err_metric(
        y_obs = truth[baseline_rows, , drop = FALSE],
        y_hat = baseline_oob$prediction[baseline_rows, , drop = FALSE],
        metric = metric,
        weight = weights[baseline_rows, , drop = FALSE]
      )
      permuted_loss <- matrix(NA_real_, nrow = repeats, ncol = length(x))

      for (r in seq_len(repeats)) {
        for (j in seq_along(x)) {
          permuted_data <- training_data
          permuted_data[, x[j]] <- permuted_data[
            sample.int(nrow(permuted_data)),
            x[j]
          ]
          permuted_oob <- .rf_aces_varimp_oob_predict(
            object = object,
            data = permuted_data,
            x = x,
            method = method
          )
          permuted_loss[r, j] <- err_metric(
            y_obs = truth[baseline_rows, , drop = FALSE],
            y_hat = permuted_oob$prediction[baseline_rows, , drop = FALSE],
            metric = metric,
            weight = weights[baseline_rows, , drop = FALSE]
          )
        }
      }
      comparison_loss <- colMeans(permuted_loss, na.rm = TRUE)
    } else {
      baseline_loss <- comparison_loss <- rep(NA_real_, length(x))
      fallback_used <- FALSE

      for (j in seq_along(x)) {
        loco <- .aces_loco_spec(object, dropped_input = x[j])
        fallback_used <- fallback_used || loco$used_fallback
        refit_args <- c(
          list(data = training_data, x = loco$x, y = y),
          loco$spec
        )
        if (verbose) {
          reduced_fit <- do.call(rf_aces, refit_args)
        } else {
          invisible(utils::capture.output(
            reduced_fit <- do.call(rf_aces, refit_args)
          ))
        }
        reduced_oob <- .rf_aces_varimp_oob_predict(
          object = reduced_fit,
          data = training_data,
          x = loco$x,
          method = method
        )
        common_rows <- intersect(baseline_oob$rows, reduced_oob$rows)
        if (length(common_rows) == 0L) next

        baseline_loss[j] <- err_metric(
          y_obs = truth[common_rows, , drop = FALSE],
          y_hat = baseline_oob$prediction[common_rows, , drop = FALSE],
          metric = metric,
          weight = weights[common_rows, , drop = FALSE]
        )
        comparison_loss[j] <- err_metric(
          y_obs = truth[common_rows, , drop = FALSE],
          y_hat = reduced_oob$prediction[common_rows, , drop = FALSE],
          metric = metric,
          weight = weights[common_rows, , drop = FALSE]
        )
      }

      if (fallback_used) {
        warning(
          "The original fitting specification was incomplete; at least one LOCO refit used an automatic knot grid.",
          call. = FALSE
        )
      }
    }

    ranking <- .aces_varimp_comparison_ranking(
      variable_names = variable_names,
      baseline_loss = baseline_loss,
      comparison_loss = comparison_loss,
      normalize = normalize
    )
    attr(ranking, "importance") <- importance
    attr(ranking, "evaluation") <- "oob"
    attr(ranking, "method") <- method
    return(ranking)
  }

  evaluation_is_training <- is.null(eval_data)
  if (evaluation_is_training) {
    evaluation_data <- training_data
    weights <- .aces_varimp_training_weights(training_data, x, y)
  } else {
    evaluation_data <- as.data.frame(eval_data, check.names = FALSE)
    if (ncol(evaluation_data) < max(c(x, y))) {
      stop("'eval_data' does not contain all model columns.", call. = FALSE)
    }
    weights <- matrix(1, nrow = nrow(evaluation_data), ncol = length(y))
  }

  truth <- as.matrix(evaluation_data[, y, drop = FALSE])
  baseline_prediction <- .aces_varimp_predict(
    object = object,
    data = evaluation_data,
    x = x,
    method = method
  )
  baseline_loss <- err_metric(
    y_obs = truth,
    y_hat = baseline_prediction,
    metric = metric,
    weight = weights
  )

  comparison_loss <- rep(NA_real_, length(x))

  if (importance == "permutation") {
    if (nrow(evaluation_data) < 2) {
      stop("Permutation importance requires at least two evaluation rows.", call. = FALSE)
    }

    permuted_loss <- matrix(NA_real_, nrow = repeats, ncol = length(x))
    for (r in seq_len(repeats)) {
      for (j in seq_along(x)) {
        permuted_data <- evaluation_data
        permuted_data[, x[j]] <- permuted_data[sample.int(nrow(permuted_data)), x[j]]
        permuted_prediction <- .aces_varimp_predict(
          object = object,
          data = permuted_data,
          x = x,
          method = method
        )
        permuted_loss[r, j] <- err_metric(
          y_obs = truth,
          y_hat = permuted_prediction,
          metric = metric,
          weight = weights
        )
      }
    }
    comparison_loss <- colMeans(permuted_loss, na.rm = TRUE)
  } else {
    fallback_used <- FALSE

    for (j in seq_along(x)) {
      loco <- .aces_loco_spec(object, dropped_input = x[j])
      fallback_used <- fallback_used || loco$used_fallback
      refit_args <- c(
        list(data = training_data, x = loco$x, y = y),
        loco$spec
      )

      fit_function <- if (model_family == "rf_aces") rf_aces else aces
      if (verbose) {
        reduced_fit <- do.call(fit_function, refit_args)
      } else {
        invisible(utils::capture.output(
          reduced_fit <- do.call(fit_function, refit_args)
        ))
      }

      reduced_prediction <- .aces_varimp_predict(
        object = reduced_fit,
        data = evaluation_data,
        x = loco$x,
        method = method
      )
      comparison_loss[j] <- err_metric(
        y_obs = truth,
        y_hat = reduced_prediction,
        metric = metric,
        weight = weights
      )
    }

    if (fallback_used) {
      warning(
        "The original fitting specification was incomplete; at least one LOCO refit used an automatic knot grid.",
        call. = FALSE
      )
    }
  }

  ranking <- .aces_varimp_comparison_ranking(
    variable_names = variable_names,
    baseline_loss = baseline_loss,
    comparison_loss = comparison_loss,
    normalize = normalize
  )
  attr(ranking, "importance") <- importance
  attr(ranking, "evaluation") <- if (evaluation_is_training) "training" else "external"
  attr(ranking, "method") <- method
  ranking
}

#' @title Run the ACES Algorithm
#'
#' @description
#' Runs the complete internal ACES procedure: data expansion, knot selection,
#' estimation of maximum attainable outputs, backward pruning, smoothing,
#' output refinement, and construction of the production technology. This
#' lower-level function expects validated and, when requested, scaled data.
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
#' each output; otherwise, use the spline model to estimate maximum attainable
#' outputs before constructing the technology.
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
#' @param kn_penalty
#' Penalty per knot used to compute GCV.
#'
#' @importFrom dplyr desc
#'
#' @return
#'
#' An \code{aces} object.
#'
#' @export

aces_algorithm <- function(
  data,
  x_vars,
  y_vars,
  quick_aces,
  max_degree,
  inter_cost,
  metric,
  shape,
  max_terms,
  err_red,
  minspan,
  endspan,
  kn_grid,
  kn_penalty
) {
  # save a copy of the original data
  DMUs <- data

  # variable names
  var_names <- colnames(DMUs[, c(x_vars, y_vars)])

  # data in [x, y] format with interaction of variables included
  data <- set_data(
    data = data,
    x = x_vars,
    y = y_vars,
    max_degree = max_degree
  )

  # samples size
  N <- nrow(data)

  # set "x" and "y" indexes in data
  x <- 1:(ncol(data) - length(y_vars))
  y <- (length(x) + 1):ncol(data)

  # set number of inputs and outputs
  nX <- length(x)
  nY <- length(y)

  # ================== #
  # VARIABLE FILTERING #
  # ================== #

  # Best candidate reduction for each prepared input in the current iteration.
  # NA means that the input has not been evaluated in that iteration.
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
    returns = "variable",
    type = "objective"
  )[, 1]

  for (out in 1:nY) {
    table_scores[, 1 + out] <- rad_out(
      tech_xmat = as.matrix(DMUs[, x_filtered]),
      tech_ymat = as.matrix(DMUs[, y_vars[out]]),
      eval_xmat = as.matrix(DMUs[, x_filtered]),
      eval_ymat = as.matrix(DMUs[, y_vars[out]]),
      convexity = TRUE,
      returns = "variable",
      type = "objective"
    )[, 1]
  }

  # weights for error metrics based on DEA
  dea_scores <- table_scores[, 2:ncol(table_scores), drop = FALSE]

  # ========== #
  # FDH SCORES #
  # ========== #

  table_scores[, 1] <- rad_out(
    tech_xmat = as.matrix(DMUs[, x_filtered]),
    tech_ymat = as.matrix(DMUs[, y_vars]),
    eval_xmat = as.matrix(DMUs[, x_filtered]),
    eval_ymat = as.matrix(DMUs[, y_vars]),
    convexity = FALSE,
    returns = "variable",
    type = "objective"
  )[, 1]

  for (out in 1:nY) {
    table_scores[, 1 + out] <- rad_out(
      tech_xmat = as.matrix(DMUs[, x_filtered]),
      tech_ymat = as.matrix(DMUs[, y_vars[out]]),
      eval_xmat = as.matrix(DMUs[, x_filtered]),
      eval_ymat = as.matrix(DMUs[, y_vars[out]]),
      convexity = FALSE,
      returns = "variable",
      type = "objective"
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
  aces_forward <- list(
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

  # list to save technologies created through ACES
  technology <- list()

  # initial error
  err_min <- err

  # iteration counter
  iter <- 0

  while (length(aces_forward[["bf_set"]]) + 2 < max_terms) {
    # add 2 new basis functions to the model:
    B_bf_knt_err <- add_basis_function(
      data = data,
      x = x,
      y = y,
      xi_degree = xi_degree,
      inter_cost = inter_cost,
      dea_scores = dea_scores,
      fdh_scores = fdh_scores,
      metric = metric,
      forward_model = aces_forward,
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

    if (!is.list(B_bf_knt_err)) break

    # new best error
    new_err <- B_bf_knt_err[[5]]

    # update model
    if (new_err[1] < err[1] * (1 - err_red[1])) {
      # update iteration counter
      iter <- iter + 1

      # compute relative reduction
      rel_reduction <- (err[1] - new_err[1]) / err[1]

      # update B
      aces_forward[["B"]] <- B_bf_knt_err[[1]]

      # update basis functions
      aces_forward[["bf_set"]] <- B_bf_knt_err[[2]]

      # update the knots list
      kn_list <- B_bf_knt_err[[3]]

      # update the Bp list
      Bp_list <- B_bf_knt_err[[4]]

      # updated error
      err <- new_err

      # updated variable importance matrix
      var_imp <- B_bf_knt_err[[6]]
      var_imp <- rbind(var_imp, rep(NA_real_, nX))

      cat(
        sprintf(
          paste0(
            "Iteration %d completed — error reduced by %.1f%% ",
            "(%.2f -> %.2f).\n"
          ),
          iter,
          100 * rel_reduction,
          err[1] / (1 - rel_reduction),
          err[1]
        )
      )
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
  # forward aces
  # ==

  variable_reduction <- if (iter > 0) {
    var_imp[seq_len(iter), , drop = FALSE]
  } else {
    var_imp[FALSE, , drop = FALSE]
  }

  rownames(variable_reduction) <- paste0(
    "iteration_",
    seq_len(nrow(variable_reduction))
  )

  aces_forward <- list(
    "basis" = aces_forward[["bf_set"]],
    "Bmatx" = aces_forward[["B"]],
    "knots" = kn_forward,
    "coefs" = rev(aces_forward[["bf_set"]])[[1]][["coefs"]],
    "variable_reduction" = variable_reduction
  )

  # generate technology
  technology[["aces_forward"]] <- generate_technology(
    var_names = var_names,
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = aces_forward[["Bmatx"]] %*% aces_forward[["coefs"]],
    table_scores = table_scores
  )

  # ====================== #
  #   BACKWARD ALGORITHM   #
  # ====================== #

  aces_submodels <- aces_pruning(
    data = data,
    x = x,
    y = y,
    xi_degree = xi_degree,
    dea_scores = dea_scores,
    fdh_scores = fdh_scores,
    metric = metric,
    forward_model = aces_forward,
    Bp_list = Bp_list,
    shape = shape,
    kn_penalty = kn_penalty
  )

  # generalized cross-validation for each model
  GCVs <- sapply(aces_submodels, function(x) x[["GCV"]])

  # model with minimum error (excluding the model without knots)
  aces_backward <- aces_submodels[[which.min(GCVs[1:(length(aces_submodels) - 1)])]]

  # set of surviving knots
  kn_backward <- do.call(rbind.data.frame, aces_backward[["t"]])

  # sort the knots by "xi", "status", "side", "t"
  knots_backward_order <- with(
    kn_backward,
    {
      order(
        xi,
        status,
        ifelse(status == "paired", t, desc(side)),
        ifelse(status == "paired", desc(side), t)
      )
    }
  )

  # update the set of knots
  kn_backward <- kn_backward[knots_backward_order, ]

  # ==
  # aces
  # ==

  aces <- list(
    "aces_submodels" = aces_submodels,
    "Bmatx" = aces_backward[["B"]],
    "knots" = kn_backward,
    "coefs" = aces_backward[["coefs"]],
    "GCV" = aces_backward[["GCV"]]
  )

  # generate technology
  technology[["aces"]] <- generate_technology(
    var_names = var_names,
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = aces[["Bmatx"]] %*% aces[["coefs"]],
    table_scores = table_scores
  )

  # ======================= #
  #   SMOOTHING PROCEDURE   #
  # ======================= #

  # sub-models sorted by gcv value
  aces_submodels_gcv <- aces_submodels[order(sapply(aces_submodels, "[[", "GCV"))]

  # initialize a list with smoothed sub-models
  aces_smoothed_submodels <- vector("list", length(aces_submodels_gcv))

  # distance between knots in smooth models
  wc <- seq(1, 2, length.out = 5)
  wq <- seq(8 / 7, 1.5, length.out = 5)

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
      kn_side_loc <- side_knot_location(
        data = data,
        nX = nX,
        knots = kn_smoothed
      )

      # ==
      # smoothing cubic aces
      # ==

      cubic_aces_models[[m]] <- cubic_aces(
        data = data,
        x = x,
        y = y,
        dea_scores = dea_scores,
        fdh_scores = fdh_scores,
        metric = metric,
        shape = shape,
        kn_grid = kn_smoothed,
        kn_side_loc = kn_side_loc,
        kn_penalty = kn_penalty,
        xi_degree = xi_degree,
        wc = wc
      )

      # ==
      # smoothing quintic aces
      # ==

      quintic_aces_models[[m]] <- quintic_aces(
        data = data,
        x = x,
        y = y,
        dea_scores = dea_scores,
        fdh_scores = fdh_scores,
        metric = metric,
        shape = shape,
        kn_grid = kn_smoothed,
        kn_side_loc = kn_side_loc,
        kn_penalty = kn_penalty,
        xi_degree = xi_degree,
        wq = wq
      )
    }
  }

  # GCVs of cubic models
  aces_cubic_gcvs <- sapply(
    cubic_aces_models,
    function(x) {
      ifelse(
        is.null(x[["GCV"]]),
        Inf,
        x[["GCV"]]
      )
    }
  )

  min_gcv <- which.min(aces_cubic_gcvs)

  # cubic aces
  aces_cubic <- cubic_aces_models[[min_gcv]]

  # generate technology
  technology[["aces_cubic"]] <- generate_technology(
    var_names = var_names,
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = aces_cubic[["Bmatx"]] %*% aces_cubic[["coefs"]],
    table_scores = table_scores
  )

  # GCVs of quintic models
  aces_quintic_gcvs <- sapply(
    quintic_aces_models,
    function(x) {
      ifelse(
        is.null(x[["GCV"]]),
        Inf,
        x[["GCV"]]
      )
    }
  )

  min_gcv <- which.min(aces_quintic_gcvs)

  # quintic aces
  aces_quintic <- quintic_aces_models[[min_gcv]]

  # generate technology
  technology[["aces_quintic"]] <- generate_technology(
    var_names = var_names,
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = aces_quintic[["Bmatx"]] %*% aces_quintic[["coefs"]],
    table_scores = table_scores
  )

  # =========== #
  # ACES OBJECT #
  # =========== #

  ACES <- aces_object(
    data = DMUs,
    x = x_vars,
    y = y_vars,
    quick_aces = quick_aces,
    max_degree = max_degree,
    inter_cost = inter_cost,
    xi_degree = xi_degree,
    metric = metric,
    shape = shape,
    max_terms = ncol(aces_forward[["Bmatx"]]),
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = kn_penalty,
    psi = 0.05,
    wc = aces_cubic[["w"]],
    wq = aces_quintic[["w"]],
    aces_forward = aces_forward,
    aces = aces,
    aces_cubic = aces_cubic,
    aces_quintic = aces_quintic,
    technology = technology
  )

  return(ACES)
}

#' @title Create an ACES Object
#'
#' @description
#'
#' Collects the fitted models, tuning controls, training data, and estimated
#' technologies in the common structure returned by \code{aces}.
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
#' @param quick_aces
#' If \code{TRUE}, Quick ACES was used.
#'
#' @param max_degree
#' Maximum interaction degree or explicit interaction list.
#'
#' @param inter_cost
#' Minimum relative improvement required to add a higher-degree basis function.
#'
#' @param xi_degree
#' A matrix that records the degree of each input in every expanded variable.
#'
#' @param metric
#' Character string specifying the lack-of-fit measure.
#'
#' @param shape
#' Logical settings used during spline estimation. When both \code{mono} and
#' \code{conc} are \code{TRUE}, the stored output-specific models are production
#' functions; otherwise, they are intermediate predictors of maximum attainable
#' output.
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
#' Knot grid used to fit the model.
#'
#' @param kn_penalty
#' Penalty per knot used to compute GCV.
#'
#' @param psi
#' Threshold used to refine multi-output predictions. Predicted components are
#' retained when the joint and output-specific radial scores differ by no more
#' than this value; otherwise, observed components are used.
#'
#' @param wc
#' Side-knot distance parameter for cubic smoothing.
#'
#' @param wq
#' Side-knot distance parameter for quintic smoothing.
#'
#' @param aces_forward
#' Results from the ACES forward step.
#'
#' @param aces
#' Results from the pruned ACES model.
#'
#' @param aces_cubic
#' Results from the cubic-smoothed ACES model.
#'
#' @param aces_quintic
#' Results from the quintic-smoothed ACES model.
#'
#' @param technology
#' Reference points that define the production technology for each fitted model.
#'
#' @return
#'
#' An \code{aces} object.

aces_object <- function(
  data,
  x,
  y,
  quick_aces,
  max_degree,
  inter_cost,
  xi_degree,
  metric,
  shape,
  max_terms,
  err_red,
  minspan,
  endspan,
  kn_grid,
  kn_penalty,
  psi,
  wc,
  wq,
  aces_forward,
  aces,
  aces_cubic,
  aces_quintic,
  technology
) {
  object <- list()

  object[["data"]] <- list(
    "df" = data,
    "x" = x,
    "y" = y,
    "xnames" = colnames(data)[x],
    "ynames" = colnames(data)[y],
    "rownames" = rownames(data)
  )

  object[["control"]] <- list(
    "quick_aces" = quick_aces,
    "max_degree" = max_degree,
    "inter_cost" = inter_cost,
    "xi_degree" = xi_degree,
    "metric" = metric,
    "shape" = shape,
    "max_terms" = max_terms,
    "err_red" = err_red,
    "minspan" = minspan,
    "endspan" = endspan,
    "kn_grid" = kn_grid,
    "kn_penalty" = kn_penalty,
    "psi" = psi,
    "wc" = wc,
    "wq" = wq
  )

  object[["methods"]] <- list(
    "aces_forward" = aces_forward,
    "aces" = aces,
    "aces_cubic" = aces_cubic,
    "aces_quintic" = aces_quintic
  )

  object[["technology"]] <- technology

  return(object)
}

#' @title Prepare Data for ACES
#'
#' @description
#' Selects the input and output columns and creates the requested multiplicative
#' input interactions. Interaction names use \code{:}, for example
#' \code{capital:labour}. The returned matrix places original inputs first,
#' followed by interaction terms and outputs.
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
#' @param max_degree
#' Maximum interaction degree, or a list of input-index vectors defining the
#' interactions to add.
#'
#' @return
#' A matrix containing the inputs, input interactions, and outputs, in that order.

set_data <- function(
  data,
  x,
  y,
  max_degree
) {
  # A data frame can be extended by a named interaction column; a matrix cannot.
  # Convert here and return a matrix at the end, as documented.
  data <- as.data.frame(data, check.names = FALSE)

  # 1. generate interaction effects
  if (is.list(max_degree) || max_degree > 1) {
    if (is.list(max_degree)) {
      degree <- lapply(max_degree, function(group) match(group, x))
    } else {
      # create a list with all the possible combinations between 1 and as much len(x) elements
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
      name <- paste(name_vars, collapse = ":")

      # create the new variable
      data[, name] <- apply(data[, vars], 1, prod)
    }
  } else {
    new_x <- x
  }

  # 2. data correctly sorted
  data <- data[, c(new_x, y)]

  return(as.matrix(data))
}

#' @title Compute a Model Error
#'
#' @description
#' Computes a weighted lack-of-fit measure from observed and predicted outputs.
#' Negative predictions receive infinite error when all observed outputs are
#' positive.
#'
#' @param y_obs
#' A numeric matrix of observed outputs.
#'
#' @param y_hat
#' A numeric matrix of predicted outputs.
#'
#' @param metric
#' Character string specifying the error measure:
#' \itemize{
#' \item{\code{"mae"}: Mean absolute error.}
#' \item{\code{"mape"}: Mean absolute percentage error.}
#' \item{\code{"mse"}: Mean squared error.}
#' \item{\code{"msle"}: Mean squared logarithmic error.}
#' \item{\code{"rmse"}: Root mean squared error.}
#' \item{\code{"nrmse1"}: Root mean squared error divided by the mean.}
#' \item{\code{"nrmse2"}: Root mean squared error divided by the range.}
#' }
#'
#' @param weight
#' A numeric matrix of observation weights.
#'
#' @return
#' A single numeric error value.

err_metric <- function(
  y_obs,
  y_hat,
  metric,
  weight
) {
  # samples in data
  N <- nrow(y_obs)

  # number of outputs
  nY <- ncol(y_obs)

  if (all(y_obs > 0) && any(y_hat < 0)) {
    # do not "allow" negative predictions (negative outputs)
    error <- Inf
  } else if (metric == "mae") {
    # mean absolute error
    devtn <- abs(y_hat - y_obs)
    error <- sum(weight * devtn) / (N * nY)
  } else if (metric == "mape") {
    # mean absolute percentage error
    devtn <- abs(y_hat - y_obs) / y_obs
    error <- sum(weight * devtn) / (N * nY) * 100
  } else if (metric == "mse") {
    # mean squared error
    devtn <- (y_hat - y_obs)^2
    error <- sum(weight * devtn) / (N * nY)
  } else if (metric == "msle") {
    # mean squared logarithmic error
    devtn <- (log(y_hat + 1) - log(y_obs + 1))^2
    error <- sum(weight * devtn) / (N * nY)
  } else if (metric == "rmse") {
    # root mean squared error
    devtn <- (y_hat - y_obs)^2
    error <- sqrt(sum(weight * devtn) / (N * nY))
  } else if (metric == "nrmse1") {
    # normalized root mean squared error by the mean
    devtn <- (y_hat - y_obs)^2
    error <- sqrt(sum(weight * devtn) / (N * nY)) / mean(y_obs)
  } else {
    # compute the mean of column-wise maximums and minimums in y
    ymax <- mean(apply(y_obs, 2, max))
    ymin <- mean(apply(y_obs, 2, min))

    # normalized root mean squared error by the range
    devtn <- (y_hat - y_obs)^2
    error <- sqrt(sum(weight * devtn) / (N * nY)) / (ymax - ymin)
  }

  return(error)
}

#' @title Compute Knot Spacing
#'
#' @description
#' Computes the minimum number of observations required between adjacent knots
#' and at both ends of each input range. The spacing can be fixed by the user or
#' obtained from the Friedman or Zhang rules.
#'
#' @param kn_grid
#' A list of knot candidates, one element per input.
#'
#' @param minspan
#' Minimum number of observations between adjacent knots. Use \code{-2} for
#' the Zhang rule, \code{-1} for the Friedman rule, or a positive integer.
#'
#' @param endspan
#' Minimum number of observations before the first knot and after the last knot.
#' Use \code{-2} for the Zhang rule, \code{-1} for the Friedman rule, or a
#' positive integer.
#'
#' @param n_input
#' Number of inputs, including contextual variables.
#'
#' @references
#'
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{zhang1994}{aces}
#'
#' @return
#' A \code{list} with two components:
#' \itemize{
#'   \item{\code{minspan}: Computed spacing between adjacent knots.}
#'   \item{\code{endspan}: Computed spacing at the ends of each input range.}
#' }

compute_span <- function(
  kn_grid,
  minspan,
  endspan,
  n_input
) {
  # data.frame
  data <- do.call(cbind, kn_grid)

  # sample size
  N <- nrow(data)

  # minimum span (L)
  if (minspan == -2) { # Zhang approach

    L <- numeric(n_input)

    # fixed log_factor
    log_factor <- log2(-(1 / N) * log(0.95))

    for (var in 1:n_input) {
      # sorted variable
      sorted_var <- sort(data[, var])

      # 3 highest values
      max3 <- tail(sorted_var, 3)

      # 3 lowest values
      min3 <- head(sorted_var, 3)

      m1 <- -(max3[1] - min3[1]) / (2.5 * (N - 1)) * log_factor
      m2 <- (1 / N) * sum(max3 - min3)

      # Lvar limited for the 10% of the DMUs
      L[var] <- floor(min(N * 0.10, max(m1, m2)))
    }
  } else if (minspan == -1) { # Friedman approach (this value must be computed later)

    L <- -1
  } else {
    L <- min(N * 0.10, minspan)
  }

  # end span (Le)
  if (endspan == -2) { # Zhang approach

    Le <- numeric(n_input)

    # fixed log_factor
    log_factor <- log2(-(1 / N) * log(0.95))

    for (var in 1:n_input) {
      # sorted variable
      sorted_var <- sort(data[, var])

      # 3 highest values
      max3 <- tail(sorted_var, 3)

      # 3 lowest values
      min3 <- head(sorted_var, 3)

      m1 <- -(max3[1] - min3[1]) / (2.5 * (N - 1)) * log_factor
      m2 <- (1 / N) * sum(max3 - min3)

      Le[var] <- floor(min(N * 0.10, max(m1, m2)))
    }
  } else if (endspan == -1) { # Friedman approach

    Le <- floor(min(N * 0.1, 3 - log2(0.05 / n_input)))
  } else {
    Le <- min(N * 0.1, endspan)
  }

  return(list(L, Le))
}

#' @title Build the Knot Grid
#'
#' @description
#' Creates the knot grid used by ACES. A user-supplied grid is expanded when
#' interaction variables are present; otherwise candidates are obtained from
#' the prepared inputs. Quick ACES uses DEA scores to retain knot neighborhoods
#' around efficient DMUs.
#'
#' @param data
#' A matrix containing the prepared model variables.
#'
#' @param n_input_1
#' Number of original inputs, including contextual variables.
#'
#' @param n_input_2
#' Number of inputs after adding interactions.
#'
#' @param kn_grid
#' A custom list of knot candidates, or \code{-1} for automatic selection.
#'
#' @param quick_aces
#' If \code{TRUE}, retain knot neighborhoods identified by the Quick ACES
#' DEA-guided reduction rule.
#'
#' @param dea_scores
#' A matrix of output-oriented DEA-VRS efficiency scores.
#'
#' @return
#' A list of knot candidates, one element per prepared input.

set_knots_grid <- function(
  data,
  n_input_1,
  n_input_2,
  kn_grid,
  quick_aces,
  dea_scores
) {
  # Case 1: kn_grid is provided (list) and new variables are created (nX > inputs):
  # expand the kn_grid list.

  # Case 2: kn_grid is provided (list) and new variables are not created (nX = inputs):
  # keep the same kn_grid list.

  # Case 3: kn_grid is not provided:
  # create the kn_grid list.

  if (is.list(kn_grid)) { # if kn_grid is provided

    if (n_input_2 > n_input_1) {
      # number of new variables (through interactions)
      new_vars <- n_input_2 - n_input_1

      for (v in seq_len(new_vars)) {
        # variable index
        var_idx <- n_input_2 - new_vars + v

        # variable name
        var_name <- colnames(data)[var_idx]

        # variable data
        var_data <- data[, var_idx]

        # length of the maximum grid
        max_len_grid <- max(sapply(kn_grid, length))

        # grid of knots for the new variable
        kn_grid[[var_name]] <- seq(
          from = min(var_data),
          to = max(var_data),
          length.out = max_len_grid
        )
      }
    }
  } else { # if kn_grid is not provided, create it

    if (quick_aces) {
      # identify efficient DMUs
      eff_dmus <- data[abs(dea_scores - 1) < 0.001, ]

      # grid of knots: inputs and interactions
      kn_grid <- lapply(1:n_input_2, function(i) eff_dmus[, i])

      for (j in 1:length(kn_grid)) {
        # sort the unique values of the dimension
        sorted_values <- sort(unique(data[, j]))

        # find the index of the current DMU value in the sorted list
        matched_indices <- which(!is.na(match(data[, j], kn_grid[[j]])))

        # add neighbourhood
        matched_indices_l <- matched_indices - 1
        matched_indices_r <- matched_indices + 1

        # remove NA values
        matched_indices <- sort(unique(c(matched_indices, matched_indices_l, matched_indices_r)))
        matched_indices <- matched_indices[matched_indices > 0 & matched_indices <= nrow(data)]

        # final knots grid
        kn_grid[[j]] <- data[matched_indices, j]
      }
    } else {
      kn_grid <- lapply(1:n_input_2, function(i) data[, i])
    }

  }

  # Keep prepared-input names (including interactions such as x1:x2) aligned
  # with the knot grids, whether the grids were automatic or user supplied.
  names(kn_grid) <- colnames(data)[seq_len(n_input_2)]

  return(kn_grid)
}
