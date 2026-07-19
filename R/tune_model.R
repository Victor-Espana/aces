#' Create Data Partitions for ACES Tuning
#'
#' Creates reusable resampling splits for [tune_model_aces()] and
#' [tune_model_rf_aces()]. The same partition can be used to compare models on
#' identical samples. The partition stores all information required for tuning,
#' so candidate configurations are evaluated under the same validation design.
#'
#' @param data A data frame or matrix containing the model variables.
#' @param x Column indexes of input variables in `data`.
#' @param y Column indexes of output variables in `data`.
#' @param validation Resampling method: `"kfold"`, `"repeated_kfold"`,
#'   `"holdout"`, `"bootstrap"`, `"loo"`, or `"group_kfold"`.
#' @param resampling_args A list of method-specific resampling parameters.
#'   Valid elements depend on `validation`:
#'   \itemize{
#'     \item `kfold`: `folds` (default 5), `strata` (default `NULL`).
#'     \item `repeated_kfold`: `folds` (default 5), `repeats` (default 1),
#'       `strata` (default `NULL`).
#'     \item `holdout`: `repeats` (default 1), `train_prop` (default 0.8),
#'       `strata` (default `NULL`).
#'     \item `bootstrap`: `repeats` (default 1).
#'     \item `loo`: no arguments needed.
#'     \item `group_kfold`: `folds` (default 5), `groups` (required).
#'   }
#'   Unrecognized elements trigger an error.
#' @param delta Maximum output-oriented DEA score included in validation. This
#'   limits the influence of observations that lie far beyond the training
#'   technology; it must be at least 1.
#' @param seed Optional integer seed used to make the resampling splits
#'   reproducible.
#'
#' @return An object of class `aces_partition` containing the data, input and
#'   output indexes, resampling indexes, `delta`, and seed.
#'
#' @examples
#' partition <- data_partition(
#'   data = cobb_douglas_XnY1(N = 30, nX = 3),
#'   x = 1:3,
#'   y = 4,
#'   resampling_args = list(folds = 3),
#'   delta = 1.5,
#'   seed = 42
#' )
#'
#' @md
#' @export
data_partition <- function(
  data,
  x,
  y,
  validation = "kfold",
  resampling_args = list(),
  delta = 1.05,
  seed = NULL
) {
  data <- as.data.frame(data)
  n <- nrow(data)

  if (n < 3L) stop("'data' must contain at least 3 observations.")

  x <- .validate_indices(data, x, "x")
  y <- .validate_indices(data, y, "y")

  if (length(intersect(x, y))) stop("'x' and 'y' must be disjoint.")
  if (!all(vapply(data[c(x, y)], is.numeric, logical(1)))) {
    stop("All input and output columns must be numeric.")
  }

  validation <- match.arg(
    tolower(validation),
    c("kfold", "repeated_kfold", "holdout", "bootstrap", "loo", "group_kfold")
  )

  if (!is.list(resampling_args)) {
    stop("'resampling_args' must be a list.")
  }

  # Get valid arguments and defaults for the chosen method
  spec <- .resampling_defaults(validation)

  unknown <- setdiff(names(resampling_args), spec$valid)
  if (length(unknown)) {
    stop(
      "Unknown resampling arguments for '", validation, "': ",
      paste(unknown, collapse = ", "), ". Valid arguments: ",
      paste(spec$valid, collapse = ", "), "."
    )
  }

  args <- modifyList(spec$defaults, resampling_args)

  # Extract and validate resampling parameters
  folds     <- args$folds
  repeats   <- args$repeats
  train_prop <- args$train_prop
  strata    <- args$strata
  groups    <- args$groups

  if (!is.null(folds)) {
    folds <- as.integer(folds)
    if (length(folds) != 1L || is.na(folds) || folds < 2L) {
      stop("'folds' must be an integer greater than 1.")
    }
  }
  if (!is.null(repeats)) {
    repeats <- as.integer(repeats)
    if (length(repeats) != 1L || is.na(repeats) || repeats < 1L) {
      stop("'repeats' must be a positive integer.")
    }
  }
  if (!is.null(train_prop)) {
    if (!is.numeric(train_prop) || length(train_prop) != 1L ||
        train_prop <= 0 || train_prop >= 1) {
      stop("'train_prop' must be between 0 and 1.")
    }
  }
  if (!is.numeric(delta) || length(delta) != 1L || delta < 1) {
    stop("'delta' must be a number greater than or equal to 1.")
  }

  strata <- .partition_variable(data, strata, "strata")
  groups <- .partition_variable(data, groups, "groups")

  if (validation %in% c("kfold", "repeated_kfold") && folds > n) {
    stop("'folds' cannot exceed the number of observations.")
  }
  if (validation == "group_kfold") {
    if (is.null(groups)) stop("'groups' is required for group_kfold.")
    if (folds > length(unique(groups))) {
      stop("'folds' cannot exceed the number of groups.")
    }
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  out <- list(
    data = data,
    x = x,
    y = y,
    splits = .make_folds(
      n, validation, folds, repeats, train_prop, strata, groups
    ),
    validation = validation,
    delta = delta,
    seed = seed
  )
  class(out) <- "aces_partition"
  out
}


.resampling_defaults <- function(validation) {
  switch(
    validation,
    kfold = list(
      valid = c("folds", "strata"),
      defaults = list(folds = 5L, strata = NULL)
    ),
    repeated_kfold = list(
      valid = c("folds", "repeats", "strata"),
      defaults = list(folds = 5L, repeats = 1L, strata = NULL)
    ),
    holdout = list(
      valid = c("repeats", "train_prop", "strata"),
      defaults = list(repeats = 1L, train_prop = 0.8, strata = NULL)
    ),
    bootstrap = list(
      valid = c("repeats"),
      defaults = list(repeats = 1L)
    ),
    loo = list(
      valid = character(0),
      defaults = list()
    ),
    group_kfold = list(
      valid = c("folds", "groups"),
      defaults = list(folds = 5L, groups = NULL)
    )
  )
}


#' Tune an ACES Model
#'
#' Builds a parameter grid and selects the ACES configuration with the smallest
#' frontier-aware validation error. Each model argument may contain one or more
#' candidate values. After comparing all configurations on the stored splits, the
#' best configuration is fitted once more using the complete data set.
#'
#' Arguments have the same form as in [aces()]. To tune elements of a list
#' argument, place vectors inside the list. For example,
#' `mul_BF = list(max_degree = c(1, 2), inter_cost = c(0.01, 0.05))`.
#'
#' @param partition An object created by [data_partition()].
#' @param scale_data Candidate logical values. If `TRUE`, divide each input and
#'   output by its mean before fitting; this can improve solver convergence.
#' @param quick_aces Candidate logical values. If `TRUE`, reduce the candidate
#'   basis functions considered during the forward step.
#' @param metric Candidate lack-of-fit measures used while fitting ACES. See
#'   [aces()] for the available values.
#' @param mul_BF A list with candidate values for `max_degree`, the maximum
#'   interaction degree, and `inter_cost`, the extra improvement required to add
#'   a higher-degree basis function.
#' @param shape A list with candidate logical values for `mono`, which enforces
#'   non-decreasing monotonicity, and `conc`, which enforces concavity. When both
#'   are `TRUE`, the output-specific fits are production functions; otherwise,
#'   they are intermediate predictors of maximum attainable output used to
#'   construct the final technology.
#' @param max_terms,err_red,kn_grid,minspan,endspan,kn_penalty Additional [aces()]
#'   parameters. Each may contain one or several values. A custom list supplied
#'   to `kn_grid` is treated as one value; use a list of custom grids to compare
#'   several of them.
#' @param validation_metric One or more validation error measures. All supplied
#'   measures are computed; the first selects the best configuration. Options:
#'   `"mae"`, `"mape"`, `"mse"`, `"msle"`, `"rmse"`, `"rmsle"`,
#'   `"nrmse_mean"`, `"nrmse_range"`.
#' @param returns Returns-to-scale assumption used by the frontier-aware
#'   validation weights: `"variable"` or `"constant"`.
#' @param method ACES model used to generate validation predictions. Available
#'   methods are `"aces_forward"`, `"aces"`, `"aces_cubic"`, and
#'   `"aces_quintic"`.
#' @param verbose If `TRUE`, print the current configuration and split while
#'   tuning.
#'
#' @return An object of class `aces_tune` with the best parameters, best
#'   scores for every computed metric, results for every configuration, and
#'   the final model fitted to all data.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' production_data <- cobb_douglas_XnY1(N = 100, nX = 6)
#' partition <- data_partition(
#'   production_data,
#'   x = 1:6,
#'   y = 7,
#'   validation = "holdout",
#'   resampling_args = list(repeats = 2),
#'   delta = 1.5,
#'   seed = 42
#' )
#' tuned <- tune_model_aces(
#'   partition,
#'   quick_aces = TRUE,
#'   mul_BF = list(
#'     max_degree = c(1, 2),
#'     inter_cost = c(0.01, 0.05)
#'   ),
#'   err_red = c(0.01, 0.05),
#'   validation_metric = "rmse"
#' )
#' }
#'
#' @md
#' @export
tune_model_aces <- function(
  partition,
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
  kn_penalty = 2,
  validation_metric = "mse",
  returns = "variable",
  method = "aces",
  verbose = TRUE
) {
  grid <- .model_grid(list(
    scale_data = scale_data,
    quick_aces = quick_aces,
    mul_BF = .nested_grid(mul_BF, c("max_degree", "inter_cost")),
    metric = metric,
    shape = .nested_grid(shape, c("mono", "conc")),
    max_terms = max_terms,
    err_red = err_red,
    kn_grid = .grid_values(kn_grid),
    minspan = minspan,
    endspan = endspan,
    kn_penalty = kn_penalty
  ))

  .tune_model(
    partition = partition,
    grid = grid,
    model = "aces",
    validation_metric = validation_metric,
    returns = returns,
    method = method,
    verbose = verbose
  )
}


#' Tune an RF-ACES Model
#'
#' Builds a parameter grid and selects the RF-ACES configuration with the
#' smallest frontier-aware validation error. After comparing all configurations
#' on the stored splits, the best configuration is fitted once more using the
#' complete data set.
#'
#' The model arguments have the same form as in [rf_aces()]. Vectors inside
#' `mul_BF`, `shape`, and `early_stopping` are expanded internally.
#'
#' @param learners Candidate positive integers giving the maximum number of ACES
#'   learners in the forest.
#' @param bag_size Candidate positive integers giving the number of observations
#'   sampled with replacement for each learner.
#' @param max_feats Candidate values for the number of inputs randomly considered
#'   by each learner.
#' @param early_stopping A list with candidate values for `ma_window`, the OOB
#'   moving-average window, and `tolerance`, the allowed consecutive iterations
#'   without improvement.
#' @param method RF-ACES model used to generate validation predictions. Available
#'   methods are `"rf_aces"`, `"rf_aces_cubic"`, and `"rf_aces_quintic"`.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' production_data <- cobb_douglas_XnY1(N = 100, nX = 6)
#' partition <- data_partition(
#'   production_data,
#'   x = 1:6,
#'   y = 7,
#'   validation = "holdout",
#'   resampling_args = list(repeats = 2),
#'   delta = 1.5,
#'   seed = 42
#' )
#' tuned <- tune_model_rf_aces(
#'   partition,
#'   quick_aces = TRUE,
#'   mul_BF = list(
#'     max_degree = c(1, 2),
#'     inter_cost = 0.05
#'   ),
#'   learners = 3,
#'   bag_size = c(60, 80),
#'   max_feats = c(2, 4),
#'   validation_metric = "rmse"
#' )
#' }
#'
#' @md
#' @rdname tune_model_aces
#' @export
tune_model_rf_aces <- function(
  partition,
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
  bag_size = nrow(partition$data),
  max_feats = length(partition$x) / 3,
  early_stopping = list(
    "ma_window" = 10,
    "tolerance" = 5
  ),
  max_terms = 50,
  err_red = 0.01,
  kn_grid = -1,
  minspan = -1,
  endspan = -1,
  validation_metric = "mse",
  returns = "variable",
  method = "rf_aces",
  verbose = TRUE
) {
  parameters <- list(
    scale_data = scale_data,
    quick_aces = quick_aces,
    mul_BF = .nested_grid(mul_BF, c("max_degree", "inter_cost")),
    metric = metric,
    shape = .nested_grid(shape, c("mono", "conc")),
    learners = learners,
    bag_size = bag_size,
    max_feats = max_feats,
    early_stopping = .nested_grid(
      early_stopping,
      c("ma_window", "tolerance")
    ),
    max_terms = max_terms,
    err_red = err_red,
    kn_grid = .grid_values(kn_grid),
    minspan = minspan,
    endspan = endspan
  )

  .tune_model(
    partition = partition,
    grid = .model_grid(parameters),
    model = "rf_aces",
    validation_metric = validation_metric,
    returns = returns,
    method = method,
    verbose = verbose
  )
}


.tune_model <- function(
  partition,
  grid,
  model,
  validation_metric,
  returns,
  method,
  verbose
) {
  if (!inherits(partition, "aces_partition")) {
    stop("'partition' must be created with data_partition().")
  }

  valid_metrics <- c(
    "mae", "mape", "mse", "msle", "rmse", "rmsle",
    "nrmse_mean", "nrmse_range"
  )
  validation_metrics <- unique(tolower(validation_metric))
  invalid <- setdiff(validation_metrics, valid_metrics)
  if (length(invalid)) {
    stop(
      "Unknown validation metric(s): ",
      paste(invalid, collapse = ", "), "."
    )
  }
  returns <- match.arg(tolower(returns), c("variable", "constant"))

  valid_methods <- if (model == "aces") {
    c("aces_forward", "aces", "aces_cubic", "aces_quintic")
  } else {
    c("rf_aces", "rf_aces_cubic", "rf_aces_quintic")
  }
  method <- match.arg(method, valid_methods)

  data <- partition$data
  x <- partition$x
  y <- partition$y
  splits <- partition$splits
  delta <- partition$delta
  seed <- partition$seed
  weights <- lapply(splits, function(split) {
    .dea_weights(data, x, y, split$analysis, split$assessment, delta, returns)
  })

  n_configs <- nrow(grid)
  n_splits  <- length(splits)
  scores <- lapply(validation_metrics, function(m) {
    matrix(NA_real_, nrow = n_configs, ncol = n_splits)
  })
  names(scores) <- validation_metrics
  fit_fun <- if (model == "aces") aces else rf_aces

  for (i in seq_len(n_configs)) {
    config <- .grid_row(grid, i)

    for (j in seq_along(splits)) {
      if (verbose) {
        message(sprintf(
          "[%s] configuration %d/%d | split %d/%d",
          paste0("tune_model_", model), i, n_configs, j, n_splits
        ))
      }

      if (!is.null(seed)) set.seed(as.integer(seed) + j)
      split <- splits[[j]]
      selected <- weights[[j]]$selected
      fit_data <- data[split$analysis, , drop = FALSE]

      fit <- do.call(
        fit_fun,
        c(list(data = fit_data, x = x, y = y), config)
      )
      estimate <- .predict_frontier(
        fit, data[selected, , drop = FALSE], x, model, method
      )
      truth <- as.matrix(data[selected, y, drop = FALSE])

      for (m in validation_metrics) {
        scores[[m]][i, j] <- .weighted_metric(
          truth, estimate, weights[[j]]$weights, m
        )
      }
    }
  }

  # Summarise each metric across splits
  score_cols <- do.call(cbind, lapply(validation_metrics, function(m) {
    df <- data.frame(
      mean = rowMeans(scores[[m]]),
      sd   = apply(scores[[m]], 1L, stats::sd)
    )
    names(df) <- paste0(c("mean_", "sd_"), m)
    df
  }))

  primary       <- validation_metrics[1L]
  primary_means <- score_cols[[paste0("mean_", primary)]]
  best          <- which.min(primary_means)

  results <- cbind(
    data.frame(config_id = seq_len(n_configs)),
    grid,
    score_cols
  )
  results <- results[
    order(results[[paste0("mean_", primary)]]), , drop = FALSE
  ]
  rownames(results) <- NULL

  if (verbose) message("Refitting the best configuration on all data")
  if (!is.null(seed)) set.seed(as.integer(seed))
  final_model <- do.call(
    fit_fun,
    c(
      list(data = data, x = x, y = y),
      .grid_row(grid, best)
    )
  )

  best_scores <- vapply(validation_metrics, function(m) {
    rowMeans(scores[[m]])[best]
  }, numeric(1))

  best_params <- .flatten_params(grid[best, , drop = FALSE])

  out <- list(
    best_params        = best_params,
    best_scores        = best_scores,
    results            = results,
    final_model        = final_model,
    model              = model,
    validation_metrics = validation_metrics
  )
  class(out) <- "aces_tune"
  out
}


#' @export
print.aces_tune <- function(x, ...) {
  primary <- x$validation_metrics[1L]
  cat("Frontier-aware tuning (", x$model, ")\n", sep = "")
  cat("Best configuration (selected by ", primary, "):\n", sep = "")
  for (m in x$validation_metrics) {
    cat("  ", m, ": ", format(x$best_scores[m], digits = 5), "\n", sep = "")
  }
  cat("Best parameters:\n")
  print(x$best_params, row.names = FALSE)
  invisible(x)
}


#' Select the Best Tuning Configuration
#'
#' Returns the evaluated configuration with the lowest value of `metric`.
#'
#' @param object An object of class `aces_tune`.
#' @param metric Validation metric used for ranking. It must have been computed
#'   during tuning.
#'
#' @return A list with:
#'   \itemize{
#'     \item `best_params`: a one-row data frame with the model parameters.
#'     \item `best_scores`: a named numeric vector with the mean score of the
#'       selected configuration for every computed metric.
#'     \item `config_id`: integer ID of the selected configuration.
#'   }
#'
#' @examples
#' \dontrun{
#' tuned <- tune_model_aces(
#'   partition,
#'   validation_metric = c("rmse", "mae", "mse")
#' )
#' select_best(tuned, metric = "mae")
#' }
#'
#' @md
#' @export
select_best <- function(object, metric) {
  if (!inherits(object, "aces_tune")) {
    stop("'object' must be an 'aces_tune' object.")
  }
  metric <- tolower(metric)
  if (!metric %in% object$validation_metrics) {
    stop(
      "'metric' must be one of: ",
      paste(object$validation_metrics, collapse = ", "), "."
    )
  }

  mean_col <- paste0("mean_", metric)
  best_idx <- which.min(object$results[[mean_col]])
  best_row <- object$results[best_idx, , drop = FALSE]

  # Extract scores for all metrics
  best_scores <- vapply(object$validation_metrics, function(m) {
    best_row[[paste0("mean_", m)]]
  }, numeric(1))

  # Extract param columns (exclude config_id and score columns)
  score_col_names <- c(
    paste0("mean_", object$validation_metrics),
    paste0("sd_", object$validation_metrics)
  )
  param_cols <- setdiff(names(best_row), c("config_id", score_col_names))
  best_params <- .flatten_params(best_row[param_cols])

  list(
    best_params = best_params,
    best_scores = best_scores,
    config_id   = best_row$config_id
  )
}


.flatten_params <- function(params) {
  cols <- lapply(names(params), function(nm) {
    val <- params[[nm]]
    # list-columns from expand.grid wrap the value in an extra list
    if (is.list(val) && length(val) == 1L) val <- val[[1L]]
    if (is.list(val) && length(val) > 0L && !is.null(names(val)) &&
        all(vapply(val, function(v) length(v) == 1L, logical(1)))) {
      # Named list with scalar elements -> expand to individual columns
      as.data.frame(val, stringsAsFactors = FALSE)
    } else {
      df <- data.frame(x = I(list(val)), stringsAsFactors = FALSE)
      names(df) <- nm
      df
    }
  })
  out <- do.call(cbind, cols)
  rownames(out) <- NULL
  out
}


.model_grid <- function(parameters) {
  if (any(lengths(parameters) == 0L)) {
    stop("Model parameters must contain at least one value.")
  }
  expand.grid(
    parameters,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}


.grid_values <- function(x) {
  if (!is.list(x)) return(x)
  if (length(x) && all(vapply(x, is.list, logical(1)))) return(x)
  list(x)
}


.nested_grid <- function(parameter, required_names) {
  if (!is.list(parameter) || !all(required_names %in% names(parameter))) {
    stop(
      "Expected a list with: ",
      paste(required_names, collapse = ", "),
      "."
    )
  }

  values <- lapply(parameter[required_names], .grid_values)
  grid <- .model_grid(values)
  lapply(seq_len(nrow(grid)), function(row) .grid_row(grid, row))
}


.grid_row <- function(grid, row) {
  lapply(as.list(grid[row, , drop = FALSE]), function(value) {
    if (is.factor(value)) value <- as.character(value)
    if (is.list(value) && length(value) == 1L) value[[1L]] else value
  })
}


.validate_indices <- function(data, cols, name) {
  if (!is.numeric(cols) || !length(cols) || anyNA(cols) ||
      any(cols != as.integer(cols)) || any(cols < 1) || any(cols > ncol(data))) {
    stop("'", name, "' must contain valid column indexes.")
  }
  as.integer(cols)
}


.partition_variable <- function(data, value, name) {
  if (is.character(value) && length(value) == 1L && value %in% names(data)) {
    value <- data[[value]]
  }
  if (!is.null(value) && length(value) != nrow(data)) {
    stop("'", name, "' must have one value per observation.")
  }
  value
}


.make_folds <- function(n, validation, folds, repeats, train_prop, strata, groups) {
  rows <- seq_len(n)

  if (validation == "kfold") {
    assignment <- .stratified_folds(n, folds, strata)
    return(lapply(seq_len(folds), function(fold) {
      list(
        analysis = rows[assignment != fold],
        assessment = rows[assignment == fold],
        id = paste0("Fold", fold)
      )
    }))
  }

  if (validation == "repeated_kfold") {
    out <- list()
    for (repeat_id in seq_len(repeats)) {
      assignment <- .stratified_folds(n, folds, strata)
      for (fold in seq_len(folds)) {
        out[[length(out) + 1L]] <- list(
          analysis = rows[assignment != fold],
          assessment = rows[assignment == fold],
          id = paste0("Repeat", repeat_id, "_Fold", fold)
        )
      }
    }
    return(out)
  }

  if (validation == "holdout") {
    return(lapply(seq_len(repeats), function(repeat_id) {
      train_n <- floor(train_prop * n)
      analysis <- if (is.null(strata)) {
        sample(rows, train_n)
      } else {
        .stratified_sample(rows, strata, train_n)
      }
      list(
        analysis = analysis,
        assessment = setdiff(rows, analysis),
        id = paste0("Holdout", repeat_id)
      )
    }))
  }

  if (validation == "bootstrap") {
    return(lapply(seq_len(repeats), function(repeat_id) {
      analysis <- sample(rows, n, replace = TRUE)
      assessment <- setdiff(rows, unique(analysis))
      if (!length(assessment)) {
        stop("A bootstrap sample has no out-of-bag observations.")
      }
      list(
        analysis = analysis,
        assessment = assessment,
        id = paste0("Bootstrap", repeat_id)
      )
    }))
  }

  if (validation == "loo") {
    return(lapply(rows, function(row) {
      list(
        analysis = setdiff(rows, row),
        assessment = row,
        id = paste0("LOO", row)
      )
    }))
  }

  unique_groups <- unique(groups)
  assignment <- sample(rep(seq_len(folds), length.out = length(unique_groups)))
  names(assignment) <- unique_groups
  lapply(seq_len(folds), function(fold) {
    assessment_groups <- unique_groups[assignment == fold]
    assessment <- rows[groups %in% assessment_groups]
    list(
      analysis = setdiff(rows, assessment),
      assessment = assessment,
      id = paste0("GroupFold", fold)
    )
  })
}


.stratified_folds <- function(n, folds, strata) {
  if (is.null(strata)) return(sample(rep(seq_len(folds), length.out = n)))

  assignment <- integer(n)
  for (level in unique(strata)) {
    rows <- which(strata == level)
    assignment[rows] <- sample(rep(seq_len(folds), length.out = length(rows)))
  }
  assignment
}


.stratified_sample <- function(rows, strata, target) {
  selected <- integer()
  for (level in unique(strata)) {
    level_rows <- rows[strata == level]
    size <- round(target * length(level_rows) / length(rows))
    selected <- c(selected, sample(level_rows, min(size, length(level_rows))))
  }

  if (length(selected) > target) selected <- sample(selected, target)
  if (length(selected) < target) {
    selected <- c(selected, sample(setdiff(rows, selected), target - length(selected)))
  }
  selected
}


.dea_weights <- function(data, x, y, analysis, assessment, delta, returns) {
  technology <- unique(c(analysis, assessment))
  phi <- as.numeric(rad_out(
    tech_xmat = as.matrix(data[technology, x, drop = FALSE]),
    tech_ymat = as.matrix(data[technology, y, drop = FALSE]),
    eval_xmat = as.matrix(data[assessment, x, drop = FALSE]),
    eval_ymat = as.matrix(data[assessment, y, drop = FALSE]),
    convexity = TRUE,
    returns = returns,
    type = "objective"
  ))

  phi[is.finite(phi)] <- pmax(1, phi[is.finite(phi)])
  keep <- is.finite(phi) & phi <= delta
  if (!any(keep)) {
    stop(
      "No validation observations satisfy phi <= ", delta,
      ". Increase 'delta' or use another partition."
    )
  }

  list(
    selected = assessment[keep],
    phi = phi[keep],
    weights = 1 / phi[keep]
  )
}


.predict_frontier <- function(object, newdata, x, model, method) {
  if (model == "rf_aces") {
    return(as.matrix(rf_aces_predict(
      object = object,
      eval_data = newdata,
      x = x,
      method = method
    )))
  }

  fitted_method <- object[["methods"]][[method]]
  if (is.null(fitted_method)) {
    stop("Fitted ACES object does not contain method '", method, "'.")
  }

  scaling <- object[["control"]][["scale"]]
  working_data <- newdata
  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    working_data[, x] <- sweep(
      as.matrix(newdata[, x, drop = FALSE]), 2L, scaling$mean_x, "/"
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
    nrow = nrow(newdata),
    ncol = length(object[["data"]][["y"]])
  )

  if (!is.null(scaling) && isTRUE(scaling$is_scaled)) {
    estimate <- sweep(estimate, 2L, scaling$mean_y, "*")
  }
  estimate
}


.weighted_metric <- function(truth, estimate, weights, metric) {
  truth <- as.matrix(truth)
  estimate <- as.matrix(estimate)
  weights <- weights / sum(weights)
  weight_matrix <- matrix(
    weights / ncol(truth),
    nrow = length(weights),
    ncol = ncol(truth)
  )
  error <- estimate - truth

  switch(
    metric,
    mae = sum(weight_matrix * abs(error)),
    mape = {
      if (any(truth == 0)) stop("MAPE is undefined when an output is zero.")
      100 * sum(weight_matrix * abs(error / truth))
    },
    mse = sum(weight_matrix * error^2),
    msle = {
      if (any(truth < 0) || any(estimate < 0)) {
        stop("MSLE requires non-negative values.")
      }
      sum(weight_matrix * (log1p(estimate) - log1p(truth))^2)
    },
    rmse = sqrt(sum(weight_matrix * error^2)),
    rmsle = {
      if (any(truth < 0) || any(estimate < 0)) {
        stop("RMSLE requires non-negative values.")
      }
      sqrt(sum(weight_matrix * (log1p(estimate) - log1p(truth))^2))
    },
    nrmse_mean = {
      rmse <- sqrt(sum(weight_matrix * error^2))
      denominator <- sum(weight_matrix * abs(truth))
      if (denominator == 0) stop("NRMSE mean denominator is zero.")
      rmse / denominator
    },
    nrmse_range = {
      rmse <- sqrt(sum(weight_matrix * error^2))
      denominator <- mean(apply(truth, 2L, function(column) diff(range(column))))
      if (denominator == 0) stop("NRMSE range denominator is zero.")
      rmse / denominator
    }
  )
}
