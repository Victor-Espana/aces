#' @title Validate ACES Arguments
#'
#' @description
#' Checks the arguments used by \code{\link{aces}} and stops with an informative
#' message when an argument is invalid.
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
#' If \code{TRUE}, divide each input and output by its mean before fitting.
#'
#' @param quick_aces
#' If \code{TRUE}, use Quick ACES to reduce the number of candidate basis
#' functions.
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
#' Character string specifying the lack-of-fit measure. Options are
#' \code{"mae"}, \code{"mape"},
#' \code{"mse"}, \code{"msle"}, \code{"rmse"}, \code{"nrmse1"}, \code{"nrmse2"}.
#'
#' @param max_terms
#' Maximum number of terms created during the forward step.
#'
#' @param err_red
#' Minimum relative error reduction required to add a new pair of first-degree
#' basis functions.
#'
#' @param kn_grid
#' Knot candidates. Use \code{-1} for automatic selection or supply a list with
#' one numeric vector for each input.
#'
#' @param minspan
#' Minimum number of observations between adjacent knots. Use \code{-2},
#' \code{-1}, or a positive integer.
#'
#' @param endspan
#' Minimum number of observations before the first knot and after the last knot.
#' Use \code{-2}, \code{-1}, or a positive integer.
#'
#' @param kn_penalty
#' Non-negative penalty per knot used to compute GCV.
#'
#' @importFrom stats na.omit
#'
#' @return
#' \code{NULL} invisibly. The function is used for its validation side effects.

display_errors_aces <- function(
  data,
  x,
  y,
  scale_data,
  quick_aces,
  max_degree,
  inter_cost,
  metric,
  max_terms,
  err_red,
  kn_grid,
  minspan,
  endspan,
  kn_penalty
  ) {

  # =========== #
  #    DATA     #
  # =========== #

  # data must be a data.frame or a matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data.frame or a matrix.")
  }

  # data must have at least one row
  if (nrow(data) == 0) {
    stop("'data' has no rows.")
  }

  # =========== #
  #   x and y   #
  # =========== #

  # x must be a non-empty numeric/integer vector
  if (!is.numeric(x) && !is.integer(x) || length(x) == 0) {
    stop("'x' must be a non-empty numeric vector of column indexes.")
  }

  # y must be a non-empty numeric/integer vector
  if (!is.numeric(y) && !is.integer(y) || length(y) == 0) {
    stop("'y' must be a non-empty numeric vector of column indexes.")
  }

  # x indexes must be valid column positions
  if (any(x < 1) || any(x > ncol(data)) || any(floor(x) != x)) {
    stop(paste0(
      "'x' contains invalid column index/indexes. ",
      "Indexes must be positive integers in [1, ", ncol(data), "]."
    ))
  }

  # y indexes must be valid column positions
  if (any(y < 1) || any(y > ncol(data)) || any(floor(y) != y)) {
    stop(paste0(
      "'y' contains invalid column index/indexes. ",
      "Indexes must be positive integers in [1, ", ncol(data), "]."
    ))
  }

  # x and y must not overlap
  if (length(intersect(x, y)) > 0) {
    stop(paste0(
      "'x' and 'y' share column index/indexes: ",
      paste(intersect(x, y), collapse = ", "), ". ",
      "Input and output variables must be distinct."
    ))
  }

  # all selected variables must be numeric
  if (!all(sapply(as.data.frame(data[, c(x, y)]), is.numeric))) {
    non_numeric <- colnames(data)[c(x, y)][
      !sapply(as.data.frame(data[, c(x, y)]), is.numeric)
    ]
    stop(paste0(
      "All input/output variables must be numeric. ",
      "Non-numeric column(s) detected: ",
      paste(non_numeric, collapse = ", "), "."
    ))
  }

  # warn about NA values (they will be dropped internally)
  if (any(is.na(data[, c(x, y)]))) {
    warning("'data' contains NA values in the selected columns. Rows with NAs will be omitted.")
  }

  # ============= #
  #   BOOLEANS    #
  # ============= #

  if (!is.logical(scale_data) || length(scale_data) != 1 || is.na(scale_data)) {
    stop("'scale_data' must be TRUE or FALSE.")
  }

  if (!is.logical(quick_aces) || length(quick_aces) != 1 || is.na(quick_aces)) {
    stop("'quick_aces' must be TRUE or FALSE.")
  }

  # ============= #
  #   mul_BF      #
  # ============= #

  # max_degree: numeric scalar or list of index vectors
  if (is.list(max_degree)) {
    if (length(max_degree) == 0) {
      stop("'mul_BF$max_degree' is an empty list. Provide at least one interaction group.")
    }
    # each element must be a numeric/integer vector of valid input indexes
    for (k in seq_along(max_degree)) {
      grp <- max_degree[[k]]
      if (!is.numeric(grp) && !is.integer(grp) || length(grp) < 2) {
        stop(paste0(
          "'mul_BF$max_degree[[", k, "]]' must be a numeric vector of at least ",
          "two input column indexes defining an interaction term."
        ))
      }
      if (any(!grp %in% x)) {
        stop(paste0(
          "'mul_BF$max_degree[[", k, "]]' contains index/indexes not present in 'x'."
        ))
      }
    }
  } else {
    if (!is.numeric(max_degree) || length(max_degree) != 1 ||
      is.na(max_degree) || floor(max_degree) != max_degree || max_degree < 1) {
      stop("'mul_BF$max_degree' must be a positive integer or a list of interaction index vectors.")
    }
    if (max_degree > length(x)) {
      stop(paste0(
        "'mul_BF$max_degree' (", max_degree, ") cannot exceed the number of inputs (", length(x), ")."
      ))
    }
  }

  # inter_cost must be a numeric value in [0, 1]
  if (!is.numeric(inter_cost) || length(inter_cost) != 1 || is.na(inter_cost) ||
    inter_cost < 0 || inter_cost > 1) {
    stop("'mul_BF$inter_cost' must be a numeric value in [0, 1].")
  }

  # ============= #
  #    METRIC     #
  # ============= #

  valid_metrics <- c("mae", "mape", "mse", "msle", "rmse", "nrmse1", "nrmse2")

  if (!is.character(metric) || length(metric) != 1 || !metric %in% valid_metrics) {
    stop(paste0(
      "'metric' must be one of: ", paste(valid_metrics, collapse = ", "), ". ",
      "See help(\"aces\") for details."
    ))
  }

  # ============= #
  #   max_terms   #
  # ============= #

  if (!is.numeric(max_terms) || length(max_terms) != 1 || is.na(max_terms) ||
    max_terms <= 0 || floor(max_terms) != max_terms) {
    stop("'max_terms' must be a positive integer.")
  }

  # ============= #
  #   err_red     #
  # ============= #

  if (!is.numeric(err_red) || length(err_red) != 1 || is.na(err_red) ||
    err_red < 0 || err_red > 1) {
    stop("'err_red' must be a numeric value in [0, 1].")
  }

  # ============= #
  #     span      #
  # ============= #

  # minspan must be -2, -1, or a non-negative integer
  if (!is.numeric(minspan) || length(minspan) != 1 || is.na(minspan) ||
    floor(minspan) != minspan || (minspan < -2)) {
    stop("'minspan' must be -2, -1, or a non-negative integer.")
  }

  # endspan must be -2, -1, or a non-negative integer
  if (!is.numeric(endspan) || length(endspan) != 1 || is.na(endspan) ||
    floor(endspan) != endspan || (endspan < -2)) {
    stop("'endspan' must be -2, -1, or a non-negative integer.")
  }

  # ============= #
  #   kn_grid     #
  # ============= #

  if (!identical(kn_grid, -1)) {
    if (!is.list(kn_grid)) {
      stop("'kn_grid' must be -1 (automatic knot selection) or a list of numeric knot vectors.")
    }
    if (length(kn_grid) != length(x)) {
      stop(paste0(
        "'kn_grid' must be a list of length ", length(x),
        " (one knot vector per input variable), but has length ", length(kn_grid), "."
      ))
    }
    for (k in seq_along(kn_grid)) {
      if (!is.numeric(kn_grid[[k]]) || length(kn_grid[[k]]) == 0) {
        stop(paste0(
          "'kn_grid[[", k, "]]' must be a non-empty numeric vector of knot candidates."
        ))
      }
    }
  }

  # ============= #
  #   kn_penalty  #
  # ============= #

  if (!is.numeric(kn_penalty) || length(kn_penalty) != 1 || is.na(kn_penalty) ||
    kn_penalty < 0) {
    stop("'kn_penalty' must be a non-negative numeric value.")
  }

  invisible(NULL)
}

#' @title Validate Efficiency Arguments
#'
#' @description
#' Checks arguments used to compute efficiency scores or targets and stops with
#' an informative message when an argument is invalid.
#'
#' @param data
#' A \code{data.frame} or \code{matrix} containing the DMUs to evaluate.
#'
#' @param x
#' Column indexes of input variables in \code{data}.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param object
#' An \code{aces} or \code{rf_aces} object.
#'
#' @param method
#' Prediction method used to define the technology.
#'
#' @param measure
#' Efficiency measure to compute.
#'
#' @param returns
#' Returns to scale: \code{"constant"} or \code{"variable"}.
#'
#' @param direction
#' Direction vectors used by the directional distance function.
#'
#' @return
#' \code{NULL} invisibly. The function is used for its validation side effects.

display_errors_scores <- function(
  data,
  x,
  y,
  object,
  method,
  measure,
  returns,
  direction
) {
  var_indexes <- sort(
    c(
      object[["data"]][["x"]],
      object[["data"]][["y"]]
    )
  )

  # check if training and test names are the same
  tr_names <- colnames(object[["data"]][["df"]])[var_indexes]
  ev_names <- colnames(data[c(x, y)])

  if (!all(ev_names %in% tr_names)) {
    stop("Different variable names in training and evaluated data.")
  }

  # object must be an aces or rf_aces object
  if (!inherits(object, c("aces", "rf_aces"))) {
    stop(paste(deparse(substitute(object)), "must be an aces or rf_aces object."))
  }

  # Valid methods per object class
  if (inherits(object, "rf_aces")) {
    allowed_methods <- c("rf_aces", "rf_aces_cubic", "rf_aces_quintic")
    if (!(method %in% allowed_methods)) {
      stop(paste0(
        "Invalid method '", method, "' for object of class 'rf_aces'. ",
        "Please choose one of: ", paste(allowed_methods, collapse = ", "), "."
      ))
    }
  }

  if (inherits(object, "aces")) {
    allowed_methods <- c("aces", "aces_cubic", "aces_quintic", "aces_forward")
    if (!(method %in% allowed_methods)) {
      stop(
        paste0(
          "Invalid method '", method, "' for object of class 'aces'. ",
          "Please choose one of: ", paste(allowed_methods, collapse = ", "), "."
        )
      )
    }
  }

  # Valid measures
  valid_measures <- c(
    "rad_out", "rad_inp", "ddf", "rsl_out", "rsl_inp",
    "wam_mip", "wam_nor", "wam_ram", "wam_bam",
    "rf_aces_rad_out"
  )

  if (!(measure %in% valid_measures)) {
    stop(
      paste0(
        "Invalid measure '", measure, "'. ",
        "Please choose one of: ", paste(valid_measures, collapse = ", "), "."
      )
    )
  }

  # Enforce compatibility of rf_aces_rad_out with class rf_aces only
  if (measure == "rf_aces_rad_out" && !inherits(object, "rf_aces")) {
    stop("Measure 'rf_aces_rad_out' can only be used with objects of class 'rf_aces'.")
  }

  # returns must be valid
  if (!is.null(returns) && !returns %in% c("constant", "variable")) {
    stop(paste(returns, "is not available. See help(\"get_scores\")."))
  }

  # direction must be a matrix with size [n x (nX+nY)]
  if (!is.null(direction)) {
    # Check if it is a valid data structure
    if (is.matrix(direction) || is.data.frame(direction)) {
      # Check rows: must match evaluated DMUs
      if (nrow(direction) != nrow(data)) {
        stop(paste("The direction matrix must have", nrow(data), "rows (one per DMU)."))
      }

      # Check columns: must be exactly nX + nY
      if (ncol(direction) != (length(x) + length(y))) {
        stop(paste("The direction matrix must have", length(x) + length(y), "columns (inputs + outputs)."))
      }

      # Optional but recommended: check if numeric
      if (!all(sapply(as.data.frame(direction), is.numeric))) {
        stop("All elements in the direction matrix must be numeric.")
      }
    } else {
      stop("Invalid 'direction' type. A numeric matrix or data.frame is required.")
    }
  } else if (measure == "ddf") {
    stop("Directional Distance Function (ddf) requires a 'direction' matrix.")
  }
}

#' @title Validate RF-ACES Arguments
#'
#' @description
#' Checks the arguments used by \code{\link{rf_aces}} and stops with an
#' informative message when an argument is invalid.
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
#' If \code{TRUE}, divide each input and output by its mean before fitting.
#'
#' @param quick_aces
#' If \code{TRUE}, use Quick ACES to reduce the number of candidate basis
#' functions in each learner.
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
#' Character string specifying the lack-of-fit measure. Options are
#' \code{"mae"}, \code{"mape"}, \code{"mse"}, \code{"msle"}, \code{"rmse"},
#' \code{"nrmse1"}, and \code{"nrmse2"}.
#'
#' @param learners
#' Number of ACES learners in the forest.
#'
#' @param bag_size
#' Number of observations sampled with replacement for each learner.
#'
#' @param max_feats
#' Number of inputs randomly considered by each learner.
#'
#' @param max_terms
#' Maximum number of terms created during the forward step.
#'
#' @param err_red
#' Minimum relative error reduction required to add a new pair of first-degree
#' basis functions.
#'
#' @param minspan
#' Minimum number of observations between adjacent knots. Use \code{-2},
#' \code{-1}, or a positive integer.
#'
#' @param endspan
#' Minimum number of observations before the first knot and after the last knot.
#' Use \code{-2}, \code{-1}, or a positive integer.
#'
#' @param kn_grid
#' Knot candidates. Use \code{-1} for automatic selection or supply a list with
#' one numeric vector for each input.
#'
#' @importFrom stats na.omit
#'
#' @return
#' \code{NULL} invisibly. The function is used for its validation side effects.

display_errors_rf_aces <- function(
  data,
  x,
  y,
  scale_data,
  quick_aces,
  max_degree,
  inter_cost,
  metric,
  learners,
  bag_size,
  max_feats,
  max_terms,
  err_red,
  minspan,
  endspan,
  kn_grid
) {
  # data is a matrix or a data.frame
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("data must be a data.frame or a matrix")
  }

  # x in data
  tryCatch(data[, x], error = function(e) {
    stop("Index values from x are not in data.")
  })

  # y in data
  tryCatch(data[, y], error = function(e) {
    stop("Index values from y are not in data.")
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

  if (quick_aces != TRUE && quick_aces != FALSE) {
    stop("quick_aces must be a boolean")
  }

  # learners must be greater than 0
  if (learners <= 0) {
    stop("learners must be greater than 0.")
  }

  # max_feats must be greater than 1 and lower than the number of inputs
  if (max_feats > length(x) || max_feats < 1) {
    stop("max_feats must be greater than 1 and lower than the number of inputs.")
  }

  # sample size must be lower training size
  if (bag_size > nrow(data)) {
    stop("sample size must be lower training size.")
  }

  # max_degree must be a valid number
  if (!is.list(max_degree) && max_degree > length(x)) {
    stop("max_degree must be lower than the number of inputs.")
  }

  # the lack-of-fit criterion must be a valid measure
  if (!metric %in% c("mae", "mape", "mse", "msle", "rmse", "nrmse1", "nrmse2")) {
    stop(paste(metric, "is not available. Please, check help(\"aces\")"))
  }

  # max_terms must be a positive integer
  if (!is.numeric(max_terms) || max_terms <= 0 || floor(max_terms) != max_terms) {
    stop("max_terms must be a positive integer.")
  }

  # err_red must be between 0 and 1
  if (!is.null(err_red) && !(err_red >= 0 && err_red <= 1)) {
    stop("err_red must be between 0 and 1.")
  }

  # inter_cost must be between 0 and 1
  if (!(inter_cost >= 0 && inter_cost <= 1)) {
    stop("inter_cost must be between 0 and 1.")
  }

  # minspan must -2, -1 or a positive integer
  if (minspan < -2 || floor(minspan) != minspan) {
    stop("minspan must be - 2, - 1 or a positive integer.")
  }

  # minspan must be -2, -1 or a positive integer
  if (!(minspan %in% c(-2, -1)) && (!is.numeric(minspan) || minspan < 0 || floor(minspan) != minspan)) {
    stop("minspan must be -2, -1 or a positive integer.")
  }

  # endspan must be -2, -1 or a positive integer
  if (!(endspan %in% c(-2, -1)) && (!is.numeric(endspan) || endspan < 0 || floor(endspan) != endspan)) {
    stop("endspan must be -2, -1 or a positive integer.")
  }

  # the kn_grid must have the same length that x
  if (is.list(kn_grid) && length(kn_grid) != length(x)) {
    stop("If kn_grid is entered, the length must be equal to length(x).")
  }
}
