#' @title Error Messaging in Adaptive Constrained Enveloping Splines (ACES) functions
#'
#' @description
#' This function displays error messages if hyperparameters are incorrectly specified in ACES functions aimed at estimating production functions.
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
#' Column indexes of contextual variables in \code{data}.
#'
#' @param quick_aces
#' A \code{logical} indicating if the fast version of ACES should be employed.
#'
#' @param model_type
#' A \code{character} string specifying the nature of the estimated production frontier. It can be either: \code{"envelopment"} or \code{"stochastic"}.
#'
#' @param error_type
#' A \code{character} string specifying the error structure when fitting the model. It can be either: \code{"add"} or \code{"mul"}.
#'
#' @param max_degree
#' A \code{list} with input indexes for interaction of variables, or a \code{numeric} specifying the maximum max_degree of interaction.
#'
#' @param compl_cost
#' A \code{numeric} value specifying the minimum percentage of improvement over the best 1-degree BF to incorporate a higher degree BF.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion employed to evaluate the model performance. Options are: \code{"mae"}, \code{"mape"}, \code{"mse"}, \code{"rmse"}, \code{"nrmse1"} or \code{"nrmse2"}.
#'
#' @param nterms
#' A positive \code{integer} specifying the maximum number of terms created during the forward step.
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree BFs.
#'
#' @param minspan
#' A \code{numeric} value specifying the minimum number of observations between two adjacent knots. Options are: \code{"-2"}, \code{"-1"} or \code{"m"}.
#'
#' @param endspan
#' A \code{numeric} value specifying the minimum number of observations before the first and after the final knot. Options are: \code{"-2"}, \code{"-1"} or \code{"m"}.
#'
#' @param kn_grid
#' Grid of knots to perform ACES. It can be: \code{-1} or a \code{list}.
#'
#' @param kn_penalty
#' A positive \code{numeric} value specifying the Generalized Cross Validation penalty per knot.
#'
#' @importFrom stats na.omit
#'
#' @return
#' This function returns error messages if hyperparameters are incorrectly specified.

display_errors_aces <- function (
    data,
    x,
    y,
    z,
    quick_aces,
    error_type,
    max_degree,
    compl_cost,
    metric,
    nterms,
    err_red,
    minspan,
    endspan,
    kn_grid,
    kn_penalty
    ) {

  # data is a matrix or a data.frame
  if (!is.data.frame(data) && !is.matrix(data)) {
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

  # z in data
  if (!is.null(z)) {
    tryCatch(data[, z], error = function(e) {
      message("Index values from z are not in data.")
    })
  }

  # variables classes are valid
  if (!all(sapply(data[, c(x, y)], is.numeric))) {
    stop("data variables must be numeric. Please, check sapply(data, is.numeric))")
  }

  # NA values
  if (any(is.na(data[, c(x, y, z)]))) {
    data <- na.omit(data[, c(x, y, z)])
    warning("Rows with NA values have been omitted .\n")
  }

  if (quick_aces != TRUE && quick_aces != FALSE) {
    stop("quick_aces must be a boolean")
  }

  # error_type must be "add" or "mul"
  if (!is.null(error_type) && !error_type %in% c("add", "mul")) {
    stop("Not available error_type. Please, check help(\"aces\")")
  }

  # max_degree must be a valid number
  if (!is.list(max_degree) && max_degree > length(x)) {
    stop("max_degree must be lower than the number of inputs.")
  }

  # the lack-of-fit criterion must be a valid measure
  if (!metric %in% c("mae", "mape", "mse", "msle", "rmse", "nrmse1", "nrmse2")) {
    stop(paste(metric, "is not available. Please, check help(\"aces\")"))
  }

  # nterms must be a positive integer
  if (!is.numeric(nterms) || nterms <= 0 || floor(nterms) != nterms) {
    stop("nterms must be a positive integer.")
  }

  # err_red must be between 0 and 1
  if (!is.null(err_red) && !(err_red >= 0 && err_red <= 1)) {
    stop("err_red must be between 0 and 1.")
  }

  # compl_cost must be between 0 and 1
  if (!(compl_cost >= 0 && compl_cost <= 1)) {
    stop("compl_cost must be between 0 and 1.")
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
  if (is.list(kn_grid) && length(kn_grid) != length(x) + length(z)) {
    stop ("If kn_grid is entered, the length must be equal to length(x + z).")
  }

  # kn_penalty must be a semi-positive number
  if (!is.null(kn_penalty) && kn_penalty < 0) {
    stop("kn_penalty must be greater than 0.")
  }

}

#' @title Error Messaging in Adaptive Constrained Enveloping Splines (ACES) functions II
#'
#' @description
#' This function displays error messages if hyperparameters are incorrectly specified in ACES functions aimed at computing efficiency scores.
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
#' @param weights
#' Weights for the additive model.
#'
#' @param digits
#' Number of digits to round efficiency scores.
#'
#' @importFrom stats na.omit
#'
#' @return
#' This function returns error messages if hyperparameters are incorrectly specified.

display_errors_scores <- function (
    data,
    x,
    y,
    object,
    measure,
    returns,
    direction,
    weights,
    digits
    ) {

  # index of variables
  # inputs are common in all the models
  # outputs and netputs are common (but interchangeable) in all the models
  # then, we can select just the first element of the ACES object.

  var_indexes <- sort (
    c (
      object[["data"]][["x"]],
      object[["data"]][["y"]]
    ))

  # check if training and test names are the same
  tr_names <- colnames(object[["data"]][["df"]])[var_indexes]
  ev_names <- colnames(data[c(x, y)])

  if (!all(ev_names %in% tr_names)) {
    stop("Different variable names in training and evaluated data.")
  }

  # object must be an aces, s_aces or rf_aces object
  if (!inherits(object, c("aces", "s_aces", "rf_aces"))) {
    stop(paste(deparse(substitute(object)), "must be an aces, s_aces or rf_aces object."))
  }

  # returns must be valid
  if (!is.null(returns) && !returns %in% c("constant", "variable")) {
    stop(paste(returns, "is not available. Please, check help(\"aces_scores\")"))
  }

  # direction must be valid
  if (!is.null(direction) && !direction %in% c("mean", "briec")) {
    stop(paste(direction, "is not available. Please, check help(\"aces_scores\")"))
  }

  # weights must be valid
  if (!is.null(weights) && !weights %in% c("WAM", "MIP", "NOR", "RAM", "BAM")) {
    stop(paste(weights, "is not available. Please, check help(\"aces_scores\")"))
  }

  # digits must be a non-negative integer
  if (!is.numeric(digits) || digits < 0 || floor(digits) != digits) {
    stop("digits must be a non-negative integer.")
  }
}

#' @title Error Messaging in Random-Forest Adaptive Constrained Enveloping Splines (RF-ACES) functions
#'
#' @description
#' This function displays error messages if hyperparameters are incorrectly specified in RF-ACES functions aimed at estimating production functions.
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
#' Column indexes of contextual variables in \code{data}.
#'
#' @param quick_aces
#' A \code{logical} indicating if the fast version of ACES should be employed.
#'
#' @param error_type
#' A \code{character} string specifying the error structure that the function will use when fitting the model. It can be either: \code{"add"} or \code{"mul"}.
#'
#' @param learners
#' An \code{integer} indicating the number of models for bagging.
#'
#' @param nvars
#' An \code{integer} indicating the number of variables randomly chosen at each split in RF-ACES.
#'
#' @param sample_size
#' An \code{integer} indicating the number of samples to draw from \code{data} to train each base estimator.
#'
#' @param max_degree
#' A \code{list} with input indexes for interaction of variables, or a \code{numeric} specifying the maximum max_degree of interaction.
#'
#' @param compl_cost
#' A \code{numeric} value specifying the minimum percentage of improvement over the best 1-max_degree basis function to incorporate a higher max_degree basis function.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance. Options are: \code{"mae"}, \code{"mape"}, \code{"mse"}, \code{"rmse"}, \code{"nrmse1"} or \code{"nrmse2"}.
#'
#' @param nterms
#' A positive \code{integer} specifying the maximum number of terms created before pruning.
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-max_degree basis functions.
#'
#' @param minspan
#' A \code{numeric} value specifying the minimum number of observations between two adjacent knots. Options are: \code{"-2"}, \code{"-1"} or \code{"m"}.
#'
#' @param endspan
#' A \code{numeric} value specifying the minimum number of observations before the first and after the final knot. Options are: \code{"-2"}, \code{"-1"} or \code{"m"}.
#'
#' @param kn_grid
#' Grid of knots to perform RF-ACES. It can be: \code{-1} or a \code{list}.
#'
#' @importFrom stats na.omit
#'
#' @return
#' This function return error messages if hyperparameters are incorrectly specified.

display_errors_rf_aces <- function (
    data,
    x,
    y,
    z,
    quick_aces,
    error_type,
    learners,
    nvars,
    sample_size,
    max_degree,
    compl_cost,
    metric,
    nterms,
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
    message("Index values from x are not in data.")
  })

  # y in data
  tryCatch(data[, y], error = function(e) {
    message("Index values from y are not in data.")
  })

  # z in data
  if (!is.null(z)) {
    tryCatch(data[, z], error = function(e) {
      message("Index values from z are not in data.")
    })
  }

  # variables classes are valid
  if (!all(sapply(data[, c(x, y)], is.numeric))) {
    stop("data variables must be numeric. Please, check sapply(data, is.numeric))")
  }

  # NA values
  if (any(is.na(data[, c(x, y, z)]))) {
    data <- na.omit(data[, c(x, y, z)])
    warning("Rows with NA values have been omitted .\n")
  }

  if (quick_aces != TRUE && quick_aces != FALSE) {
    stop("quick_aces must be a boolean")
  }

  # error_type must be "add" or "mul"
  if (!is.null(error_type) && !error_type %in% c("add", "mul")) {
    stop("Not available error_type. Please, check help(\"aces\")")
  }

  # learners must be greater than 0
  if (learners <= 0) {
    stop("learners must be greater than 0.")
  }

  # nvars must be lower than the number of inputs
  if (nvars > length(x)) {
    stop("nvars must be lower than the number of inputs.")
  }

  # sample size must be lower training size
  if (sample_size > nrow(data)) {
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

  # nterms must be a positive integer
  if (!is.numeric(nterms) || nterms <= 0 || floor(nterms) != nterms) {
    stop("nterms must be a positive integer.")
  }

  # err_red must be between 0 and 1
  if (!is.null(err_red) && !(err_red >= 0 && err_red <= 1)) {
    stop("err_red must be between 0 and 1.")
  }

  # compl_cost must be between 0 and 1
  if (!(compl_cost >= 0 && compl_cost <= 1)) {
    stop("compl_cost must be between 0 and 1.")
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
  if (is.list(kn_grid) && length(kn_grid) != length(x) + length(z)) {
    stop ("If kn_grid is entered, the length must be equal to length(x + z).")
  }
}

