#' @title Error Messaging in Adaptive Constrained Enveloping Splines (ACES) functions.
#'
#' @description
#' This function displays error messages if hyperparameters are incorrectly specified in ACES functions.
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
#' A \code{character} string specifying the nature of the production frontier that the function estimates. It can be either: \code{"env"} or \code{"sto"}.
#'
#' @param error_type
#'  A \code{character} string specifying the error structure that the function will use when fitting the model. It can be either: \code{"add"} or \code{"mul"}.
#'
#' @param nvars
#' An \code{integer} indicating the number of variables randomly chosen at each split in RF-ACES.
#'
#' @param degree
#' A \code{list} with input indexes for interaction of variables, or a \code{numeric} specifying the maximum degree of interaction.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance. Options are: \code{"mae"}, \code{"mape"}, \code{"mse"}, \code{"rmse"}, \code{"nrmse1"} or \code{"nrmse2"}.
#'
#' @param nterms
#' A positive \code{integer} specifying the maximum number of terms created before pruning.
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree basis functions.
#'
#' @param hd_cost
#' A \code{numeric} value specifying the minimum percentage of improvement over the best 1-degree basis function to incorporate a higher degree basis function.
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
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
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
    convexity,
    returns,
    direction,
    digits
    ) {

  if (caller == "aces") {

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
    nets <- ifelse(y_type == "ind", length(y) - 1, 0)

    if (!is.null(nvars) && nvars > (length(x) + nets)) {
      stop("nvars must be lower than the number of inputs.")
    }

    # degree must be a valid number
    if (!is.list(degree) && degree > length(x)) {
      stop("degree must be lower than the number of inputs.")
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

    # hd_cost must be between 0 and 1
    if (!(hd_cost >= 0 && hd_cost <= 1)) {
      stop("hd_cost must be between 0 and 1.")
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

    # object must be an aces or rf_aces object
    if (!inherits(object, c("aces", "rf_aces"))) {
      stop(paste(deparse(substitute(object)), "must be an aces or rf_aces object."))
    }

    # convexity must be assumed
    if (!is.null(convexity) && !convexity) {
      stop(paste(convexity, "is not available. Please, check help(\"aces_scores\")"))
    }

    # returns must be valid
    if (!is.null(returns) && !returns %in% c("constant", "variable")) {
      stop(paste(returns, "is not available. Please, check help(\"aces_scores\")"))
    }

    # direction must be valid
    if (!is.null(direction) && !direction %in% c("mean", "briec")) {
      stop(paste(direction, "is not available. Please, check help(\"aces_scores\")"))
    }

    # digits must be a non-negative integer
    if (!is.numeric(digits) || digits < 0 || floor(digits) != digits) {
      stop("digits must be a non-negative integer.")
    }

  }
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
