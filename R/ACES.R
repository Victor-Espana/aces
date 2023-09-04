#' @title Fit Adaptive Constrained Enveloping Splines (ACES) model
#'
#' @description This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity, which is based upon the adaptation of the machine learning technique known as Multivariate Adaptive Regression Splines (MARS) developed by \insertCite{friedman1991;textual}{aces}.
#'
#' @name aces
#'
#' @param data A \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column indexes of input variables in \code{data}.
#' @param y Column indexes of output variables in \code{data}.
#' @param y_type Determines the prediction approach for \code{y}:
#' \itemize{
#' \item{\code{"individual"}} to predict one output at a time and use the rest of outputs as additional inputs (netputs).
#' \item{\code{"all"}} to predict all the outputs at the same time with the original set of inputs.
#' }
#' @param addi_x Column indexes of additional covariables in \code{data}.
#' @param degree A \code{list} with the input indexes for interaction of variables or a \code{numeric} value that determines the maximum degree of interaction between variables. Basis functions products are constrained to contain factors involving distinct variables to ensure interpretability and avoid multicollinearity.
#' @param metric The lack-of-fit criterion to evaluate the model performance:
#' \itemize{
#' \item{\code{"mae"}}: Mean Absolute Error.
#' \item{\code{"mape"}}: Mean Absolute Percentage Error.
#' \item{\code{"mse"}}: Mean Squared Error.
#' \item{\code{"rmse"}}: Root Mean Squared Error.
#' \item{\code{"nrmse1"}}: Normalized Root Mean Squared Error (using mean).
#' \item{\code{"nrmse2"}}: Normalized Root Mean Squared Error (using range).
#' \item{\code{"wape"}}: Weighted Absolute Percentage Error.
#' }
#' @param turbo Hyperparameter aimed at selecting the DMU near the DEA-VRS frontier.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#' @param nterms Maximum number of terms created before pruning. Default is \code{50}.
#' @param err_red Minimum reduced error rate for the addition of a new pair of 1-degree basis functions. Default is \code{0.01}.
#' @param hd_cost Minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.
#' @param minspan Minimum number of observations between two adjacent knots.
#' \itemize {
#' \item{\code{minspan = -2}}: computed as in \insertCite{zhang1994;textual}{aces}.
#' \item{\code{minspan = -1}}: computed as in \insertCite{friedman1991;textual}{aces}.
#' \item{\code{minspan = +m}}
#' }
#' @param endspan Minimum number of observations before the first and after the final knot.
#' \itemize {
#' \item{\code{endspan = -2}} Computed as in \insertCite{zhang1994;textual}{aces}.
#' \item{\code{endspan = -1}} Computed as in \insertCite{friedman1991;textual}{aces}.
#' \item{\code{endspan = +m}}
#' }
#' @param knots_grid Design of the grid of knots. If \code{knotsGrid == -1} (default), the original approach of \insertCite{friedman1991;textual}{aces} based on the observed data is used. Alternatively, a \code{list} with the grid of knots to perform ACES can be entered. Each element of the \code{list} contains a vector with the knots of one variable (e.g., the second element of the list contains the knot values of the second variable and so on).
#'
#' @param d Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}.
#'
#' @param wc For the cubic smoothing procedure \insertCite{friedman1991}{aces}. Let `d` be the distance between the central knot and the left side knot, and let be `e` be the distance between the central knot and the right side knot. If the condition `1 < d / e < 2` is not satisfied, then `e = wc * d`. It must be set between 1 and 2. If a \code{vector} is entered, the \code{wc} value that most reduce the lack-of-fit criterion is selected.
#' @param wq For the quintic smoothing procedure \insertCite{chen1999}{aces}. Let `d` be the distance between the central knot and the left side knot, and let be `e` be the distance between the central knot and the right side knot. If the condition `8/7 < e / d < 1.5` is not satisfied, then `d = wq * e`. It must be set between 8/7 and 1.5. If a \code{vector} is entered, the \code{wc} value that most reduce the lack-of-fit criterion is selected.
#'
#' @details
#' Arguments for the forward algorithm:
#' \itemize{
#' \item{\code{nterms}}
#' \item{\code{err_red}}
#' \item{\code{minspan}}
#' \item{\code{endspan}}
#' \item{\code{knotsGrid}}
#' }
#'
#' Arguments for the backward algorithm:
#' \itemize{
#' \item{\code{d}}
#' }
#'
#' Arguments for the smoothing procedure:
#' \itemize{
#' \item{\code{wc}}
#' \item{\code{wq}}
#' }
#'
#' @references
#' \insertRef{zhang1994}{aces} \cr
#' \cr
#' \insertRef{friedman1991}{aces} \cr
#' \cr
#' \insertRef{chen1999}{aces}
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc
#'
#' @return An \code{aces} object.
#'
#' @export

aces <- function (
    data, x, y, y_type = "individual", addi_x = NULL, degree = 1, metric = "mse",
    turbo = Inf, monotonicity = TRUE, concavity = TRUE, x0_y0 = FALSE,
    nterms = 50, err_red = 0.01, hd_cost = 0.5, minspan = 1, endspan = 1, knots_grid = - 1,
    d = 2, wc = seq(1, 2, length.out = 5), wq = seq(8 / 7, 1.5, length.out = 5)
    ) {

  print("prueba")

  # Possible error messages:
  display_errors (
    data = data,
    x = x,
    y = y,
    y_type = y_type,
    addi_x = addi_x,
    degree = degree,
    metric = metric,
    knots_grid = knots_grid,
    d = d,
    wc = wc,
    wq = wq,
    object = NULL,
    measure = NULL,
    returns = NULL,
    direction = NULL,
    digits = NULL
  )

  # discard observations far from the VRS-frontier
  ddf_scores <- ddf (
    tech_xmat = as.matrix(data[, x]),
    tech_ymat = as.matrix(data[, y]),
    eval_xmat = as.matrix(data[, x]),
    eval_ymat = as.matrix(data[, y]),
    direction = "briec",
    convexity = TRUE,
    returns = "variable"
  )

  # lower distance than 1 + turbo
  check1 <- ddf_scores <= turbo

  # we save TRUE indexes that satisfies check1
  turbo_index <- c(1:nrow(data))[check1]

  # we discard the observations far from the VRS-frontier
  data <- data[turbo_index, ]

  # List with individual ACES models
  if (y_type == "all") {
    ACES <- vector("list", 1)
  } else {
    ACES <- vector("list", length(y))
  }

  if (y_type == "all") {

    ACES[[1]] <- aces_algorithm (
      data = data,
      x = x,
      y = y,
      z = NULL,
      y_type = y_type,
      addi_x = addi_x,
      degree = degree,
      metric = metric,
      monotonicity = monotonicity,
      concavity = concavity,
      x0_y0 = x0_y0,
      nterms = nterms,
      err_red = err_red,
      hd_cost = hd_cost,
      minspan = minspan,
      endspan = endspan,
      knots_grid = knots_grid,
      d = d,
      wc = wc,
      wq = wq
      )

    # Name
    names(ACES)[1] <- "y_all"

  } else {

    for (out in 1:length(y)) {

      # Output indicator
      indx <- rep(FALSE, length(y))
      # Select the output
      indx[out] <- TRUE

      # Outputs as inputs
      nets <- y[!indx]
      outs <- y[indx]

      ACES[[out]] <- aces_algorithm (
        data = data,
        x = x,
        y = outs,
        z = nets,
        y_type = y_type,
        addi_x = addi_x,
        degree = degree,
        metric = metric,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0,
        nterms = nterms,
        err_red = err_red,
        hd_cost = hd_cost,
        minspan = minspan,
        endspan = endspan,
        knots_grid = knots_grid,
        d = d,
        wc = wc,
        wq = wq
      )

      # Name
      names(ACES)[out] <- names(data)[y[out]]
    }
  }

  class(ACES) <- "aces"

  return(ACES)
}

#' @title Algorithm of Adaptive Constrained Enveloping Splines (ACES)
#'
#' @description This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity, which is based upon the adaptation of the machine learning technique known as Multivariate Adaptive Regression Splines (MARS) developed by \insertCite{friedman1991;textual}{aces}.
#'
#' @param data A \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column indexes of input variables in \code{data}.
#' @param y Column indexes of output variables in \code{data}.
#' @param z Column indexes of netput variables in \code{data} (outputs evaluated as inputs). These variables must be considered as output when computing predictions or efficiency scores.
#' @param y_type Output prediction strategy: whether it should be done all at once or one at a time.
#' @param addi_x Column indexes of additional covariables in \code{data}.
#' @param degree Maximum degree of interaction between variables.
#' @param metric The lack-of-fit criterion to evaluate the model performance.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#' @param nterms For the Forward algorithm. Maximum number of terms created before pruning.
#' @param err_red For the Forward algorithm. Minimum reduced error rate for the addition of a new pair of 1-degree basis functions.
#' @param hd_cost Minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.
#' @param minspan For the Forward algorithm. Minimum number of observations between two adjacent knots.
#' @param endspan For the Forward algorithm. Minimum number of observations before the first and after the final knot.
#' @param knots_grid For the Forward algorithm. Design of the grid of knots.
#' @param d For the Backward algorithm. Generalized Cross Validation (GCV) penalty per knot.
#' @param wc Hyperparameter for the side knot distances in the cubic smoothing procedure \insertCite{friedman1991}{aces}.
#' @param wq Hyperparameter for the side knot distances in the quintic smoothing procedure \insertCite{chen1999}{aces}.
#'
#' @references
#' \insertRef{zhang1994}{aces} \cr
#' \cr
#' \insertRef{friedman1991}{aces} \cr
#' \cr
#' \insertRef{chen1999}{aces}
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc
#'
#' @return An \code{aces} object.
#'
#' @export

aces_algorithm <- function (
    data, x, y, z, y_type, addi_x, degree, metric, monotonicity, concavity,
    x0_y0, nterms, err_red, hd_cost, minspan, endspan, knots_grid, d, wc, wq
    ) {

  # Original indexes: inputs, outputs, netputs, environmental --> for aces_object()
  inps <- x
  outs <- y
  nets <- z
  envi <- addi_x

  # Original data --> for aces_object()
  dmus <- data

  # Data in [x, z, y] format with interaction of variables included
  data <- set_interactions (
    data = data,
    x = c(x, addi_x),
    y = y,
    z = z,
    degree = degree
    )

  # Samples in data
  N <- nrow(data)

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  # Number of inputs / outputs as inputs and number of outputs
  nX <- length(x)
  nY <- length(y)

  # matrix with index variables and degree
  xi_degree <- matrix(0, nrow = 2, ncol = length(x))

  # 1st row of xi_degree: variable index
  xi_degree[1, ] <- x

  # 2nd row of xi_degree: degree of variable index
  v <- 0
  for (i in 1:degree) {
    combs <- combn(1:(length(inps) + length(nets) + length(envi)), i)
    for (k in 1:ncol(combs)) {
      v <- v + 1
      xi_degree[2, v] <- i
    }
  }

  # ================= #
  # FORWARD ALGORITHM #
  # ================= #

  y_hat <- matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, max)

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
    "id"     = 1,
    "status" = "intercept",
    "side"   = "E",
    "Bp"     = rep(1, N),
    "xi"     = c(- 1),
    "t"      = c(- 1),
    "R"      = err_metric(data[, y, drop = F], y_hat, metric),
    "GCV"    = compute_gcv(data[, y, drop = F], y_hat, metric, matrix(rep(1, 50)), 1, 0, xi_degree),
    "GRSq"   = 0,
    "coefs"  = unname(apply(data[, y, drop = FALSE], 2, max))
  )

  # Set of knots. It saves indexes of data used as knots.
  kn_list <- vector("list", nX)

  # Set of basis functions by variable.
  Bp_list <- vector("list", nX)

  for (xi in 1:nX) {
    Bp_list[[xi]] <- list (
      "paired" = NULL,
      "right" = NULL,
      "left" = NULL
      )
  }

  # Set of basis functions (bf_set) and the matrix of basis functions (B)
  aces_forward <- list (
    "bf_set" = list(bf),
    "B" = matrix(rep(1, N))
    )

  # Error of the first basis function
  err <- bf[["R"]]

  # Minimum span (minspan) and end span (endspan)
  L_Le <- compute_span (
    data = data,
    minspan = minspan,
    endspan = endspan,
    nX = nX
    )

  # Minimum span
  L  <- L_Le[[1]]
  # End span
  Le <- L_Le[[2]]

  # Set the grid of knots
  kn_grid <- set_knots_grid (
    data = data,
    nX = nX,
    inps = length(inps),
    kn_grid = knots_grid
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
      xi_degree = xi_degree,
      y = y,
      metric = metric,
      monotonicity = monotonicity,
      concavity = concavity,
      x0_y0 = x0_y0,
      forward_model = aces_forward,
      kn_list = kn_list,
      Bp_list = Bp_list,
      L = L,
      Le = Le,
      kn_grid = kn_grid,
      err_min = err,
      hd_cost = hd_cost
      )

    if (!is.list(B_bf_knt_err)) {
      break
    } else {
      new_err <- B_bf_knt_err[[5]]
    }

    # New minimum error
    if (new_err[1] < err[1] * (1 - err_red[1])) {

      # Update B
      aces_forward[["B"]] <- B_bf_knt_err[[1]]
      # Update basis functions
      aces_forward[["bf_set"]] <- B_bf_knt_err[[2]]
      # Update the knots list
      kn_list <- B_bf_knt_err[[3]]
      # Update the Bp list
      Bp_list <-  B_bf_knt_err[[4]]
      # Update error
      err <- new_err

    } else {
      break
    }
  }

  # set of knots from forward algorithm
  var <- c(); knt <- c(); sts <- c()

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

  knots_forward <- data.frame(
    xi = var,
    t = knt,
    status = sts
    )

  # set of coefficients from forward algorithm
  coefs <- rev(aces_forward[["bf_set"]])[[1]][["coefs"]]

  # ==
  # forward aces
  # ==
  aces_forward = list (
    "basis" = aces_forward[["bf_set"]],
    "Bmatx" = aces_forward[["B"]],
    "knots" = knots_forward,
    "coefs" = coefs
  )

  # ================== #
  # BACKWARD ALGORITHM #
  # ================== #
  aces_submodels <- aces_pruning (
      data = data,
      xi_degree = xi_degree,
      y = y,
      metric = metric,
      monotonicity = monotonicity,
      concavity = concavity,
      x0_y0 = x0_y0,
      forward_model = aces_forward,
      Bp_list = Bp_list,
      d = d
    )

  # Generalized cross-validation for each model
  GCVs <- sapply(aces_submodels, function(x) x[["GCV"]])

  # Model with minimum error (excluding the model without knots)
  aces_backward <- aces_submodels[[which.min(GCVs[1:(length(aces_submodels) - 1)])]]

  # Set of surviving knots
  knots_backward <- do.call(rbind.data.frame, aces_backward[["t"]])

  # Sort the knots by "xi", "status", "side", "t"
  knots_backward_order <- order (
    knots_backward$xi,
    knots_backward$status,
    ifelse(knots_backward$status == "paired", knots_backward$t, desc(knots_backward$side)),
    ifelse(knots_backward$status == "paired", desc(knots_backward$side), knots_backward$t)
    )

  # Update the set of knots
  knots_backward <- knots_backward[knots_backward_order, ]

  # ==
  # aces
  # ==

  aces <- list (
    "aces_submodels" = aces_submodels,
    "Bmatx" = aces_backward[["B"]],
    "knots" = knots_backward,
    "coefs" = aces_backward[["coefs"]],
    "GCV" = aces_backward[["GCV"]]
  )

  # =================== #
  # SMOOTHING PROCEDURE #
  # =================== #

  # For the smoothing procedure: more restrictive approach for estimate coefficients

  # submodels sorted by gcv value
  aces_submodels_gcv <- aces_submodels[order(sapply(aces_submodels, "[[", "GCV"))]

  # initialize a list with smoothed submodels
  aces_smoothed_submodels <- vector("list", length(aces_submodels_gcv))

  for (s in 1:length(aces_submodels_gcv)) {

    # select a model to be smoothed
    aces_smoothed <- aces_submodels_gcv[[s]]

    if (is.null(aces_smoothed[["t"]])) next

    knots_smoothed <- do.call(rbind.data.frame, aces_smoothed[["t"]])

    # If monotonicity is required:
    # 1- wc in (1, 2) & wq in (8/7, 1.5)

    # If concavity is required:
    # 1- wc in (1, 2) & wq in (8/7, 1.5)
    # 2- unpaired right basis functions are not allowed

    check1 <- knots_smoothed$side == "R"
    check2 <- knots_smoothed$status == "unpaired"

    if (concavity && max(check1 + check2) == 2) {
      next
    } else {
      aces_smoothed_submodels[[s]][["Model"]] <- aces_smoothed
      aces_smoothed_submodels[[s]][["Knots"]] <- knots_smoothed
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
      knots_smoothed <- aces_smoothed_submodels[[m]][["Knots"]]

      # generate the space for side knots location
      x_space <- side_knots_location (
        data = data,
        nX = nX,
        knots = knots_smoothed
      )

      # ==
      # smoothing cubic aces
      # ==

      cubic_aces_models[[m]] <- cubic_aces (
        data = data,
        y = y,
        nX = nX,
        metric = metric,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0,
        knots = knots_smoothed,
        x_space = x_space,
        wc = wc,
        d = 0.25
      )

      # ==
      # smoothing quintic aces
      # ==

      quintic_aces_models[[m]] <- quintinc_aces (
        data = data,
        y = y,
        nX = nX,
        metric = metric,
        monotonicity = monotonicity,
        concavity = concavity,
        x0_y0 = x0_y0,
        knots = knots_smoothed,
        x_space = x_space,
        wq = wq,
        d = 0.25
      )
    }
  }

  # gcvs of cubic models
  aces_cubic_gcvs <- sapply(cubic_aces_models, function(x) ifelse(is.null(x[["GCV"]]), Inf, x[["GCV"]]))
  min_gcv <- which.min(aces_cubic_gcvs)

  # cubic aces
  aces_cubic <- cubic_aces_models[[min_gcv]]

  # gcvs of quintic models
  aces_quintic_gcvs <- sapply(quintic_aces_models, function(x) ifelse(is.null(x[["GCV"]]), Inf, x[["GCV"]]))
  min_gcv <- which.min(aces_quintic_gcvs)

  # quintic aces
  aces_quintic <- quintic_aces_models[[min_gcv]]

  # =========== #
  # aces object #
  # =========== #

  ACES <- aces_object (
    data = dmus,
    x = inps,
    y = outs,
    z = nets,
    y_type = y_type,
    addi_x = envi,
    degree = degree,
    metric = metric,
    monotonicity = monotonicity,
    concavity = concavity,
    x0_y0 = x0_y0,
    nterms = nterms,
    err_red = err_red,
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
#' @description This function saves information about the Adaptive Constrained Enveloping Splines (ACES) model.
#'
#' @param data A \code{matrix} containing the variables in the model.
#' @param x Column indexes of input variables in \code{data}.
#' @param y Column indexes of output variables in \code{data}.
#' @param z Column indexes of netput variables in \code{data} (outputs evaluated as inputs).
#' @param y_type Determines the prediction approach for \code{y}: whether it should be done all at once or one at a time.
#' @param addi_x Column indexes of additional covariables in \code{data}.
#' @param degree Maximum degree of interaction between variables.
#' @param metric The lack-of-fit criterion to evaluate the model performance.
#' @param monotonicity \code{logical} indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.
#' @param concavity \code{logical} indicating whether to enforce the constraint of concavity in the estimator.
#' @param x0_y0 \code{logical} indicating if f(0) = 0.
#' @param nterms Maximum number of terms created by the forward algorithm (before pruning).
#' @param err_red Minimum reduced error rate for the addition of two new basis functions.
#' @param minspan Minimum number of observations between two adjacent knots.
#' @param endspan Minimum number of observations before the first and after the final knot.
#' @param kn_grid Grid of knots to perform ACES.
#' @param d Generalized Cross Validation (GCV) penalty per knot.
#' @param wc Hyperparameter for the side knot distances in the cubic smoothing procedure \insertCite{friedman1991}{aces}.
#' @param wq Hyperparameter for the side knot distances in the quintic smoothing procedure \insertCite{chen1999}{aces}.
#' @param aces_forward A \code{list} containing the forward Adaptive Constrained Enveloping Splines model.
#' @param aces A \code{list} containing the Adaptive Constrained Enveloping Splines model.
#' @param aces_cubic A \code{list} containing the smoothed version of Adaptive Constrained Enveloping Splines through cubic basis functions.
#' @param aces_quintic A \code{list} containing the smoothed version of Adaptive Constrained Enveloping Splines through quintic basis functions.
#'
#' @references
#' \insertRef{friedman1991}{aces} \cr
#' \cr
#' \insertRef{chen1999}{aces}
#'
#' @return An \code{aces} object.

aces_object <- function (
    data, x, y, z, y_type, addi_x, degree, metric,
    monotonicity, concavity, x0_y0,
    nterms, err_red, minspan, endspan, kn_grid,
    d, wc, wq,
    aces_forward, aces, aces_cubic, aces_quintic
    ) {

  object <- list()

  object[["data"]] <- list (
    "df" = data,
    "x" = x,
    "z" = z,
    "y" = y,
    "addi_x" = addi_x,
    "xnames" = colnames(data)[c(x, z, addi_x)],
    "ynames" = colnames(data)[y],
    "y_type" = y_type,
    "rownames" = rownames(data)
  )

  object[["control"]] <- list (
    "degree" = degree,
    "metric" = metric,
    "monotonicity" = monotonicity,
    "concavity" = concavity,
    "x0_y0" = x0_y0,
    "nterms" = nterms,
    "err_red" = err_red,
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

#' @title Error Messaging in ACES functions.
#'
#' @description This function displays error messages if hyperparameters or data are bad introduced.
#'
#' @param data A \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column indexes of input variables in \code{data}.
#' @param y Column indexes of output variables in \code{data}.
#' @param y_type The mode of predicting the output: whether it should be done all at once or one at a time.
#' @param addi_x Column indexes of additional covariables in \code{data}.
#' @param degree Maximum degree of interaction between variables.
#' @param metric The lack-of-fit criterion to evaluate the model performance.
#' @param knots_grid Grid of knots to perform ACES.
#' @param d Generalized Cross Validation (GCV) penalty per knot.
#' @param wc Hyperparameter for the side knot distances in the cubic smoothing procedure \insertCite{friedman1991}{aces}.
#' @param wq Hyperparameter for the side knot distances in the quintic smoothing procedure \insertCite{chen1999}{aces}.
#' @param object An \code{aces} object.
#' @param measure Mathematical programming model to compute the efficiency scores.
#' @param returns Type of returns to scale for computing the efficiency scores.
#' @param direction Vector direction for DMU projection in Directional Distance Function when computing the efficiency scores.
#' @param digits Number of digits to round efficiency scores.
#'
#' @references
#' \insertRef{friedman1991}{aces} \cr
#' \cr
#' \insertRef{chen1999}{aces}
#'
#' @importFrom stats na.omit
#'
#' @return Possible error messages.

display_errors <- function (
    data, x, y, y_type, addi_x, degree, metric, knots_grid, d, wc, wq,
    object, measure, returns, direction, digits
    ) {

  # Function that calls display_errors
  caller <- as.character(sys.calls()[[1]])[[1]]

  if (caller == "aces") {

    # x in data
    tryCatch(data[, x], error = function(e) {
      message("Index values from x are not in data.")
    })

    # y in data
    tryCatch(data[, y], error = function(e) {
      message("Index values from y are not in data.")
    })

    # data is a matrix or a data.frame
    if (is.list(data) && !is.data.frame(data)) {
      stop("data must be a data.frame or a matrix")
    }

    # variables classes are valid
    if (!all(sapply(data, is.numeric))) {
      stop("data variables must be numeric. Please, check sapply(data, is.numeric))")
    }

    # NA values
    if (any(is.na(data))) {
      data <- na.omit(data)
      warning("Rows with NA values have been omitted .\n")
    }

    # y_type must be "individual" or "all"
    if (!y_type %in% c("individual", "all")) {
      stop("Not available y_type. Please, check help(\"aces\")")
    }

    # degree must be a valid number
    if (!is.list(degree)) {
      case1 <- y_type == "individual" && degree > length(x) + length(addi_x) + length(y) - 1
      case2 <- y_type == "all" && degree > length(x) + length(addi_x)

      if (case1 || case2) {
        stop("degree must be lower than the number of inputs.")
      }
    }

    # the lack-of-fit criterion must be a valid measure
    if (!metric %in% c("mae", "mape", "mse", "rmse", "nrmse1", "nrmse2", "wape")) {
      stop(paste(metric, "is not available. Please, check help(\"aces\")"))
    }

    # the knots_grid must have the same length that x
    if (is.list(knots_grid) && length(knots_grid) != length(x)) {
      stop ("If knots_grid is entered, it must have the same length that x.")
    }

    # d must be a semi-positive number
    if (any(d < 0)) {
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

    if (!class(object) == "aces") {
      stop(paste(deparse(substitute(object)), "must be an aces object."))
    }

    if (!measure %in% c("rad_out", "rad_inp", "ddf", "rsl_out", "rsl_inp", "wam")) {
      stop(paste(measure, "is not available. Please, check help(\"aces_scores\")"))
    }

    if (!returns %in% c("constant", "variable")) {
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
    data, x, y, z, degree
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
      name <- paste("IVar", p, sep = "")
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
#' @description This function computes an error metric for model evaluation.
#'
#' @param y_obs Vector of observed data.
#' @param y_hat Vector of predicted values.
#' @param metric Lack-of-fit criterion to evaluate the model performance:
#' \itemize{
#' \item{\code{mae}}: Mean Absolute Error
#' \item{\code{mape}}: Mean Absolute Percentage Error
#' \item{\code{mse}}: Mean Squared Error
#' \item{\code{rmse}}: Root Mean Squared Error
#' \item{\code{nrmse1}}: Normalized Root Mean Squared Error (using mean)
#' \item{\code{nrmse2}}: Normalized Root Mean Squared Error (using range)
#' \item{\code{wape}}: Weighted Absolute Percentage Error
#' }
#'
#' @return The calculated error metric.

err_metric <- function(y_obs, y_hat, metric){

  # samples in data
  N <- nrow(y_obs)

  # number of outputs
  nY <- ncol(y_obs)

  if (metric == "mae") {
    # mean absolute error
    error <- sum(abs(y_hat - y_obs)) / (N * nY)

  } else if (metric == "mape") {
    # mean absolute percentage error
    error <- sum(abs(y_hat - y_obs) / y_obs) / (N * nY) * 100

  } else if (metric == "mse") {
    # mean squared error
    error <- sum((y_hat - y_obs) ^ 2) / (N * nY)

  } else if (metric == "rmse") {
    # root mean squared error
    error <- sqrt(sum((y_hat - y_obs) ^ 2) / (N * nY))

  } else if (metric == "nrmse1") {
    # normalized root mean squared error by the mean
    error <- sqrt(sum((y_hat - y_obs) ^ 2) / (N * nY)) / mean(y_obs)

  } else if (metric == "nrmse2") {
    # calculate the mean of column-wise maximums and minimums in y
    ymax <- mean(apply(y_obs, 2, max))
    ymin <- mean(apply(y_obs, 2, min))

    # normalized root mean squared error by the range
    error <- sqrt(sum((y_hat - y_obs) ^ 2) / (N * nY)) / (ymax - ymin)

  } else {
    # weighted absolute percentage error
    error <- sum(abs(y_hat - y_obs)) / (sum(abs(y_obs)))
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

set_knots_grid <- function(data, nX, inps, kn_grid) {

  # Case 1: knots_grid is provided (list) and new variables are created (nX > inputs):
    # expand the knots_grid list.
  # Case 2: knots_grid is provided (list) and new variables are not created (nX = inputs):
    # keep the same knots_grid list.
  # Case 3: knots_grid is not provided:
    # create the knots_grid list.

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
      kn_grid <- knots_grid
    }

  } else {
    kn_grid <- lapply(1:nX, function(i) data[, i])
  }

  return(kn_grid)
}
