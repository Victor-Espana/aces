#' @title Fit a Stochastic Adaptive Constrained Enveloping Splines (S-ACES) model
#'
#' @description
#'
#' This function estimates a stochastic production frontier that adheres to classical production theory axioms, such as monotonicity and concavity. The estimation is based on adaptations of the Multivariate Adaptive Regression Splines (MARS) technique, developed by \insertCite{friedman1991;textual}{aces}. For comprehensive details on the methodology and implementation, please refer to \insertCite{espana2024;textual}{aces}.
#'
#' @name s_aces
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
#' A \code{character} string specifying the error structure when fitting the model. Options are:
#' \itemize{
#'   \item{\code{"add"}}: Additive error structure.
#'   \item{\code{"mul"}}: Multiplicative error structure.
#' }
#'
#' @param mul_BF
#' A \code{list} specifying the maximum degree of basis functions (BFs) and the cost of introducing a higher-degree BF. Items include:
#' \itemize{
#' \item{\code{max_degree}}: A \code{list} with input indexes for interaction of variables, or a \code{numeric} specifying the maximum degree of interaction. BFs products are constrained to contain factors involving distinct variables to ensure interpretability and avoid multicollinearity.
#' \item{\code{compl_cost}}: A \code{numeric} specifying the minimum percentage improvement over the best 1-degree BF to incorporate a higher degree BF. Default is \code{0.05}.
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
#'   \item{\code{ptto}}: A \code{logical} indicating if (0, 0) should be included in the technology (only for piece-wise linear version).
#' }
#'
#' @param nterms
#' A positive \code{integer} specifying the maximum number of terms created during the forward step. Default is \code{50}.
#'
#' @param err_red
#' A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree BFs. Default is \code{0.01}.
#'
#' @param kn_grid
#' Specifies the grid of knots to be used in the ACES algorithm. Accepts two options:
#' \itemize{
#'    \item{\code{-1}} (default): Uses the original approach by \insertCite{friedman1991;textual}{aces}, where the knots are automatically selected based on the observed data.
#'    \item{\code{list}}: A user-defined grid of knots for each variable. This must be a \code{list} where each element corresponds to a vector of knots for a specific variable (e.g., the first element of the \code{list} contains the knot values for the first variable, the second element contains the knots for the second variable, and so on). This option allows for greater control over the knot placement when customizing the estimation process.
#'
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
#' @param kn_penalty
#' A positive \code{numeric} value specifying the Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}.
#'
#' @details
#'
#' This function generates a production frontier that adheres to classical production theory axioms, such as monotonicity and concavity. The algorithm comprises three main procedures:
#'
#' Forward Selection Algorithm. This initial step constructs a set of linear BFs to model the production frontier. While it may initially overfit the training data by capturing complex patterns and interactions, it provides a comprehensive foundation for the model.
#'
#' Backward Elimination Algorithm. To enhance model generalization and prevent overfitting, this procedure systematically removes BFs that do not significantly contribute to the model's performance. It retains only those functions that have a meaningful impact on predicting the output.
#'
#' Smoothing Procedures. After refining the set of BFs, the algorithm offers two smoothing options to produce a smooth and continuous production frontier:
#'
#' * Cubic Smoothing: Applies cubic functions to achieve a balance between flexibility and smoothness, suitable for capturing moderate curvature in the data.
#'
#' * Quintic Smoothing: Utilizes quintic functions for a higher degree of smoothness and continuity, ideal for modelling more complex relationships.
#'
#' The frontier's shape is estimated without imposing enveloping constraints on the observations, allowing for random noise and statistical variability in the data. In the second stage, the expected value of inefficiency is estimated using the residuals obtained from the first stage. This provides insights into the deviation of actual production from the optimal frontier, accounting for inefficiency and random error.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{friedman1991}{aces} \cr \cr
#' \insertRef{zhang1994}{aces} \cr
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc group_by across all_of summarise
#' @importFrom extraDistr rsign
#'
#' @return
#'
#' A \code{s_aces} object.
#'
#' @export

s_aces <- function (
    data,
    x,
    y,
    z = NULL,
    quick_aces = TRUE,
    error_type = "add",
    mul_BF = list (
      "max_degree" = 1,
      "compl_cost" = 0.05
    ),
    metric = "mse",
    shape = list (
      "mono" = TRUE,
      "conc" = TRUE,
      "ptto" = FALSE,
    ),
    nterms = 50,
    err_red = 0.01,
    kn_grid = - 1,
    minspan = - 1,
    endspan = - 1,
    kn_penalty = 2
    ) {

  # possible error messages:
  display_errors_aces (
    data = data,
    x = x,
    y = y,
    z = z,
    quick_aces = quick_aces,
    error_type = error_type,
    max_degree = mul_BF[["max_degree"]],
    compl_cost = mul_BF[["compl_cost"]],
    metric = metric,
    nterms = nterms,
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = kn_penalty
  )

  if (shape[["ptto"]]) {
    data <- rbind(data, rep(1e-10, ncol(data)))
  }

  S_ACES <- s_aces_algorithm (
    data = data,
    x_vars = x,
    y_vars = y,
    z_vars = z,
    quick_aces = quick_aces,
    error_type = error_type,
    max_degree = mul_BF[["max_degree"]],
    compl_cost = mul_BF[["compl_cost"]],
    metric = metric,
    shape = list (
      "mono" = shape[["mono"]],
      "conc" = shape[["conc"]],
      "ptto" = shape[["ptto"]]
    ),
    nterms = nterms,
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = kn_penalty
  )

  # type of object
  class(S_ACES) <- "s_aces"

  return(S_ACES)

}

#' @title Algorithm of Stochastic Adaptive Constrained Enveloping Splines (S-ACES).
#'
#' @description
#' This function implements the Stochastic Adaptive Constrained Enveloping Splines (S-ACES) algorithm.
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
#' @param z_vars
#' Column indexes of contextual variables in \code{data}.
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
#' @param compl_cost
#' Minimum percentage of improvement over the best 1 degree BF to incorporate a higher degree BF.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion employed to evaluate the model performance.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.
#'
#' @param nterms
#' Maximum number of terms created during the forward step.
#'
#' @param err_red
#' Minimum reduced error rate for the addition of a new pair of 1-degree BFs.
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
#' @param kn_penalty
#' Penalty per knot for computing Generalized Cross Validation.
#'
#' @importFrom Rdpack reprompt
#' @importFrom dplyr desc
#'
#' @return
#'
#' A \code{s-aces} object.
#'
#' @export

s_aces_algorithm <- function (
    data,
    x_vars,
    y_vars,
    z_vars,
    quick_aces,
    error_type,
    max_degree,
    compl_cost,
    metric,
    shape,
    nterms,
    err_red,
    minspan,
    endspan,
    kn_grid,
    kn_penalty
    ) {

  # save a copy of the original data
  DMUs <- data

  # data in [x, z, y] format with interaction of variables included
  data <- prepare_data (
    data = data,
    x = x_vars,
    y = y_vars,
    z = z_vars,
    max_degree = max_degree,
    error_type = error_type
  )

  # samples size
  N <- nrow(data)

  # set "x", "z" and "y" indexes in data
  x <- 1:(ncol(data) - length(z_vars) - length(y_vars))

  if (!is.null(z_vars)) {
    z <- (length(x) + 1):(ncol(data) - length(y_vars))
  } else {
    z <- integer(0)
  }

  y <- (length(x) + length(z) + 1):ncol(data)

  # set number of inputs, contextual variables and outputs
  nX <- length(x)
  nZ <- length(z)
  nY <- length(y)

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

  y_hat <- matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, mean)

  # lack-of-fit
  LOF <- err_metric (
    y_obs = data[, y, drop = F],
    y_hat = y_hat,
    metric = metric,
    weight = rep(1, N)
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
  kn_list <- vector("list", nX + nZ)

  # set of basis functions by variable
  Bp_list <- vector("list", nX + nZ)

  for (xi in 1:(nX + nZ)) {
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
  err <- bf[["LOF"]]

  # set the grid of knots
  kn_grid <- set_knots_grid (
    data = data,
    n_input_1 = length(x_vars) + length(z_vars),
    n_input_2 = nX + nZ,
    nZ = nZ,
    kn_grid = kn_grid,
    quick_aces = quick_aces,
    de = table_scores[, 1]
  )

  # minimum span (minspan) and end span (endspan)
  L_Le <- compute_span (
    kn_grid = kn_grid,
    minspan = minspan,
    endspan = endspan,
    n_input = nX + nZ
  )

  # list to save technologies created through S_ACES
  technology <- list()

  # initial error
  err_min <- err

  # variable importance
  var_imp <- matrix (
    rep(0, nX + nZ),
    nrow = 1
  )
  colnames(var_imp) <- colnames(data)[1:(nX + nZ)]

  # remove variables with low correlation
  if (quick_aces) {

    # Spearman’s Rank Correlation
    spearman_corr <- cor(data, method = "spearman")[1:length(x), (length(x) + 1):ncol(data)]

    # Kendall’s Tau
    kendall_corr <- cor(data, method = "kendall")[1:length(x), (length(x) + 1):ncol(data)]

    # compute threshold for removing variables
    threshold_1 <- min(0.1, quantile(spearman_corr[spearman_corr > 0], probs = 0.2))
    threshold_2 <- min(0.1, quantile(kendall_corr[kendall_corr > 0], probs = 0.2))

    # iterate over the variables
    for (j in 1:(nX + nZ)) {

      # check for both Spearman and Kendall correlations
      if (spearman_corr[j] < threshold_1 & kendall_corr[j] < threshold_2) {

        var_imp[1, j] <- - 1

      }
    }
  }

  while(length(aces_forward[["bf_set"]]) + 2 < nterms) {

    # add 2 new basis functions to the model:
    B_bf_knt_err <- add_basis_function (
      data = data,
      x = x,
      y = y,
      z = z,
      xi_degree = xi_degree,
      compl_cost = compl_cost,
      model_type = "stochastic",
      dea_scores = rep(1, N),
      metric = metric,
      forward_model = aces_forward,
      Bp_list = Bp_list,
      shape = shape,
      kn_list = kn_list,
      kn_grid = kn_grid,
      span = c(L_Le[[1]], L_Le[[2]]),
      err_min = err,
      var_imp = var_imp,
      quick_aces = quick_aces
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

      # updated error
      err <- new_err

      # updated variable importance matrix
      var_imp <- B_bf_knt_err[[6]]
      var_imp <- rbind(var_imp, rep(0, nX + nZ))

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
    x = x,
    y = y,
    z = z,
    xi_degree = xi_degree,
    model_type = "stochastic",
    dea_scores = rep(1, N),
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
  knots_backward_order <- with (
    kn_backward, {
      order (
        xi,
        status,
        ifelse(status == "paired", t, desc(side)),
        ifelse(status == "paired", desc(side), t))
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
        dea_scores = rep(1, N),
        model_type = "stochastic",
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

      quintic_aces_models[[m]] <- quintic_aces (
        data = data,
        x = x,
        y = y,
        dea_scores = rep(1, N),
        model_type = "stochastic",
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

  model_list <- list (
    aces_forward = aces_forward,
    aces = aces,
    aces_cubic = aces_cubic,
    aces_quintic = aces_quintic
  )

  for (model in names(model_list)) {

    # =========
    # Residuals
    # =========

    # model's prediction
    y_hat <- model_list[[model]][["Bmatx"]] %*% model_list[[model]][["coefs"]]

    # change scale if necessary
    mean_pred <- if (error_type == "mul") exp(y_hat) else y_hat

    # DEA model to ensure properties
    mean_pred <- rad_out (
      tech_xmat = as.matrix(DMUs[, x_vars]),
      tech_ymat = as.matrix(mean_pred),
      eval_xmat = as.matrix(data[, x_vars]),
      eval_ymat = as.matrix(mean_pred),
      convexity = TRUE,
      returns = "variable"
    ) * mean_pred

    # change scale (again) if necessary
    mean_pred <- if (error_type == "mul") log(mean_pred) else mean_pred

    # compute residuals
    residuals <- data[, y] - mean_pred

    # H0: residuals are symmetrically distributed
    # H1: residuals are negatively skewed

    sqrt_b1_btp <- c()
    for (m in 1:100) {
      # number of residuals
      nresid <- length(residuals)

      # stage 1: re-centred residuals
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

    # recompute M3 if necessary
    M3 <- ifelse(M3 > 0, - 0.0001, M3)

    # compute statistic
    sqrt_b1 <- M3 / M2 ^ (3 / 2)

    # shape of residuals
    residual_shape <- ifelse(sqrt_b1 > critical_value, "symmetric", "left-skewed")

    # warning (
    #   "The statistical evidence is not sufficient to reject the hypothesis
    # that the residuals are symmetrically distributed."
    # )

    # ==
    # Efficiency estimation by the method of moments
    # ==

    # standard deviation for inefficiency term
    std_u_mm <- (M3 / (sqrt(2 / pi) * (1 - 4 / pi))) ^ (1 / 3)

    # standard deviation for error term
    std_v_mm <- sqrt(M2 - ((pi - 2) / pi) * std_u_mm ^ 2)
    std_v_mm <- ifelse(is.nan(std_v_mm), 0.0001, std_v_mm)

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
    std_u_pl <- ifelse(std_v_pl < 0, 0.0001, std_v_pl)

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

    model_list[[model]][["stochastic"]] <- stochastic_analysis

    # generate technology
    technology[[model]][["xmat"]] <- as.matrix(DMUs[, x_vars])
    technology[[model]][["ymat"]] <- as.matrix(y_hat)

    if (ptto) {
      tech_ymat[nrow(tech_ymat), ] <- if (error_type == "add") 1e-10 else log(1 + 1e-10)
    }

    names(technology[[model]][["xmat"]]) <- names(DMUs)[x_vars]
    names(technology[[model]][["ymat"]]) <- names(DMUs)[y_vars]

  }

  # ============= #
  # S_ACES OBJECT #
  # ============= #

  S_ACES <- aces_object (
    data = DMUs,
    x = x_vars,
    y = y_vars,
    z = z_vars,
    error_type = error_type,
    max_degree = max_degree,
    compl_cost = compl_cost,
    xi_degree = xi_degree,
    metric = metric,
    shape = shape,
    nterms = ncol(aces_forward[["Bmatx"]]),
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = kn_penalty,
    wc = aces_cubic[["w"]],
    wq = aces_quintic[["w"]],
    aces_forward = model_list[["aces_forward"]],
    aces = model_list[["aces"]],
    aces_cubic = model_list[["aces_cubic"]],
    aces_quintic = model_list[["aces_quintic"]],
    technology = technology
  )

  return(S_ACES)

}

#' @title Conditional Efficiency Estimation
#'
#' @description This function computes efficiency scores assuming a zero-truncated normal distribution.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param mean_pred
#' A \code{vector} of mean predicted outputs.
#'
#' @param std_ineff
#' A \code{numeric} indicating the Standard Deviation of the Inefficiency Term.
#'
#' @param std_error
#' A \code{numeric} indicating the Standard Deviation of the Error Term.
#'
#' @param error_type
#' A \code{character} string specifying the error structure when fitting the model.
#'
#' @return
#'
#' A \code{vector} of \code{"numeric"} scores computed through the Conditional Estimation.

cond_eff <- function (
    data,
    y,
    mean_pred,
    std_ineff,
    std_error,
    error_type
) {

  # number of DMUs in the technology
  scores <- rep(NaN, nrow(data))

  # variance
  sigma2 <- std_ineff ^ 2 + std_error ^ 2

  # standard deviation
  sigma <- sqrt(sigma2)

  # variance_star
  sigma2_star <- std_ineff ^ 2 * std_error ^ 2 / sigma2

  # dev_star
  sigma_star <- sqrt(sigma2_star)

  # lambda
  lambda <- std_ineff / std_error

  # frontier estimation
  y_front <- mean_pred + std_ineff * sqrt(2 / pi)

  # composite error term
  comp_err <- data[, y] - y_front

  # mu_star
  mu_star <- - comp_err * std_ineff ^ 2 / sigma2

  # random variable for normal distribution
  argnorm <- comp_err * lambda / sigma

  # conditional mean (2) = (3) in Jondrow et al. (1981)
  cond_mean <- mu_star + sigma_star * (dnorm(argnorm) / (1 - pnorm(argnorm)))

  if (error_type == "mul") {
    scores <- exp( - cond_mean)

  } else {
    scores <- 1 - cond_mean / y_front

  }

  # If sigma_u == 0: all firms are diagnosed as efficient
  if (std_ineff == 0) {
    scores <- rep(1, nrow(data))
  }

  return(scores)

}



