#' @title Fit an Adaptive Constrained Enveloping Splines (ACES) model
#'
#' @description
#'
#' This function estimates a deterministic production frontier that adheres to classical production theory axioms, such as monotonicity and concavity. The estimation is based on adaptations of the Multivariate Adaptive Regression Splines (MARS) technique, developed by \insertCite{friedman1991;textual}{aces}. For comprehensive details on the methodology and implementation, please refer to \insertCite{espana2024;textual}{aces} and \insertCite{espana2025;textual}{aces}.
#'
#' @name aces
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
#' @param max_terms
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
#' An \code{aces} object.
#'
#' @export

aces <- function (
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
    max_terms = 50,
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
    quick_aces = quick_aces,
    error_type = error_type,
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

  if (shape[["ptto"]]) {
    data <- rbind(data, rep(1e-10, ncol(data)))
  }

  ACES <- aces_algorithm (
    data = data,
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
    max_terms = max_terms,
    err_red = err_red,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    kn_penalty = kn_penalty
    )

  # type of object
  class(ACES) <- "aces"

  return(ACES)

}

#' @title Algorithm of Adaptive Constrained Enveloping Splines (ACES).
#'
#' @description
#' This function implements the Adaptive Constrained Enveloping Splines (ACES) algorithm.
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
#' @param max_terms
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
#' @importFrom dplyr desc
#'
#' @return
#'
#' An \code{aces} object.
#'
#' @export

aces_algorithm <- function (
    data,
    x_vars,
    y_vars,
    quick_aces,
    error_type,
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

  # data in [x, y] format with interaction of variables included
  data <- set_data (
    data = data,
    x = x_vars,
    y = y_vars,
    max_degree = max_degree,
    error_type = error_type
    )

  # samples size
  N <- nrow(data)

  # set "x" and "y" indexes in data
  x <- 1:(ncol(data) - length(y_vars))
  y <- (length(x) + 1):ncol(data)

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
  #    LOF: mean error between true data and predicted data (B %*% coefs)
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
  aces_forward <- list (
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

  # list to save technologies created through ACES
  technology <- list()

  # initial error
  err_min <- err

  # iteration counter
  iter <- 0

  while(length(aces_forward[["bf_set"]]) + 2 < max_terms) {

    # add 2 new basis functions to the model:
    B_bf_knt_err <- add_basis_function (
      data = data,
      x = x,
      y = y,
      xi_degree = xi_degree,
      inter_cost = inter_cost,
      model_type = "envelopment",
      dea_scores = dea_scores,
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
      var_imp <- rbind(var_imp, rep(0, nX))

      cat(
        sprintf(
          paste0(
            "Iteration %d completed — error reduced by %.1f%% ",
            "(%.1f -> %.1f).\n"
          ),
          iter,
          100 * rel_reduction,
          err[1] / (1 - rel_reduction),  # old err reconstructed
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

  # generate technology
  technology[["aces_forward"]] <- generate_technology (
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = aces_forward[["Bmatx"]] %*% aces_forward[["coefs"]],
    error_type = error_type,
    table_scores = table_scores,
    ptto = shape[["ptto"]]
  )

  # ====================== #
  #   BACKWARD ALGORITHM   #
  # ====================== #

  aces_submodels <- aces_pruning (
      data = data,
      x = x,
      y = y,
      xi_degree = xi_degree,
      model_type = "envelopment",
      dea_scores = dea_scores,
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

  # generate technology
  technology[["aces"]] <- generate_technology (
    tech_xmat = DMUs[, x_vars],
    tech_ymat1 = DMUs[, y_vars],
    tech_ymat2 = aces[["Bmatx"]] %*% aces[["coefs"]],
    error_type = error_type,
    table_scores = table_scores,
    ptto = shape[["ptto"]]
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
          dea_scores = dea_scores,
          model_type = "envelopment",
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
          dea_scores = dea_scores,
          model_type = "envelopment",
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

    # generate technology
    technology[["aces_cubic"]] <- generate_technology (
      tech_xmat = DMUs[, x_vars],
      tech_ymat1 = DMUs[, y_vars],
      tech_ymat2 = aces_cubic[["Bmatx"]] %*% aces_cubic[["coefs"]],
      error_type = error_type,
      table_scores = table_scores,
      ptto = shape[["ptto"]]
    )

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

    # generate technology
    technology[["aces_quintic"]] <- generate_technology (
      tech_xmat = DMUs[, x_vars],
      tech_ymat1 = DMUs[, y_vars],
      tech_ymat2 = aces_quintic[["Bmatx"]] %*% aces_quintic[["coefs"]],
      error_type = error_type,
      table_scores = table_scores,
      ptto = shape[["ptto"]]
    )

    # =========== #
    # ACES OBJECT #
    # =========== #

    ACES <- aces_object (
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
      max_terms = ncol(aces_forward[["Bmatx"]]),
      err_red = err_red,
      minspan = minspan,
      endspan = endspan,
      kn_grid = kn_grid,
      kn_penalty = kn_penalty,
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

#' @title Create an aces object
#'
#' @description
#'
#' This function saves information about the Adaptive Constrained Enveloping Splines (ACES) model.
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
#' A \code{character} string specifying the error structure when fitting the model.
#'
#' @param max_degree
#' Maximum degree of interaction between variables.
#'
#' @param inter_cost
#' Minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.
#'
#' @param xi_degree
#' A \code{matrix} specifying the degree of each variable.
#'
#' @param metric
#' A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity.
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
#' Design of the grid of knots to perform ACES
#'
#' @param kn_penalty
#' Generalized Cross Validation (GCV) penalty per knot.
#'
#' @param wc
#' Hyperparameter for the side knot distances in the cubic smoothing procedure.
#'
#' @param wq
#' Hyperparameter for the side knot distances in the quintic smoothing procedure.
#'
#' @param aces_forward
#' A \code{list} containing the forward step of the Adaptive Constrained Enveloping Splines model.
#'
#' @param aces
#' A \code{list} containing the Adaptive Constrained Enveloping Splines model.
#'
#' @param aces_cubic
#' A \code{list} containing the Smooth Adaptive Constrained Enveloping Splines through cubic basis functions.
#'
#' @param aces_quintic
#' A \code{list} containing the Smooth Adaptive Constrained Enveloping Splines through quintic basis functions.
#'
#' @param technology
#' A \code{list} with the points that make up the technology set for each model.
#'
#' @return
#'
#' An \code{aces} object.

aces_object <- function (
    data,
    x,
    y,
    quick_aces,
    error_type,
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
    wc,
    wq,
    aces_forward,
    aces,
    aces_cubic,
    aces_quintic,
    technology
    ) {

  object <- list()

  object[["data"]] <- list (
    "df" = data,
    "x" = x,
    "y" = y,
    "xnames" = colnames(data)[x],
    "ynames" = colnames(data)[y],
    "rownames" = rownames(data)
  )

  object[["control"]] <- list (
    "quick_aces" = quick_aces,
    "error_type" = error_type,
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
    "wc" = wc,
    "wq" = wq
  )

  object[["methods"]] <- list (
    "aces_forward" = aces_forward,
    "aces" = aces,
    "aces_cubic" = aces_cubic,
    "aces_quintic" = aces_quintic
  )

  object[["technology"]] <- technology

  return(object)

}

#' @title Arrange Data for Fitting Model
#'
#' @description
#' This function prepares the data for model fitting by generating additional input variables through interactions between variables. It also performs any necessary transformations, such as changing to a logarithmic scale if the error type is multiplicative. It returns a matrix in [x, y] format, where x represents input variables and y represents output variables.
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
#' @param max_degree
#'  Maximum degree of interaction between variables. It can be a \code{list} of input indexes for interactions or a \code{numeric} value determining the maximum degree of interaction.
#'
#' @param error_type
#' A \code{character} string specifying the error structure when fitting the model.
#'
#' @return
#' A \code{matrix} in a [x, y] format with variable interactions and / or transformations included.

set_data <- function (
    data,
    x,
    y,
    max_degree,
    error_type
    ) {

  # 1. change to logarithmic scale if the error_type is multiplicative
  if (error_type == "mul") {
    data[, c(y)] <- log(data[, c(y)])
  }

  # 2. generate interaction effects
  if (is.list(max_degree) || max_degree > 1) {

    if (!is.list(max_degree)) {

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
      name <- paste(name_vars, collapse = "_")

      # create the new variable
      data[, name] <- apply(data[, vars], 1, prod)

    }

  } else {

    new_x <- x

  }

  # 3. data correctly sorted
  data <- data[, c(new_x, y)]

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
#' @param weight
#' Weights to compute error metrics.
#'
#' @return
#' The calculated error metric.

err_metric <- function (
    y_obs,
    y_hat,
    metric,
    weight
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
    error <- sum(weight * devtn) / (N * nY)

  } else if (metric == "mape") {

    # mean absolute percentage error
    devtn <- abs(y_hat - y_obs) / y_obs
    error <- sum(weight * devtn) / (N * nY) * 100

  } else if (metric == "mse") {

    # mean squared error
    devtn <- (y_hat - y_obs) ^ 2
    error <- sum(weight * devtn) / (N * nY)

  } else if (metric == "msle") {

    # mean squared logarithmic error
    devtn <- (log(y_hat + 1) - log(y_obs + 1)) ^ 2
    error <- sum(weight * devtn) / (N * nY)

  } else if (metric == "rmse") {

    # root mean squared error
    devtn <- (y_hat - y_obs) ^ 2
    error <- sqrt(sum(weight * devtn) / (N * nY))

  } else if (metric == "nrmse1") {

    # normalized root mean squared error by the mean
    devtn <- (y_hat - y_obs) ^ 2
    error <- sqrt(sum(weight * devtn) / (N * nY)) / mean(y_obs)

  } else {

    # compute the mean of column-wise maximums and minimums in y
    ymax <- mean(apply(y_obs, 2, max))
    ymin <- mean(apply(y_obs, 2, min))

    # normalized root mean squared error by the range
    devtn <- (y_hat - y_obs) ^ 2
    error <- sqrt(sum(weight * devtn) / (N * nY)) / (ymax - ymin)
  }

  return(error)

}

#' @title Compute Minimum and End Span
#'
#' @description
#' This function computes the minimum span, which is the minimum number of observations between two adjacent knots, and the end span, which is the minimum number of observations before the first knot and after the final knot.
#'
#' @param kn_grid
#' A \code{list} with the available set of knots to perform ACES.
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
#' @param n_input
#' Number of input variables and contextual variables
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
  if (minspan == - 2) { # Zhang approach

    L <- numeric(n_input)

    # fixed log_factor
    log_factor <- log2(- (1 / N) * log(0.95))

    for (var in 1:n_input) {

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

    Le <- numeric(n_input)

    # fixed log_factor
    log_factor <- log2(- (1 / N) * log(0.95))

    for (var in 1:n_input) {

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

    Le <- floor(min(N * 0.1, 3 - log2(0.05 / n_input)))

  } else {

    Le <- min(N * 0.1, endspan)

  }

  return(list(L, Le))

}

#' @title Define the Grid of Knots
#'
#' @description
#' This function defines the grid of knots to perform Adaptive Constrained Enveloping Splines (ACES).
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param n_input_1
#' Number of inputs (excluding interactions) and contextual variables.
#'
#' @param n_input_2
#' Number of inputs (including interactions) and contextual variables.
#'
#' @param kn_grid
#' A \code{list} providing a custom grid of knots for each variable. If not supplied, the function automatically generates a grid of knots for each variable based on the data.
#'
#' @param quick_aces
#' A \code{logical} indicating whether to use the fast version of ACES.
#
#' @param dea_scores
#' A \code{matrix} containing DEA-VRS efficiency scores, calculated using an output-oriented radial model.
#'
#' @return
#' A \code{list} where each element contains the generated or provided vector of knots for the corresponding variable in the model.

set_knots_grid <- function (
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
        kn_grid[[varName]] <- seq (
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

    # names
    names(kn_grid) <- colnames(data)[1:n_input_2]

  }

  return(kn_grid)

}

#' @title Generate a Technology Set
#'
#' @description
#' This function generates a technology set based on input data and an ACES model output.
#' The technology set represents the feasible combinations of inputs and outputs.
#'
#' @param tech_xmat
#' A \code{matrix} representing the input data.
#'
#' @param tech_ymat1
#' A \code{matrix} representing the output data.
#'
#' @param tech_ymat2
#' A \code{matrix} representing the output data generated from an ACES model.
#'
#' @param error_type
#' A \code{character} string specifying the error structure when fitting the model.
#'
#' @param table_scores
#' A \code{matrix} containing radial output scores using all outputs and using just each individual output.
#'
#' @param ptto
#' A \code{logical} indicating if (0, 0) should be included in the technology.
#'
#' @return
#' A \code{matrix} representing the technology set.

generate_technology <- function (
    tech_xmat,
    tech_ymat1,
    tech_ymat2,
    error_type,
    table_scores,
    ptto
    ) {

  # define technology
  tecno <- list()

  # matrix of inputs for the technology
  tecno[["xmat"]] <- as.matrix(tech_xmat)

  # matrix of outputs for the technology
  tech_ymat <- as.data.frame (
    matrix (
      NA,
      nrow = nrow(table_scores),
      ncol = ncol(table_scores) - 1
      )
    )

  colnames(tech_ymat) <- colnames(tech_ymat1)

  for (i in 1:nrow(table_scores)) {

    for (j in 2:ncol(table_scores)) {

      if (abs(table_scores[i, 1] - table_scores[i, j]) <= 0.05) {

        tech_ymat[i, j - 1] <- tech_ymat2[i, j - 1]

      } else {

        tech_ymat[i, j - 1] <- tech_ymat1[i, j - 1]

      }
    }
  }

  # modify the last row of tech_ymat based on the error type
  if (ptto) {
    tech_ymat[nrow(tech_ymat), ] <- if (error_type == "add") 1e-10 else log(1 + 1e-10)
  }

  # assign tech_ymat to tecno[["ymat"]] with the appropriate transformation
  tecno[["ymat"]] <- if (error_type == "add") as.matrix(tech_ymat) else exp(as.matrix(tech_ymat))

  return(tecno)

}

#' @title Generate a New Technology Set
#'
#' @description
#' This function constructs a new production technology by modifying the output matrix based on a refinement
#' procedure that replaces overestimated outputs. For each output, the predicted values from an ACES model are compared against the original observed values. When the predicted value exceeds the observed one by more than a predefined threshold, the observed value is retained instead.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param tech_xmat
#' A \code{matrix} representing the input data.
#'
#' @param tech_ymat
#' A \code{matrix} representing the output data.
#'
#' @param psi
#' A \code{numeric} threshold controlling the refinement of predicted outputs. If the predicted value for a given output exceeds the observed value by more than \code{psi}, the observed value is used instead.
#'
#' @return
#' An \code{aces} object with updated technologies.
#'
#' @export

change_technology <- function (
    object,
    tech_xmat,
    tech_ymat,
    psi
    ) {

  # sample size
  N <- nrow(tech_xmat)

  # number of outputs
  nY <- ncol(tech_ymat)

  # table of scores
  table_scores <- matrix (
    ncol = nY + 1,
    nrow = N,
    dimnames = list(NULL, c("y_all", paste("y", 1:nY, sep = "")))
  ) %>% as.data.frame()

  table_scores[, 1] <- rad_out (
    tech_xmat = as.matrix(tech_xmat),
    tech_ymat = as.matrix(tech_ymat),
    eval_xmat = as.matrix(tech_xmat),
    eval_ymat = as.matrix(tech_ymat),
    convexity = TRUE,
    returns = "variable"
  )[, 1]

  for (out in 1:nY) {

    table_scores[, 1 + out] <- rad_out (
      tech_xmat = as.matrix(tech_xmat),
      tech_ymat = as.matrix(tech_ymat[, out]),
      eval_xmat = as.matrix(tech_xmat),
      eval_ymat = as.matrix(tech_ymat[, out]),
      convexity = TRUE,
      returns = "variable"
    )[, 1]

  }

  for (m in names(object[["technology"]])) {

    # define technology
    tecno <- list()

    # matrix of inputs for the technology
    tecno[["xmat"]] <- as.matrix(tech_xmat)

    # prediction
    Bmat <- object[["methods"]][[m]][["Bmatx"]]
    coef <- object[["methods"]][[m]][["coefs"]]

    tech_ymat2 <- Bmat %*% coef

    colnames(tech_ymat2) <- colnames(tech_ymat)

    # matrix of outputs for the technology
    tech_out <- as.data.frame (
      matrix (
        NA,
        nrow = nrow(table_scores),
        ncol = ncol(table_scores) - 1
      )
    )

    for (i in 1:nrow(table_scores)) {

      for (j in 2:ncol(table_scores)) {

        if (abs(table_scores[i, 1] - table_scores[i, j]) <= psi) {

          tech_out[i, j - 1] <- tech_ymat2[i, j - 1]

        } else {

          tech_out[i, j - 1] <- tech_ymat[i, j - 1]

        }
      }
    }

    ptto <- object[["control"]][["shape"]][["ptto"]]

    # modify the last row of tech_ymat based on the error type
    if (ptto) {
      tech_out[nrow(tech_ymat), ] <- if (error_type == "add") 1e-10 else log(1 + 1e-10)
    }

    error_type <- object[["control"]][["error_type"]]

    # assign tech_ymat to tecno[["ymat"]] with the appropriate transformation
    tecno[["ymat"]] <- if (error_type == "add") as.matrix(tech_out) else exp(as.matrix(tech_out))

    # update technology set
    object[["technology"]][[m]][["ymat"]] <- as.matrix(tecno[["ymat"]])

  }

  return(object)

}
