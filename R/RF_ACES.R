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
#' @param z
#' Column indexes of contextual variables in \code{data}.
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
#'   \item{\code{ptto}}: A \code{logical} indicating if the estimator should satisfy \code{f(0) = 0}.
#' }
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
#' @param nterms
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
#' @details
#'
#' This function generates a production frontier that adheres to classical production theory axioms, such as monotonicity and concavity. The algorithm comprises two main procedures:
#'
#' Forward Selection Algorithm. This initial step constructs a set of linear basis functions to model the production frontier based on the Random Forest foundations.
#'
#' Smoothing Procedures. After the forward step, the algorithm offers two smoothing options to produce a smooth and continuous production frontier:
#'
#' * Cubic Smoothing: Applies cubic functions to achieve a balance between flexibility and smoothness, suitable for capturing moderate curvature in the data.
#'
#' * Quintic Smoothing: Utilizes quintic functions for a higher degree of smoothness and continuity, ideal for modelling more complex relationships.
#'
#' @references
#'
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
#' An \code{aces} object.
#'
#' @export

rf_aces <- function (
    data,
    x,
    y,
    z = NULL,
    quick_aces = TRUE,
    model_type = "env",
    error_type = "add",
    mul_BF = list (
      "max_degree" = 1,
      "compl_cost" = 0.05
    ),
    metric = "mse",
    shape = list (
      "mono" = TRUE,
      "conc" = TRUE,
      "ptto" = FALSE
      ),
    learners = 100,
    nvars = length(x) / 3,
    sample_size = nrow(data),
    nterms = 50,
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
    z = z,
    quick_aces = quick_aces,
    error_type = error_type,
    learners = learners,
    nvars = nvars,
    sample_size = sample_size,
    max_degree = mul_BF[["max_degree"]],
    compl_cost = mul_BF[["compl_cost"]],
    metric = metric,
    nterms = nterms,
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
  oob_vec <- c()
  mov_avg_oob_05 <- c()
  mov_avg_oob_10 <- c()
  mov_avg_oob_25 <- c()

  stopping_condition_counter <- 0

  for (m in 1:length(RF_ACES)) {

    # update Progress Bar
    setTxtProgressBar(pb, m)

    # index of samples in the model
    sample_bag <- sample (
      1:nrow(data),
      size = sample_size,
      replace = TRUE
    )

    # data for the m-model
    data_bag <- data[sample_bag, ]

    # out of bag data for the m-model
    inb_idxs <- 1:nrow(data) %in% sample_bag
    oob_idxs <- which(!(1:nrow(data) %in% sample_bag))
    data_oob <- data[!inb_idxs, ]

    RF_ACES[[m]] <- rf_aces_algorithm (
      data = data_bag,
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
      nvars = nvars,
      err_red = err_red,
      minspan = minspan,
      endspan = endspan,
      kn_grid = kn_grid
      )

      # predictions
      Bmatx <- RF_ACES[[m]][["methods"]][["rf_aces"]][["Bmatx"]]
      coefs <- RF_ACES[[m]][["methods"]][["rf_aces"]][["coefs"]]
      y_hat <- Bmatx %*% coefs

      # select out-of-bag indexes
      y_hat <- y_hat[oob_idxs, ]

      oob_pred[[m]] <- matrix (
        c (
          "y_hat" = y_hat,
          "idx" = oob_idxs
          ),
        ncol = nY + 1
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
      oob_vec <- c(oob_vec, RF_ACES[[m]][["OOB"]])
      mov_avg_oob_05 <- c(mov_avg_oob_05, mean(tail(oob_vec, 05)))
      mov_avg_oob_10 <- c(mov_avg_oob_10, mean(tail(oob_vec, 10)))
      mov_avg_oob_25 <- c(mov_avg_oob_25, mean(tail(oob_vec, 25)))

      if (m > min(nrow(data), 200) && oob_vec[m] >= oob_vec[m - 1]) {

        stopping_condition_counter <- ifelse (
          oob_vec[m] >= mov_avg_oob_05[m] * (1 - RF[["oob_red"]]) &
          oob_vec[m] >= mov_avg_oob_10[m] * (1 - RF[["oob_red"]]) &
          oob_vec[m] >= mov_avg_oob_25[m] * (1 - RF[["oob_red"]]),
          stopping_condition_counter + 1,
          0
        )

        if (stopping_condition_counter == 5) {

          cat("\n")
          print(paste("Out-of-bag error stabilized after", m, "models"))
          break

        }
      }
    }

    # random-forest hyperparameters
    rf_list = list (
      "apply" = RF[["apply"]],
      "sample" = RF[["sample"]],
      "models" = m,
      "nvars" = RF[["nvars"]],
      "oob_red" = RF[["oob_red"]]
    )

    ACES <- RF_ACES[1:m]

    for (j in seq_along(ACES)) {

      ACES[[j]][["control"]] <- append (
        ACES[[j]][["control"]],
        list("RF" = rf_list)
      )

    }

    models <- list()

    for (i in 1:length(ACES)) {

      models[[i]] <- ACES[[i]]

    }

    ACES <- list(models = models)

    # free up memory
    rm(RF_ACES, models)

    xmat_rf_aces <- xmat_rf_aces_cubic <- xmat_rf_aces_quintic <- NULL
    ymat_rf_aces <- ymat_rf_aces_cubic <- ymat_rf_aces_quintic <- NULL

    # save technology for Random Forest
    for (model in ACES[["models"]]) {

      xmat_rf_aces <- rbind (
        xmat_rf_aces,
        model[["technology"]][["rf_aces"]][["xmat"]]
        )

      ymat_rf_aces <- rbind (
        ymat_rf_aces,
        model[["technology"]][["rf_aces"]][["ymat"]]
      )

      xmat_rf_aces_cubic <- rbind (
        xmat_rf_aces_cubic,
        model[["technology"]][["rf_aces_cubic"]][["xmat"]]
      )

      ymat_rf_aces_cubic <- rbind (
        ymat_rf_aces_cubic,
        model[["technology"]][["rf_aces_cubic"]][["ymat"]]
      )

      xmat_rf_aces_quintic <- rbind (
        xmat_rf_aces_quintic,
        model[["technology"]][["rf_aces_quintic"]][["xmat"]]
      )

      ymat_rf_aces_quintic <- rbind (
        ymat_rf_aces_quintic,
        model[["technology"]][["rf_aces_quintic"]][["ymat"]]
      )

    }

    input_cols <- names(data)[1:length(x)]
    output_cols <- names(data)[(length(x) + 1):(length(x) + length(y))]

    colnames(xmat_rf_aces) <- colnames(xmat_rf_aces_cubic) <- colnames(xmat_rf_aces_quintic) <- input_cols
    colnames(ymat_rf_aces) <- colnames(ymat_rf_aces_cubic) <- colnames(ymat_rf_aces_quintic) <- output_cols

    rf_aces_technology <- data.frame(xmat_rf_aces, ymat_rf_aces) %>%
      group_by(across(all_of(input_cols))) %>%
      summarise(across(all_of(output_cols), mean, .names = "{col}"), .groups = 'drop')

    rf_aces_cubic_technology <- data.frame(xmat_rf_aces_cubic, ymat_rf_aces_cubic) %>%
      group_by(across(all_of(input_cols))) %>%
      summarise(across(all_of(output_cols), mean, .names = "{col}"), .groups = 'drop')

    rf_aces_quintic_technology <- data.frame(xmat_rf_aces_quintic, ymat_rf_aces_quintic) %>%
      group_by(across(all_of(input_cols))) %>%
      summarise(across(all_of(output_cols), mean, .names = "{col}"), .groups = 'drop')

    ACES[["technology"]] <- list (
      "rf_aces" = rf_aces_technology,
      "rf_aces_cubic" = rf_aces_cubic_technology,
      "rf_aces_quintic" = rf_aces_quintic_technology
    )



  # type of object
  class(ACES) <- ifelse(RF[["apply"]], "rf_aces", "aces")

  return(ACES)

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
#' A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin (only for piece-wise linear version).
#'
#' @param nterms
#' Maximum number of terms created before pruning.
#'
#' @param nvars
#' An \code{integer} indicating the number of variables randomly chosen at each split.
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
    z_vars,
    quick_aces,
    error_type,
    max_degree,
    compl_cost,
    metric,
    shape,
    nterms,
    nvars,
    err_red,
    minspan,
    endspan,
    kn_grid
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

  # table of scores
  table_scores <- matrix (
    ncol = nY + 1,
    nrow = nrow(data),
    dimnames = list(NULL, c("y_all", paste("y", 1:nY, sep = "")))
  ) %>% as.data.frame()

  dea_returns <- ifelse(shape[["conc"]], "variable", "constant")

  table_scores[, 1] <- rad_out (
    tech_xmat = as.matrix(DMUs[, x_vars]),
    tech_ymat = as.matrix(DMUs[, y_vars]),
    eval_xmat = as.matrix(DMUs[, x_vars]),
    eval_ymat = as.matrix(DMUs[, y_vars]),
    convexity = TRUE,
    returns = dea_returns
  )[, 1]

  for (out in 1:nY) {

    table_scores[, 1 + out] <- rad_out (
      tech_xmat = as.matrix(DMUs[, x_vars]),
      tech_ymat = as.matrix(DMUs[, y_vars[out]]),
      eval_xmat = as.matrix(DMUs[, x_vars]),
      eval_ymat = as.matrix(DMUs[, y_vars[out]]),
      convexity = TRUE,
      returns = dea_returns
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
    dea_scores = table_scores[, 1]
  )

  # minimum span (minspan) and end span (endspan)
  L_Le <- compute_span (
    kn_grid = kn_grid,
    minspan = minspan,
    endspan = endspan,
    n_input = nX + nZ
  )

  # list to save technologies created through ACES
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

  browser()

  while(length(rf_aces[["bf_set"]]) + 2 < nterms) {

    # negative GRSq
    last_bf <- rf_aces[["bf_set"]][[length(rf_aces[["bf_set"]])]]
    if (last_bf[["GRSq"]] < 0) break

    while (last_bf[["id"]] == 1) {

      # random input for bagging
      random_x <- sample(x, nvars)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function (
        data = data,
        x = random_x,
        y = y,
        xi_degree = xi_degree,
        dea_scores = dea_scores,
        model_type = model_type,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        L = L,
        Le = Le,
        err_min = err,
        hd_cost = hd_cost
      )

      if (is.list(B_bf_knt_err)) break

    }

    if (last_bf[["id"]] != 1) {

      # random input for bagging
      random_x <- sample(x, nvars)

      # add 2 new basis functions to the model:
      B_bf_knt_err <- add_basis_function (
        data = data,
        x = random_x,
        y = y,
        xi_degree = xi_degree,
        dea_scores = dea_scores,
        model_type = model_type,
        metric = metric,
        forward_model = rf_aces,
        Bp_list = Bp_list,
        shape = shape,
        kn_list = kn_list,
        kn_grid = kn_grid,
        L = L,
        Le = Le,
        err_min = err,
        hd_cost = hd_cost
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
    tech_xmat = dmus[, inps],
    tech_ymat1 = dmus[, outs],
    tech_ymat2 = rf_aces[["Bmatx"]] %*% rf_aces[["coefs"]],
    error_type = error_type,
    table_scores = table_scores
  )

  # ======================= #
  #   SMOOTHING PROCEDURE   #
  # ======================= #

  # select a model to be smoothed
  rf_aces_smoothed <- rf_aces

  # data.frame of knots
  kn_smoothed <- rf_aces_smoothed[["knots"]]
  kn_smoothed$side <- "R"

  kn_smoothed_left <- kn_smoothed
  kn_smoothed_left$side <- "L"

  kn_smoothed <- rbind (
    kn_smoothed,
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
    model_type = model_type,
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    d = 0,
    wc = wc
  )

  # generate technology
  technology[["rf_aces_cubic"]] <- generate_technology (
    tech_xmat = dmus[, inps],
    tech_ymat1 = dmus[, outs],
    tech_ymat2 = rf_aces_cubic[["Bmatx"]] %*% rf_aces_cubic[["coefs"]],
    error_type = error_type,
    table_scores = table_scores
  )

  rf_aces_quintic <- quintic_aces (
    data = data,
    x = x,
    y = y,
    dea_scores = dea_scores,
    model_type = model_type,
    metric = metric,
    shape = shape,
    kn_grid = kn_smoothed,
    kn_side_loc = kn_side_loc,
    d = 0,
    wq = wq
  )

  # generate technology
  technology[["rf_aces_quintic"]] <- generate_technology (
    tech_xmat = dmus[, inps],
    tech_ymat1 = dmus[, outs],
    tech_ymat2 = rf_aces_quintic[["Bmatx"]] %*% rf_aces_quintic[["coefs"]],
    error_type = error_type,
    table_scores = table_scores
  )

  # =========== #
  # ACES OBJECT #
  # =========== #

  RF_ACES <- aces_object (
    data = DMUs,
    x = inps,
    y = outs,
    model_type = model_type,
    error_type = error_type,
    degree = degree,
    metric = metric,
    shape = shape,
    nterms = ncol(rf_aces[["Bmatx"]]),
    err_red = err_red,
    hd_cost = hd_cost,
    minspan = minspan,
    endspan = endspan,
    kn_grid = kn_grid,
    d = NULL,
    wc = wc,
    wq = wq,
    aces_forward = rf_aces,
    aces = NULL,
    aces_cubic = rf_aces_cubic,
    aces_quintic = rf_aces_quintic,
    technology = technology
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

  # calculate the OOB mean squared error
  oob_idxs <- as.numeric(rownames(model_pred))
  oob_mse <- sum((data[oob_idxs, y] - model_pred[, 1:nY, drop = F]) ^ 2) / (length(oob_idxs) * nY)

  return(oob_mse)

}

#' @title Compute Variable Importance in Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).
#'
#' @description
#'
#' This function computes a metric of variable importance for a given Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) with a certain number of base learners.
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
#' @return
#'
#' This function returns a metric of variable importance for each input variable

compute_variable_importance <- function (
    data,
    x,
    y,
    object,
    repeats = 1
) {

  # ACES meta-data
  control_features <- object[[1]][[1]][["control"]]

  model_type <- control_features[["model_type"]]
  error_type <- control_features[["error_type"]]
  y_type <- object[[1]][[1]][["data"]][["y_type"]]
  degree <- control_features[["degree"]]
  metric <- control_features[["metric"]]
  mon <- control_features[["monotonicity"]]
  con <- control_features[["concavity"]]
  ori <- as.logical(control_features[["origin"]])
  err_red <- control_features[["err_red"]]
  hd_cost <- control_features[["hd_cost"]]
  minspan <- control_features[["minspan"]]
  endspan <- control_features[["endspan"]]
  wc <- control_features[["wc"]]
  wq <- control_features[["wq"]]

  # RF-ACES meta-data
  RF_features <- object[[1]][[1]][["control"]][["RF"]]

  RF_sample <- RF_features[["sample"]]
  RF_models <- RF_features[["models"]]
  RF_nvars <- RF_features[["nvars"]]
  RF_oob_red <- RF_features[["oob_red"]]

  # number of estimated technologies
  nrankings <- length(object[[1]])

  # ranking
  rankings <- vector("list", nrankings)

  # sample size
  N <- nrow(data)

  for (k in 1:nrankings) {

    # select a model
    data_info <- object[[1]][[k]][["data"]]

    # look for netputs if y_type = "individual"
    z <- data_info[["z"]]

    # output index
    y <- data_info[["y"]]

    # data in [x, z, y] format with interaction of variables included
    data_model <- set_interactions (
      data = data,
      x = x,
      y = NULL,
      z = z,
      degree = degree
    )

    # add output data to input data
    data_model <- cbind(data_model, data[, y, drop = FALSE])

    # reorder index 'x' and 'y' in data
    x <- 1:(ncol(data_model) - length(y))
    y <- (length(x) + 1):ncol(data_model)

    # number of inputs
    nX <- length(x)

    # initialize matrix of variable importance
    mat_varimp <- matrix (
      0,
      nrow = nX
    )
    rownames(mat_varimp) <- colnames(data_model)[x]
    colnames(mat_varimp) <- "importance"

    # out-of-bag error of the complete model
    oob <- object[[length(object)]][[1]][["OOB"]]

    # out-of-bag excluding the j-th input variable
    oob_j <- c()

    # compute oob shuffling each variable "j"
    for (j in 1:nX) {
      for (n in 1:repeats) {
        data_shuffled <- data_model
        xvar_shuffled <- data_model[sample(1:N), x[j]]
        data_shuffled[, x[j]] <- xvar_shuffled

        model_j <- aces (
          data = data_shuffled,
          x = x,
          y = y,
          y_type = y_type,
          model_type = model_type,
          error_type = error_type,
          RF = list (
            "apply" = TRUE,
            "sample" = RF_sample,
            "models" = RF_models,
            "nvars" = RF_nvars,
            "oob_red" = RF_oob_red
          ),
          mul_BF = list (
            "degree" = degree,
            "hd_cost" = 0
          ),
          metric = "mse",
          shape = list (
            "mon" = mon,
            "con" = con,
            "ori" = ori
          ),
          nterms = nrow(data),
          err_red = err_red,
          kn_grid = - 1,
          minspan = minspan,
          endspan = endspan,
          kn_penalty = NULL,
          smoothing = list (
            "wc" = wc,
            "wq" = wq
          )
        )

        # out-of-bag error for the model with the jth variable shuffled
        oob_j <- c(oob_j, model_j[[length(model_j)]][[1]][["OOB"]])
      }

      # add oob_j to the matrix of variable importance
      mat_varimp[j, "importance"] <- round(100 * ((mean(oob_j) - oob) / mean(oob_j)), 2)
    }

    rankings[[k]] <- mat_varimp

    return(rankings)
  }

}

