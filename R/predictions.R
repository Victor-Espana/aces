#' @title Compute Efficient Targets
#'
#' @description
#' Projects each decision-making unit (DMU) onto the production technology
#' constructed by ACES or RF-ACES. The stored input-output reference points are
#' enveloped under the selected returns-to-scale assumption, and the efficiency
#' measure determines which inputs, outputs, or both are adjusted. If
#' \code{method = NULL}, the standard model for the supplied object is used.
#'
#' @param eval_data
#' A \code{data.frame} or \code{matrix} containing the DMUs to evaluate.
#'
#' @param x
#' Column indexes of input variables in \code{eval_data}.
#'
#' @param y
#' Column indexes of output variables in \code{eval_data}.
#'
#' @param relevant
#' If \code{TRUE}, compute the efficiency result with only the original inputs
#' represented in at least one selected basis function. Inputs used through an
#' interaction are also retained.
#'
#' @param object
#' An \code{aces} or \code{rf_aces} object.
#'
#' @param method
#' Fitted model used to define the technology. When \code{NULL}, use \code{"aces"}
#' for an \code{aces} object and \code{"rf_aces"} for an \code{rf_aces} object:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_cubic"}}: Cubic Smooth Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_quintic"}}: Quintic Smooth Adaptive Constrained Enveloping Splines.
#' \item{\code{"rf_aces"}}: Random Forest Adaptive Constrained Enveloping Splines.
#' \item{\code{"rf_aces_cubic"}}: Cubic-smoothed RF-ACES.
#' \item{\code{"rf_aces_quintic"}}: Quintic-smoothed RF-ACES.
#' }
#'
#' @param measure
#' Efficiency measure used for the projection:
#' \itemize{
#' \item{\code{"rad_out"}: Proportional expansion of all outputs while inputs
#' remain fixed.}
#' \item{\code{"rad_inp"}: Proportional contraction of all inputs while outputs
#' remain fixed.}
#' \item{\code{"ddf"}: Simultaneous input contraction and output expansion along
#' the supplied direction.}
#' \item{\code{"rsl_out"}: Separate proportional adjustments for each output.}
#' \item{\code{"rsl_inp"}: Separate proportional adjustments for each input.}
#' \item{\code{"wam_mip"}: Slack-based Measure of Inefficiency Proportions.}
#' \item{\code{"wam_nor"}: Slack-based normalized weighted additive measure.}
#' \item{\code{"wam_ram"}: Slack-based range-adjusted measure.}
#' \item{\code{"wam_bam"}: Slack-based bounded adjusted measure.}
#' \item{\code{"rf_aces_rad_out"}: Targets based directly on RF-ACES
#' predictions. Available only for \code{rf_aces} objects.}
#' }
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param direction
#' Direction vectors used when \code{measure = "ddf"}. Supply a matrix or data
#' frame with one row per evaluated DMU. Input directions must come first,
#' followed by output directions. The direction determines the relative input
#' contractions and output expansions represented by the projection.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{espana2024rf}{aces} \cr \cr
#' \insertRef{banker1984}{aces} \cr \cr
#' \insertRef{chambers1998}{aces} \cr \cr
#' \insertRef{fare1978}{aces} \cr \cr
#' \insertRef{cooper1999}{aces} \cr \cr
#' \insertRef{lovell1995}{aces} \cr \cr
#' \insertRef{cooper2011}{aces}
#'
#' @importFrom dplyr summarise %>% mutate_if
#' @importFrom stats median quantile sd
#'
#' @return
#' A data frame with one row per evaluated DMU. Target input columns are followed
#' by target output columns and use the suffix \code{_hat}. Values are returned on
#' the original data scale.
#'
#' @details
#' The \code{relevant} argument applies a structural rule, not an importance
#' threshold. For an ACES object, the retained inputs come from the selected
#' basis functions of \code{method}. For an RF-ACES object, the function uses
#' the union of the inputs selected by any learner. The reduced input set is
#' used only in the envelopment problem; the model is not refitted. In
#' particular, \code{aces_varimp()} scores do not determine this set. With
#' \code{measure = "rf_aces_rad_out"}, the argument changes only the input
#' columns returned with the targets; it does not change the predicted outputs.
#'
#' @export
#'
get_targets <- function (
    eval_data,
    x,
    y,
    relevant = FALSE,
    object,
    method = NULL,
    measure = "rad_out",
    returns = "variable",
    direction = NULL
    ) {

  # default method based on object class
  if (is.null(method)) {
    method <- if (inherits(object, "rf_aces")) "rf_aces" else "aces"
  }

  # handle errors:
  display_errors_scores (
    data = eval_data,
    x = x,
    y = y,
    object = object,
    method = method,
    measure = measure,
    returns = returns,
    direction = direction
  )

  # number of inputs
  nX <- length(x)

  # number of outputs
  nY <- length(y)

  # =================== #
  # Data for technology #
  # =================== #

  # matrix of inputs
  tech_xmat <- as.matrix(object[["technology"]][[method]][["xmat"]])

  # auxiliar variable for using only relevant variables
  rel_x <- NULL

  if (relevant) {
    if (inherits(object, "rf_aces")) {
      # aggregate knots across all trees
      all_xi <- unique(unlist(lapply(object[["forest"]], function(tree) {
        tree[["methods"]][[method]][["knots"]]$xi
      })))
    } else {
      all_xi <- unique(object[["methods"]][[method]][["knots"]]$xi)
    }

    # variable degree
    xi_degree <- object[["control"]][["xi_degree"]]
    colnames(xi_degree) <- names(object[["control"]][["kn_grid"]])

    # check participating variables
    participating_vars <- intersect(xi_degree[1, ], all_xi)

    # names of participant variables
    participating_vars_names <- colnames(xi_degree)[participating_vars]

    # extract the original variables from interaction names such as "x1:x2"
    split_vars <- unique(unlist(strsplit(
      participating_vars_names,
      ":",
      fixed = TRUE
    )))

    # get the column indices for the unique variables
    rel_x <- sort(match(split_vars, colnames(xi_degree)))

    # update technology
    tech_xmat <- as.matrix(tech_xmat[, rel_x, drop = FALSE])

  } else {

    rel_x <- x

  }

  # matrix of outputs
  tech_ymat <- as.matrix(object[["technology"]][[method]][["ymat"]])

  # ======================= #
  # Data for evaluated DMUs #
  # ======================= #

  # matrix of inputs
  if (relevant && !is.null(rel_x)) {
    eval_xmat <- as.matrix(eval_data[, x[rel_x]])
  } else {
    eval_xmat <- as.matrix(eval_data[, x])
  }

  # matrix of outputs
  eval_ymat <- as.matrix(eval_data[, y])

  # scaling setup
  scaling <- object[["control"]][["scale"]]

  if (!is.null(scaling) && scaling$is_scaled) {

    if (relevant) {
      sx <- scaling$mean_x[rel_x]
    } else {
      sx <- scaling$mean_x
    }

    sy <- scaling$mean_y

    # apply scaling
    eval_xmat <- sweep(eval_xmat, 2, sx, "/")
    eval_ymat <- sweep(eval_ymat, 2, sy, "/")

    if (!is.null(direction) && (is.matrix(direction) || is.data.frame(direction))) {

      dir_x <- as.matrix(direction[, 1:nX])
      dir_y <- as.matrix(direction[, (nX + 1):(nX + nY)])

      if (relevant) {
        if (ncol(dir_x) == length(scaling$mean_x)) {
          dir_x <- dir_x[, rel_x, drop = FALSE]
        }
      }

      # apply scaling
      dir_x_scaled <- sweep(dir_x, 2, sx, "/")
      dir_y_scaled <- sweep(dir_y, 2, sy, "/")

      # rebuild direction vectors
      direction <- cbind(dir_x_scaled, dir_y_scaled)

    }

  }

  # ========== #
  # Get Scores #
  # ========== #

  if (measure == "rf_aces_rad_out") {
    # RF-ACES specific: targets are the forest predictions (already in original scale)
    y_hat_point <- rf_aces_predict(
      object = object,
      eval_data = eval_data,
      x = x,
      method = method
    )

    # inputs remain unchanged; outputs are the predicted frontier values
    if (!is.null(scaling) && scaling$is_scaled) {
      x_hat <- sweep(eval_xmat, 2, sx, "*")
    } else {
      x_hat <- eval_xmat
    }
    y_hat <- as.matrix(y_hat_point)

    colnames(x_hat) <- paste(names(eval_data)[rel_x], "_hat", sep = "")
    colnames(y_hat) <- paste(names(eval_data)[y], "_hat", sep = "")

    targets <- data.frame(x_hat, y_hat)
    return(targets)
  }

  if (measure == "rad_out") {

    scores <- rad_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1]

  } else if (measure == "rad_inp") {

    scores <- rad_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1]

  } else if (measure == "ddf") {

    scores <- ddf (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      direction = direction,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1]

  } else if (measure == "rsl_out") {

    scores <- rsl_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1:nY]

  } else if (measure == "rsl_inp") {

    scores <- rsl_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1:nX]

  } else if (grepl("wam", measure)) {

    scores <- wam (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      weights = measure,
      convexity = TRUE,
      returns = returns,
      type = "variables"
    )[, 1:(nX + nY)]

  }

  # =========== #
  # Get Targets #
  # =========== #

  if (measure == "rad_out") {
    x_hat <- eval_xmat
    y_hat <- eval_ymat * scores

  } else if (measure == "rad_inp") {
    x_hat <- eval_xmat * scores
    y_hat <- eval_ymat

  } else if (measure == "ddf") {
    x_hat <- eval_xmat - scores * direction[, 1:nX]
    y_hat <- eval_ymat + scores * direction[, (nX + 1):(nX + nY)]

  } else if (measure == "rsl_out") {
    x_hat <- eval_xmat
    y_hat <- eval_ymat * scores

  } else if (measure == "rsl_inp") {
    x_hat <- eval_xmat * scores
    y_hat <- eval_ymat

  } else if (grepl("wam", measure)) {

    # extract scaled slacks directly from the linear programming solution
    slacks_x <- scores[, 1:nX, drop = FALSE]
    slacks_y <- scores[, (nX + 1):(nX + nY), drop = FALSE]

    x_hat <- eval_xmat - slacks_x
    y_hat <- eval_ymat + slacks_y

  }

  if (!is.null(scaling) && scaling$is_scaled) {

    # restore the targets to their original scale
    x_hat <- sweep(x_hat, 2, sx, "*")
    y_hat <- sweep(y_hat, 2, sy, "*")

  }

  # save names
  colnames(x_hat) <- paste(names(eval_data)[rel_x], "_hat", sep = "")
  colnames(y_hat) <- paste(names(eval_data)[y], "_hat", sep = "")

  # return predictions
  targets <- data.frame (
    x_hat,
    y_hat
  )

  return(targets)

}

#' @title Predict with an RF-ACES Model
#'
#' @description
#'
#' Predicts attainable outputs by averaging the learners in an \code{rf_aces}
#' object. The averaged predictions are then projected onto the corresponding
#' aggregate technology to preserve the production assumptions used by RF-ACES.
#'
#' @param object
#' A \code{rf_aces} object.
#'
#' @param eval_data
#' A \code{data.frame} or \code{matrix} containing the new observations.
#'
#' @param x
#' Column indexes of input variables in \code{eval_data}.
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
#'
#' A data frame with one row per observation and one column per output. Predictions
#' are returned on the original output scale.

rf_aces_predict <- function (
    object,
    eval_data,
    x,
    method = "rf_aces"
    ) {

  # number of outputs
  nY <- length(object[["data"]][["y"]])

  # =========================== #
  # 1. SCALING INPUTS           #
  # =========================== #

  # scaling parameters
  scaling <- object[["control"]][["scale"]]

  # data to work with
  data_work <- eval_data

  if (!is.null(scaling) && scaling$is_scaled) {
    data_work[, x] <- sweep(as.matrix(eval_data[, x]), 2, scaling$mean_x, "/")
  }

  # =========================== #
  # 2. PREDICTION               #
  # =========================== #

  # data in [x, y] format with interaction of variables included
  data <- set_data (
    data = data_work,
    x = x,
    y = NULL,
    max_degree = object[["control"]][["max_degree"]]
  )

  # technology
  tecno <- object[["technology"]][[method]]

  # number of models in Random Forest
  RF_models <- length(object[["forest"]])

  # list of predictions for each model
  y_hat_RF <- vector("list", RF_models)

  for (t in 1:RF_models) {

    # output predictions
    y_hat <- as.data.frame(matrix(NA, nrow = nrow(eval_data), ncol = nY))

    # select an element from the forest
    model <- object[["forest"]][[t]]

    # method
    aces_model <- model[["methods"]][[method]]
    knots <- aces_model[["knots"]]

    # matrix of basis function
    B <- set_Bmat (
      newdata = data,
      model = aces_model,
      knots = knots,
      method = method
    )

    # prediction
    for (out in 1:nY) {
      y_hat[, out] <- pmax(0, B %*% aces_model[["coefs"]][, out, drop = F])
    }

    y_hat_RF[[t]] <- as.data.frame(y_hat)

  }

  # point estimation
  y_hat_aux <- as.data.frame(matrix(NA, nrow = nrow(eval_data), ncol = nY))

  # mean prediction
  for (var in 1:nY) {
    rf_estimation <- lapply(y_hat_RF, function(df) df[, var])
    matrix_var <- do.call(cbind, rf_estimation)
    y_hat_aux[, var] <- rowMeans(matrix_var, na.rm = TRUE)
  }

  # compute DEA scores to project onto the technology
  scores <- rad_out (
    tech_xmat = tecno[["xmat"]],
    tech_ymat = tecno[["ymat"]],
    eval_xmat = as.matrix(data_work[, x, drop = FALSE]),
    eval_ymat = as.matrix(y_hat_aux),
    convexity = TRUE,
    returns = "variable"
  )

  # remove unfeasibilities
  scores[scores < -10000 | scores > 10000] <- NA

  # predictions
  y_hat <- as.data.frame(y_hat_aux * scores)

  # =========================== #
  # 3. RESTORE ORIGINAL SCALE   #
  # =========================== #

  if (!is.null(scaling) && scaling$is_scaled) {
    y_hat <- sweep(y_hat, 2, scaling$mean_y, "*")
  }

  names(y_hat) <- paste(object[["data"]][["ynames"]], "_pred", sep = "")

  return(y_hat)

}

#' @title Predict Method for RF-ACES Models
#'
#' @description
#' S3 prediction method for \code{rf_aces} objects. It delegates to
#' \code{rf_aces_predict()} and returns the ensemble output predictions.
#'
#' @param object A \code{rf_aces} object.
#' @param newdata A data frame or matrix containing the new observations.
#' @param x Column indexes of input variables in \code{newdata}.
#' @param method RF-ACES model to use.
#' @param ... Additional arguments. They are currently ignored.
#'
#' @return A data frame of predicted outputs.
#'
#' @export

predict.rf_aces <- function (object, newdata, x, method = "rf_aces", ...) {
  rf_aces_predict(
    object = object,
    eval_data = newdata,
    x = x,
    method = method
  )
}

#' @title Compute RF-ACES Ensemble Intervals
#'
#' @description
#' Computes empirical intervals from the variation among the learners of an
#' \code{rf_aces} forest. Output intervals use the prediction from each learner
#' projected onto that learner's technology. Efficiency-score intervals
#' recompute the requested efficiency measure against each learner-specific
#' technology.
#'
#' @param object
#' A fitted \code{rf_aces} object.
#'
#' @param newdata
#' A \code{data.frame} or \code{matrix} containing the observations to predict
#' or evaluate.
#'
#' @param x
#' Column indexes of input variables in \code{newdata}.
#'
#' @param y
#' Column indexes of output variables in \code{newdata}. Required when
#' \code{type} is \code{"score"} or \code{"both"}.
#'
#' @param level
#' Interval level, as a number strictly between zero and one. For example,
#' \code{0.95} returns the 2.5 and 97.5 percent empirical quantiles.
#'
#' @param type
#' Result to compute: \code{"output"}, \code{"score"}, or \code{"both"}.
#'
#' @param method
#' RF-ACES model to use: \code{"rf_aces"}, \code{"rf_aces_cubic"}, or
#' \code{"rf_aces_quintic"}.
#'
#' @param measure
#' Efficiency measure passed to \code{get_scores()}. The default,
#' \code{"rf_aces_rad_out"}, compares each learner's predicted outputs with the
#' observed outputs. The other measures supported by \code{get_scores()} use
#' each learner-specific technology.
#'
#' @param returns
#' Returns-to-scale assumption passed to \code{get_scores()}.
#'
#' @param relevant
#' If \code{TRUE}, apply the input-selection rule described in
#' \code{get_scores()} separately to every learner.
#'
#' @param direction
#' Direction matrix used when \code{measure = "ddf"}.
#'
#' @param calibration
#' Output-interval calibration. Use \code{"none"} for learner quantiles or
#' \code{"oob"} for an approximate conformal prediction interval calibrated
#' with absolute out-of-bag residuals from the training outputs. OOB calibration
#' applies only to the \code{output} component; score intervals remain ensemble
#' intervals because the true efficiency score is not observed.
#'
#' @param min_oob
#' Minimum number of out-of-bag learner predictions required for a training DMU
#' to contribute a calibration residual. Used only when
#' \code{calibration = "oob"}.
#'
#' @details
#' These are ensemble-dispersion intervals: their limits are empirical
#' quantiles across the fitted learners. They describe sensitivity to the
#' bootstrap samples and randomized input selection used by RF-ACES. They are
#' not calibrated confidence or prediction intervals and do not guarantee a
#' stated frequentist coverage probability.
#'
#' The point estimate is computed with the complete forest using
#' \code{predict()} or \code{get_scores()}. Because RF-ACES aggregation and the
#' technology projection are nonlinear, a point estimate need not equal the
#' midpoint of its learner interval and can occasionally fall outside it.
#'
#' With \code{calibration = "oob"}, output limits are replaced by marginal
#' prediction intervals based on absolute OOB residuals. Their target is a new
#' observed output conditional on its inputs, not the unobserved true production
#' frontier. The calibration is approximate because the same forest supplies
#' both the OOB residuals and the final predictor. With multiple outputs,
#' coverage is marginal for each output rather than simultaneous for the full
#' output vector.
#'
#' @return
#' An object of class \code{rf_aces_intervals}. The \code{output} component is a
#' data frame containing each forest prediction followed by its lower and upper
#' limits. The \code{score} component has the analogous efficiency-score
#' columns. A component not requested in \code{type} is \code{NULL}. The result
#' also records \code{level}, \code{method}, \code{measure}, and the number of
#' learners used. The \code{interval_type} and \code{calibration} components
#' record the interpretation and diagnostics of the returned limits.
#'
#' @export

rf_aces_intervals <- function(
  object,
  newdata,
  x,
  y = NULL,
  level = 0.95,
  type = c("both", "output", "score"),
  method = "rf_aces",
  measure = "rf_aces_rad_out",
  returns = "variable",
  relevant = FALSE,
  direction = NULL,
  calibration = c("none", "oob"),
  min_oob = 5L
) {
  if (!inherits(object, "rf_aces")) {
    stop("'object' must be an rf_aces object.")
  }

  type <- match.arg(type)
  calibration <- match.arg(calibration)

  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be one finite number strictly between 0 and 1.")
  }

  allowed_methods <- c("rf_aces", "rf_aces_cubic", "rf_aces_quintic")
  if (!method %in% allowed_methods) {
    stop(
      paste0(
        "Invalid method '", method, "'. Please choose one of: ",
        paste(allowed_methods, collapse = ", "), "."
      )
    )
  }

  compute_output <- type %in% c("both", "output")
  compute_score <- type %in% c("both", "score")

  if (compute_score && is.null(y)) {
    stop("'y' is required when computing score intervals.")
  }

  if (calibration == "oob" && !compute_output) {
    stop("OOB calibration requires type = 'output' or type = 'both'.")
  }

  if (!is.numeric(min_oob) || length(min_oob) != 1L ||
      !is.finite(min_oob) || min_oob < 1 || min_oob %% 1 != 0) {
    stop("'min_oob' must be one positive integer.")
  }
  min_oob <- as.integer(min_oob)

  n_learners <- length(object[["forest"]])
  if (n_learners < 2L) {
    stop("At least two fitted learners are required to compute an interval.")
  }

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)

  calibration_info <- list(
    method = calibration,
    target = if (calibration == "oob") "observed_output" else "learner_dispersion",
    min_oob = if (calibration == "oob") min_oob else NULL,
    n = NULL,
    adjustment = NULL
  )

  finite_quantile <- function(z, prob) {
    z <- z[is.finite(z)]
    if (length(z) == 0L) return(NA_real_)
    unname(stats::quantile(z, probs = prob, names = FALSE, na.rm = TRUE))
  }

  learner_object <- function(index) {
    learner <- object[["forest"]][[index]]
    ans <- object
    ans[["forest"]] <- list(learner)
    ans[["technology"]] <- learner[["technology"]]
    ans[["control"]][["learners"]] <- 1L
    ans
  }

  output_result <- NULL
  if (compute_output) {
    output_point <- rf_aces_predict(
      object = object,
      eval_data = newdata,
      x = x,
      method = method
    )

    output_members <- array(
      NA_real_,
      dim = c(nrow(newdata), ncol(output_point), n_learners)
    )

    for (tree in seq_len(n_learners)) {
      output_members[, , tree] <- as.matrix(
        rf_aces_predict(
          object = learner_object(tree),
          eval_data = newdata,
          x = x,
          method = method
        )
      )
    }

    output_lower <- apply(
      output_members,
      c(1, 2),
      finite_quantile,
      prob = probs[1]
    )
    output_upper <- apply(
      output_members,
      c(1, 2),
      finite_quantile,
      prob = probs[2]
    )

    if (calibration == "oob") {
      has_bag_indices <- vapply(
        object[["forest"]],
        function(learner) !is.null(learner[["sample_bag"]]),
        logical(1)
      )
      if (!all(has_bag_indices)) {
        stop(
          "OOB calibration requires 'sample_bag' indexes in every learner. ",
          "Refit the rf_aces object with the current package version."
        )
      }

      training_data <- object[["data"]][["df"]]
      training_x <- object[["data"]][["x"]]
      training_y <- object[["data"]][["y"]]
      n_training <- nrow(training_data)

      oob_members <- array(
        NA_real_,
        dim = c(n_training, ncol(output_point), n_learners)
      )

      for (tree in seq_len(n_learners)) {
        bag_indices <- unique(object[["forest"]][[tree]][["sample_bag"]])
        oob_indices <- setdiff(seq_len(n_training), bag_indices)
        if (length(oob_indices) == 0L) next

        learner_predictions <- as.matrix(
          rf_aces_predict(
            object = learner_object(tree),
            eval_data = training_data,
            x = training_x,
            method = method
          )
        )
        oob_members[oob_indices, , tree] <- learner_predictions[oob_indices, , drop = FALSE]
      }

      oob_counts <- apply(is.finite(oob_members), c(1, 2), sum)
      oob_point <- apply(oob_members, c(1, 2), mean, na.rm = TRUE)
      observed_outputs <- as.matrix(
        training_data[, training_y, drop = FALSE]
      )

      adjustment <- numeric(ncol(output_point))
      calibration_n <- integer(ncol(output_point))

      for (out in seq_len(ncol(output_point))) {
        usable <- oob_counts[, out] >= min_oob &
          is.finite(oob_point[, out]) &
          is.finite(observed_outputs[, out])
        residuals <- abs(observed_outputs[usable, out] - oob_point[usable, out])

        if (length(residuals) == 0L) {
          stop(
            "No training observations have at least ", min_oob,
            " OOB predictions for output '", names(output_point)[out], "'."
          )
        }

        conformal_prob <- min(
          1,
          ceiling((length(residuals) + 1) * level) / length(residuals)
        )
        adjustment[out] <- unname(
          stats::quantile(
            residuals,
            probs = conformal_prob,
            type = 1,
            names = FALSE
          )
        )
        calibration_n[out] <- length(residuals)
      }

      output_lower <- matrix(
        pmax(
          0,
          sweep(as.matrix(output_point), 2, adjustment, "-")
        ),
        nrow = nrow(output_point),
        ncol = ncol(output_point)
      )
      output_upper <- sweep(
        as.matrix(output_point),
        2,
        adjustment,
        "+"
      )

      names(adjustment) <- names(output_point)
      names(calibration_n) <- names(output_point)
      calibration_info$n <- calibration_n
      calibration_info$adjustment <- adjustment
    }

    output_result <- data.frame(row.names = row.names(newdata))
    for (out in seq_len(ncol(output_point))) {
      point_name <- names(output_point)[out]
      output_result[[point_name]] <- output_point[[out]]
      output_result[[paste0(point_name, "_lower")]] <- output_lower[, out]
      output_result[[paste0(point_name, "_upper")]] <- output_upper[, out]
    }
  }

  score_result <- NULL
  if (compute_score) {
    score_point <- get_scores(
      eval_data = newdata,
      x = x,
      y = y,
      relevant = relevant,
      object = object,
      method = method,
      measure = measure,
      returns = returns,
      direction = direction
    )

    score_members <- matrix(
      NA_real_,
      nrow = nrow(newdata),
      ncol = n_learners
    )

    for (tree in seq_len(n_learners)) {
      score_members[, tree] <- get_scores(
        eval_data = newdata,
        x = x,
        y = y,
        relevant = relevant,
        object = learner_object(tree),
        method = method,
        measure = measure,
        returns = returns,
        direction = direction
      )[[1]]
    }

    score_lower <- apply(
      score_members,
      1,
      finite_quantile,
      prob = probs[1]
    )
    score_upper <- apply(
      score_members,
      1,
      finite_quantile,
      prob = probs[2]
    )

    score_name <- names(score_point)[1]
    score_result <- data.frame(row.names = row.names(newdata))
    score_result[[score_name]] <- score_point[[1]]
    score_result[[paste0(score_name, "_lower")]] <- score_lower
    score_result[[paste0(score_name, "_upper")]] <- score_upper
  }

  structure(
    list(
      output = output_result,
      score = score_result,
      level = level,
      method = method,
      measure = if (compute_score) measure else NULL,
      learners = n_learners,
      interval_type = list(
        output = if (compute_output) {
          if (calibration == "oob") "oob_prediction" else "ensemble"
        } else {
          NULL
        },
        score = if (compute_score) "ensemble" else NULL
      ),
      calibration = calibration_info
    ),
    class = "rf_aces_intervals"
  )
}

#' @title Build a Basis Matrix for Prediction
#'
#' @description
#' Evaluates a fitted model's basis functions on new data.
#'
#' @param newdata
#' A data frame or matrix containing the prepared new inputs.
#'
#' @param model
#' Fitted model information.
#'
#' @param knots
#' Knots used by \code{model}.
#'
#' @param method
#' Model type used to choose the basis functions.
#'
#' @return
#' A matrix of evaluated basis functions.

set_Bmat <- function (
    newdata,
    model,
    knots,
    method
    ) {

  # sample size
  N <- nrow(newdata)

  # initialize B matrix
  B <- matrix(rep(1, N), nrow = N)

  if (method %in% c("aces_forward", "rf_aces")) {

    for (i in 1:nrow(knots)) {

      # variable
      v <- knots[i, 1]

      # knot
      t <- knots[i, 2]

      # basis functions for the new data
      hinge1 <- pmax(0, newdata[, v] - t)
      hinge2 <- pmax(0, t - newdata[, v])

      # update B matrix
      B <- cbind(B, hinge1, hinge2)
    }

  } else if (method == "aces") {

    for (i in 1:nrow(knots)) {

      # variable
      v <- knots[i, 1]

      # knot
      t <- knots[i, 2]

      # basis functions for the new data
      if (knots[i, 3] == "R") {

        # right-hand basis function
        hinge1 <- pmax(0, newdata[, v] - t)

        # update B matrix
        B <- cbind(B, hinge1)

      } else {

        # left-hand basis function
        hinge2 <- pmax(0, t - newdata[, v])

        # update B matrix
        B <- cbind(B, hinge2)

      }
    }

  } else if (method %in% c("aces_cubic", "rf_aces_cubic")) {

    for (status in c("paired", "unpaired")) {

      # select a variable
      for (v in 1:length(knots)) {

        # the variable is not used
        if (is.null(knots[[v]])) next

        # select a knot
        for (t in 1:length(knots[[v]])) {

          # keep the order: first paired and then unpaired
          if (knots[[v]][[t]][["status"]] != status) next

          # knot side
          side <- knots[[v]][[t]][["side"]]

          # cubic knot
          Ct <- knots[[v]][[t]][["t"]]

          # t-
          t0 <- Ct[1]

          # t
          t1 <- Ct[2]

          # t+
          t2 <- Ct[3]

          # t  - t-
          d <- t1 - t0

          # t+ - t
          e <- t2 - t1

          p1 <- (2 * e - d) / (e + d) ^ 2
          r1 <- (d - e) / (e + d) ^ 3

          p2 <- (2 * d - e) / (- e - d) ^ 2
          r2 <- (e - d) / (- e - d) ^ 3

          if (status == "paired" || (status == "unpaired" && side == "R")) {

            term1 <- p1 * (newdata[, v] - t0) ^ 2
            term2 <- r1 * (newdata[, v] - t0) ^ 3

            # C1
            C1 <- ifelse(newdata[, v] <= t0,
                         0,
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2,
                                 newdata[, v] - t1)))
            B  <- cbind(B, C1)

          }

          if (status == "paired" || (status == "unpaired" && side == "L")) {

            term1 <- p2 * (newdata[, v] - t2) ^ 2
            term2 <- r2 * (newdata[, v] - t2) ^ 3

            # C2
            C2 <- ifelse(newdata[, v] <= t0,
                         t1 - newdata[, v],
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2,
                                 0)))
            B  <- cbind(B, C2)

          }
        }
      }
    }

  } else if (method %in% c("aces_quintic", "rf_aces_quintic")) {

    for (status in c("paired", "unpaired")) {

      # select a variable
      for (v in 1:length(knots)) {

        # the variable is not used
        if (is.null(knots[[v]])) next

        # select a knot
        for (t in 1:length(knots[[v]])) {

          # keep the order: first paired and then unpaired
          if (knots[[v]][[t]][["status"]] != status) next

          # knot side
          side <- knots[[v]][[t]][["side"]]

          # quintic knot
          Qt <- knots[[v]][[t]][["t"]]

          # t-
          t0 <- Qt[1]

          # t
          t1 <- Qt[2]

          # t+
          t2 <- Qt[3]

          # t+ - t-
          d  <- t2 - t0

          # t+ - t
          d1 <- t2 - t1

          # t  - t-
          d2 <- t1 - t0

          alpha1 <- (6 * d1 - 4 * d2) / d ^ 3
          alpha2 <- (4 * d1 - 6 * d2) / d ^ 3

          beta1 <- (7 * d2 - 8 * d1) / d ^ 4
          beta2 <- (7 * d1 - 8 * d2) / d ^ 4

          gamma1 <- gamma2 <- (3 * d1 - 3 * d2) / d ^ 5

          if (status == "paired" || (status == "unpaired" && side == "R")) {

            term1 <- alpha1 * (newdata[, v] - t0) ^ 3
            term2 <- beta1  * (newdata[, v] - t0) ^ 4
            term3 <- gamma1 * (newdata[, v] - t0) ^ 5

            # Q1
            Q1 <- ifelse(newdata[, v] <= t0,
                         0,
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2 + term3,
                                 newdata[, v] - t1)))

            B  <- cbind(B, Q1)

          }

          if (status == "paired" || (status == "unpaired" && side == "L")) {

            term1 <- alpha2 * (newdata[, v] - t2) ^ 3
            term2 <- beta2  * (newdata[, v] - t2) ^ 4
            term3 <- gamma2 * (newdata[, v] - t2) ^ 5

            # Q2
            Q2 <- ifelse(newdata[, v] <= t0,
                         t1 - newdata[, v],
                         (ifelse((newdata[, v] > t0) & (newdata[, v] < t2),
                                 term1 + term2 + term3,
                                 0)))
            B  <- cbind(B, Q2)

          }
        }
      }
    }

  } else {

    stop("Not available method. Please, check help(\"predict\")")

  }

  return(B)

}
