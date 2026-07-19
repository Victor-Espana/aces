#' @title Rank DMUs with ACES or RF-ACES
#'
#' @description
#' Ranks decision-making units (DMUs) from their efficiency scores. For an
#' RF-ACES object, the score is computed separately against every eligible
#' learner-specific technology and then averaged. For a standard ACES object,
#' the ranking uses its single fitted technology.
#'
#' @param object
#' A fitted \code{aces} or \code{rf_aces} object.
#'
#' @param eval_data
#' Evaluation specification. Use \code{"oob"} to rank the RF-ACES training DMUs
#' only with learners for which each DMU was out of bag. To rank external DMUs,
#' supply a named list with elements \code{data}, \code{x}, and \code{y}, where
#' \code{data} is a data frame or matrix and \code{x} and \code{y} are its input
#' and output column indexes.
#'
#' @param relevant
#' If \code{TRUE}, compute each efficiency score using only original inputs
#' represented in the fitted basis functions, following \code{get_scores()}.
#'
#' @param method
#' Fitted representation used to define each technology. When \code{NULL}, use
#' \code{"aces"} for an \code{aces} object and \code{"rf_aces"} for an
#' \code{rf_aces} object. Smoothed representations supported by
#' \code{get_scores()} are also available.
#'
#' @param measure
#' Efficiency measure passed to \code{get_scores()}: \code{"rad_out"},
#' \code{"rad_inp"}, \code{"ddf"}, \code{"rsl_out"}, \code{"rsl_inp"},
#' \code{"wam_mip"}, \code{"wam_nor"}, \code{"wam_ram"},
#' \code{"wam_bam"}, or \code{"rf_aces_rad_out"} for RF-ACES.
#'
#' @param returns
#' Returns-to-scale assumption passed to \code{get_scores()}: \code{"constant"}
#' or \code{"variable"}.
#'
#' @param direction
#' Direction matrix used when \code{measure = "ddf"}, with one row per evaluated
#' DMU and the input directions followed by the output directions.
#'
#' @param digits
#' Non-negative integer giving the number of decimal places displayed for the
#' mean and standard deviation. Ranks are computed before rounding.
#'
#' @details
#' For external data, every RF-ACES learner is eligible. If \eqn{s_{ib}} is the
#' score of DMU \eqn{i} against learner \eqn{b}, the reported score is
#' \eqn{B^{-1}\sum_b s_{ib}}. With \code{eval_data = "oob"}, learner \eqn{b}
#' contributes to DMU \eqn{i} only when \eqn{i} was absent from that learner's
#' bootstrap sample. Consequently, different training DMUs may use different
#' numbers of learners.
#'
#' The learner mean is generally not equal to the score obtained from the
#' aggregated RF-ACES technology: averaging technologies before solving the
#' efficiency problem and averaging scores after solving it are distinct
#' operations. The standard deviation describes dispersion across learner
#' scores; it is not a calibrated standard error or confidence interval.
#'
#' Ranks follow the economic interpretation of the selected measure. Larger
#' values are better for input-oriented radial and Russell scores. Smaller
#' values are better for output-oriented, directional-distance, weighted
#' additive, and direct RF-ACES output scores. Ties receive the same minimum
#' rank.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{espana2024rf}{aces} \cr \cr
#' \insertRef{banker1984}{aces} \cr \cr
#' \insertRef{chambers1998}{aces} \cr \cr
#' \insertRef{fare1978}{aces}
#'
#' @return
#' A data frame sorted by rank with columns \code{unit}, \code{rank},
#' \code{score_mean}, \code{score_sd}, and \code{learners_used}. The result has
#' attributes identifying the evaluation mode, method, measure, and ranking
#' direction. For ACES, \code{learners_used} is one and \code{score_sd} is
#' unavailable.
#'
#' @export
ranking_aces <- function(
  object,
  eval_data = "oob",
  relevant = FALSE,
  method = NULL,
  measure = "rad_out",
  returns = "variable",
  direction = NULL,
  digits = 3
) {
  if (!inherits(object, "aces") && !inherits(object, "rf_aces")) {
    stop("'object' must be an aces or rf_aces object.", call. = FALSE)
  }
  if (!is.logical(relevant) || length(relevant) != 1L || is.na(relevant)) {
    stop("'relevant' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) ||
      digits < 0 || digits != floor(digits)) {
    stop("'digits' must be a non-negative integer.", call. = FALSE)
  }

  is_rf <- inherits(object, "rf_aces")
  if (is.null(method)) method <- if (is_rf) "rf_aces" else "aces"

  evaluation <- .ranking_aces_evaluation(object, eval_data)
  data <- evaluation$data
  x <- evaluation$x
  y <- evaluation$y

  # Reuse the public score validation so method, measure, direction, and
  # variable names follow exactly the same contract as get_scores().
  display_errors_scores(
    data = data,
    x = x,
    y = y,
    object = object,
    method = method,
    measure = measure,
    returns = returns,
    direction = direction
  )

  if (is_rf) {
    learners <- object[["forest"]]
    if (length(learners) == 0L) {
      stop("The RF-ACES object contains no learners.", call. = FALSE)
    }

    score_members <- matrix(
      NA_real_,
      nrow = nrow(data),
      ncol = length(learners)
    )

    for (tree in seq_along(learners)) {
      learner_object <- .ranking_aces_learner_object(object, tree)
      learner_scores <- .ranking_aces_get_scores(
        eval_data = data,
        x = x,
        y = y,
        relevant = relevant,
        object = learner_object,
        method = method,
        measure = measure,
        returns = returns,
        direction = direction
      )
      if (length(learner_scores) != nrow(data)) {
        stop("A learner returned an invalid number of efficiency scores.", call. = FALSE)
      }

      if (evaluation$mode == "oob") {
        sample_bag <- learners[[tree]][["sample_bag"]]
        if (is.null(sample_bag)) {
          stop(
            "OOB ranking is unavailable because bag membership is not stored in every learner.",
            call. = FALSE
          )
        }
        if (anyNA(sample_bag) || any(sample_bag < 1L | sample_bag > nrow(data))) {
          stop("Stored bootstrap indexes are invalid for the training data.", call. = FALSE)
        }
        eligible <- which(!seq_len(nrow(data)) %in% unique(sample_bag))
      } else {
        eligible <- seq_len(nrow(data))
      }

      score_members[eligible, tree] <- learner_scores[eligible]
    }
  } else {
    if (evaluation$mode == "oob") {
      stop("OOB ranking is available only for RF-ACES objects.", call. = FALSE)
    }
    single_score <- .ranking_aces_get_scores(
      eval_data = data,
      x = x,
      y = y,
      relevant = relevant,
      object = object,
      method = method,
      measure = measure,
      returns = returns,
      direction = direction
    )
    if (length(single_score) != nrow(data)) {
      stop("The ACES model returned an invalid number of efficiency scores.", call. = FALSE)
    }
    score_members <- matrix(single_score, ncol = 1L)
  }

  result <- .ranking_aces_summarize(
    score_members = score_members,
    unit_names = row.names(data),
    measure = measure,
    digits = digits
  )
  attr(result, "evaluation") <- evaluation$mode
  attr(result, "method") <- method
  attr(result, "measure") <- measure
  attr(result, "ranking_direction") <- if (
    .ranking_aces_higher_is_better(measure)
  ) "descending" else "ascending"
  class(result) <- c("aces_ranking", "data.frame")
  result
}

.ranking_aces_evaluation <- function(object, eval_data) {
  if (is.character(eval_data)) {
    if (length(eval_data) != 1L || is.na(eval_data) ||
        !identical(tolower(eval_data), "oob")) {
      stop("Character 'eval_data' must be exactly 'oob'.", call. = FALSE)
    }
    if (!inherits(object, "rf_aces")) {
      stop("OOB ranking is available only for RF-ACES objects.", call. = FALSE)
    }
    data <- as.data.frame(object[["data"]][["df"]], check.names = FALSE)
    return(list(
      data = data,
      x = object[["data"]][["x"]],
      y = object[["data"]][["y"]],
      mode = "oob"
    ))
  }

  if (!is.list(eval_data)) {
    stop(
      "'eval_data' must be 'oob' or a named list with data, x, and y.",
      call. = FALSE
    )
  }
  list_names <- names(eval_data)
  if (is.null(list_names) || any(list_names == "") || anyDuplicated(list_names)) {
    stop("Every element of 'eval_data' must have a unique name.", call. = FALSE)
  }
  required <- c("data", "x", "y")
  missing <- setdiff(required, list_names)
  unknown <- setdiff(list_names, required)
  if (length(missing) > 0L) {
    stop(
      "'eval_data' is missing: ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (length(unknown) > 0L) {
    stop(
      "Unknown element in 'eval_data': ", paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.data.frame(eval_data$data) && !is.matrix(eval_data$data)) {
    stop("'eval_data$data' must be a data frame or matrix.", call. = FALSE)
  }

  data <- as.data.frame(eval_data$data, check.names = FALSE)
  validate_indexes <- function(index, label) {
    if (!is.numeric(index) || length(index) == 0L || anyNA(index) ||
        any(!is.finite(index)) ||
        any(index != floor(index)) || any(index < 1L | index > ncol(data))) {
      stop("'eval_data$", label, "' must contain valid column indexes.", call. = FALSE)
    }
    as.integer(index)
  }
  x <- validate_indexes(eval_data$x, "x")
  y <- validate_indexes(eval_data$y, "y")
  if (anyDuplicated(x) || anyDuplicated(y) || length(intersect(x, y)) > 0L) {
    stop("'eval_data$x' and 'eval_data$y' must contain distinct columns.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("'eval_data$data' must contain at least one DMU.", call. = FALSE)
  }
  if (!all(vapply(data[, c(x, y), drop = FALSE], is.numeric, logical(1)))) {
    stop("All input and output columns in 'eval_data$data' must be numeric.", call. = FALSE)
  }
  if (anyNA(data[, c(x, y), drop = FALSE])) {
    stop("Input and output columns in 'eval_data$data' cannot contain missing values.", call. = FALSE)
  }
  expected_x <- object[["data"]][["xnames"]]
  expected_y <- object[["data"]][["ynames"]]
  if (!identical(colnames(data)[x], expected_x) ||
      !identical(colnames(data)[y], expected_y)) {
    stop(
      "Input and output names in 'eval_data$data' must match the fitted model in the same order.",
      call. = FALSE
    )
  }

  list(data = data, x = x, y = y, mode = "external")
}

.ranking_aces_learner_object <- function(object, index) {
  learner <- object[["forest"]][[index]]
  answer <- object
  answer[["forest"]] <- list(learner)
  answer[["technology"]] <- learner[["technology"]]
  answer[["control"]][["learners"]] <- 1L
  answer
}

.ranking_aces_get_scores <- function(...) {
  as.numeric(get_scores(...)[[1]])
}

.ranking_aces_higher_is_better <- function(measure) {
  measure %in% c("rad_inp", "rsl_inp")
}

.ranking_aces_summarize <- function(score_members, unit_names, measure, digits) {
  finite <- is.finite(score_members)
  learners_used <- rowSums(finite)
  score_mean <- rowSums(replace(score_members, !finite, 0))
  score_mean[learners_used > 0L] <- score_mean[learners_used > 0L] /
    learners_used[learners_used > 0L]
  score_mean[learners_used == 0L] <- NA_real_

  score_sd <- vapply(seq_len(nrow(score_members)), function(i) {
    values <- score_members[i, finite[i, ]]
    if (length(values) > 1L) stats::sd(values) else NA_real_
  }, numeric(1))

  if (all(learners_used == 0L)) {
    stop("No finite efficiency score is available for any DMU.", call. = FALSE)
  }
  if (any(learners_used == 0L)) {
    warning(
      sum(learners_used == 0L),
      " DMU(s) have no eligible finite learner score and receive no rank.",
      call. = FALSE
    )
  }

  rank_value <- if (.ranking_aces_higher_is_better(measure)) {
    -score_mean
  } else {
    score_mean
  }
  unit_rank <- rank(rank_value, ties.method = "min", na.last = "keep")

  result <- data.frame(
    unit = as.character(unit_names),
    rank = as.integer(unit_rank),
    score_mean = round(score_mean, digits),
    score_sd = round(score_sd, digits),
    learners_used = as.integer(learners_used),
    check.names = FALSE,
    row.names = NULL
  )
  result <- result[order(result$rank, na.last = TRUE), , drop = FALSE]
  rownames(result) <- NULL
  result
}
