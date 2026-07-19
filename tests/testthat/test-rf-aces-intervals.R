make_interval_tree <- function(frontier) {
  model <- list(
    knots = data.frame(
      xi = 1,
      t = 1,
      status = "paired"
    ),
    coefs = matrix(c(frontier, 0, 0), ncol = 1)
  )

  list(
    methods = list(rf_aces = model),
    technology = list(
      rf_aces = list(
        xmat = matrix(c(1, 2), ncol = 1),
        ymat = matrix(rep(frontier, 2), ncol = 1)
      )
    )
  )
}

make_interval_forest <- function() {
  trees <- lapply(c(2, 4, 6), make_interval_tree)
  trees[[1]]$sample_bag <- 2
  trees[[2]]$sample_bag <- 1
  trees[[3]]$sample_bag <- c(1, 2)

  structure(
    list(
      data = list(
        df = data.frame(x = c(1, 2), y = c(1, 3)),
        x = 1,
        y = 2,
        xnames = "x",
        ynames = "y",
        rownames = c("a", "b")
      ),
      control = list(
        max_degree = 1,
        learners = 3,
        scale = list(
          is_scaled = FALSE,
          mean_x = 1,
          mean_y = 1
        )
      ),
      forest = trees,
      technology = list(
        rf_aces = list(
          xmat = matrix(c(1, 2), ncol = 1),
          ymat = matrix(c(4, 4), ncol = 1)
        )
      )
    ),
    class = "rf_aces"
  )
}

test_that("rf_aces_intervals returns output and score learner quantiles", {
  forest <- make_interval_forest()
  newdata <- data.frame(x = 1.5, y = 2, row.names = "dmu_1")

  intervals <- rf_aces_intervals(
    object = forest,
    newdata = newdata,
    x = 1,
    y = 2,
    level = 0.5
  )

  expect_s3_class(intervals, "rf_aces_intervals")
  expect_identical(intervals$learners, 3L)
  expect_named(
    intervals$output,
    c("y_pred", "y_pred_lower", "y_pred_upper")
  )
  expect_equal(as.numeric(intervals$output[1, ]), c(4, 3, 5))

  score_name <- "aces_vrt_rf_aces_rad_out"
  expect_named(
    intervals$score,
    c(score_name, paste0(score_name, "_lower"), paste0(score_name, "_upper"))
  )
  expect_equal(as.numeric(intervals$score[1, ]), c(2, 1.5, 2.5))
  expect_identical(row.names(intervals$output), "dmu_1")
  expect_identical(row.names(intervals$score), "dmu_1")

  technology_intervals <- rf_aces_intervals(
    object = forest,
    newdata = newdata,
    x = 1,
    y = 2,
    level = 0.5,
    type = "score",
    measure = "rad_out"
  )
  expect_equal(
    as.numeric(technology_intervals$score[1, ]),
    c(2, 1.5, 2.5)
  )
})

test_that("rf_aces_intervals validates interval-specific arguments", {
  forest <- make_interval_forest()
  newdata <- data.frame(x = 1.5, y = 2)

  expect_error(
    rf_aces_intervals(forest, newdata, x = 1, type = "score"),
    "'y' is required"
  )
  expect_error(
    rf_aces_intervals(forest, newdata, x = 1, level = 1),
    "strictly between"
  )

  output_only <- rf_aces_intervals(
    forest,
    newdata,
    x = 1,
    type = "output",
    level = 0.8
  )
  expect_null(output_only$score)
  expect_equal(output_only$level, 0.8)
})

test_that("rf_aces_intervals calibrates output limits with OOB residuals", {
  forest <- make_interval_forest()
  newdata <- data.frame(x = 1.5, y = 2)

  calibrated <- rf_aces_intervals(
    forest,
    newdata,
    x = 1,
    type = "output",
    level = 0.5,
    calibration = "oob",
    min_oob = 1
  )

  expect_equal(as.numeric(calibrated$output[1, ]), c(4, 3, 5))
  expect_identical(calibrated$interval_type$output, "oob_prediction")
  expect_identical(calibrated$interval_type$score, NULL)
  expect_identical(calibrated$calibration$target, "observed_output")
  expect_equal(unname(calibrated$calibration$adjustment), 1)
  expect_equal(unname(calibrated$calibration$n), 2)

  expect_error(
    rf_aces_intervals(
      forest,
      newdata,
      x = 1,
      y = 2,
      type = "score",
      calibration = "oob"
    ),
    "requires type"
  )
})
