simulate_production_data <- function(seed) {
  set.seed(seed)
  cobb_douglas_XnY1(N = 100, nX = 3)
}


test_that("tune_model_aces tunes a complex simulated frontier", {
  data <- simulate_production_data(2026)
  partition <- data_partition(
    data = data,
    x = 1:3,
    y = 4,
    validation = "holdout",
    resampling_args = list(repeats = 1, train_prop = 0.8),
    delta = 1.25,
    seed = 2026
  )

  invisible(capture.output(
    tuned <- tune_model_aces(
      partition = partition,
      mul_BF = list(
        max_degree = c(2),
        inter_cost = c(0.01, 0.05)
      ),
      metric = "mae",
      err_red = c(0.01, 0.05),
      kn_penalty = c(1, 2),
      validation_metric = c("rmse", "mae"),
      verbose = TRUE
    )
  ))

  expect_equal(nrow(data), 100)
  expect_s3_class(partition, "aces_partition")
  expect_equal(partition$delta, 1.5)
  expect_equal(partition$seed, 2026)
  expect_s3_class(tuned, "aces_tune")
  expect_s3_class(tuned$final_model, "aces")
  expect_equal(nrow(tuned$results), 8)
  expect_length(tuned$validation_metrics, 2)
  expect_named(tuned$best_scores, c("rmse", "mae"))
  expect_true(all(is.finite(tuned$best_scores)))
  expect_true(all(is.finite(tuned$results$mean_rmse)))
  expect_true(all(is.finite(tuned$results$mean_mae)))
  model_args <- names(formals(aces))[-(1:3)]
  expect_identical(
    names(formals(tune_model_aces))[seq_along(model_args) + 1L],
    model_args
  )
  expect_false(any(c("delta", "seed") %in% names(formals(tune_model_aces))))
  expect_named(
    tuned$best_params,
    c(
      "scale_data", "quick_aces", "max_degree", "inter_cost",
      "metric", "mono", "conc",
      "max_terms", "err_red", "kn_grid", "minspan", "endspan",
      "kn_penalty"
    )
  )
  expect_true(tuned$best_params$max_degree %in% 2)
  expect_true(tuned$best_params$inter_cost %in% c(0.01, 0.05))

  # select_best
  alt <- select_best(tuned, metric = "mae")
  expect_named(alt, c("best_params", "best_scores", "config_id"))
  expect_named(alt$best_scores, c("rmse", "mae"))
  expect_true(alt$config_id %in% tuned$results$config_id)

  expect_error(data_partition(data, x = "x1", y = 7), "column indexes")
})


test_that("tune_model_rf_aces tunes a complex simulated forest", {
  data <- simulate_production_data(2027)
  partition <- data_partition(
    data = data,
    x = 1:6,
    y = 7,
    validation = "holdout",
    resampling_args = list(repeats = 2, train_prop = 0.8),
    delta = 1.5,
    seed = 2027
  )

  invisible(capture.output(
    tuned <- tune_model_rf_aces(
      partition = partition,
      quick_aces = TRUE,
      mul_BF = list(
        max_degree = c(1, 2),
        inter_cost = 0.05
      ),
      learners = 3,
      bag_size = c(60, 80),
      max_feats = c(2, 4),
      early_stopping = list(
        ma_window = 2,
        tolerance = 2
      ),
      max_terms = 15,
      err_red = 0.05,
      validation_metric = "rmse",
      verbose = FALSE
    )
  ))

  expect_equal(nrow(data), 100)
  expect_s3_class(partition, "aces_partition")
  expect_equal(partition$delta, 1.5)
  expect_equal(partition$seed, 2027)
  expect_s3_class(tuned, "aces_tune")
  expect_s3_class(tuned$final_model, "rf_aces")
  expect_equal(nrow(tuned$results), 8)
  expect_length(tuned$validation_metrics, 1)
  expect_named(tuned$best_scores, "rmse")
  expect_true(all(is.finite(tuned$best_scores)))
  expect_true(all(is.finite(tuned$results$mean_rmse)))
  model_args <- names(formals(rf_aces))[-(1:3)]
  expect_identical(
    names(formals(tune_model_rf_aces))[seq_along(model_args) + 1L],
    model_args
  )
  expect_false(any(c("delta", "seed") %in% names(formals(tune_model_rf_aces))))
  expect_true("ma_window" %in% names(tuned$best_params))
  expect_true("tolerance" %in% names(tuned$best_params))
  expect_true(tuned$best_params$bag_size %in% c(60, 80))
  expect_true(tuned$best_params$max_feats %in% c(2, 4))
})
