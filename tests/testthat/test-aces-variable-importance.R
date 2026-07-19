test_that("prepared interactions use colon notation", {
  data <- data.frame(
    capital = c(1, 2, 3),
    output = c(2, 4, 8),
    labour = c(3, 4, 5),
    materials = c(2, 3, 4),
    check.names = FALSE
  )

  pairwise <- set_data(
    data = data,
    x = c(1, 3),
    y = 2,
    max_degree = 2
  )

  expect_equal(
    colnames(pairwise),
    c("capital", "labour", "capital:labour", "output")
  )
  expect_equal(
    unname(pairwise[, "capital:labour"]),
    data$capital * data$labour
  )

  pairwise_matrix <- set_data(
    data = as.matrix(data),
    x = c(1, 3),
    y = 2,
    max_degree = 2
  )
  expect_equal(colnames(pairwise_matrix), colnames(pairwise))

  explicit <- set_data(
    data = data,
    x = c(1, 3, 4),
    y = 2,
    max_degree = list(c(1, 4))
  )

  expect_equal(
    colnames(explicit),
    c("capital", "labour", "materials", "capital:materials", "output")
  )

  grids <- set_knots_grid(
    data = pairwise,
    n_input_1 = 2,
    n_input_2 = 3,
    kn_grid = list(1:3, 3:5),
    quick_aces = FALSE,
    dea_scores = rep(1, nrow(pairwise))
  )

  expect_equal(names(grids), c("capital", "labour", "capital:labour"))
})

test_that("ACES importance credits the best candidate for every input", {
  reduction <- matrix(
    c(
      0.20, 0.10, NA,
      0.05, 0.30, NA
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("iteration_1", "iteration_2"),
      c("capital", "labour", "capital:labour")
    )
  )

  object <- structure(
    list(
      methods = list(
        aces_forward = list(variable_reduction = reduction)
      )
    ),
    class = "aces"
  )

  raw <- aces_varimp(object, normalize = FALSE)
  expect_equal(raw$variable, c("labour", "capital", "capital:labour"))
  expect_equal(raw$importance, c(0.4, 0.25, 0))
  expect_equal(raw$iterations_evaluated, c(2L, 2L, 0L))

  relative <- aces_varimp(object)
  expect_equal(relative$importance, c(100, 62.5, 0))

  history <- aces_varimp(object, control = list(type = "iterations"))
  expect_equal(history$iteration, 1:2)
  expect_equal(names(history)[-1], colnames(reduction))
  expect_equal(as.matrix(history[, -1]), reduction, ignore_attr = TRUE)
})

test_that("ACES importance validates its inputs", {
  expect_error(aces_varimp(list()), "must be an aces or rf_aces object", fixed = TRUE)

  old_object <- structure(
    list(methods = list(aces_forward = list())),
    class = "aces"
  )

  expect_error(
    aces_varimp(old_object),
    "refit it with the current version",
    fixed = TRUE
  )
})

test_that("RF-ACES supports forward, external permutation, and OOB permutation", {
  data <- data.frame(x = 1:6, y = 2 * (1:6))
  learner_model <- list(
    knots = data.frame(xi = 1, t = 0),
    coefs = matrix(c(0, 2, 0), ncol = 1),
    variable_reduction = matrix(
      c(0.2, 0.1),
      ncol = 1,
      dimnames = list(c("iteration_1", "iteration_2"), "x")
    )
  )
  learner_1 <- list(
    methods = list(rf_aces = learner_model),
    sample_bag = c(1, 2, 3, 3)
  )
  learner_2 <- list(
    methods = list(rf_aces = learner_model),
    sample_bag = c(4, 5, 6, 6)
  )
  object <- structure(
    list(
      data = list(
        df = data,
        x = 1,
        y = 2,
        xnames = "x",
        ynames = "y"
      ),
      control = list(
        scale = list(is_scaled = FALSE, mean_x = 1, mean_y = 1),
        max_degree = 1,
        metric = "mse"
      ),
      forest = list(learner_1, learner_2)
    ),
    class = "rf_aces"
  )

  forward <- aces_varimp(object, normalize = FALSE)
  expect_equal(forward$importance, 0.3)
  expect_equal(forward$learners_evaluated, 2L)
  expect_equal(forward$iterations_evaluated, 4L)

  history <- aces_varimp(object, control = list(type = "iterations"))
  expect_equal(history$learner, c(1L, 1L, 2L, 2L))
  expect_equal(history$iteration, c(1L, 2L, 1L, 2L))

  external <- aces_varimp(
    object,
    importance = "permutation",
    control = list(eval_data = data, repeats = 2, seed = 1)
  )
  expect_equal(external$baseline_loss, 0)
  expect_equal(external$importance, 100)
  expect_equal(attr(external, "evaluation"), "external")
  expect_equal(attr(external, "method"), "rf_aces")

  oob <- aces_varimp(
    object,
    importance = "permutation",
    control = list(repeats = 2, seed = 1)
  )
  expect_equal(oob$baseline_loss, 0)
  expect_equal(oob$importance, 100)
  expect_equal(attr(oob, "evaluation"), "oob")

  expect_error(
    aces_varimp(
      object,
      importance = "permutation",
      control = list(eval_data = NULL)
    ),
    "must be 'oob' or a data frame or matrix",
    fixed = TRUE
  )
})

test_that("permutation importance evaluates the fitted ACES model", {
  data <- data.frame(x = 1:6, y = 2 * (1:6))
  model <- list(
    knots = data.frame(
      xi = 1,
      t = 0,
      side = "R",
      status = "unpaired"
    ),
    coefs = matrix(c(0, 2), ncol = 1)
  )
  object <- structure(
    list(
      data = list(
        df = data,
        x = 1,
        y = 2,
        xnames = "x",
        ynames = "y"
      ),
      control = list(
        scale = list(is_scaled = FALSE, mean_x = 1, mean_y = 1),
        max_degree = 1,
        metric = "mse"
      ),
      methods = list(aces = model)
    ),
    class = "aces"
  )

  result <- aces_varimp(
    object,
    importance = "permutation",
    control = list(
      eval_data = data,
      repeats = 5,
      seed = 42
    )
  )

  expect_equal(result$variable, "x")
  expect_equal(result$importance, 100)
  expect_equal(result$baseline_loss, 0)
  expect_gt(result$comparison_loss, 0)
  expect_equal(attr(result, "importance"), "permutation")
  expect_equal(attr(result, "evaluation"), "external")
})

test_that("importance controls are specific to each approach", {
  object <- structure(
    list(methods = list(aces_forward = list(
      variable_reduction = matrix(0.1, nrow = 1, dimnames = list(NULL, "x"))
    ))),
    class = "aces"
  )

  expect_error(
    aces_varimp(object, control = list(repeats = 2)),
    "Unknown control option for importance = 'forward'",
    fixed = TRUE
  )
  expect_error(
    aces_varimp(
      object,
      importance = "permutation",
      control = list(verbose = TRUE)
    ),
    "Unknown control option for importance = 'permutation'",
    fixed = TRUE
  )
  expect_error(
    aces_varimp(object, importance = "loco", control = list(seed = 1)),
    "Unknown control option for importance = 'loco'",
    fixed = TRUE
  )
})

test_that("LOCO removes an input from interactions and custom knot grids", {
  object <- list(
    data = list(x = c(1, 3, 4)),
    control = list(
      fit_spec = list(
        mul_BF = list(
          max_degree = list(c(1, 3), c(1, 4)),
          inter_cost = 0.05
        ),
        kn_grid = list(1:3, 4:6, 7:9)
      )
    )
  )

  spec <- .aces_loco_spec(object, dropped_input = 3)

  expect_equal(spec$x, c(1, 4))
  expect_equal(spec$spec$mul_BF$max_degree, list(c(1, 4)))
  expect_equal(spec$spec$kn_grid, list(1:3, 7:9))
  expect_false(spec$used_fallback)
})

test_that("RF-ACES LOCO retains forest controls and caps max_feats", {
  object <- structure(
    list(
      data = list(x = c(1, 3, 4)),
      control = list(
        fit_spec = list(
          mul_BF = list(
            max_degree = list(c(1, 3), c(1, 4)),
            inter_cost = 0.05
          ),
          kn_grid = list(1:3, 4:6, 7:9),
          max_feats = 3
        )
      )
    ),
    class = "rf_aces"
  )

  spec <- .aces_loco_spec(object, dropped_input = 3)

  expect_equal(spec$x, c(1, 4))
  expect_equal(spec$spec$mul_BF$max_degree, list(c(1, 4)))
  expect_equal(spec$spec$kn_grid, list(1:3, 7:9))
  expect_equal(spec$spec$max_feats, 2)
})
