make_ranking_tree <- function(frontier, sample_bag) {
  list(
    sample_bag = sample_bag,
    technology = list(
      rf_aces = list(
        xmat = matrix(c(1, 2), ncol = 1),
        ymat = matrix(rep(frontier, 2), ncol = 1)
      )
    )
  )
}

make_ranking_forest <- function() {
  trees <- list(
    make_ranking_tree(2, 2),
    make_ranking_tree(4, 1),
    make_ranking_tree(6, c(1, 2))
  )

  structure(
    list(
      data = list(
        df = data.frame(
          x = c(1, 2),
          y = c(1, 3),
          row.names = c("a", "b")
        ),
        x = 1,
        y = 2,
        xnames = "x",
        ynames = "y",
        rownames = c("a", "b")
      ),
      control = list(
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

test_that("ranking_aces averages external scores across learners", {
  forest <- make_ranking_forest()
  new_data <- data.frame(
    x = c(1.5, 1.5),
    y = c(1, 2),
    row.names = c("weak", "strong")
  )

  ranking <- ranking_aces(
    object = forest,
    eval_data = list(data = new_data, x = 1, y = 2),
    method = "rf_aces",
    measure = "rad_out"
  )

  expect_s3_class(ranking, "aces_ranking")
  expect_equal(ranking$unit, c("strong", "weak"))
  expect_equal(ranking$rank, c(1L, 2L))
  expect_equal(ranking$score_mean, c(2, 4))
  expect_equal(ranking$score_sd, c(1, 2))
  expect_equal(ranking$learners_used, c(3L, 3L))
  expect_equal(attr(ranking, "evaluation"), "external")
  expect_equal(attr(ranking, "ranking_direction"), "ascending")
})

test_that("ranking_aces uses only eligible OOB learners", {
  forest <- make_ranking_forest()

  ranking <- ranking_aces(
    object = forest,
    eval_data = "oob",
    method = "rf_aces",
    measure = "rad_out"
  )

  expect_equal(ranking$unit, c("b", "a"))
  expect_equal(ranking$rank, c(1L, 2L))
  expect_equal(ranking$score_mean, c(1.333, 2))
  expect_true(all(is.na(ranking$score_sd)))
  expect_equal(ranking$learners_used, c(1L, 1L))
  expect_equal(attr(ranking, "evaluation"), "oob")
})

test_that("ranking direction follows the efficiency measure", {
  members <- matrix(
    c(
      0.8, 0.9,
      0.6, 0.7
    ),
    nrow = 2,
    byrow = TRUE
  )

  ranking <- .ranking_aces_summarize(
    score_members = members,
    unit_names = c("more_efficient", "less_efficient"),
    measure = "rad_inp",
    digits = 3
  )

  expect_equal(ranking$unit, c("more_efficient", "less_efficient"))
  expect_equal(ranking$rank, c(1L, 2L))
})

test_that("ranking_aces validates the evaluation specification", {
  forest <- make_ranking_forest()

  expect_error(
    ranking_aces(
      object = forest,
      eval_data = list(data = forest$data$df, x = 1)
    ),
    "missing: y",
    fixed = TRUE
  )
  expect_error(
    ranking_aces(object = structure(list(), class = "aces")),
    "available only for RF-ACES",
    fixed = TRUE
  )
  expect_error(
    ranking_aces(
      object = forest,
      eval_data = list(
        data = data.frame(y = 1:2, x = 2:3),
        x = 1,
        y = 2
      )
    ),
    "must match the fitted model in the same order",
    fixed = TRUE
  )
})
