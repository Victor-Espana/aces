test_that("multi-output filtering checks every output", {
  input_1 <- 1:20
  input_2 <- 20:1

  data <- cbind(
    input_1 = input_1,
    input_2 = input_2,
    output_1 = input_1,
    output_2 = input_2
  )

  result <- quick_aces_correlation_filter(
    data = data,
    x = 1:2,
    y = 3:4
  )

  expect_equal(result$keep, c(input_1 = TRUE, input_2 = TRUE))
  expect_equal(result$spearman["input_2", "output_1"], -1)
  expect_equal(result$spearman["input_2", "output_2"], 1)
})

test_that("an input is removed only when it fails for every output", {
  output <- 1:20

  data <- cbind(
    increasing = output,
    decreasing = rev(output),
    constant = rep(1, length(output)),
    output = output
  )

  result <- quick_aces_correlation_filter(
    data = data,
    x = 1:3,
    y = 4
  )

  expect_equal(
    result$keep,
    c(increasing = TRUE, decreasing = FALSE, constant = FALSE)
  )
})

test_that("the filter handles data without positive correlations", {
  data <- cbind(
    input = 20:1,
    output = 1:20
  )

  result <- quick_aces_correlation_filter(
    data = data,
    x = 1,
    y = 2
  )

  expect_false(unname(result$keep))
  expect_equal(result$thresholds, c(spearman = 0.1, kendall = 0.1))
})

test_that("Quick ACES maps retained terms to original input columns", {
  expect_equal(
    quick_aces_retained_inputs(
      keep = c(TRUE, FALSE, TRUE),
      x_vars = c(2, 4, 7),
      max_degree = 1
    ),
    c(2, 7)
  )
})

test_that("a retained interaction keeps each of its original inputs", {
  expect_equal(
    quick_aces_retained_inputs(
      keep = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
      x_vars = c(2, 4, 7),
      max_degree = 2
    ),
    c(2, 4)
  )

  expect_equal(
    quick_aces_retained_inputs(
      keep = c(FALSE, FALSE, FALSE, TRUE),
      x_vars = c(2, 4, 7),
      max_degree = list(c(2, 7))
    ),
    c(2, 7)
  )
})
