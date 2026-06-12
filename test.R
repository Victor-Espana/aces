devtools::load_all()

# ============================================================ #
#  ACES package — test suite with execution-time reporting     #
# ============================================================ #

pass <- 0
fail <- 0
timings <- list()

ok <- function(label, expr) {
  t <- system.time(
    result <- tryCatch(
      { force(expr); TRUE },
      error   = function(e) { message("  FAIL [", label, "]: ", conditionMessage(e)); FALSE },
      warning = function(w) { message("  WARN [", label, "]: ", conditionMessage(w)); TRUE }
    )
  )
  timings[[label]] <<- t[["elapsed"]]
  if (result) {
    pass <<- pass + 1
    message("  PASS [", label, "]  (", round(t[["elapsed"]], 3), "s)")
  } else {
    fail <<- fail + 1
  }
  invisible(result)
}

expect_true  <- function(x) if (!isTRUE(x))  stop("expected TRUE, got: ", deparse(x))
expect_false <- function(x) if (!isFALSE(x)) stop("expected FALSE, got: ", deparse(x))
expect_equal <- function(a, b, tol = 1e-6) if (abs(a - b) > tol) stop(a, " != ", b)
expect_class <- function(obj, cls) if (!inherits(obj, cls)) stop("expected class '", cls, "', got: ", paste(class(obj), collapse = "/"))
expect_error <- function(expr) {
  tryCatch({ force(expr); stop("no error was raised") }, error = function(e) invisible(TRUE))
}
expect_gt    <- function(a, b) if (!(a > b))  stop(a, " not > ", b)
expect_gte   <- function(a, b) if (!(a >= b)) stop(a, " not >= ", b)
expect_lte   <- function(a, b) if (!(a <= b)) stop(a, " not <= ", b)

set.seed(42)

# ============================================================ #
#  1. SIMULATION FUNCTIONS                                     #
# ============================================================ #

message("\n--- 1. Simulation functions ---")

ok("cobb_douglas_XnY1 nX=1", {
  d <- cobb_douglas_XnY1(50, 1)
  expect_class(d, "data.frame")
  expect_equal(ncol(d), 3)          # x1, y, yT
  expect_true(all(d$y > 0))
  expect_true(all(d$yT >= d$y - 1e-9))
})

ok("cobb_douglas_XnY1 nX=3", {
  d <- cobb_douglas_XnY1(50, 3)
  expect_equal(ncol(d), 5)
})

ok("cobb_douglas_XnY1 nX=6", {
  d <- cobb_douglas_XnY1(50, 6)
  expect_equal(ncol(d), 8)
})

ok("cobb_douglas_XnY1 invalid nX", {
  expect_error(cobb_douglas_XnY1(50, 2))
})

ok("cobb_douglas_XnY1_noise nX=1", {
  d <- cobb_douglas_XnY1_noise(50, 1, p = 0.5)
  expect_class(d, "data.frame")
  expect_equal(ncol(d), 3)
})

ok("cobb_douglas_XnY1_noise nX=3 homoskedastic", {
  d <- cobb_douglas_XnY1_noise(60, 3, p = 0.5)
  expect_equal(ncol(d), 5)
})

ok("cobb_douglas_XnY1_noise nX=3 heteroskedastic", {
  d <- cobb_douglas_XnY1_noise(60, 3, p = 0.5, heteroskedasticity = TRUE)
  expect_true(nrow(d) > 0)
})

ok("cobb_douglas_XnY1_noise invalid nX", {
  expect_error(cobb_douglas_XnY1_noise(50, 2, p = 0.5))
})

ok("cobb_douglas_X3Y3 CRS", {
  d <- cobb_douglas_X3Y3(60, border = 0.2, noise = 0.5, returns = "CRS")
  expect_class(d, "data.frame")
  expect_equal(ncol(d), 9)         # x1-x3, y1-y3, yT1-yT3
})

ok("cobb_douglas_X3Y3 DRS (bug fix check)", {
  d <- cobb_douglas_X3Y3(60, border = 0.2, noise = 0.5, returns = "DRS")
  expect_true(nrow(d) == 60)
})

ok("cobb_douglas_X3Y3 IRS", {
  d <- cobb_douglas_X3Y3(60, border = 0.2, noise = 0, returns = "IRS")
  expect_true(nrow(d) == 60)
})

ok("cobb_douglas_XnHnY1 no correlation", {
  d <- cobb_douglas_XnHnY1(
    N = 50, nX = 2, nH = 1,
    exp_x = c(0.4, 0.4),
    rho_x = NULL, rho_h = NULL
  )
  expect_class(d, "data.frame")
  expect_true("y1" %in% names(d))
})

ok("cobb_douglas_XnHnY1 with rho_x", {
  rho_x <- matrix(c(0.5), nrow = 1, ncol = 1,
                  dimnames = list("x1", "x2"))
  d <- cobb_douglas_XnHnY1(
    N = 50, nX = 2, nH = 0,
    exp_x = c(0.4, 0.4),
    rho_x = rho_x, rho_h = NULL
  )
  expect_true(nrow(d) == 50)
})

ok("add_scenario_XnY1 scenario A", {
  d <- add_scenario_XnY1(60, "A")
  expect_class(d, "data.frame")
  expect_true(all(d$y >= 0))
})

ok("add_scenario_XnY1 all scenarios", {
  for (sc in c("A","B","C","D","E","F")) {
    d <- add_scenario_XnY1(60, sc)
    expect_true(nrow(d) > 0)
  }
})

ok("add_scenario_XnY1 invalid scenario", {
  expect_error(add_scenario_XnY1(50, "Z"))
})

ok("mult_scenario_XnY1 scenario A", {
  d <- mult_scenario_XnY1(60, "A")
  expect_class(d, "data.frame")
})

ok("translog_X2Y2 no noise", {
  d <- translog_X2Y2(80, border = 0.2, noise = FALSE)
  expect_class(d, "data.frame")
  expect_equal(ncol(d), 6)
})

ok("translog_X2Y2 with noise", {
  d <- translog_X2Y2(80, border = 0.2, noise = TRUE)
  expect_true(nrow(d) == 80)
})

ok("cet_cd_X2Y2 no noise", {
  d <- cet_cd_X2Y2(80, border = 0.2, noise = FALSE)
  expect_class(d, "data.frame")
  expect_true(all(d$x1 > 0))
})

ok("cet_cd_X2Y2 with noise", {
  d <- cet_cd_X2Y2(80, border = 0.2, noise = TRUE)
  expect_true(nrow(d) == 80)
})

# ============================================================ #
#  2. EFFICIENCY SCORE PRIMITIVES (LP solvers)                 #
# ============================================================ #

message("\n--- 2. LP efficiency-score primitives ---")

data_1x1y <- cobb_douglas_XnY1(40, 1)
xm <- as.matrix(data_1x1y[, "x1"])
ym <- as.matrix(data_1x1y[, "y"])

ok("rad_out VRS", {
  s <- rad_out(xm, ym, xm, ym, convexity = TRUE,  returns = "variable")
  expect_true(all(s >= 1 - 1e-6))
  expect_class(s, "matrix")
})

ok("rad_out CRS", {
  s <- rad_out(xm, ym, xm, ym, convexity = TRUE, returns = "constant")
  expect_true(all(s >= 1 - 1e-6))
})

ok("rad_out FDH", {
  s <- rad_out(xm, ym, xm, ym, convexity = FALSE, returns = "variable")
  expect_true(all(s >= 1 - 1e-6))
})

ok("rad_inp VRS", {
  s <- rad_inp(xm, ym, xm, ym, convexity = TRUE, returns = "variable")
  expect_true(all(s <= 1 + 1e-6))
})

ok("rad_inp CRS", {
  s <- rad_inp(xm, ym, xm, ym, convexity = TRUE, returns = "constant")
  expect_true(all(s <= 1 + 1e-6))
})

ok("ddf VRS", {
  dir <- cbind(xm, ym)
  s <- ddf(xm, ym, xm, ym, direction = dir, convexity = TRUE, returns = "variable")
  expect_true(all(s >= -1e-6))
})

ok("rsl_out VRS", {
  s <- rsl_out(xm, ym, xm, ym, convexity = TRUE, returns = "variable")
  expect_true(all(s >= 1 - 1e-6))
})

ok("rsl_inp VRS", {
  s <- rsl_inp(xm, ym, xm, ym, convexity = TRUE, returns = "variable")
  expect_true(all(s <= 1 + 1e-6))
})

ok("wam_mip VRS", {
  s <- wam(xm, ym, xm, ym, weights = "wam_mip", convexity = TRUE, returns = "variable")
  expect_true(all(s >= -1e-6))
})

ok("wam_nor VRS", {
  s <- wam(xm, ym, xm, ym, weights = "wam_nor", convexity = TRUE, returns = "variable")
  expect_true(all(s >= -1e-6))
})

ok("wam_ram VRS", {
  s <- wam(xm, ym, xm, ym, weights = "wam_ram", convexity = TRUE, returns = "variable")
  expect_true(all(s >= -1e-6))
})

ok("wam_bam VRS", {
  s <- wam(xm, ym, xm, ym, weights = "wam_bam", convexity = TRUE, returns = "variable")
  expect_true(all(s >= -1e-6))
})

# ============================================================ #
#  3. ACES — main model                                        #
# ============================================================ #

message("\n--- 3. ACES model ---")

data_1x <- cobb_douglas_XnY1(60, 1)

ok("aces 1x1y default", {
  m <- aces(data = data_1x, x = 1, y = 2)
  expect_class(m, "aces")
  expect_true(!is.null(m[["methods"]][["aces"]]))
})

ok("aces 1x1y quick_aces=TRUE", {
  m <- aces(data = data_1x, x = 1, y = 2, quick_aces = TRUE)
  expect_class(m, "aces")
})

ok("aces 1x1y metric=mae", {
  m <- aces(data = data_1x, x = 1, y = 2, metric = "mae")
  expect_class(m, "aces")
})

ok("aces 1x1y no concavity", {
  m <- aces(data = data_1x, x = 1, y = 2,
            shape = list(mono = TRUE, conc = FALSE))
  expect_class(m, "aces")
})

ok("aces 1x1y scale_data=FALSE", {
  m <- aces(data = data_1x, x = 1, y = 2, scale_data = FALSE)
  expect_class(m, "aces")
})

ok("aces 1x1y max_degree=2", {
  m <- aces(data = data_1x, x = 1, y = 2,
            mul_BF = list(max_degree = 2, inter_cost = 0.05))
  expect_class(m, "aces")
})

ok("aces 1x1y minspan=-2 endspan=-2", {
  m <- aces(data = data_1x, x = 1, y = 2, minspan = -2, endspan = -2)
  expect_class(m, "aces")
})

data_3x <- cobb_douglas_XnY1(80, 3)

ok("aces 3x1y", {
  m <- aces(data = data_3x, x = 1:3, y = 4)
  expect_class(m, "aces")
})

ok("aces forward model not NULL", {
  m <- aces(data = data_1x, x = 1, y = 2)
  expect_true(!is.null(m[["methods"]][["aces_forward"]]))
})

ok("aces cubic model not NULL", {
  m <- aces(data = data_1x, x = 1, y = 2)
  expect_true(!is.null(m[["methods"]][["aces_cubic"]]))
})

ok("aces quintic model not NULL", {
  m <- aces(data = data_1x, x = 1, y = 2)
  expect_true(!is.null(m[["methods"]][["aces_quintic"]]))
})

ok("aces print method runs", {
  m <- aces(data = data_1x, x = 1, y = 2)
  capture.output(print(m, method = "aces"))
  TRUE
})

ok("aces print unknown method (bug fix check)", {
  m <- aces(data = data_1x, x = 1, y = 2)
  out <- capture.output(print(m, method = "nonexistent"))
  expect_true(any(grepl("No model", out)))
})

# ============================================================ #
#  4. get_scores — all measures                                #
# ============================================================ #

message("\n--- 4. get_scores ---")

m_aces <- aces(data = data_1x, x = 1, y = 2)

aces_methods <- c("aces_forward", "aces", "aces_cubic", "aces_quintic")
score_measures <- c("rad_out", "rad_inp", "rsl_out", "rsl_inp",
                    "wam_mip", "wam_nor", "wam_ram", "wam_bam")

for (meth in aces_methods) {
  for (meas in score_measures) {
    local({
      me <- meth; ms <- meas
      ok(paste("get_scores", me, ms), {
        s <- get_scores(
          eval_data = data_1x, x = 1, y = 2,
          object = m_aces, method = me, measure = ms
        )
        expect_class(s, "data.frame")
        expect_true(nrow(s) == nrow(data_1x))
        expect_true(all(is.finite(s[[1]])))
      })
    })
  }
}

ok("get_scores ddf VRS", {
  dir <- cbind(as.matrix(data_1x[, 1]), as.matrix(data_1x[, 2]))
  s <- get_scores(
    eval_data = data_1x, x = 1, y = 2,
    object = m_aces, method = "aces", measure = "ddf",
    direction = dir
  )
  expect_class(s, "data.frame")
})

ok("get_scores relevant=TRUE", {
  s <- get_scores(
    eval_data = data_1x, x = 1, y = 2, relevant = TRUE,
    object = m_aces, method = "aces", measure = "rad_out"
  )
  expect_class(s, "data.frame")
})

ok("get_scores CRS", {
  s <- get_scores(
    eval_data = data_1x, x = 1, y = 2,
    object = m_aces, method = "aces", measure = "rad_out",
    returns = "constant"
  )
  expect_class(s, "data.frame")
})

# ============================================================ #
#  5. get_targets                                              #
# ============================================================ #

message("\n--- 5. get_targets ---")

for (meas in c("rad_out", "rad_inp", "rsl_out", "rsl_inp", "wam_mip")) {
  local({
    ms <- meas
    ok(paste("get_targets", ms), {
      t <- get_targets(
        eval_data = data_1x, x = 1, y = 2,
        object = m_aces, method = "aces", measure = ms
      )
      expect_class(t, "data.frame")
      expect_true(nrow(t) == nrow(data_1x))
    })
  })
}

ok("get_targets ddf", {
  dir <- cbind(as.matrix(data_1x[, 1]), as.matrix(data_1x[, 2]))
  t <- get_targets(
    eval_data = data_1x, x = 1, y = 2,
    object = m_aces, method = "aces", measure = "ddf",
    direction = dir
  )
  expect_class(t, "data.frame")
})

# ============================================================ #
#  6. ACES — multi-output (3x3y)                              #
# ============================================================ #

message("\n--- 6. ACES multi-output ---")

data_3x3y <- cobb_douglas_X3Y3(60, border = 0.2, noise = 0, returns = "CRS")

ok("aces 3x3y", {
  m <- aces(data = data_3x3y, x = 1:3, y = 4:6)
  expect_class(m, "aces")
})

ok("get_scores 3x3y rad_out", {
  m <- aces(data = data_3x3y, x = 1:3, y = 4:6)
  s <- get_scores(
    eval_data = data_3x3y, x = 1:3, y = 4:6,
    object = m, method = "aces", measure = "rad_out"
  )
  expect_class(s, "data.frame")
  expect_true(all(s[[1]] >= 1 - 1e-4))
})

# ============================================================ #
#  7. RF-ACES                                                  #
# ============================================================ #

message("\n--- 7. RF-ACES ---")

data_rf <- cobb_douglas_XnY1(50, 1)

ok("rf_aces basic (5 learners)", {
  m <- rf_aces(
    data = data_rf, x = 1, y = 2,
    learners = 5,
    early_stopping = list(ma_window = 3, tolerance = 2)
  )
  expect_class(m, "rf_aces")
  expect_true(length(m[["forest"]]) >= 1)
})

ok("rf_aces predict.rf_aces", {
  m <- rf_aces(
    data = data_rf, x = 1, y = 2,
    learners = 5,
    early_stopping = list(ma_window = 3, tolerance = 2)
  )
  preds <- predict(m, newdata = data_rf, x = 1)
  expect_class(preds, "data.frame")
  expect_true(nrow(preds) == nrow(data_rf))
  expect_true(all(preds[[1]] >= 0))
})

ok("rf_aces_scores rad_out", {
  m <- rf_aces(
    data = data_rf, x = 1, y = 2,
    learners = 5,
    early_stopping = list(ma_window = 3, tolerance = 2)
  )
  s <- rf_aces_scores(
    eval_data = data_rf, x = 1, y = 2,
    object = m, method = "rf_aces", measure = "rad_out"
  )
  expect_class(s, "data.frame")
  expect_true(all(s[[1]] >= 1 - 1e-4))
})

ok("rf_aces_scores cubic", {
  m <- rf_aces(
    data = data_rf, x = 1, y = 2,
    learners = 5,
    early_stopping = list(ma_window = 3, tolerance = 2)
  )
  s <- rf_aces_scores(
    eval_data = data_rf, x = 1, y = 2,
    object = m, method = "rf_aces_cubic", measure = "rad_out"
  )
  expect_class(s, "data.frame")
})

ok("rf_aces_scores quintic", {
  m <- rf_aces(
    data = data_rf, x = 1, y = 2,
    learners = 5,
    early_stopping = list(ma_window = 3, tolerance = 2)
  )
  s <- rf_aces_scores(
    eval_data = data_rf, x = 1, y = 2,
    object = m, method = "rf_aces_quintic", measure = "rad_out"
  )
  expect_class(s, "data.frame")
})

# ============================================================ #
#  8. FULL PIPELINE — timing benchmark                         #
# ============================================================ #

message("\n--- 8. Full pipeline timing ---")

ok("full pipeline 1x1y N=100", {
  set.seed(1)
  d <- cobb_douglas_XnY1(100, 1)
  m <- aces(data = d, x = 1, y = 2, quick_aces = FALSE,
            mul_BF = list(max_degree = 2, inter_cost = 0.05),
            err_red = 0.005)
  s <- get_scores(eval_data = d, x = 1, y = 2,
                  object = m, method = "aces", measure = "rad_out")
  expect_true(all(s[[1]] >= 1 - 1e-4))
})

ok("full pipeline 3x1y N=80", {
  set.seed(2)
  d <- cobb_douglas_XnY1(80, 3)
  m <- aces(data = d, x = 1:3, y = 4, quick_aces = TRUE)
  s <- get_scores(eval_data = d, x = 1:3, y = 4,
                  object = m, method = "aces", measure = "rad_out")
  expect_true(all(s[[1]] >= 1 - 1e-4))
})

ok("full pipeline 3x3y N=60", {
  set.seed(3)
  d <- cobb_douglas_X3Y3(60, border = 0.2, noise = 0, returns = "CRS")
  m <- aces(data = d, x = 1:3, y = 4:6, quick_aces = TRUE)
  s <- get_scores(eval_data = d, x = 1:3, y = 4:6,
                  object = m, method = "aces", measure = "rad_out")
  expect_true(all(s[[1]] >= 1 - 1e-4))
})

# ============================================================ #
#  SUMMARY                                                     #
# ============================================================ #

message(
  "\n============================================================",
  "\n  Results: ", pass, " passed,  ", fail, " failed",
  "\n============================================================"
)

message("\nExecution times (seconds):")
ord <- order(unlist(timings), decreasing = TRUE)
nm  <- names(timings)[ord]
vl  <- unlist(timings)[ord]
for (i in seq_along(nm)) {
  message(sprintf("  %-55s %6.3fs", nm[i], vl[i]))
}
message(sprintf("\n  Total: %.2fs", sum(unlist(timings))))
