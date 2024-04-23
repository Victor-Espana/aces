library("plotly")
library("latex2exp")
library("ggplot2")
library("tictoc")
library("Benchmarking")

devtools::document()
devtools::load_all()

# ============================================================================ #
#                           Evaluación escenarios                              #
# ============================================================================ #

set.seed(123)

for (i in 1:20) {

  print(i)

  data <- reffcy (
    DGP = "cobb_douglas_XnY1_noise",
    parms = list (
      N = 100,
      nX = 1,
      p = 1,
      heteroskedasticity = FALSE
    )
  )

  data$eff <- data[, "y"] / data[, "yD"]

  x <- 1
  y <- 2

  # modelo
  ACES1 <- aces (
    data = data,
    x = x,
    y = y,
    y_type = "individual",
    model_type = "stochastic",
    error_type = "mul",
    bagging = list (
      apply = FALSE,
      sample = nrow(data),
      models = 10,
      nvars = ceiling(length(x) / 3)
      ),
    degree = 1,
    hd_cost = 0.05,
    metric = "mse",
    shape = list (
      monotonicity = T,
      concavity = T,
      origin = F
    ),
    nterms = 100,
    err_red = 0,
    minspan = 1,
    endspan = 1,
    kn_grid = - 1,
    d = 1,
    wc = seq(1, 2, length.out = 5),
    wq = seq(8 / 7, 1.5, length.out = 5)
  )

  data$f_aces_moments <- predict (
    object = ACES1,
    newdata = data,
    x = x,
    method = "aces_forward",
    stochastic_pred = "moments"
  )$y_pred

  data$f_aces_pse <- predict (
    object = ACES1,
    newdata = data,
    x = x,
    method = "aces_forward",
    stochastic_pred = "pseudolikelihood"
  )$y_pred

  mean_aces <- predict (
    object = ACES1,
    newdata = data,
    x = x,
    method = "aces_forward",
    stochastic_pred = "mean"
  )$y_pred

  data$aces_moments <- predict (
    object = ACES1,
    newdata = data,
    x = x,
    method = "aces",
    stochastic_pred = "moments"
  )$y_pred

  data$aces_pse <- predict (
    object = ACES1,
    newdata = data,
    x = x,
    method = "aces",
    stochastic_pred = "pseudolikelihood"
  )$y_pred

  stoned_fitting_mom <- stoned (
    X = as.matrix(data[, x]),
    Y = as.matrix(data[, y]),
    RTS = "vrs",
    MULT = 1,
    METHOD = "MM"
  )

  mean_stoned <- stoned_fitting_mom[["yhat"]]

  stoned_fitting_pse <- stoned (
    X = as.matrix(data[, x]),
    Y = as.matrix(data[, y]),
    RTS = "vrs",
    MULT = 1,
    METHOD = "PSL"
  )

  data$stoned_mom <- stoned_fitting_mom[["front"]][, 1]
  data$stoned_pse <- stoned_fitting_pse[["front"]][, 1]

  # mean error aces
  print(paste("Mean err ACES:", mean((mean_aces - data[, y]) ^ 2)))

  # mean error stoned
  print(paste("Mean err StoNED:", mean((mean_stoned - data[, y]) ^ 2)))

  err_f_aces_mom <- mean((data[, 3] - data[, 5]) ^ 2)
  err_f_aces_pse <- mean((data[, 3] - data[, 6]) ^ 2)
  err_b_aces_mom <- mean((data[, 3] - data[, 7]) ^ 2)
  err_b_aces_pse <- mean((data[, 3] - data[, 8]) ^ 2)
  err_stoned_mom <- mean((data[, 3] - data[, 9]) ^ 2)
  err_stoned_pse <- mean((data[, 3] - data[, 10]) ^ 2)

  print (
    paste (
      "aces_f_mom:", round(err_f_aces_mom, 4),
      "aces_f_pse:", round(err_f_aces_pse, 4),
      "aces_mom:", round(err_b_aces_mom, 4),
      "aces_pse:", round(err_b_aces_pse, 4),
      "stoned_mom:", round(err_stoned_mom, 4),
      "stoned_pse:", round(err_stoned_pse, 4)
    )
  )
}

data$aces_mean <- predict (
  object = ACES1,
  newdata = data,
  x = x,
  method = "aces",
  stochastic_pred = "moments"
)$y_pred

ggplot(data) +
  # geom_point(aes(x = x1, y = y)) +
  geom_line(aes(x = 1:100, y = aces_mean), color = "red") +
  geom_point(aes(x = 1:100, y = aces_mean), color = "red") +
  # geom_line(aes(x = 1:100, y = stoned_mom), color = "blue") +
  # geom_point(aes(x = 1:100, y = stoned_mom), color = "blue") +
  geom_line(aes(x = 1:100, y = yD), color = "green") +
  geom_point(aes(x = 1:100, y = yD), color = "green") +
  theme_bw()

ggplot(data) +
  geom_point(aes(x = x1, y = y)) +
  geom_line(aes(x = x1, y = yD, color = "DGP"), linewidth = 1) +
  geom_line(aes(x = x1, y = aces_mean, color = "aces_mom"), linewidth = 1) +
  # geom_line(aes(x = x1, y = stoned, color = "stoned"), linewidth = 1) +
  theme_bw() +
  guides(color = guide_legend(title = "model")) +
  theme(legend.position = c(0.8, 0.2))

scores <- aces_scores (
  tech_data = data,
  eval_data = data,
  x = x,
  y = y,
  object = ACES1,
  method = "aces",
  proximity = ifelse(length(y) == 1, Inf, 0.05),
  measure = "rad_out",
  convexity = TRUE,
  returns = "variable",
  direction = NULL,
  weights = NULL,
  digits = 3
)

data$y_hat1 <- scores$dea_vrt_rad_out * data[, 2]

# ============================================================================ #
#                                 Concavidad                                   #
# ============================================================================ #

# la concavidad se rompe si queda un intervalo vacío distinto al último por la derecha

# ejemplo 1: primer intervalo
set.seed(123)

data <- reffcy(
  DGP = "add_scenario_XnY1",
  parms = list(
    N = 50,
    scenario = "A"
  ))

a1 <- 0.8
a2 <- - 0.4

k1 <- 3.084632
k2 <- 5.080007

hinge1 <- pmax(0, data[, 1] - k1)
hinge2 <- pmax(0, data[, 1] - k2)

data$y_hat <- 0.5 + a1 * hinge1 + a2 * hinge2

ggplot(data) +
  geom_line(aes(x = x1, y = y_hat)) +
  geom_vline(xintercept = k1) +
  geom_vline(xintercept = k2) +
  theme_bw()

# ejemplo 2: segundo intervalo
set.seed(123)

data <- reffcy(
  DGP = "add_scenario_XnY1",
  parms = list(
    N = 50,
    scenario = "A"
  ))

b1 <- - 0.8
a2 <- 0.4

k1 <- 3.084632
k2 <- 5.080007

hinge1 <- pmax(0, k1 - data[, 1])
hinge2 <- pmax(0, data[, 1] - k2)

data$y_hat <- 0.5 + b1 * hinge1 + a2 * hinge2

ggplot(data) +
  geom_line(aes(x = x1, y = y_hat)) +
  geom_vline(xintercept = k1) +
  geom_vline(xintercept = k2) +
  theme_bw()

# ============================================================================ #
#                                   X2 ~ Y2                                    #
# ============================================================================ #

data <- reffcy (
  DGP = "translog_X2Y2",
  parms = list (
    N = 100,
    border = 0.1,
    noise = FALSE
  ))

x <- 1:2
y <- 3:4

# modelo
ACES <- aces (
  data = data,
  x = x,
  y = y,
  y_type = "all",
  model_type = "env",
  error_type = "add",
  bagging = list (
    apply = FALSE,
    sample = nrow(data),
    models = 50,
    nvars = ceiling(length(x) / 3),
    oob_red = 0.01,
  ),
  degree = 2,
  hd_cost = 0,
  metric = "mse",
  shape = list (
    monotonicity = T,
    concavity = T,
    origin = F
  ),
  nterms = nrow(data),
  err_red = 0,
  minspan = - 1,
  endspan = - 1,
  kn_grid = - 1,
  d = 1,
  wc = seq(1, 2, length.out = 5),
  wq = seq(8 / 7, 1.5, length.out = 5)
)

data[, c("y1_pred2", "y2_pred2")] <- predict (
  object = ACES,
  newdata = data,
  x = x,
  method = "aces_forward"
)

plot_ly() %>%
  add_trace (
    data = data, x = data$x1, y = data$x2, z = data$y1,
    type = "mesh3d"
  ) %>%
  add_trace (
    data = data, x = data$x1, y = data$x2, z = data$y1_pred,
    type = "mesh3d"
  ) %>%
  add_markers (
    data = data, x = data$x1, y = data$x2, z = data$y1,
    marker = list(size = 5, color = "red", opacity = 0.8)
  )

# ============================================================================ #
#                                 RandomForest                                 #
# ============================================================================ #

devtools::load_all()

data <- reffcy (
  DGP = "translog_X2Y2",
  parms = list (
    N = 25,
    border = 0.1,
    noise = FALSE
  ))

x <- 1:2
y <- 3:4

ACES <- aces (
  data = data,
  x = x,
  y = y,
  y_type = "all",
  model_type = "env",
  error_type = "add",
  RF = list (
    "apply" = TRUE,
    "sample" = nrow(data),
    "models" = 200,
    "nvars" = 1,
    "oob_red" = 0.001
  ),
  mul_BF = list (
    "degree" = 1,
    "hd_cost" = 0
  ),
  metric = "mse",
  shape = list (
    "mon" = T,
    "con" = T,
    "ori" = F
  ),
  nterms = nrow(data),
  err_red = 0.005,
  kn_grid = - 1,
  minspan = - 1,
  endspan = 0,
  kn_penalty = 1,
  smoothing = list (
    "wc" = seq(1, 2, length.out = 5),
    "wq" = seq(8 / 7, 1.5, length.out = 5)
  )
)

data[, c("y1_pred2RF", "y2_pred2RF")] <- predict (
  object = ACES,
  newdata = data,
  x = x,
  method = "aces_forward"
)

compute_variable_importance (
  data = data,
  x = x,
  y = y,
  object = ACES,
  repeats = 1
)

# OOB
oob <- c()
for (j in 1:60) {
  oob <- c(oob, ACES[[1]][[j]][[1]][["OOB"]])
}

oob_data <- data.frame (
  model = 1:60,
  oob_err = oob
)

ggplot(oob_data) +
  geom_point(aes(x = model, y = oob_err), color = "#FFABAB") +
  geom_line(aes(x = model, y = oob_err), color = "#85E3FF") +
  theme_bw() +
  theme (
    axis.title.x = element_text(
      size = 12, face = "bold", color = "#921F30",
      margin = margin(t = 10)),
    axis.title.y = element_text(
      size = 12, face = "bold", color = "#921F30",
      margin = margin(r = 10)),
    axis.text = element_text(
      size = 12, color = "black"),
    plot.margin = unit(c(1.25, 1.25, 1.25, 1.25), "lines"),
    plot.title = element_text(
      size = 12, face = "bold", color = "#921F30",
      margin = margin(b = 10)
    )
  )
