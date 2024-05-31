devtools::load_all()
library("ggplot2")

data <- reffcy (
  DGP = "translog_X2Y2",
  parms = list (
    N = 100,
    border = 0.1,
    noise = FALSE
  ))

x <- 1:2
y <- 3:4

model <- aces (
  data = data,
  x = x,
  y = y,
  y_type = "all",
  model_type = "env",
  error_type = "add",
  RF = list (
    "apply" = FALSE,
    "sample" = nrow(data),
    "models" = 10,
    "nvars" = length(x),
    "oob_red" = 0.001
  ),
  mul_BF = list (
    "degree" = 2,
    "hd_cost" = 0.05
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
  ),
  nterms = nrow(data),
  err_red = 0.005,
  kn_grid = - 1,
  minspan = - 1,
  endspan = - 1,
  kn_penalty = 2,
  smoothing = list(
    "wc" = seq(1, 2, length.out = 5),
    "wq" = seq(8 / 7, 1.5, length.out = 5)
  )
)

data$score <- data[, 5] / data[, 3]

data$aces <- aces_scores (
  tech_data = data,
  eval_data = data,
  x = x,
  y = y,
  object = model,
  method = "aces",
  convexity = TRUE,
  returns = "variable",
  direction = NULL,
  weights = NULL,
  digits = 3
  )$aces_vrt_rad_out

data$dea <- rad_out (
  tech_xmat = as.matrix(data[, x]),
  tech_ymat = as.matrix(data[, y]),
  eval_xmat = as.matrix(data[, x]),
  eval_ymat = as.matrix(data[, y]),
  convexity = TRUE,
  returns = "variable"
)[, 1]

sum((data[, 8] - data[, 7]) ^ 2) / 100
sum((data[, 9] - data[, 7]) ^ 2) / 100

ggplot (data) +
  geom_point(aes(x = x1, y = y)) +
  geom_line(aes(x = x1, y = yD)) +
  geom_line(aes(x = x1, y = y_hat), color = "red", linewidth = 1) +
  theme_bw()

