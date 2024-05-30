devtools::load_all()
library("ggplot2")

data <- reffcy (
  DGP = "add_scenario_XnY1",
  parms = list (
    N = 100,
    scenario = "A"
  ))

x <- 1
y <- 2

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
    "degree" = 1,
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

data$y_hat <- predict (
  object = model,
  newdata = data,
  x = 1,
  method = "aces_cubic"
)$y_pred

ggplot (data) +
  geom_point(aes(x = x1, y = y)) +
  geom_line(aes(x = x1, y = yD)) +
  geom_line(aes(x = x1, y = y_hat), color = "red", linewidth = 1) +
  theme_bw()

