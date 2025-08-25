devtools::load_all()
devtools::document()
library("ggplot2")

set.seed(314)

data <- cobb_douglas_XnY1 (
  N = 50,
  nX = 3
)

x <- 1:3
y <- 4

# ACES
model <- aces (
  data = data, #
  x = x, #
  y = y, #
  quick_aces = FALSE,
  error_type = "add",
  mul_BF = list ( #
    "max_degree" = 1,
    "inter_cost" = 0.05
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
  ),
  max_terms = nrow(data), #
  err_red = 0.0001, #
  kn_grid = - 1,
  minspan = - 1,
  endspan = - 1
)

y_hat <- predict (
  object = model,
  newdata = rbind(data, rep(0, 5)),
  x = x,
  method = "aces"
)$y_pred

round(tail(y_hat), 3)

aces_scores (
  eval_data = data,
  x = x,
  y = y,
  relevant = TRUE,
  object = model,
  method = "aces",
  measure = "wam",
  returns = "variable",
  direction = NULL,
  weights = "NOR",
  digits = 3
)

# RF-ACES

model <- rf_aces (
  data = data,
  x = x,
  y = y,
  quick_aces = FALSE,
  error_type = "add",
  mul_BF = list (
    "max_degree" = 2,
    "inter_cost" = 0.05
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = T,
    "ptto" = F
  ),
  learners = 100,
  bag_size = nrow(data),
  max_feats = length(x) / 3,
  early_stopping = list (
    "ma_window" = 20,
    "tolerance" = 10
  ),
  max_terms = nrow(data),
  err_red = 0.01,
  kn_grid = - 1,
  minspan = - 1,
  endspan = 0
)

importance <- rf_varimp (
    data = data,
    x = x,
    y = y,
    object = model,
    repeats = 3
    )

p_predict2 <- rf_aces_p_predict (
  object = model,
  newdata = data,
  x = x,
  y = y,
  p = 1,
  method = "rf_aces"
)

scores <- rf_aces_scores (
  eval_data = data[, c(x, y)],
  x = x,
  y = y,
  object = model,
  method = "rf_aces",
  measure = "rf_aces_rad_out",
  returns = "variable",
  direction = NULL,
  weights = NULL,
  digits = 3
)

toc()

data$y_hat <- predict (
  object = model,
  newdata = data,
  x = x,
  method = "rf_aces"
)$y_pred

aces_scores (
  eval_data = data,
  x = x,
  y = y,
  relevant = FALSE,
  object = model,
  method = "aces",
  measure = "wam",
  returns = "variable",
  direction = NULL,
  weights = "NOR",
  digits = 3
)

# rf_scores
# interval rf
