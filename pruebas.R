devtools::load_all()
devtools::document()
library("ggplot2")

data <- cobb_douglas_XnY1 (
  N = 100,
  nX = 3
)

x <- 1:3
y <- 4

# ACES
model <- aces (
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
  max_terms = nrow(data),
  err_red = 0.01,
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
  relevant = FALSE,
  object = model,
  method = "aces",
  measure = "wam",
  returns = "variable",
  direction = NULL,
  weights = "NOR",
  digits = 3
)

# RF-ACES
tic()
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
  learners = 50,
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
toc()

data$y_hat <- predict (
  object = model,
  newdata = data,
  x = x,
  method = "rf_aces"
)$y_pred

# aces_scores
# rf_scores
# interval rf

# S-ACES
model <- s_aces (
  data = data,
  x = x,
  y = y,
  z = z,
  quick_aces = TRUE,
  error_type = "add",
  mul_BF = list (
    "max_degree" = 1,
    "compl_cost" = 0.05
  ),
  metric = "mse",
  shape = list (
    "mono" = T,
    "conc" = F,
    "ptto" = F
  ),
  nterms = nrow(data),
  err_red = 0.01,
  kn_grid = - 1,
  minspan = - 1,
  endspan = 0,
  kn_penalty = 1
)

data$aces <- predict (
  object = model,
  newdata = data,
  x = x,
  method = "aces"
)$y_pred

ggplot(data) +
  geom_point(aes(x = x1, y = y)) +
  geom_line(aes(x = x1, y = aces), color = "red") +
  theme_bw()
