library("plotly")
library("latex2exp")
library("ggplot2")

devtools::document()
devtools::load_all()

# ===================== #
# Evaluación escenarios #
# ===================== #

set.seed(314)     # error: 0.019

data <- reffcy (
  DGP = "translog_X2Y2",
  parms = list (
    border = 0.20,
    N = 100,
    noise = FALSE
  )
)

data[, "x1x2"] <- data[, 1] * data[, 2]

# data <- reffcy (
#   DGP = "add_scenario_XnY1",
#   parms = list (
#     scenario = "F",
#     N = 100
#   )
# )

library("psych")

pairs.panels(data)

x <- 1:2
y <- 3:4

# modelo
ACES <- aces (
  data = data,
  x = x,
  y = y,
  y_type = "individual",
  addi_x = NULL,
  degree = 2,
  metric = "mse",
  turbo = Inf,
  monotonicity = T,
  concavity = T,
  x0_y0 = T,
  nterms = 50,
  err_red = 0.01,
  minspan = - 1,
  endspan = - 1,
  knots_grid = - 1,
  d = 2,
  wc = seq(1, 2, length.out = 5),
  wq = seq(8 / 7, 1.5, length.out = 5)
)

y_hat <- predict (
  object = ACES,
  newdata = data,
  x = c(x),
  method = "aces_forward"
)[, c("y1_pred", "y2_pred")]

pruebas <- data.frame (
  y_hat1 = y_hat[, 1],
  y_hat2 = y_hat[, 2],
  y1 = data[, 3],
  y2 = data[, 4],
  yD1 = data[, 5],
  yD2 = data[, 6],
  esti_dif1 = round(y_hat[, 1] - data[, 3], 3),
  estm_dif2 = round(y_hat[, 2] - data[, 4], 3),
  real_dif1 = round(y_hat[, 1] - data[, 5], 3),
  real_dif2 = round(y_hat[, 2] - data[, 6], 3)
)

ratio_data <- data.frame (
  r_yhat1_x1 = y_hat[, 1] / data[, 1],
  r_yhat1_x2 = y_hat[, 1] / data[, 2],
  r_yhat2_x1 = y_hat[, 2] / data[, 1],
  r_yhat2_x2 = y_hat[, 2] / data[, 2],
  r_y1_x1 = data[, 3] / data[, 1],
  r_y1_x2 = data[, 3] / data[, 2],
  r_y2_x1 = data[, 4] / data[, 1],
  r_y2_x2 = data[, 4] / data[, 2]
)

library("MASS")
library("ggdist")

mds1 <- data.frame(cmdscale(dist(ratio_data)))

ggplot(mds1, aes(x = mds1[, 1], y = mds1[, 2])) +
  geom_point() +
  geom_text(aes(label = row.names(mds1)), vjust = -0.5) +
  labs(title = "Co-plot - Datos 1", x = "Dimensión 1", y = "Dimensión 2")



# predicción
devtools::load_all()

scores <- aces_scores (
  tech_data = data,
  eval_data = data,
  x = x,
  y = y,
  addi_x = 7,
  object = ACES,
  method = "aces",
  proximity = Inf,
  measure = "rad_out",
  convexity = TRUE,
  returns = "variable",
  direction = NULL,
  weights = NULL,
  digits = 10
)







ggplot(data) +
  geom_point(aes(x = x1, y = y)) +
  geom_line(aes(x = x1, y = yD, color = "Data Generation Process")) +
  geom_line(aes(x = x1, y = y_hat, color = "aces")) +
  theme_bw() +
  guides(color = guide_legend(title = "model")) +
  theme(legend.position = c(0.8, 0.2)) +
  expand_limits(x = 0, y = 0)

# ============================= #
# Unidades cerca de la frontera #
# ============================= #

x <- 1:2
y <- 3:4

data$scores <- ddf (
  tech_xmat = as.matrix(data[, x]),
  tech_ymat = as.matrix(data[, y]),
  eval_xmat = as.matrix(data[, x]),
  eval_ymat = as.matrix(data[, y]),
  direction = "briec",
  convexity = TRUE,
  returns = "variable"
)[, 1]

data$eff <- factor(ifelse(data$scores < 0.05, 1, 0))

ggplot(data) +
  geom_line(aes(x = x1, y = yD)) +
  geom_point(aes(x = x1, y = y, color = eff)) +
  theme_bw()

# ========== #
# Concavidad #
# ========== #

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







  geom_line(aes(x = x1, y = DEA.OUT, color = "Best Practices Frontier")) +
  # geom_point(aes(x = x1, y = DEA.OUT, color = "Output projection")) +
  # geom_point(aes(x = DEA.INP, y = y, color = "Input projection")) +
  geom_segment(aes(x = 3.678058, y = 5.615856 + 0.05, xend = 3.678058, yend = 5.615856 * 1.094396 - 0.05),
               color = "#c77cff", size = 1,
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(aes(x = 3.678058 + 0.35, y = 5.615856 + 0.25, label = TeX("$\\phi = 1.09$", output = "character")),
            parse = TRUE, size = 4) +
  geom_segment(aes(x = 3.678058 - 0.05, y = 5.615856, xend = 3.678058 * 0.7936781 + 0.05, yend = 5.615856),
               color = "#00bfc4", size = 1,
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(aes(x = 3.678058 - 0.35, y = 5.615856 - 0.15, label = TeX("$\\theta = 0.79$", output = "character")),
            parse = TRUE, size = 4) +








pred <- predict(model[[1]], data, x, y1, 4)

dmu  <- nrow(data)
xmat <- as.matrix(data[, x])
ymat <- as.matrix(data[, y2])
nX   <- length(x)
nY   <- length(y2)

data$BCC.OUT <- (AAFS_BCC.OUT(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu))[, 1]
data$BCC.INP <- (AAFS_BCC.INP(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu))[, 1]
data$DEA.OUT <- (AAFS_BCC.OUT(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu) * data[, y2])[, 1]
data$DEA.INP <- (AAFS_BCC.INP(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu) * data[, x])[, 1]

ggplot(data) +
  geom_line(aes(x = x1, y = yD, color = "Data Generation Process")) +
  geom_point(aes(x = x1, y = y)) +
  geom_line(aes(x = x1, y = DEA.OUT, color = "Best Practices Frontier")) +
  # geom_point(aes(x = x1, y = DEA.OUT, color = "Output projection")) +
  # geom_point(aes(x = DEA.INP, y = y, color = "Input projection")) +
  geom_segment(aes(x = 3.678058, y = 5.615856 + 0.05, xend = 3.678058, yend = 5.615856 * 1.094396 - 0.05),
               color = "#c77cff", size = 1,
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(aes(x = 3.678058 + 0.35, y = 5.615856 + 0.25, label = TeX("$\\phi = 1.09$", output = "character")),
            parse = TRUE, size = 4) +
  geom_segment(aes(x = 3.678058 - 0.05, y = 5.615856, xend = 3.678058 * 0.7936781 + 0.05, yend = 5.615856),
               color = "#00bfc4", size = 1,
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(aes(x = 3.678058 - 0.35, y = 5.615856 - 0.15, label = TeX("$\\theta = 0.79$", output = "character")),
            parse = TRUE, size = 4) +
  theme_bw()


data[, c("By1", "By2")] <- predict(model, data, x, NULL, 2)[, c("y1_pred", "y2_pred")]
data[, c("Cy1", "Cy2")] <- predict(model, data, x, NULL, 3)[, c("y1_pred", "y2_pred")]
data[, c("Qy1", "Qy2")] <- predict(model, data, x, NULL, 4)[, c("y1_pred", "y2_pred")]

# Hinge function + logarithm
set.seed(314)
data <- X2Y2.sim(1000, 1)

t1 <- 37
data$BF1 <- pmax(0, t1 - data$x1)
data$BF2 <- log(pmax(1, data$x1 - t1 + 1))

t2 <- 21
data$BF3 <- pmax(0, t2 - data$x1)
data$BF4 <- log(pmax(1, data$x1 - t2 + 1))

t3 <- 16
data$BF5 <- pmax(0, t3 - data$x1)
data$BF6 <- log(pmax(1, data$x1 - t3 + 1))

data$ypred <- 0.8 - 2 * data$BF1 + 1 * data$BF2 - 3 * data$BF3 + 2 * data$BF4 - 1.8 * data$BF5 + 1.2 * data$BF6

ggplot(data) +
  geom_line(aes(x = x1, y = ypred, color = "Prediction"))

for (nx in 1:3) {
  for (N in c(25, 50, 75, 100, 150)) {
    print("Minimum span")
    print(paste("nX: ", nx, "N: ", N, "log: ", - log2((- 1 / (nx * N)) * log(1 - 0.05)) / 2.5))

    print("End span")
    print(paste("nX: ", nx, "N: ", N, "log: ", 3 - log2(0.05 / nx)))
  }
}
