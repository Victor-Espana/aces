devtools::load_all()
library("ggplot2")

set.seed(123)

data <- reffcy(
  DGP = "add_scenario_XnY1",
  parms = list(
    N = 200,
    scenario = "A"
  ))

x <- 1
y <- 2

library("cluster")
library("factoextra")

sil <- numeric(15)
for (k in 2:15) {
  km_res <- kmeans(data[, 1:2], centers = k, nstart = 25)
  ss <- silhouette(km_res$cluster, dist(data[, 1:2]))
  sil[k] <- mean(ss[, 3])
}

clusters <- which.max(sil)

km_res <- kmeans(data[, 1:2], centers = 2, nstart = 25)

data$cluster <- as.factor(km_res$cluster)

ACES <- vector("list", 2)

for (i in 1:2) {

  data_k <- data[data[, "cluster"] == i, ]

  ACES[[i]] <- aces (
    data = data_k,
    x = x,
    y = y,
    y_type = "individual",
    model_type = "envelopment",
    error_type = "additive",
    degree = 1,
    hd_cost = 0.5,
    metric = "mse",
    shape = list (
      monotonicity = T,
      concavity = T,
      origin = F
    ),
    nterms = 50,
    err_red = 0.005,
    minspan = - 1,
    endspan = - 1,
    kn_grid = - 1,
    d = 1,
    wc = seq(1, 2, length.out = 5),
    wq = seq(8 / 7, 1.5, length.out = 5)
  )
}

data[data[, "cluster"] == 1, "y_hat"] <- predict (
  object = ACES[[1]],
  newdata = data[data[, "cluster"] == 1, ],
  x = x,
  method = "aces"
)$y_pred

data[data[, "cluster"] == 2, "y_hat"] <- predict (
  object = ACES[[2]],
  newdata = data[data[, "cluster"] == 2, ],
  x = x,
  method = "aces"
)$y_pred

data[, "dea"] <- rad_out (
  tech_xmat = as.matrix(data[, x]),
  tech_ymat = as.matrix(data[, 5]),
  eval_xmat = as.matrix(data[, x]),
  eval_ymat = as.matrix(data[, 5]),
  convexity = TRUE,
  returns = "variable"
  )[, 1] * data$y_hat

ggplot(data) +
  geom_point(aes(x = x1, y = y, color = cluster)) +
  geom_line(aes(x = x1, y = dea), color = "red") +
  geom_line(aes(x = x1, y = yD)) +
  theme_bw()
