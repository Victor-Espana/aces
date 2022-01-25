#' @title Single Output Data Generation
#'
#' @description This function is used to simulate the data in a single output scenario.
#'
#' @param N Sample size.
#' @param nX Number of inputs. Possible values: \code{1}, \code{3}, \code{6}, \code{9}, \code{12} and \code{15}.
#'
#' @importFrom dplyr %>%
#' @importFrom stats runif rnorm
#'
#' @return \code{data.frame} with the simulated data.
#'
#' @export
Y1.sim <- function(N, nX) {
  if(!(nX %in% c(1, 3, 6, 9, 12, 15))){
    stop(paste(nX, "is not allowed"))
  }

  colnames <- c(paste("x", 1:nX, sep = ""), "y")

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  for (x in 1:nX){
    data[, x] <- runif(n = N, min = 1, max = 10)
  }

  u <- abs(rnorm(n = N, mean = 0, sd = 0.4))

  if (nX == 1){
    y <- 3 * log(data[, "x1"])
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y

  } else if (nX == 3){
    y <- 3 * (data[, "x1"] ** 0.05) * (data[, "x2"] ** 0.15) * (data[, "x3"] ** 0.3)
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y

  } else if (nX == 6){
    y <- 3 * (data[, "x1"] ** 0.05) * (data[, "x2"] ** 0.001) * (data[, "x3"] ** 0.004) *
      (data[, "x4"] ** 0.045) * (data[, "x5"] ** 0.1) * (data[, "x6"] ** 0.3)
    data[, "y"] <- y * exp(-u)
    data[, "yD"] <- y

  } else if (nX == 9){
    y <- 3 * (data[, "x1"] ** 0.005) * (data[, "x2"] ** 0.001) * (data[, "x3"] ** 0.004) *
      (data[, "x4"] ** 0.005) * (data[, "x5"] ** 0.001) * (data[, "x6"] ** 0.004) *
      (data[, "x7"] ** 0.08) * (data[, "x8"] ** 0.1) * (data[, "x9"] ** 0.3)
    data["y"] <- y * exp(-u)
    data["yD"] <- y

  } else if (nX == 12){
    y <- 3 * (data[, "x1"] ** 0.005) * (data[, "x2"] ** 0.001) * (data[, "x3"] ** 0.004) *
      (data[, "x4"] ** 0.005) * (data[, "x5"] ** 0.001) * (data[, "x6"] ** 0.004) *
      (data[, "x7"] ** 0.08) * (data[, "x8"] ** 0.05) * (data[, "x9"] ** 0.05) *
      (data[, "x10"] ** 0.075) * (data[, "x11"] ** 0.025) * (data[, "x12"] ** 0.2)
    data["y"] <- y * exp(-u)
    data["yD"] <- y

  } else {
    y <- 3 * (data[, "x1"] ** 0.005) * (data[, "x2"] ** 0.001) * (data[, "x3"] ** 0.004) *
      (data[, "x4"] ** 0.005) * (data[, "x5"] ** 0.001) * (data[, "x6"] ** 0.004) *
      (data[, "x7"] ** 0.08) * (data[, "x8"] ** 0.05) * (data[, "x9"] ** 0.05) *
      (data[, "x10"] ** 0.05) * (data[, "x11"] ** 0.025) * (data[, "x12"] ** 0.025) *
      (data[, "x13"] ** 0.025) * (data[, "x14"] ** 0.025) * (data[, "x15"] ** 0.15)
    data["y"] <- y * exp(-u)
    data["yD"] = y
  }

  return(data)
}

#' @title 2 Inputs & 2 Outputs Data Generation
#'
#' @description This function is used to simulate the data in a scenario with 2 inputs and 2 outputs.
#'
#' @param N Sample size.
#' @param border Percentage of DMUs in the frontier.
#' @param noise Random noise.
#'
#' @importFrom dplyr %>%
#' @importFrom stats runif rnorm
#'
#' @return \code{data.frame} with the simulated data.
#'
#' @export
X2Y2.sim <- function(N, border, noise = NULL) {
  nX <- 2
  nY <- 2

  colnames <- c(paste("x", 1:nX, sep = ""), paste("y", 1:nY, sep = ""))

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  data[, 1:nX] <- runif(N, 5, 50)

  z <- runif(N, -1.5, 1.5)

  ln_x1 <- log(data[, "x1"])
  ln_x2 <- log(data[, "x2"])

  op1 <- -1 + 0.5 * z + 0.25 * (z**2) - 1.5 * ln_x1

  op2 <- -0.6 * ln_x2 + 0.2 * (ln_x1**2) + 0.05 * (ln_x2**2) - 0.1 * ln_x1 * ln_x2

  op3 <- 0.05 * ln_x1 * z - 0.05 * ln_x2 * z

  ln_y1_ast <- -(op1 + op2 + op3)

  data[, "y1"] <- exp(ln_y1_ast)

  data[, "y2"] <- exp(ln_y1_ast + z)

  if (border > 0) {
    index <- sample(1:N, N * border)

    N_sample <- length(index)

    half_normal <- rnorm(N_sample, 0, 0.3**(1 / 2)) %>%
      abs()

    half_normal <- exp(half_normal)

    if (!is.null(noise)) {
      normal1 <- rnorm(N_sample, 0, 0.01**(1 / 2))
      normal1 <- exp(normal1)

      normal2 <- rnorm(N_sample, 0, 0.01**(1 / 2))
      normal2 <- exp(normal2)

      data[index, "y1"] <- data[index, "y1"] / (half_normal * normal1)
      data[index, "y2"] <- data[index, "y2"] / (half_normal * normal2)
    } else {
      data[index, "y1"] <- data[index, "y1"] / half_normal
      data[index, "y2"] <- data[index, "y2"] / half_normal
    }
  }

  return(data)
}

#' @title Additive Scenarios for 1 Output
#'
#' @description This function is used to simulate the data in different additive scenarios as in \insertCite{kuosmanen2010;textual}{mafs}.
#'
#' @param N Sample size.
#' @param scenario \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, \code{"E"} or \code{"F"}. For details, check Table 2.
#'
#' @importFrom dplyr %>%
#' @importFrom stats runif rnorm
#' @importFrom Rdpack reprompt
#'
#' @return \code{data.frame} with the simulated data.
#'
#' \insertRef{kuosmanen2010}{mafs}
#'
#' @export
AddScenario <- function(N, scenario) {
  if(!(scenario %in% c("A", "B", "C", "D", "E", "F"))){
    stop(paste(scenario, "is not allowed"))
  }

  if (scenario %in% c("A", "B")) {
    nX <- 1

  } else if (scenario %in% c("C", "E")) {
    nX <- 2

  } else {
    nX <- 3
  }

  colnames <- c(paste("x", 1:nX, sep = ""), "y")

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  for (x in 1:nX){
    data[, x] <- runif(n = N, min = 1, max = 10)
  }

  u <- abs(rnorm(n = N, mean = 0, sd = 0.4))

  if (scenario == "A"){
    x1 <- data[, "x1"]
    y <- log(x1) + 3
    data[, "y"] <- y - u
    data[, "yD"] <- y

  } else if (scenario == "B"){
    x1 <- data[, "x1"]
    y <- 3 + sqrt(x1) + log(x1)
    data[, "y"] <- y - u
    data[, "yD"] <- y

  } else if (scenario == "C"){
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.3 * sqrt(x1 * x2)
    data[, "y"] <- y - u
    data[, "yD"] <- y

  } else if (scenario == "D"){
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 3)
    data["y"] <- y - u
    data["yD"] <- y

  } else if (scenario == "E"){
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.3 * (x1 * x2) ^ (1 / 3)
    data["y"] <- y - u
    data["yD"] <- y

  } else {
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 4)
    data["y"] <- y - u
    data["yD"] <- y
  }

  return(data)
}

#' @title Multiplicative Scenarios for 1 Output
#'
#' @description This function is used to simulate the data in different multiplicative scenarios as in \insertCite{kuosmanen2010;textual}{mafs}.
#'
#' @param N Sample size.
#' @param scenario \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, \code{"E"} or \code{"F"}. For details, check Table 4.
#'
#' @importFrom dplyr %>%
#' @importFrom stats runif rnorm
#' @importFrom Rdpack reprompt
#'
#' @return \code{data.frame} with the simulated data.
#'
#' \insertRef{kuosmanen2010}{mafs}
#'
#' @export
MultScenario <- function(N, scenario) {
  if(!(scenario %in% c("A", "B", "C", "D", "E", "F"))){
    stop(paste(scenario, "is not allowed"))
  }

  if (scenario %in% c("A", "B")) {
    nX <- 1

  } else if (scenario %in% c("C", "E")) {
    nX <- 2

  } else {
    nX <- 3
  }

  colnames <- c(paste("x", 1:nX, sep = ""), "y")

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  for (x in 1:nX){
    data[, x] <- runif(n = N, min = 1, max = 10)
  }

  u <- abs(rnorm(n = N, mean = 0, sd = 0.4))

  if (scenario == "A"){
    x1 <- data[, "x1"]
    y <- log(x1) + 3
    data[, "y"] <- y / (1 + u)
    data[, "yD"] <- y

  } else if (scenario == "B"){
    x1 <- data[, "x1"]
    y <- 3 + sqrt(x1) + log(x1)
    data[, "y"] <- y / (1 + u)
    data[, "yD"] <- y

  } else if (scenario == "C"){
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.3 * sqrt(x1 * x2)
    data[, "y"] <- y / (1 + u)
    data[, "yD"] <- y

  } else if (scenario == "D"){
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 3)
    data["y"] <- y / (1 + u)
    data["yD"] <- y

  } else if (scenario == "E"){
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.3 * (x1 * x2) ^ (1 / 3)
    data["y"] <- y / (1 + u)
    data["yD"] <- y

  } else {
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]
    y <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 4)
    data["y"] <- y / (1 + u)
    data["yD"] <- y
  }

  return(data)
}


#' @title 3 Inputs & 3 Outputs Data Generation
#'
#' @description This function is used to simulate the data in a scenario with 3 inputs and 3 outputs.
#'
#' @param N Sample size.
#' @param noise Random noise. 4 possible levels: \code{0}, \code{0.5}, \code{1} or \code{2}.
#' @param inefficiency Percentage of inneficient DMUs: \code{0.9}, \code{0.8} or \code{0.7}.
#' @param return_scale Returns to scale. \code{"CRS"} for Constant Return to Scale, \code{"DRS"} for Decreasing Return to Scale and \code{"IRS"} for Increasing Return to Scale.
#'
#' @importFrom stats runif rnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom dplyr %>%
#'
#' @return \code{data.frame} with the simulated data.
#'
#' @export
X3Y3.sim <- function(N, noise, inefficiency, return_scale) {
  nX <- 3
  nY <- 3

  colnames <- c(paste("x", 1:nX, sep = ""), paste("y", 1:nY, sep = ""), paste("yD", 1:nY, sep = ""))

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # Inputs
  data[, 1:nX] <- runif(N, 5, 15)

  b0 <- 0

  if (return_scale == "CRS") {
    b1 <- b2 <- b3 <- 1 / 3

  } else if (return_scale == "DRS") {
    b1 <- b2 <- b3 <- 0.267

  } else {
    b1 <- b2 <- b3 <- 0.4

  }

  Yd_aux <- exp(b0) * data[, 1] ^ b1 * data[, 2] ^ b2 * data[, 3] ^ b3

  for (i in 1:nY) {
    alpha <- rtruncnorm(N, 0, 1, 1 / nY, 1 / nY ^ 2)
    a <- alpha / sum(alpha)

    # Yd
    data[, nX + nY + i] <- sqrt(a) * Yd_aux

    # Yobs
    # Random Noise
    v <- rnorm(N, 0, (noise * 0.136) ^ 2)

    # Inneficiency
    if (inefficiency == 0.9) {
      mu <- rtruncnorm(N, 0, Inf, 0, 0.136 ^ 2)

    } else if (inefficiency == 0.8) {
      mu <- rtruncnorm(N, 0, Inf, 0, 0.299 ^ 2)

    } else if (inefficiency == 0.7) {
      mu <- rtruncnorm(N, 0, Inf, 0, 0.488 ^ 2)

    }

    Yobs_aux <- Yd_aux * exp(- mu) * exp(v)

    data[, nY + i] <- sqrt(a) * Yobs_aux

  }
  return(data)
}
