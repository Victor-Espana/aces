#' @title 1 output ~ nX inputs Cobb-Douglas Data Generation Process
#'
#' @description
#' This function simulates a 1 output ~ nX inputs \code{data.frame} with a Cobb-Douglas Data Generation Process.
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param nX
#' An integer representing the number of inputs. Possible values: \code{1}, \code{3}, \code{6}, \code{9}, \code{12} and \code{15}.
#'
#' @importFrom dplyr %>%
#' @importFrom stats runif rnorm
#'
#' @return
#' A \code{data.frame} with the simulated data: nX inputs, 1 output (y) and the theoretical frontier (yT).
#'
#' @export

cobb_douglas_XnY1 <- function (
    N,
    nX
    ) {

  if(!(nX %in% c(1, 3, 6, 9, 12, 15))){
    stop(paste(nX, "is not allowed"))
  }

  colnames <- c(paste("x", 1:nX, sep = ""), "y")

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # input generation
  for (x in 1:nX){
    data[, x] <- runif(n = N, min = 0, max = 10)
  }

  # inefficiency generation
  u <- abs(rnorm(n = N, mean = 0, sd = 0.4))

  if (nX == 1) {

    # theoretical frontier
    data[, "yT"] <- 3 * data[, "x1"] ** 0.5

    # output generation
    data[, "y"]  <- data[, "yT"] * exp(- u)

  } else if (nX == 3) {

    # theoretical frontier
    data[, "yT"] <- 3 * data[, "x1"] ** 0.05 * data[, "x2"] ** 0.15 * data[, "x3"] ** 0.3

    # output generation
    data[, "y"]  <- data[, "yT"] * exp(- u)

  } else if (nX == 6) {

    # theoretical frontier
    data[, "yT"]  <- 3 * data[, "x1"] ** 0.05 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.045 * data[, "x5"] ** 0.1 * data[, "x6"] ** 0.3

    # output generation
    data[, "y"]  <- data[, "yT"] * exp(- u)

  } else if (nX == 9) {

    # theoretical frontier
    data["yT"] <- 3 * data[, "x1"] ** 0.005 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.005 * data[, "x5"] ** 0.001 * data[, "x6"] ** 0.004 *
      data[, "x7"] ** 0.08 * data[, "x8"] ** 0.1 * data[, "x9"] ** 0.3

    # output generation
    data["y"]  <- data["yT"] * exp(- u)

  } else if (nX == 12) {

    # theoretical frontier
    data["yT"] <- 3 * data[, "x1"] ** 0.005 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.005 * data[, "x5"] ** 0.001 * data[, "x6"] ** 0.004 *
      data[, "x7"] ** 0.08 * data[, "x8"] ** 0.05 * data[, "x9"] ** 0.05 *
      data[, "x10"] ** 0.075 * data[, "x11"] ** 0.025 * data[, "x12"] ** 0.2

    # output generation
    data["y"]  <- data["yT"]  * exp(- u)

  } else {

    # theoretical frontier
    data["yT"] <- 3 * data[, "x1"] ** 0.005 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.005 * data[, "x5"] ** 0.001 * data[, "x6"] ** 0.004 *
      data[, "x7"] ** 0.08 * data[, "x8"] ** 0.05 * data[, "x9"] ** 0.05 *
      data[, "x10"] ** 0.05 * data[, "x11"] ** 0.025 * data[, "x12"] ** 0.025 *
      data[, "x13"] ** 0.025 * data[, "x14"] ** 0.025 * data[, "x15"] ** 0.15

    # output generation
    data["y"]  <- data["yT"] * exp(- u)
  }

  return(data)

}

#' @title 1 output ~ nX input Cobb-Douglas Data Generation Process with Noise
#'
#' @description
#' This function simulates a 1 output ~ nX input \code{data.frame} with a Cobb-Douglas Data Generation Process including statistical noise as in \insertCite{simar2011;textual}{aces}.
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param nX
#' An integer representing the number of inputs.
#'
#' @param p
#' Noise-to-signal ratio.
#'
#' @param heteroskedasticity
#' A \code{"logical"} value indicating if heteroskedasticity should be introduced in the inefficiency term. Only available for 3 inputs (\code{nX = 3}).
#'
#' @importFrom dplyr %>%
#' @importFrom stats runif rnorm
#'
#' @return
#'
#' A \code{data.frame} with the simulated data: 1 inputs (x1), 1 output (y) and the theoretical frontier (yT).
#'
#' @references
#'
#' \insertRef{simar2011}{aces}
#'
#' @export

cobb_douglas_XnY1_noise <- function (
    N,
    nX,
    p,
    heteroskedasticity = FALSE
    ) {

  if(!(nX %in% c(1, 3))){
    stop(paste(nX, "is not allowed"))
  }

  if (nX == 1 && heteroskedasticity == TRUE) {
    stop("heteroskedasticity is only availble for 3 inputs")
  }

  colnames <- c(paste("x", 1:nX, sep = ""), "y")

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # input generation
  for (x in 1:nX){
    data[, x] <- runif(n = N, min = 0, max = 1)
  }

  # inefficiency generation
  if (heteroskedasticity) {

    # standard deviations
    sds <- 0.3 * (data[, 1] + data[, 2]) / 2

    # inefficiency term
    u <- sapply(sds, function(sd) abs(rnorm(1, mean = 0, sd = sd)))

  } else {

    u <- rexp(n = N, rate = 6)

  }

  # error generation
  if (heteroskedasticity) {

    v <- rnorm(n = N, mean = 0, sd = p * 0.3 * sqrt((pi - 2) / pi))

  } else {

    v <- rnorm(n = N, mean = 0, sd = p * (1 / 6))

  }

  # theoretical frontier

  if (nX == 1) {

    data[, "yT"] <- data[, "x1"] ^ 0.5

  } else {

    data[, "yT"] <- data[, "x1"] ^ 0.4 * data[, "x2"] ^ 0.3 * data[, "x3"] ^ 0.2

  }

  # output generation
  data[, "y"]  <- data[, "yT"] * exp(- u) * exp(v)

  return(data)

}

#' @title 3 outputs ~ 3 inputs Cobb-Douglas Data Generation Process
#'
#' @description
#' This function simulates a 3 outputs ~ 3 inputs \code{data.frame}  with a Cobb-Douglas Data Generation Process as in \insertCite{ahn2023;textual}{aces}.
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param border
#' Proportion of DMUs in the frontier: \code{0.1}, \code{0.2} or \code{0.3}.
#'
#' @param noise
#' Random noise. 4 possible levels: \code{0}, \code{0.5}, \code{1} or \code{2}.
#'
#' @param returns
#' Returns to scale. \code{"CRS"} for Constant Return to Scale, \code{"DRS"} for Decreasing Return to Scale and \code{"IRS"} for Increasing Return to Scale.
#'
#' @importFrom stats runif rnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom dplyr %>%
#'
#' @return
#' A \code{data.frame} with the simulated data: 3 inputs (x1, x2, x3), 3 outputs (y1, y2, y3) and the theoretical frontier (yT1, yT2, yT3)
#'
#' @references
#' \insertRef{ahn2023}{aces}
#'
#' @export

cobb_douglas_X3Y3 <- function (
    N,
    border,
    noise,
    returns
) {

  # number of inputs
  nX <- 3

  # number of outputs
  nY <- 3

  colnames <- c(
    paste("x", 1:nX, sep = ""),
    paste("y", 1:nY, sep = ""),
    paste("yT", 1:nY, sep = "")
  )

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # input generation
  for (i in 1:nX) {
    data[, i] <- runif(N, 5, 15)
  }

  b0 <- 0

  if (returns == "CRS") {

    b1 <- b2 <- b3 <- 1 / 3

  } else if (return_scale == "DRS") {

    b1 <- b2 <- b3 <- 0.267

  } else {

    b1 <- b2 <- b3 <- 0.4

  }

  yT_aux <- exp(b0) * data[, 1] ^ b1 * data[, 2] ^ b2 * data[, 3] ^ b3

  for (i in 1:nY) {

    alpha <- rtruncnorm(N, 0, 1, 1 / nY, 1 / nY ^ 2)
    a <- alpha / sum(alpha)

    # theoretical frontier
    data[, nX + nY + i] <- sqrt(a) * yT_aux

    # random Noise
    v <- rnorm(N, 0, (noise * 0.136) ^ 2)

    # inefficiency
    if (border == 0.1) {

      mu <- rtruncnorm(N, 0, Inf, 0, 0.136 ^ 2)

    } else if (border == 0.2) {

      mu <- rtruncnorm(N, 0, Inf, 0, 0.299 ^ 2)

    } else if (border == 0.3) {

      mu <- rtruncnorm(N, 0, Inf, 0, 0.488 ^ 2)

    }

    Yobs_aux <- yT_aux * exp(- mu) * exp(v)

    # output generation
    data[, nY + i] <- sqrt(a) * Yobs_aux
  }

  return(data)

}

#' @title 1 output ~ nX relevant inputs and nH irrelevant inputs Cobb-Douglas Data Generation Process
#'
#' @description
#' This function simulates a dataset with one output, \code{nX} relevant inputs, and \code{nH} irrelevant inputs using a Cobb-Douglas production function. Relevant inputs can have specified correlations as in \insertCite{nataraja2011;textual}{aces}, and irrelevant inputs are constructed as linear combinations of relevant inputs and random noise as in \insertCite{peyrache2020;textual}{aces}.
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param nX
#' An integer representing the number of relevant inputs.
#'
#' @param nH
#' An integer representing the number of irrelevant inputs.
#'
#' @param exp_x
#' A numeric vector of exponents (elasticities) for relevant inputs in the Cobb-Douglas function. Should be of length \code{nX}.
#'
#' @param rho_x
#' A numeric indicating correlations between relevant inputs. Must have variable names as row and column names (e.g., 'x1', 'x2', etc.). Each column represents a variable that is correlated with at most one variable in the rows.
#'
#' @param rho_h
#' A numeric value between 0 and 0.2. Weight for relevant inputs in constructing irrelevant inputs.
#'
#' @importFrom stats runif rnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom dplyr %>%
#'
#' @return
#' A \code{data.frame} containing the simulated data: relevant inputs (\code{x1}, \code{x2}, ...), irrelevant inputs (\code{h1}, \code{h2}, ...), the theoretical output (\code{yT1}), and the observed output (\code{y1}).
#'
#' @details
#' - **Relevant Inputs**: Generated with specified correlations. Each variable in \code{rho_x} columns is correlated with at most one variable in \code{rho_x} rows.
#'
#' - **Irrelevant Inputs**: Constructed as linear combinations of a random relevant input and random input generation, controlled by \code{rho_h}.
#'
#' - **Output Variables**:
#'   - \code{yT1}: Theoretical output calculated using the Cobb-Douglas function.
#'   - \code{y1}: Observed output including inefficiency, modelled as a random variable.
#'
#' @references
#' \insertRef{nataraja2011}{aces} \cr \cr
#' \insertRef{peyrache2020}{aces}
#'
#' @export

cobb_douglas_XnHnY1 <- function (
    N,
    nX,
    nH,
    exp_x,
    rho_x,
    rho_h
) {

  # input validation
  if (!is.numeric(N) || length(N) != 1 || N <= 0 || N != as.integer(N)) {
    stop("N must be a positive integer.")
  }

  if (!is.numeric(nX) || length(nX) != 1 || nX <= 0 || nX != as.integer(nX)) {
    stop("nX must be a positive integer.")
  }

  if (!is.numeric(nH) || length(nH) != 1 || nH < 0 || nH != as.integer(nH)) {
    stop("nH must be a non-negative integer.")
  }

  if (!is.numeric(exp_x) || length(exp_x) != nX) {
    stop("exp_x must be a numeric vector of length nX.")
  }

  if (!is.null(rho_x)) {

    if (!is.matrix(rho_x)) {
      stop("rho_x must be a matrix or NULL.")
    }

    if (any(abs(rho_x) > 1)) {
      stop("Correlation coefficients in rho_x must be between -1 and 1.")
    }

    if (!all(rownames(rho_x) %in% paste0("x", 1:nX))) {
      stop("Row names of rho_x must be among 'x1', 'x2', ..., 'x{nX}'.")
    }

    if (!all(colnames(rho_x) %in% paste0("x", 1:nX))) {
      stop("Column names of rho_x must be among 'x1', 'x2', ..., 'x{nX}'.")
    }

    # check that each column in rho_x has at most one non-zero value
    corr_check <- colSums(rho_x != 0)
    if (any(corr_check > 1)) {
      stop("Each column in rho_x must have at most one non-zero value.")
    }

  }

  if (!is.null(rho_h)) {

    if (!is.matrix(rho_h)) {
      stop("rho_h must be a matrix or NULL.")
    }

    if (any(abs(rho_h) > 1)) {
      stop("Correlation coefficients in rho_h must be between -1 and 1.")
    }

    if (!all(rownames(rho_h) %in% paste0("x", 1:nX))) {
      stop("Row names of rho_h must be among 'x1', 'x2', ..., 'x{nX}'.")
    }

    if (!all(colnames(rho_h) %in% paste0("h", 1:nH))) {
      stop("Column names of rho_h must be among 'h1', 'h2', ..., 'h{nH}'.")
    }

    # check that each column in rho_h has at most one non-zero value
    corr_check <- colSums(rho_h != 0)
    if (any(corr_check > 1)) {
      stop("Each column in rho_h must have at most one non-zero value.")
    }

  }

  # initialize data frame
  colnames <- c (
    paste("x", 1:nX, sep = ""),
    if (nH > 0) paste0("h", 1:nH, sep = ""),
    "y1",
    "yT1"
  )

  data <- matrix (
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # relevant input generation
  if (!is.null(rho_x)) {

    for (v in rownames(rho_x)) {
      data[[v]] <- runif(N, 10, 20)
    }

    # relevant input generation: correlated variables
    for (j in 1:ncol(rho_x)) {

      for (i in 1:nrow(rho_x)) {

        if (rho_x[i, j] == 0) next

        var1 <- colnames(rho_x)[j]
        var2 <- rownames(rho_x)[i]
        rvar <- runif(N, 10, 20)

        data[[var1]] <- rho_x[i, j] * data[[var2]] + rvar * sqrt(1 - rho_x[i, j] ^ 2)

      }
    }

  } else {

    for (v in 1:nX) {
      data[, v] <- runif(N, 10, 20)
    }

  }

  # irrelevant input generation
  if (!is.null(rho_h)) {

    for (v in colnames(rho_h)) {
      data[[v]] <- runif(N, 10, 20)
    }

    # relevant input generation: correlated variables
    for (j in 1:ncol(rho_h)) {

      for (i in 1:nrow(rho_h)) {

        if (rho_h[i, j] == 0) next

        var1 <- colnames(rho_h)[j]
        var2 <- rownames(rho_h)[i]
        rvar <- runif(N, 10, 20)

        data[[var1]] <- rho_h[i, j] * data[[var2]] + rvar * sqrt(1 - rho_h[i, j] ^ 2)

      }
    }

  } else {

    for (v in 1:nH) {
      data[, nX + v] <- runif(N, 10, 20)
    }

  }

  # theoretical frontier
  data[, "yT1"] <- apply(data[, 1:nX, drop = FALSE], 1, function(row) prod(row ^ exp_x))

  # inefficient term
  u <- abs(rnorm(n = N, mean = 0, sd = 0.4))

  # observed output
  data[, "y1"] <- data[, "yT1"] * exp(- u)

  return(data)

}

#' @title 1 output ~ nX inputs Additive Data Generation Process
#'
#' @description
#' This function simulates a 1 output ~ nX inputs \code{data.frame} with an Additive Data Generation Process as in \insertCite{kuosmanen2010;textual}{aces}.
#'
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param scenario
#' A character string specifying the scenario. Must be one of \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, \code{"E"} or \code{"F"}. For details, see Table 2 in \insertCite{kuosmanen2010;textual}{aces}.
#'
#' @importFrom dplyr %>% filter
#' @importFrom stats runif rnorm
#' @importFrom Rdpack reprompt
#'
#' @return
#' A \code{data.frame} with the simulated data: 1-3 inputs, 1 output (y) and the theoretical frontier (yT).
#'
#' @references
#' \insertRef{kuosmanen2010}{aces}
#'
#' @export

add_scenario_XnY1 <- function (
    N,
    scenario
) {

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

  data <- matrix (
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # input generation
  for (x in 1:nX){
    data[, x] <- runif(n = N, min = 1, max = 10)
  }

  # inefficiency generation
  u <- abs(rnorm(n = N, mean = 0, sd = 0.4))

  if (scenario == "A") {

    # input
    x1 <- data[, "x1"]

    # theoretical frontier
    data[, "yT"] <- log(x1) + 3

    # output
    data[, "y"] <- data[, "yT"] - u

  } else if (scenario == "B") {

    # input
    x1 <- data[, "x1"]

    # theoretical frontier
    data[, "yT"] <- 3 + sqrt(x1) + log(x1)

    # output
    data[, "y"] <- data[, "yT"] - u

  } else if (scenario == "C") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data[, "yT"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * sqrt(x1 * x2)

    # output
    data[, "y"] <- data[, "yT"] - u

  } else if (scenario == "D") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yT"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 3)

    # output
    data["y"] <- data["yT"] - u

  } else if (scenario == "E") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data["yT"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * (x1 * x2) ^ (1 / 3)

    # output
    data["y"] <- data["yT"] - u

  } else {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yT"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 4)

    # output
    data["y"] <- data["yT"]  - u
  }

  data <- data %>% filter(y >= 0)

  return(data)

}

#' @title 1 output ~ nX inputs Multiplicative Data Generation Process
#'
#' @description
#' This function simulates a 1 output ~ nX inputs \code{data.frame} with a Multiplicative Data Generation Process as in \insertCite{kuosmanen2010;textual}{aces}.
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param scenario
#' A character string specifying the scenario. Must be one of \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, \code{"E"} or \code{"F"}. For details, see Table 2 in \insertCite{kuosmanen2010;textual}{aces}.
#'
#' @importFrom dplyr %>% filter
#' @importFrom stats runif rnorm
#' @importFrom Rdpack reprompt
#'
#' @return
#' A \code{data.frame} with the simulated data: 1-3 inputs, 1 output (y) and the theoretical frontier (yT).
#'
#' @references
#' \insertRef{kuosmanen2010}{aces}
#'
#' @export

mult_scenario_XnY1 <- function (
    N,
    scenario
) {

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

  # input generation
  for (x in 1:nX){
    data[, x] <- runif(n = N, min = 1, max = 10)
  }

  # inefficiency generation
  u <- abs(rnorm(n = N, mean = 0, sd = 0.4))

  if (scenario == "A") {

    # input
    x1 <- data[, "x1"]

    # theoretical frontier
    data[, "yT"] <- log(x1) + 3

    # output
    data[, "y"]  <- y / (1 + u)

  } else if (scenario == "B") {

    # input
    x1 <- data[, "x1"]

    # theoretical frontier
    data[, "yT"] <- 3 + sqrt(x1) + log(x1)

    # output
    data[, "y"]  <- y / (1 + u)

  } else if (scenario == "C") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data[, "yT"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * sqrt(x1 * x2)

    # output
    data[, "y"]  <- y / (1 + u)

  } else if (scenario == "D") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yT"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 3)

    # output
    data["y"]  <- y / (1 + u)

  } else if (scenario == "E") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data["yT"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * (x1 * x2) ^ (1 / 3)

    # output
    data["y"]  <- y / (1 + u)

  } else {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yT"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 4)

    # output
    data["y"]  <- y / (1 + u)
  }

  return(data)

}

#' @title 2 outputs ~ 2 inputs Translog Data Generation Process
#'
#' @description
#' This function simulates a 2 outputs ~ 2 inputs \code{data.frame} with a Translog Data Generation Process as in \insertCite{santin2009;textual}{aces}.
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param border
#' Proportion of DMUs in the frontier.
#'
#' @param noise
#' A \code{logical} indicating whether to add random noise.
#'
#' @importFrom stats runif rnorm
#' @importFrom Rdpack reprompt
#'
#' @return
#' A \code{data.frame} with the simulated data: 2 inputs (x1, x2), 2 outputs (y1, y2) and the theoretical frontier (yT1, yT2).
#'
#' @references
#' \insertRef{santin2009}{aces}
#'
#' @export

translog_X2Y2 <- function (
    N,
    border,
    noise
    ) {

  nX <- 2
  nY <- 2

  colnames <- c(paste("x", 1:nX, sep = ""), paste("y", 1:nY, sep = ""))

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # Input generation
  data[, 1] <- runif(N, 5, 50)
  data[, 2] <- runif(N, 5, 50)

  z <- runif(N, - 1.5, 1.5)

  ln_x1 <- log(data[, "x1"])
  ln_x2 <- log(data[, "x2"])

  op1 <- - 1 + 0.5 * z + 0.25 * (z ^ 2) - 1.5 * ln_x1

  op2 <- - 0.6 * ln_x2 + 0.20 * (ln_x1 ^ 2) + 0.05 * (ln_x2 ^ 2) - 0.1 * ln_x1 * ln_x2

  op3 <- + 0.05 * ln_x1 * z - 0.05 * ln_x2 * z

  ln_y1_ast <- - (op1 + op2 + op3)

  # theoretical frontier
  data[, c("y1", "yT1")] <- exp(ln_y1_ast)
  data[, c("y2", "yT2")] <- exp(ln_y1_ast + z)

  if (border < 1) {

    # output generation
    index <- sample(1:N, N * (1 - border))
    N_sample <- length(index)
    half_normal <- exp(abs(rnorm(N_sample, 0, 0.3 ** (1 / 2))))

    if (noise) {

      normal1 <- exp(rnorm(N_sample, 0, 0.01 ** (1 / 2)))
      normal2 <- exp(rnorm(N_sample, 0, 0.01 ** (1 / 2)))

      data[index, "y1"] <- data[index, "yT1"] / (half_normal * normal1)
      data[index, "y2"] <- data[index, "yT2"] / (half_normal * normal2)

    } else {

      data[index, "y1"] <- data[index, "yT1"] / half_normal
      data[index, "y2"] <- data[index, "yT2"] / half_normal

    }
  }

  return(data)

}

#' @title 2 outputs ~ 2 inputs Constant Elasticity of Transformation - Cobb-Douglas Data Generation Process
#'
#' @description
#'
#' This function simulates a 2 outputs ~ 2 inputs \code{data.frame} using a Data Generation Process (DGP) based on
#' Constant Elasticity of Transformation (CET) on the output side and Cobb-Douglas (CD) on the input side as in \insertCite{fare1994;textual}{aces}.
#'
#' The dataset consists of \code{N} observations divided into four groups with equal proportions:
#' \itemize{
#'   \item \strong{Small-size producers}: \eqn{h(u) = 25}, \eqn{f(x) = 25}, \eqn{\delta = 0.898}
#'   \item \strong{Medium-size producers (Lower bound)}: \eqn{h(u) = 50}, \eqn{f(x) = 50}, \eqn{\delta = 1.000}
#'   \item \strong{Medium-size producers (Upper bound)}: \eqn{h(u) = 75}, \eqn{f(x) = 75}, \eqn{\delta = 1.000}
#'   \item \strong{Large-size producers}: \eqn{h(u) = 100}, \eqn{f(x) = 100}, \eqn{\delta = 0.927}
#' }
#'
#' This structure ensures the artificial technology exhibits regions of increasing, constant, and decreasing returns to scale.
#'
#' @param N
#' An integer representing the sample size.
#'
#' @param border
#' Proportion of DMUs in the frontier.
#'
#' @param noise
#' A \code{logical} indicating whether to add random noise.
#'
#' @references
#' \insertRef{fare1994}{aces}
#'
#' @return
#' A \code{data.frame} with the simulated data: 2 inputs (x1, x2), 2 outputs (y1, y2), the theoretical frontier (yT1, yT2).
#'
#' @export

cet_cd_X2Y2 <- function (
    N,
    border,
    noise
    ) {

  # define producer groups and their parameters
  groups <- data.frame (
    f = c(25, 50, 75, 100),
    h = c(25, 50, 75, 100),
    min_value_x = c(20, 30, 50, 90),
    max_value_x = c(60, 80, 100, 230),
    min_value_y = c(10, 15, 25, 45),
    max_value_y = c(35, 70, 100, 135),
    delta = c(0.898, 1.000, 1.000, 0.927),
    proportion = c(0.25, 0.25, 0.25, 0.25)
  )

  # assign each observation to a group based on proportions
  group_assignments <- sample (
    1:nrow(groups),
    size = N,
    replace = TRUE,
    prob = groups$proportion
    )

  # initialize the data frame
  data <- data.frame(
    x1 = numeric(N),
    x2 = numeric(N),
    y1 = numeric(N),
    y2 = numeric(N),
    yT1 = numeric(N),
    yT2 = numeric(N),
    fx = numeric(N),
    hy = numeric(N)
  )

  # generate inputs and outputs for each group
  for (i in 1:nrow(groups)) {

    # input parameters for group
    delta <- groups[i, "delta"]
    min_value_x <- groups[i, "min_value_x"]
    max_value_x <- groups[i, "max_value_x"]
    fx <- groups[i, "f"]

    # output parameters for group
    hy <- groups[i, "h"]
    min_value_y <- groups[i, "min_value_y"]
    max_value_y <- groups[i, "max_value_y"]

    # indices of observations belonging to the group
    idx_group <- which(group_assignments == i)

    # generate x1
    data$x1[idx_group] <- runif (
      length(idx_group),
      min = min_value_x,
      max = max_value_x
      )

    # generate x2
    data$x2[idx_group] <- (fx / data$x1[idx_group] ^ (0.5 * delta)) ^ (1 / (0.5 * delta))

    # generate yT1
    data$yT1[idx_group] <- runif (
      length(idx_group),
      min = min_value_y,
      max = max_value_y
    )

    # generate yT2
    data$yT2[idx_group] <- sqrt(2 * hy ^ 2 - data$yT1[idx_group] ^ 2)

    # generate y1
    data$y1[idx_group] <- data$yT1[idx_group]

    # generate y2
    data$y2[idx_group] <- data$yT2[idx_group]

    # Log-linear Cobb-Douglas
    data$fx[idx_group] <- (sqrt(data$x1[idx_group]) * sqrt(data$x2[idx_group])) ^ delta

  }

  # generate inefficiency
  if (border < 1) {

    # output generation: inefficiency
    idx_ineff <- sample(1:N, N * (1 - border))

    # generate inefficiency
    half_normal <- exp(abs(rnorm(length(idx_ineff), 0, 0.3 ** (1 / 2))))

    if (noise) {

      normal1 <- exp(rnorm(length(idx_ineff), 0, 0.01 ** (1 / 2)))
      normal2 <- exp(rnorm(length(idx_ineff), 0, 0.01 ** (1 / 2)))

      data[idx_ineff, "y1"] <- data[idx_ineff, "yT1"] / (half_normal * normal1)
      data[idx_ineff, "y2"] <- data[idx_ineff, "yT2"] / (half_normal * normal2)

    } else {

      data[idx_ineff, "y1"] <- data[idx_ineff, "yT1"] / half_normal
      data[idx_ineff, "y2"] <- data[idx_ineff, "yT2"] / half_normal

    }
  }

  # Constant Elasticity of Transformation
  data$hy <- sqrt((1 / 2) * data$y1 ^ 2 + (1 / 2) * data$y2 ^ 2)

  data <- round(data, 4)

  return(data)

}
