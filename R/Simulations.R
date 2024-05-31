#' @title Random Sample for Efficiency Analysis
#'
#' @description
#'
#' This function is used to simulate a \code{data.frame} with some given Data Generation Process.
#'
#' @param DGP
#' Data Generation Process:
#' \itemize {
#'    \item{\code{"cobb_douglas_XnY1"}}: Cobb-Douglas Data Generation Process for a XnY1 scenario.
#'    See \code{help("cobb_douglas_XnY1")}.
#'    \item{\code{"cobb_douglas_X1Y1_noise"}}: Cobb-Douglas Data Generation Process for a X1Y1 scenario with noise.
#'    See \code{help("cobb_douglas_X1Y1")}.
#'    \item{\code{"translog_X2Y2"}}: Translog Data Generation Process for a X2Y2 scenario.
#'    See \code{help("translog_X2Y2")}.
#'    \item{\code{"add_scenario_XnY1"}}: Additive Data Generation Process for a XnY1 scenario.
#'    See \code{help("add_scenario_XnY1")}.
#'    \item{\code{"mult_scenario_XnY1"}}: Multiplicative Data Generation Process for a XnY1 scenario.
#'    See \code{help("mult_scenario_XnY1")}.
#'    \item{\code{cobb_douglas_X3Y3}}: Cobb-Douglas Data Generation Process for a X3Y3 scenario.
#'    See \code{help("cobb_douglas_X3Y3")}.
#'  }
#'
#' @param parms
#' A \code{list} containing the parameters for the simulation functions.
#' \itemize{
#'    \item{\code{"N"}}: An integer representing the sample size.
#'    \item{\code{"nX"}}: An integer representing the number of inputs.
#'    Applicable  for the \code{cobb_douglas_XnY1(N, nX)} function.
#'    \item{\code{"border"}}: Proportion of DMUs in the frontier.
#'    Applicable  for the \code{translog_X2Y2(N, border, noise)} and the \code{cobb_douglas_X3Y3(N, border, noise, returns)} functions.
#'    \item{\code{"noise"}}: \code{logical} indicating presence of random noise.
#'    Applicable  for the \code{translog_X2Y2(N, border, noise)} and the \code{cobb_douglas_X3Y3(N, border, noise, returns)} functions.
#'    \item{\code{"scenario"}}: Determine the type of Data Generation Process.
#'    Applicable  for the \code{add_scenario_XnY1(N, scenario)} and the \code{mult_scenario_XnY1(N, scenario)} functions.
#'    \item{\code{"returns"}}: Returns to scale.
#'    Applicable  for the \code{cobb_douglas_X3Y3(N, border, noise, returns)} function.
#'    \item{\code{"p"}}: Signal-to-noise ratio.
#'    Applicable  for the \code{cobb_douglas_X1Y1_noise(N, p)} function.
#'  }
#'
#' @details
#' Please refer to the help manuals of the mentioned functions for usage examples and more detailed information on the parameters.
#'
#'
#' @return
#'
#' A \code{data.frame} simulated with the selected Data Generation Process.
#'
#' @export

reffcy <- function (
    DGP, parms
    ) {

  if (DGP == "cobb_douglas_XnY1") {

    data <- cobb_douglas_XnY1 (
      N = parms[["N"]],
      nX = parms[["nX"]]
      )

  } else if (DGP == "cobb_douglas_XnY1_noise") {

    data <- cobb_douglas_XnY1_noise (
      N = parms[["N"]],
      nX = parms[["nX"]],
      p = parms[["p"]]
    )

  } else if (DGP == "cobb_douglas_X3Y3") {

    data <- cobb_douglas_X3Y3 (
      N = parms[["N"]],
      border = parms[["border"]],
      noise = parms[["noise"]],
      returns = parms[["returns"]]
    )

  } else if (DGP == "translog_X2Y2") {

    data <- translog_X2Y2 (
      N = parms[["N"]],
      border = parms[["border"]],
      noise = parms[["noise"]]
    )

  } else if (DGP == "add_scenario_XnY1") {

    data <- add_scenario_XnY1 (
      N = parms[["N"]],
      scenario = parms[["scenario"]]
    )

  } else if (DGP == "mult_scenario_XnY1") {

    data <- mult_scenario_XnY1 (
      N = parms[["N"]],
      scenario = parms[["scenario"]]
    )

  } else (

    stop(paste(DGP, "not availble."))

  )

  return(data)

}

#' @title 1 output ~ nX inputs Cobb-Douglas Data Generation Process
#'
#' @description
#' This function is used to simulate a 1 output ~ nX inputs \code{data.frame} with a Cobb-Douglas Data Generation Process.
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
#' A \code{data.frame} with the simulated data: nX inputs, 1 output (y) and the theoretical frontier (yD).
#'
#' @keywords internal

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
    data[, "yD"] <- 3 * data[, "x1"] ** 0.5

    # output generation
    data[, "y"]  <- data[, "yD"] * exp(- u)

  } else if (nX == 3) {

    # theoretical frontier
    data[, "yD"] <- 3 * data[, "x1"] ** 0.05 * data[, "x2"] ** 0.15 * data[, "x3"] ** 0.3

    # output generation
    data[, "y"]  <- data[, "yD"] * exp(- u)

  } else if (nX == 6) {

    # theoretical frontier
    data[, "yD"]  <- 3 * data[, "x1"] ** 0.05 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.045 * data[, "x5"] ** 0.1 * data[, "x6"] ** 0.3

    # output generation
    data[, "y"]  <- data[, "yD"] * exp(- u)

  } else if (nX == 9) {

    # theoretical frontier
    data["yD"] <- 3 * data[, "x1"] ** 0.005 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.005 * data[, "x5"] ** 0.001 * data[, "x6"] ** 0.004 *
      data[, "x7"] ** 0.08 * data[, "x8"] ** 0.1 * data[, "x9"] ** 0.3

    # output generation
    data["y"]  <- data["yD"] * exp(- u)

  } else if (nX == 12) {

    # theoretical frontier
    data["yD"] <- 3 * data[, "x1"] ** 0.005 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.005 * data[, "x5"] ** 0.001 * data[, "x6"] ** 0.004 *
      data[, "x7"] ** 0.08 * data[, "x8"] ** 0.05 * data[, "x9"] ** 0.05 *
      data[, "x10"] ** 0.075 * data[, "x11"] ** 0.025 * data[, "x12"] ** 0.2

    # output generation
    data["y"]  <- data["yD"]  * exp(- u)

  } else {

    # theoretical frontier
    data["yD"] <- 3 * data[, "x1"] ** 0.005 * data[, "x2"] ** 0.001 * data[, "x3"] ** 0.004 *
      data[, "x4"] ** 0.005 * data[, "x5"] ** 0.001 * data[, "x6"] ** 0.004 *
      data[, "x7"] ** 0.08 * data[, "x8"] ** 0.05 * data[, "x9"] ** 0.05 *
      data[, "x10"] ** 0.05 * data[, "x11"] ** 0.025 * data[, "x12"] ** 0.025 *
      data[, "x13"] ** 0.025 * data[, "x14"] ** 0.025 * data[, "x15"] ** 0.15

    # output generation
    data["y"]  <- data["yD"] * exp(- u)
  }

  return(data)
}

#' @title 1 output ~ nX input Cobb-Douglas Data Generation Process with Noise
#'
#' @description
#' This function is used to simulate a 1 output ~ nX input \code{data.frame} with a Cobb-Douglas Data Generation Process including statistical noise as in \insertCite{simar2011;textual}{aces}.
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
#' A \code{data.frame} with the simulated data: 1 inputs (x1), 1 output (y) and the theoretical frontier (yD).
#'
#' @references
#'
#' \insertRef{simar2011}{aces}
#'
#' @keywords internal

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

    data[, "yD"] <- data[, "x1"] ^ 0.5

  } else {

    data[, "yD"] <- data[, "x1"] ^ 0.4 * data[, "x2"] ^ 0.3 * data[, "x3"] ^ 0.2

  }

  # output generation
  data[, "y"]  <- data[, "yD"] * exp(- u) * exp(v)

  return(data)
}

#' @title 2 outputs ~ 2 inputs Translog Data Generation Process
#'
#' @description
#' This function is used to simulate a 2 outputs ~ 2 inputs \code{data.frame} with a Translog Data Generation Process as in \insertCite{santin2009;textual}{aces}.
#'
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
#' A \code{data.frame} with the simulated data: 2 inputs (x1, x2), 2 outputs (y1, y2) and the theoretical frontier (yD1, yD2).
#'
#' @references
#' \insertRef{santin2009}{aces}
#'
#' @keywords internal

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
  data[, c("y1", "yD1")] <- exp(ln_y1_ast)
  data[, c("y2", "yD2")] <- exp(ln_y1_ast + z)

  if (border < 1) {

    # output generation
    index <- sample(1:N, N * (1 - border))
    N_sample <- length(index)
    half_normal <- exp(abs(rnorm(N_sample, 0, 0.3 ** (1 / 2))))

    if (noise) {

      normal1 <- exp(rnorm(N_sample, 0, 0.01 ** (1 / 2)))
      normal2 <- exp(rnorm(N_sample, 0, 0.01 ** (1 / 2)))

      data[index, "y1"] <- data[index, "yD1"] / (half_normal * normal1)
      data[index, "y2"] <- data[index, "yD2"] / (half_normal * normal2)

    } else {

      data[index, "y1"] <- data[index, "yD1"] / half_normal
      data[index, "y2"] <- data[index, "yD2"] / half_normal

    }
  }

  return(data)
}

#' @title 1 output ~ nX inputs Additive Data Generation Process
#'
#' @description
#' This function is used to simulate a 1 output ~ nX inputs \code{data.frame} with an Additive Data Generation Process as in \insertCite{kuosmanen2010;textual}{aces}.
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
#' A \code{data.frame} with the simulated data: 1-3 inputs, 1 output (y) and the theoretical frontier (yD).
#'
#' @references
#' \insertRef{kuosmanen2010}{aces}
#'
#' @keywords internal

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
    data[, "yD"] <- log(x1) + 3

    # output
    data[, "y"] <- data[, "yD"] - u

  } else if (scenario == "B") {

    # input
    x1 <- data[, "x1"]

    # theoretical frontier
    data[, "yD"] <- 3 + sqrt(x1) + log(x1)

    # output
    data[, "y"] <- data[, "yD"] - u

  } else if (scenario == "C") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data[, "yD"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * sqrt(x1 * x2)

    # output
    data[, "y"] <- data[, "yD"] - u

  } else if (scenario == "D") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yD"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 3)

    # output
    data["y"] <- data["yD"] - u

  } else if (scenario == "E") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data["yD"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * (x1 * x2) ^ (1 / 3)

    # output
    data["y"] <- data["yD"] - u

  } else {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yD"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 4)

    # output
    data["y"] <- data["yD"]  - u
  }

  data <- data %>% filter(y >= 0)

  return(data)

}

#' @title 1 output ~ nX inputs Multiplicative Data Generation Process
#'
#' @description
#' This function is used to simulate a 1 output ~ nX inputs \code{data.frame} with a Multiplicative Data Generation Process as in \insertCite{kuosmanen2010;textual}{aces}.
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
#' A \code{data.frame} with the simulated data: 1-3 inputs, 1 output (y) and the theoretical frontier (yD).
#'
#' @references
#' \insertRef{kuosmanen2010}{aces}
#'
#' @keywords internal

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
    data[, "yD"] <- log(x1) + 3

    # output
    data[, "y"]  <- y / (1 + u)

  } else if (scenario == "B") {

    # input
    x1 <- data[, "x1"]

    # theoretical frontier
    data[, "yD"] <- 3 + sqrt(x1) + log(x1)

    # output
    data[, "y"]  <- y / (1 + u)

  } else if (scenario == "C") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data[, "yD"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * sqrt(x1 * x2)

    # output
    data[, "y"]  <- y / (1 + u)

  } else if (scenario == "D") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yD"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 3)

    # output
    data["y"]  <- y / (1 + u)

  } else if (scenario == "E") {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]

    # theoretical frontier
    data["yD"] <- 0.1 * x1 + 0.1 * x2 + 0.3 * (x1 * x2) ^ (1 / 3)

    # output
    data["y"]  <- y / (1 + u)

  } else {

    # inputs
    x1 <- data[, "x1"]
    x2 <- data[, "x2"]
    x3 <- data[, "x3"]

    # theoretical frontier
    data["yD"] <- 0.1 * x1 + 0.1 * x2 + 0.1 * x3 + 0.3 * (x1 * x2 * x3) ^ (1 / 4)

    # output
    data["y"]  <- y / (1 + u)
  }

  return(data)

}

#' @title 3 outputs ~ 3 inputs Cobb-Douglas Data Generation Process
#'
#' @description
#' This function is used to simulate a 3 outputs ~ 3 inputs \code{data.frame}  with a Cobb-Douglas Data Generation Process as in \insertCite{ahn2023;textual}{aces}.
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
#' A \code{data.frame} with the simulated data: 3 inputs (x1, x2, x3), 3 outputs (y1, y2, y3) and the theoretical frontier (yD1, yD2, yD3)
#'
#' @references
#' \insertRef{ahn2023}{aces}
#'
#' @keywords internal

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
    paste("yD", 1:nY, sep = "")
    )

  data <- matrix(
    ncol = length(colnames),
    nrow = N,
    dimnames = list(NULL, colnames)
  ) %>% as.data.frame()

  # input generation
  data[, 1:nX] <- runif(N, 5, 15)

  b0 <- 0

  if (returns == "CRS") {

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

    # theoretical frontier
    data[, nX + nY + i] <- sqrt(a) * Yd_aux

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

    Yobs_aux <- Yd_aux * exp(- mu) * exp(v)

    # output generation
    data[, nY + i] <- sqrt(a) * Yobs_aux
  }

  return(data)

}
