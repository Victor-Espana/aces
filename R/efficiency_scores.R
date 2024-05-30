#' @title The output-oriented radial model in the envelopment format
#'
#' @description
#' This function computes the efficiency scores through the output-oriented radial model in the envelopment format.
#'
#' @param tech_xmat
#' A \code{data.frame} or \code{matrix} containing the observed inputs to determine the technology.
#'
#' @param tech_ymat
#' A \code{data.frame} or \code{matrix} containing the observed outputs to determine the technology.
#'
#' @param eval_xmat
#' A \code{data.frame} or \code{matrix} containing the containing the input data of the DMUs to be evaluated.
#'
#' @param eval_ymat
#' A \code{data.frame} or \code{matrix} containing the containing the output data of the DMUs to be evaluated.
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
#'
#' @param returns
#' Type of returns to scale.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective get.variables
#'
#' @return
#' A \code{vector} of \code{"numeric"} scores computed through the output-oriented radial model in the envelopment format.

rad_out <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    convexity,
    returns
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  for (d in 1:eval_dmu) {

    objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
    objVal[1] <- 1

    lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      add.constraint(lps, xt = c(0, tech_xmat[, xi]), "<=",  rhs = eval_xmat[d, xi])
    }

    # outputs
    for (yi in 1:nY) {
      add.constraint(lps, xt = c(- eval_ymat[d, yi], tech_ymat[, yi]), ">=", rhs = 0)
    }

    # technology
    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
      }
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The input-oriented radial model in the envelopment format
#'
#' @description This function computes efficiency scores through the input-oriented radial model in the envelopment format.
#'
#' @param tech_xmat A \code{data.frame} or \code{matrix} containing the observed inputs to determine the technology.
#' @param tech_ymat A \code{data.frame} or \code{matrix} containing the observed outputs to determine the technology.
#' @param eval_xmat A \code{data.frame} or \code{matrix} containing the containing the input data of the DMUs to be evaluated.
#' @param eval_ymat A \code{data.frame} or \code{matrix} containing the containing the output data of the DMUs to be evaluated.
#' @param convexity A \code{logical} value indicating if a convex technology is assumed.
#' @param returns Type of returns to scale.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A \code{vector} of \code{"numeric"} scores computed through the input-oriented radial model in the envelopment format.

rad_inp <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    convexity,
    returns
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  for (d in 1:eval_dmu) {

    objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
    objVal[1] <- 1

    lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
    lp.control(lps, sense = 'min')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      add.constraint(lps, xt = c(- eval_xmat[d, xi], tech_xmat[, xi]), "<=",  rhs = 0)
    }

    # outputs
    for (yi in 1:nY) {
      add.constraint(lps, xt = c(0, tech_ymat[, yi]), ">=", rhs = eval_ymat[d, yi])
    }

    # technology
    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
      }
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The Directional Distance Function
#'
#' @description This function computes efficiency scores through the directional distance function in the envelopment format.
#'
#' @param tech_xmat A \code{data.frame} or \code{matrix} containing the observed inputs to determine the technology.
#' @param tech_ymat A \code{data.frame} or \code{matrix} containing the observed outputs to determine the technology.
#' @param eval_xmat A \code{data.frame} or \code{matrix} containing the containing the input data of the DMUs to be evaluated.
#' @param eval_ymat A \code{data.frame} or \code{matrix} containing the containing the output data of the DMUs to be evaluated.
#' @param direction Direction of the vector to project on the frontier. Two possibilities:
#' \itemize{
#' \item{\code{"mean"}} Projection vector given by the average value of inputs and outputs of all DMUs.
#' \item{\code{"briec"}} Projection vector given by the value of inputs and outputs of the evaluated DMU.
#' }
#' @param convexity A \code{logical} value indicating if a convex technology is assumed.
#' @param returns Type of returns to scale.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A \code{vector} of \code{"numeric"} scores computed through the input-oriented radial model in the envelopment format.

ddf <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    direction,
    convexity,
    returns
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  for (d in 1:eval_dmu) {
    objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
    objVal[1] <- 1

    # structure for lpSolve
    lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    if (direction == "mean") {
      G_x <- matrix(colMeans(tech_xmat), nrow = 1)
      G_y <- matrix(colMeans(tech_ymat), nrow = 1)

    } else {
      G_x <- matrix(eval_xmat[d, ], nrow = 1)
      G_y <- matrix(eval_ymat[d, ], nrow = 1)
    }

    # inputs
    for (xi in 1:nX) {
      add.constraint(lps, xt = c(G_x[, xi], tech_xmat[, xi]), "<=",  rhs = eval_xmat[d, xi])
    }

    # outputs
    for (yi in 1:nY) {
      add.constraint(lps, xt = c(- G_y[, yi], tech_ymat[, yi]), ">=", rhs =  eval_ymat[d, yi])
    }

    # technology
    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
      }
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The output-oriented Russell model in the envelopment format
#'
#' @description This function computes efficiency scores through the output-oriented Russell model in the envelopment format.
#'
#' @param tech_xmat A \code{data.frame} or \code{matrix} containing the observed inputs to determine the technology.
#' @param tech_ymat A \code{data.frame} or \code{matrix} containing the observed outputs to determine the technology.
#' @param eval_xmat A \code{data.frame} or \code{matrix} containing the containing the input data of the DMUs to be evaluated.
#' @param eval_ymat A \code{data.frame} or \code{matrix} containing the containing the output data of the DMUs to be evaluated.
#' @param convexity A \code{logical} value indicating if a convex technology is assumed.
#' @param returns Type of returns to scale.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A \code{vector} of \code{"numeric"} scores computed through the output-oriented Russell model.

rsl_out <- function (
    tech_xmat, tech_ymat, eval_xmat, eval_ymat, convexity, returns
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  for (d in 1:eval_dmu) {

    objVal <- matrix(ncol = nY + tech_dmu, nrow = 1)
    objVal[1:nY] <- 1 / nY

    # structure for lpSolve
    lps <- make.lp(nrow = 0, ncol = tech_dmu + nY)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      add.constraint(lps, xt = c(rep(0, nY), tech_xmat[, xi]), "<=",  rhs = eval_xmat[d, xi])
    }

    # outputs
    for (yi in 1:nY) {
      phi <- rep(0, nY)
      phi[yi] <- - eval_ymat[d, yi]
      add.constraint(lps, xt = c(phi, tech_ymat[, yi]), ">=", rhs = 0)
    }

    # lower bounds: phi >= 1
    phi.idx <- 1:nY
    set.bounds(lps, lower = rep(1, nY), upper = rep(Inf, nY), columns = phi.idx)

    # technology
    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
      }
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The input-oriented Russell model in the envelopment format
#'
#' @description This function computes efficiency scores through the input-oriented Russell model in the envelopment format.
#'
#' @param tech_xmat A \code{data.frame} or \code{matrix} containing the observed inputs to determine the technology.
#' @param tech_ymat A \code{data.frame} or \code{matrix} containing the observed outputs to determine the technology.
#' @param eval_xmat A \code{data.frame} or \code{matrix} containing the containing the input data of the DMUs to be evaluated.
#' @param eval_ymat A \code{data.frame} or \code{matrix} containing the containing the output data of the DMUs to be evaluated.
#' @param convexity A \code{logical} value indicating if a convex technology is assumed.
#' @param returns Type of returns to scale.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A \code{vector} of \code{"numeric"} scores computed through the input-oriented Russell model.

rsl_inp <- function (
    tech_xmat, tech_ymat, eval_xmat, eval_ymat, convexity, returns
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  for (d in 1:eval_dmu) {

    objVal <- matrix(ncol = nX + tech_dmu, nrow = 1)
    objVal[1:nX] <- 1 / nX

    # structure for lpSolve
    lps <- make.lp(nrow = 0, ncol = nX + tech_dmu)
    lp.control(lps, sense = 'min')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      the <- rep(0, nX)
      the[xi] <- - eval_xmat[d, xi]
      add.constraint(lps, xt = c(the, tech_xmat[, yi]), "<=",  rhs = 0)
    }

    # outputs
    for (yi in 1:nY) {
      add.constraint(lps, xt = c(rep(0, nX), tech_xmat[, yi]), ">=", rhs = ymat[d, yi])
    }

    # upper bounds: theta <= 1
    the.idx <- 1:nX
    set.bounds(lps, lower = rep(0, nX), upper = rep(1, nX), columns = the.idx)

    # technology
    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
      }
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The Weighted Additive Model
#'
#' @description This function computes efficiency scores through a Weighted Additive Model.
#'
#' @param tech_xmat A \code{data.frame} or \code{matrix} containing the observed inputs to determine the technology.
#' @param tech_ymat A \code{data.frame} or \code{matrix} containing the observed outputs to determine the technology.
#' @param eval_xmat A \code{data.frame} or \code{matrix} containing the containing the input data of the DMUs to be evaluated.
#' @param eval_ymat A \code{data.frame} or \code{matrix} containing the containing the output data of the DMUs to be evaluated.
#' @param weights Weights for the additive model:
#' \itemize{
#' \item{\code{"WAM"}} Weighted Additive Model.
#' \item{\code{"MIP"}} Measure of Inefficiency Proportions.
#' \item{\code{"NOR"}} Normalized Weighted Additive Model.
#' \item{\code{"RAM"}} Range Adjusted Measure.
#' \item{\code{"BAM"}} Bounded Adjusted Measure.
#' }
#' @param convexity \code{logical} value indicating if a convex technology is assumed.
#' @param returns Type of returns to scale.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A \code{vector} of \code{"numeric"} scores computed through the Weighted Additive Model.

wam <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    weights,
    convexity,
    returns
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  for (d in 1:eval_dmu) {

    # objective function
    objVal <- matrix(ncol = nX + nY + tech_dmu, nrow = 1)

    # Weights
    if (weights == "WAM") {
      # Weighted Additive Model
      objVal[1:(nX + nY)] <- 1

    } else if (weights == "MIP") {
      # Measure of Inefficiency Proportions
      objVal[1:(nX + nY)] <- c(1 / eval_xmat[d, ], 1 / eval_ymat[d, ])

    } else if (weights == "NOR") {
      # Normalized Weighted Additive Model
      objVal[1:(nX + nY)] <- c(1 / apply(eval_xmat, 2, sd), 1 / apply(eval_ymat, 2, sd))

    } else if (weights == "RAM") {
      # Range Adjusted Measure
      xranges <- apply(eval_xmat, 2, max) - apply(eval_xmat, 2, min)
      yranges <- apply(eval_ymat, 2, max) - apply(eval_ymat, 2, min)
      objVal[1:(nX + nY)] <- c(1 / ((nX + nY) * xranges), 1 / ((nX + nY) * yranges))

    } else if (weights == "BAM") {
      # Bounded Adjusted Measure
      p1 <- eval_xmat[d, ] - apply(eval_xmat, 2, min)
      p2 <- apply(eval_ymat, 2, max) - eval_ymat[d, ]
      objVal[1:(nX + nY)] <- c(1 / ((nX + nY) * p1), 1 / ((nX + nY) * p2))

    } else {
      stop(print(paste(weights, "no disponibles")))
    }

    # structure for lpSolve
    lps <- make.lp(nrow = 0, ncol = nX + nY + tech_dmu)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      x_slack <- rep(0, nX)
      x_slack[xi] <- 1
      slacks <- c(x_slack, rep(0, nY))

      add.constraint(lps, xt = c(slacks, tech_xmat[, xi]), "=", rhs = eval_xmat[d, xi])
    }

    # outputs
    for (yi in 1:nY) {
      y_slack <- rep(0, nY)
      y_slack[yi] <- - 1
      slacks <- c(rep(0, nX), y_slack)

      add.constraint(lps, xt = c(slacks, tech_ymat[, yi]), "=", rhs = eval_ymat[d, yi])
    }

    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(rep(0, nX + nY), rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(rep(0, nX + nY), rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + (nX + nY), type = c("binary"))
      }
    }

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title Conditional Efficiency Estimation
#'
#' @description This function computes efficiency scores assuming a zero-truncated normal distribution.
#'
#' @param data
#' A \code{matrix} containing the variables in the model.
#'
#' @param y
#' Column indexes of output variables in \code{data}.
#'
#' @param mean_pred
#' A \code{vector} of mean predicted outputs.
#'
#' @param std_ineff
#' A \code{numeric} indicating the Standard Deviation of the Inefficiency Term.
#'
#' @param std_error
#' A \code{numeric} indicating the Standard Deviation of the Error Term.
#'
#' @param error_type
#' A \code{character} string specifying the error structure that the function will use when fitting the model.
#'
#' @return
#'
#' A \code{vector} of \code{"numeric"} scores computed through the Conditional Estimation.


cond_eff <- function (
    data,
    y,
    mean_pred,
    std_ineff,
    std_error,
    error_type
    ) {

  # number of DMUs in the technology
  scores <- rep(NaN, nrow(data))

  # variance
  sigma2 <- std_ineff ^ 2 + std_error ^ 2

  # standard deviation
  sigma <- sqrt(sigma2)

  # variance_star
  sigma2_star <- std_ineff ^ 2 * std_error ^ 2 / sigma2

  # dev_star
  sigma_star <- sqrt(sigma2_star)

  # lambda
  lambda <- std_ineff / std_error

  # frontier estimation
  y_front <- mean_pred + std_ineff * sqrt(2 / pi)

  # composite error term
  comp_err <- data[, y] - y_front

  # mu_star
  mu_star <- - comp_err * std_ineff ^ 2 / sigma2

  # random variable for normal distribution
  argnorm <- comp_err * lambda / sigma

  # conditional mean (2) = (3) in Jondrow et al. (1981)
  cond_mean <- mu_star + sigma_star * (dnorm(argnorm) / (1 - pnorm(argnorm)))

  if (error_type == "mul") {
    scores <- exp( - cond_mean)

  } else {
    scores <- 1 - cond_mean / y_front

  }

  # If sigma_u == 0: all firms are diagnosed as efficient
  if (std_ineff == 0) {
    scores <- rep(1, nrow(data))
  }

  return(scores)
}

#' @title Compute Efficiency Scores using an Adaptive Constrained Enveloping Splines model.
#'
#' @description
#'
#' This function calculates the efficiency scores for each Decision-Making-Unit through an Adaptive Constrained Enveloping Splines model and a specific technology determined by the user.
#'
#' @param tech_data
#'
#' A \code{data.frame} or a \code{matrix} containing the observed DMUs to determine the technology.
#'
#' @param eval_data
#' A \code{data.frame} or a \code{matrix} containing the DMUs to be evaluated.
#'
#' @param x
#' Column indexes of input variables in \code{tech_data} and \code{eval_data}.
#'
#' @param y
#' Column indexes of output and netput variables in \code{tech_data} and \code{eval_data}.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param method
#' Model prediction used to compute predictions of inputs and obtain a new vector of outputs for \code{eval_data}:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_cubic"}}: Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_quintic"}}: Quintic Smoothed Adaptive Constrained Enveloping Splines model.
#' }
#'
#' @param measure
#' Mathematical programming model to calculate scores:
#' \itemize{
#' \item{\code{rad_out}} The output-oriented radial measure. Efficiency level at 1. (default)
#' \item{\code{rad_inp}} The input-oriented radial measure. Efficiency level at 1.
#' \item{\code{ddf}}     The directional distance function. Efficiency level at 0.
#' \item{\code{rsl_out}} The output-oriented Russell measure. Efficiency level at 1.
#' \item{\code{rsl_inp}} The input-oriented Russell measure. Efficiency level at 1.
#' \item{\code{wam_mip}} The weighted additive model: measure of inefficiency proportions. Efficiency level at 0.
#' \item{\code{wam_ram}} The weighted additive model: range adjusted measure of inefficiency. Efficiency level at 0.
#' }
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed. Only \code{TRUE} is available.
#'
#' @param returns
#' Type of returns to scale:
#' \itemize{
#' \item{\code{"variable"}} Variable returns to scale (default).
#' \item{\code{"constant"}} Constant returns to scale.
#' }
#'
#' @param direction
#' Only applied if \code{measure = "ddf"}. Direction of the vector to project on the frontier:
#' \itemize{
#' \item{\code{"mean"}} Projection vector given by the average value of inputs and outputs of all DMUs. Applied if \code{direction = NULL}.
#' \item{\code{"briec"}} Projection vector given by the value of inputs and outputs of the evaluated DMU.
#' }
#'
#' @param weights
#' Weights for the additive model:
#' \itemize{
#' \item{\code{"WAM"}} Weighted Additive Model.
#' \item{\code{"MIP"}} Measure of Inefficiency Proportions.
#' \item{\code{"NOR"}} Normalized Weighted Additive Model.
#' \item{\code{"RAM"}} Range Adjusted Measure.
#' \item{\code{"BAM"}} Bounded Adjusted Measure.
#' }
#'
#' @param digits
#' Number of decimal places to which efficiency scores are rounded.
#'
#' @importFrom dplyr summarise %>% mutate_if
#' @importFrom stats median quantile sd
#'
#' @export
#'
#' @return
#' A \code{data.frame} with the efficiency scores computed through an Adaptive Constrained Enveloping Splines model.

aces_scores <- function (
    tech_data,
    eval_data,
    x,
    y,
    object,
    method,
    measure = "rad_out",
    convexity = TRUE,
    returns = "variable",
    direction = NULL,
    weights = NULL,
    digits = 3
    ) {

  if (object[[1]][["control"]][["error_type"]] != "add") {
    stop("This function must be used to estimate efficiency scores in additive modeling.")
  }

  # Possible error messages:
  display_errors (
    caller = "aces_scores",
    data = NULL,
    x = NULL,
    y = NULL,
    y_type = NULL,
    model_type = NULL,
    error_type = NULL,
    degree = NULL,
    metric = NULL,
    nterms = NULL,
    err_red = NULL,
    hd_cost = NULL,
    minspan = NULL,
    endspan = NULL,
    kn_grid = NULL,
    d = NULL,
    wc = NULL,
    wq = NULL,
    object = object,
    measure = measure,
    convexity = convexity,
    returns = returns,
    direction = direction,
    digits = digits
  )

  # index of variables
  # inputs are common in all the models
  # outputs and netputs are common (but interchangeable) in all the models
  # then, we can select just the first element of the ACES object.

  var_indexes <- sort (
    c (
    object[[1]][["data"]][["x"]],
    object[[1]][["data"]][["z"]],
    object[[1]][["data"]][["y"]]
    ))

  # check if training and test names are the same
  tr_names <- colnames(object[[1]][["data"]][["df"]])[var_indexes]
  te_names <- colnames(tech_data[c(x, y)])
  ev_names <- colnames(eval_data[c(x, y)])

  # training and test names are different
  check1 <- all(te_names %in% tr_names)

  # training and evaluated names are different
  check2 <- all(ev_names %in% tr_names)

  if (!(check1 & check2)) {
    stop("Different variable names in training, test or evaluated data.")
  }

  # number of inputs
  nX <- length(x)

  # number of outputs
  nY <- length(y)

  # =================== #
  # Data for technology #
  # =================== #

  # matrix of inputs
  tech_xmat <- as.matrix(tech_data[, x])

  # matrix of outputs
  y_hat <- predict (
    object = object,
    newdata = tech_data,
    x = c(x),
    method = method
  )

  tech_ymat <- as.data.frame(tech_ymat)

  # ======================= #
  # Data for evaluated DMUs #
  # ======================= #

  # matrix of inputs
  eval_xmat <- as.matrix(eval_data[, x])

  # matrix of outputs
  eval_ymat <- as.matrix(eval_data[, y])

  if (measure == "rad_out") {
    scores <- rad_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = convexity,
      returns = returns
      )

  } else if (measure == "rad_inp") {
    scores <- rad_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = convexity,
      returns = returns
    )

  } else if (measure == "ddf") {

    if (is.null(direction)) {
      direction <- "mean"
    }

    scores <- ddf (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      direction = direction,
      convexity = convexity,
      returns = returns
    )

  } else if (measure == "rsl_out") {
    scores <- rsl_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = convexity,
      returns = returns
      )

  } else if (measure == "rsl_inp") {
    scores <- rsl_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = convexity,
      returns = returns
    )

  } else if (measure == "wam") {
    scores <- wam (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      weights = weights,
      convexity = convexity,
      returns = returns
    )
  }

  # model name
  if (returns == "constant") {
    model <- "dea_crt"

  } else {
    if (convexity) {
      model <- "dea_vrt"

    } else {
      model <- "aces_fdh"

    }
  }

  model <- paste(model, measure, sep = "_")

  # scores as data.frame
  scores <- as.data.frame(scores)

  # column names
  names(scores) <- model

  # row names
  rownames(scores) <- row.names(tech_data)

  return(round(scores, digits))
}

#' @title Compute Efficiency Scores using a Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) model.
#'
#' @description
#'
#' This function calculates the efficiency scores for each Decision-Making-Unit through a Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) model.
#'
#' @param eval_data
#' A \code{data.frame} or a \code{matrix} containing the DMUs to be evaluated.
#'
#' @param x
#' Column indexes of input variables in \code{tech_data} and \code{eval_data}.
#'
#' @param y
#' Column indexes of output and netput variables in \code{tech_data} and \code{eval_data}.
#'
#' @param object
#' An \code{rf_aces} object.
#'
#' @param method
#' Model prediction used to compute predictions of inputs and obtain a new vector of outputs for \code{eval_data}:
#' \itemize{
#' \item{\code{"aces"}}: Random Forest Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_cubic"}}: Random Forest Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"aces_quintic"}}: Random Forest Quintic Smoothed Adaptive Constrained Enveloping Splines model.
#' }
#'
#' @param digits
#' Number of decimal places to which efficiency scores are rounded.
#'
#' @importFrom dplyr summarise %>% mutate_if
#' @importFrom stats median quantile sd
#'
#' @export
#'
#' @return
#' A \code{data.frame} with the efficiency scores computed through a Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) model.

rf_aces_scores <- function (
    eval_data,
    x,
    y,
    object,
    method,
    digits = 3
    ) {

  if (object[[1]][[1]][["control"]][["error_type"]] != "add") {
    stop("This function must be used to estimate efficiency scores in additive modeling.")
  }

  # Possible error messages:
  display_errors (
    caller = "rf_aces_scores",
    data = NULL,
    x = NULL,
    y = NULL,
    y_type = NULL,
    model_type = NULL,
    error_type = NULL,
    degree = NULL,
    metric = NULL,
    nterms = NULL,
    err_red = NULL,
    hd_cost = NULL,
    minspan = NULL,
    endspan = NULL,
    kn_grid = NULL,
    d = NULL,
    wc = NULL,
    wq = NULL,
    object = object,
    measure = NULL,
    returns = NULL,
    direction = NULL,
    digits = digits
  )

  # index of variables
  # inputs are common in all the models
  # outputs and netputs are common (but interchangeable) in all the models
  # then, we can select just the first element of the ACES object.

  var_indexes <- sort (
    c (
      object[[1]][[1]][["data"]][["x"]],
      object[[1]][[1]][["data"]][["z"]],
      object[[1]][[1]][["data"]][["y"]]
    ))

  # Check if training and test names are the same
  tr_names <- colnames(object[[1]][[1]][["data"]][["df"]])[var_indexes]
  ev_names <- colnames(eval_data[c(x, y)])

  if (!all(ev_names %in% tr_names)) {
    stop("Different variable names in training, test or evaluated data.")
  }

  # number of inputs
  nX <- length(x)

  # number of outputs
  nY <- length(y)

  # =============================== #
  # Radial output efficiency scores #
  # =============================== #

  # point estimation for each profile of inputs
  y_hat_point <- predict (
    object = object,
    newdata = eval_data,
    x = x,
    method = method
  )[["point_estimation"]]

  # ratio y_hat / y_obs
  ratios <- as.data.frame(y_hat_point / eval_data[, y])

  # score: minimum ratio by row
  scores <- apply(ratios, 1, min)

  # scores as data.frame
  scores <- as.data.frame(scores);
  names(scores) <- "rf_aces_scores"
  rownames(scores) <- row.names(eval_data)

  return(round(scores, digits))

}
