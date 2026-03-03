#' @title The output-oriented radial model
#'
#' @description
#' This function computes the efficiency scores through the output-oriented radial model in the envelopment format.
#'
#' @param tech_xmat
#' A \code{matrix} containing the observed input variables used to define the technology.
#'
#' @param tech_ymat
#' A \code{matrix} containing the observed output variables used to define the technology.
#'
#' @param eval_xmat
#' A \code{matrix} containing the input data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param eval_ymat
#' A \code{matrix} containing the output data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
#'
#' @param returns
#' Type of returns to scale. Either \code{"constant"} (CRS) or \code{"variable"} (VRS).
#'
#' @param type
#' A \code{character} string specifying the return value. Options are \code{"objective"} to return the optimal value of the objective function, and \code{"variables"} to return the optimal values of the decision variables related to efficiency measure.
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
    returns,
    type
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

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

    # solve model
    solve(lps)

    # get scores
    if (type == "objective") {
      scores[d, 1] <- get.objective(lps)
    } else if (type == "variables") {
      scores[d, 1] <- get.variables(lps)[1]
    }

  }

  return(scores)

}

#' @title The input-oriented radial model
#'
#' @description
#' This function computes efficiency scores through the input-oriented radial model in the envelopment format.
#'
#' @param tech_xmat
#' A \code{matrix} containing the observed input variables used to define the technology.
#'
#' @param tech_ymat
#' A \code{matrix} containing the observed output variables used to define the technology.
#'
#' @param eval_xmat
#' A \code{matrix} containing the input data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param eval_ymat
#' A \code{matrix} containing the output data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
#'
#' @param returns
#' Type of returns to scale. Either \code{"constant"} (CRS) or \code{"variable"} (VRS).
#'
#' @param type
#' A \code{character} string specifying the return value. Options are \code{"objective"} to return the optimal value of the objective function, and \code{"variables"} to return the optimal values of the decision variables related to efficiency measure.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return
#' A \code{vector} of \code{"numeric"} scores computed through the input-oriented radial model in the envelopment format.

rad_inp <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    convexity,
    returns,
    type
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

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

    # solve model
    solve(lps)

    # get scores
    if (type == "objective") {
      scores[d, 1] <- get.objective(lps)
    } else if (type == "variables") {
      scores[d, 1] <- get.variables(lps)[1]
    }

  }

  return(scores)

}


#' @title The Directional Distance Function
#'
#' @description
#' This function computes efficiency scores through the Directional Distance Function in the envelopment format.
#'
#' @param tech_xmat
#' A \code{matrix} containing the observed input variables used to define the technology.
#'
#' @param tech_ymat
#' A \code{matrix} containing the observed output variables used to define the technology.
#'
#' @param eval_xmat
#' A \code{matrix} containing the input data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param eval_ymat
#' A \code{matrix} containing the output data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param direction
#' Direction of the vector to project on the frontier.
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
#'
#' @param returns
#' Type of returns to scale. Either \code{"constant"} (CRS) or \code{"variable"} (VRS).
#'
#' @param type
#' A \code{character} string specifying the return value. Options are \code{"objective"} to return the optimal value of the objective function, and \code{"variables"} to return the optimal values of the decision variables related to efficiency measure.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return
#' A \code{vector} of \code{"numeric"} scores computed through the directional distance function in the envelopment format.

ddf <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    direction,
    convexity,
    returns,
    type
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # direction of projection
  G_x <- as.matrix(direction[, 1:nX])
  G_y <- as.matrix(direction[, (nX + 1):(nX + nY)])

  for (d in 1:eval_dmu) {

    objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
    objVal[1] <- 1

    # structure for lpSolve
    lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      add.constraint(lps, xt = c(G_x[d, xi], tech_xmat[, xi]), "<=",  rhs = eval_xmat[d, xi])
    }

    # outputs
    for (yi in 1:nY) {
      add.constraint(lps, xt = c(- G_y[d, yi], tech_ymat[, yi]), ">=", rhs =  eval_ymat[d, yi])
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

    set.bounds(lps, lower = c(- Inf, rep(0, tech_dmu)))

    # solve model
    solve(lps)

    # get scores
    if (type == "objective") {
      scores[d, 1] <- get.objective(lps)
    } else if (type == "variables") {
      scores[d, 1] <- get.variables(lps)[1]
    }

  }

  return(scores)

}

#' @title The output-oriented Russell model
#'
#' @description
#' This function computes efficiency scores through the output-oriented Russell model in the envelopment format.
#'
#' @param tech_xmat
#' A \code{matrix} containing the observed input variables used to define the technology.
#'
#' @param tech_ymat
#' A \code{matrix} containing the observed output variables used to define the technology.
#'
#' @param eval_xmat
#' A \code{matrix} containing the input data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param eval_ymat
#' A \code{matrix} containing the output data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
#'
#' @param returns
#' Type of returns to scale. Either \code{"constant"} (CRS) or \code{"variable"} (VRS).
#'
#' @param type
#' A \code{character} string specifying the return value. Options are \code{"objective"} to return the optimal value of the objective function, and \code{"variables"} to return the optimal values of the decision variables related to efficiency measure.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return
#' A \code{vector} of \code{"numeric"} scores computed through the output-oriented Russell model.

rsl_out <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    convexity,
    returns,
    type
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # initialize vector of scores
  if (type == "objective") {
    scores <- matrix(nrow = eval_dmu, ncol = 1)
  } else if (type == "variables") {
    scores <- matrix(nrow = eval_dmu, ncol = nY)
  }

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
      phi[yi] <- eval_ymat[d, yi]
      add.constraint(lps, xt = c(- phi, tech_ymat[, yi]), ">=", rhs = 0)
    }

    # lower bounds: phi >= 1
    set.bounds(lps, lower = rep(1, nY), columns = 1:nY)

    # technology
    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(rep(0, nY), rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(rep(0, nY), rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + nY, type = c("binary"))
      }
    }

    # solve model
    solve(lps)

    # get scores
    if (type == "objective") {
      scores[d, 1] <- get.objective(lps)
    } else if (type == "variables") {
      scores[d, 1:nY] <- get.variables(lps)[1:nY]
    }

  }

  return(scores)

}

#' @title The input-oriented Russell model
#'
#' @description
#' This function computes efficiency scores through the input-oriented Russell model in the envelopment format.
#'
#' @param tech_xmat
#' A \code{matrix} containing the observed input variables used to define the technology.
#'
#' @param tech_ymat
#' A \code{matrix} containing the observed output variables used to define the technology.
#'
#' @param eval_xmat
#' A \code{matrix} containing the input data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param eval_ymat
#' A \code{matrix} containing the output data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
#'
#' @param returns
#' Type of returns to scale. Either \code{"constant"} (CRS) or \code{"variable"} (VRS).
#'
#' @param type
#' A \code{character} string specifying the return value. Options are \code{"objective"} to return the optimal value of the objective function, and \code{"variables"} to return the optimal values of the decision variables related to efficiency measure.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return
#' A \code{vector} of \code{"numeric"} scores computed through the input-oriented Russell model.

rsl_inp <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    convexity,
    returns,
    type
    ) {

  # number of DMUs in theta technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # initialize vector of scores
  if (type == "objective") {
    scores <- matrix(nrow = eval_dmu, ncol = 1)
  } else if (type == "variables") {
    scores <- matrix(nrow = eval_dmu, ncol = nX)
  }

  for (d in 1:eval_dmu) {

    objVal <- matrix(ncol = nX + tech_dmu, nrow = 1)
    objVal[1:nX] <- 1 / nX

    # structure for lpSolve
    lps <- make.lp(nrow = 0, ncol = nX + tech_dmu)
    lp.control(lps, sense = 'min')
    set.objfn(lps, objVal)

    # inputs
    for (xi in 1:nX) {
      theta <- rep(0, nX)
      theta[xi] <- eval_xmat[d, xi]
      add.constraint(lps, xt = c(- theta, tech_xmat[, xi]), "<=",  rhs = 0)
    }

    # outputs
    for (yi in 1:nY) {
      add.constraint(lps, xt = c(rep(0, nX), tech_ymat[, yi]), ">=", rhs = eval_ymat[d, yi])
    }

    # upper bounds: theta <= 1
    set.bounds(lps, upper = rep(1, nX), columns = 1:nX)

    # technology
    if (returns == "variable") {
      if (convexity) {
        add.constraint(lprec = lps, xt = c(rep(0, nX), rep(1, tech_dmu)), type = "=", rhs = 1)
      } else {
        add.constraint(lprec = lps, xt = c(rep(0, nX), rep(1, tech_dmu)), type = "=", rhs = 1)
        set.type(lps, columns = 1:tech_dmu + nX, type = c("binary"))
      }
    }

    # solve model
    solve(lps)

    # get scores
    if (type == "objective") {
      scores[d, 1] <- get.objective(lps)
    } else if (type == "variables") {
      scores[d, 1:nX] <- get.variables(lps)[1:nX]
    }

  }

  return(scores)

}

#' @title The Weighted Additive Model
#'
#' @description
#' This function computes efficiency scores through a Weighted Additive Model.
#'
#' @param tech_xmat
#' A \code{matrix} containing the observed input variables used to define the technology.
#'
#' @param tech_ymat
#' A \code{matrix} containing the observed output variables used to define the technology.
#'
#' @param eval_xmat
#' A \code{matrix} containing the input data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param eval_ymat
#' A \code{matrix} containing the output data of the Decision-Making Units (DMUs) to be evaluated.
#'
#' @param weights Weights for the additive model:
#' \itemize{
#' \item{\code{"wam.one"}} Weighted Additive Model.
#' \item{\code{"wam.mip"}} Measure of Inefficiency Proportions.
#' \item{\code{"wam.nor"}} Normalized Weighted Additive Model.
#' \item{\code{"wam.ram"}} Range Adjusted Measure.
#' \item{\code{"wam.bam"}} Bounded Adjusted Measure.
#' }
#'
#' @param convexity
#' A \code{logical} value indicating if a convex technology is assumed.
#'
#' @param returns
#' Type of returns to scale. Either \code{"constant"} (CRS) or \code{"variable"} (VRS).
#'
#' @param type
#' A \code{character} string specifying the return value. Options are \code{"objective"} to return the optimal value of the objective function, and \code{"variables"} to return the optimal values of the decision variables related to efficiency measure.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return
#' A \code{vector} of \code{"numeric"} scores computed through the Weighted Additive Model.

wam <- function (
    tech_xmat,
    tech_ymat,
    eval_xmat,
    eval_ymat,
    weights,
    convexity,
    returns,
    type
    ) {

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # initialize vector of scores
  if (type == "objective") {
    scores <- matrix(nrow = eval_dmu, ncol = 1)
  } else if (type == "variables") {
    scores <- matrix(nrow = eval_dmu, ncol = nX + nY)
  }

  for (d in 1:eval_dmu) {

    # objective function
    objVal <- matrix(ncol = nX + nY + tech_dmu, nrow = 1)

    # Weights
    if (weights == "wam_mip") {

      # Measure of Inefficiency Proportions
      objVal[1:(nX + nY)] <- c(1 / eval_xmat[d, ], 1 / eval_ymat[d, ])

    } else if (weights == "wam_nor") {

      # Normalized Weighted Additive Model
      objVal[1:(nX + nY)] <- c(1 / apply(eval_xmat, 2, sd), 1 / apply(eval_ymat, 2, sd))

    } else if (weights == "wam_ram") {

      # Range Adjusted Measure
      xranges <- apply(eval_xmat, 2, max) - apply(eval_xmat, 2, min)
      yranges <- apply(eval_ymat, 2, max) - apply(eval_ymat, 2, min)
      objVal[1:(nX + nY)] <- c(1 / ((nX + nY) * xranges), 1 / ((nX + nY) * yranges))

    } else if (weights == "wam_bam") {

      # Bounded Adjusted Measure
      p1 <- eval_xmat[d, ] - apply(eval_xmat, 2, min)
      p2 <- apply(eval_ymat, 2, max) - eval_ymat[d, ]
      objVal[1:(nX + nY)] <- c(1 / ((nX + nY) * p1), 1 / ((nX + nY) * p2))

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

    # solve model
    solve(lps)

    # get scores
    if (type == "objective") {
      scores[d, 1] <- get.objective(lps)
    } else if (type == "variables") {
      scores[d, 1:(nX + nY)] <- get.variables(lps)[1:(nX + nY)]
    }

  }

  return(scores)

}

#' @title Compute Efficiency Scores using an Adaptive Constrained Enveloping Splines model.
#'
#' @description
#' This function computes the efficiency scores for each Decision-Making-Unit (DMU) using an Adaptive Constrained Enveloping Splines (ACES) model.
#'
#' @param eval_data
#' A \code{data.frame} or a \code{matrix} containing the DMUs to be evaluated.
#'
#' @param x
#' Column indexes of input variables in \code{eval_data}.
#'
#' @param y
#' Column indexes of output variables in \code{eval_data}.
#'
#' @param relevant
#' A \code{logical} indicating if only relevant variables should be included in the technology definition.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param method
#' Model prediction method used to compute predictions of inputs and obtain a new vector of outputs for \code{eval_data}:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_cubic"}}: Cubic Smooth Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_quintic"}}: Quintic Smooth Adaptive Constrained Enveloping Splines.
#' }
#'
#' @param measure
#' Mathematical programming model to calculate scores:
#' \itemize{
#' \item{\code{rad_out}} Output-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{rad_inp}} Input-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{ddf}}     Directional distance function proposed by \insertCite{chambers1998;textual}{aces}.
#' \item{\code{rsl_out}} Output-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{rsl_inp}} Input-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{wam_mip}} Measure of Inefficiency Proportions proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{wam_nor}} Normalized Weighted Additive Model proposed by \insertCite{lovell1995;textual}{aces}.
#' \item{\code{wam_ram}} Range Adjusted Measure proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{wam_bam}} Bounded Adjusted Measure proposed by \insertCite{cooper2011;textual}{aces}.
#' }
#'
#' @param returns
#' Type of returns to scale:
#' \itemize{
#' \item{\code{"constant"}} Constant Returns to Scale.
#' \item{\code{"variable"}} Variable Returns to Scale (default).
#' }
#'
#' @param direction
#' Direction of the vector to project on the frontier. Only applied if \code{measure = "ddf"}. A \code{matrix} or \code{data.frame} with \code{n} rows (number of DMUs to be evaluated) and \code{nX + nY} columns, containing the direction of the input variables followed by the direction of the output variables in the same order as they appear in the data.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{banker1984}{aces} \cr \cr
#' \insertRef{chambers1998}{aces} \cr \cr
#' \insertRef{fare1978}{aces} \cr \cr
#' \insertRef{cooper1999}{aces} \cr \cr
#' \insertRef{lovell1995}{aces} \cr \cr
#' \insertRef{cooper2011}{aces}
#'
#' @importFrom dplyr summarise %>% mutate_if
#' @importFrom stats median quantile sd
#'
#' @return
#' A \code{data.frame} with the efficiency scores computed through an Adaptive Constrained Enveloping Splines model.
#'
#' @export

get_scores <- function (
    eval_data,
    x,
    y,
    relevant = FALSE,
    object,
    method = "aces",
    measure = "rad_out",
    returns = "variable",
    direction = NULL
    ) {

  # handle errors:
  display_errors_scores (
    data = eval_data,
    x = x,
    y = y,
    object = object,
    method = method,
    measure = measure,
    returns = returns,
    direction = direction
    )

  # number of inputs
  nX <- length(x)

  # number of outputs
  nY <- length(y)

  # =================== #
  # Data for technology #
  # =================== #

  # matrix of inputs
  tech_xmat <- as.matrix(object[["technology"]][[method]][["xmat"]])

  # auxiliar variable for using only relevant variables
  rel_x <- NULL

  if (relevant) {

    # set of knots
    knots <- object[["methods"]][[method]][["knots"]]

    # variable degree
    xi_degree <- object[["control"]][["xi_degree"]]
    colnames(xi_degree) <- names(object[["control"]][["kn_grid"]])

    # check participating variables
    participating_vars <- intersect(xi_degree[1, ], unique(knots$xi))

    # names of participant variables
    participating_vars_names <- colnames(xi_degree)[participating_vars]

    # extract individual variables assuming interaction variable names use "_"
    split_vars <- unique(unlist(strsplit(participating_vars_names, "_")))

    # get the column indices for the unique variables
    rel_x <- sort(match(split_vars, colnames(xi_degree)))

    # update technology
    tech_xmat <- as.matrix(tech_xmat[, rel_x, drop = FALSE])

  }

  # matrix of outputs
  tech_ymat <- as.matrix(object[["technology"]][[method]][["ymat"]])

  # ======================= #
  # Data for evaluated DMUs #
  # ======================= #

  # matrix of inputs
  if (relevant && !is.null(rel_x)) {
    eval_xmat <- as.matrix(eval_data[, x[rel_x]])
  } else {
    eval_xmat <- as.matrix(eval_data[, x])
  }

  # matrix of outputs
  eval_ymat <- as.matrix(eval_data[, y])

  # scaling setup
  scaling <- object[["control"]][["scale"]]

  if (!is.null(scaling) && scaling$is_scaled) {

    if (relevant) {
      sx <- scaling$mean_x[rel_x]
    } else {
      sx <- scaling$mean_x
    }

    sy <- scaling$mean_y

    # apply scaling
    eval_xmat <- sweep(eval_xmat, 2, sx, "/")
    eval_ymat <- sweep(eval_ymat, 2, sy, "/")

    if (!is.null(direction) && (is.matrix(direction) || is.data.frame(direction))) {

      dir_x <- as.matrix(direction[, 1:nX])
      dir_y <- as.matrix(direction[, (nX + 1):(nX + nY)])

      if (relevant) {
        if (ncol(dir_x) == length(scaling$mean_x)) {
          dir_x <- dir_x[, rel_x, drop = FALSE]
        }
      }

      # apply scaling
      dir_x_scaled <- sweep(dir_x, 2, sx, "/")
      dir_y_scaled <- sweep(dir_y, 2, sy, "/")

      # rebuild direction vectors
      direction <- cbind(dir_x_scaled, dir_y_scaled)

    }

  }

  if (measure == "rad_out") {

    scores <- rad_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
      )

  } else if (measure == "rad_inp") {

    scores <- rad_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )

  } else if (measure == "ddf") {

    scores <- ddf (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      direction = direction,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )

  } else if (measure == "rsl_out") {

    scores <- rsl_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
      )

  } else if (measure == "rsl_inp") {

    scores <- rsl_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )

  } else if (grepl("wam", measure)) {

    scores <- wam (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      weights = measure,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )

  }

  # model name
  model <- ifelse(returns == "constant", "aces_crt", "aces_vrt")
  model <- paste(model, measure, sep = "_")

  # transform scores as a data.frame format
  scores <- as.data.frame(scores)

  # column names of the new data.frame
  colnames(scores) <- model

  # row names of the new data.frame
  rownames(scores) <- row.names(eval_data)

  return(scores)

}

#' @title Compute Efficiency Scores using a Random Forest-Adaptive Constrained Enveloping Splines (RF-ACES) model.
#'
#' @description
#'
#' This function computes the efficiency scores for each Decision-Making-Unit (DMU) by first generating a virtual dataset using a Random Forest-Adaptive Constrained Enveloping Splines (ACES) model. It then constructs a standard DEA technology based on this virtual dataset and evaluates the efficiency scores according to a user-defined measure.
#'
#' @param eval_data
#' A \code{data.frame} or a \code{matrix} containing the DMUs to be evaluated.
#'
#' @param x
#' Column indexes of input variables in \code{eval_data}. The user may choose to include only the variables considered relevant based on the results of the \code{rf_aces_varimp} function.
#'
#' @param y
#' Column indexes of output variables in \code{eval_data}.
#'
#' @param object
#' A \code{rf_aces} object.
#'
#' @param method
#' Model prediction used to compute predictions of inputs and obtain a new vector of outputs for \code{eval_data}:
#' \itemize{
#' \item{\code{"rf_aces"}}: Random Forest Adaptive Constrained Enveloping Splines model.
#' \item{\code{"rf_aces_cubic"}}: Random Forest Cubic Smoothed Adaptive Constrained Enveloping Splines model.
#' \item{\code{"rf_aces_quintic"}}: Random Forest Quintic Smoothed Adaptive Constrained Enveloping Splines model.
#' }
#'
#' #' @param measure
#' Mathematical programming model to calculate scores:
#' \itemize{
#' \item{\code{rad_out}} The output-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{rad_inp}} The input-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{ddf}}     The directional distance function proposed by \insertCite{chambers1998;textual}{aces}.
#' \item{\code{rsl_out}} The output-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{rsl_inp}} The input-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{wam}}     A weighted additive model.
#' \item{\code{rf_aces_rad_out}} The output-oriented radial efficiency measure derived from the Random Forest-Adaptive Constrained Enveloping Splines model proposed by \insertCite{espana2024rf;textual}{aces}.
#' }
#'
#' @param returns
#' Type of returns to scale:
#' \itemize{
#' \item{\code{"constant"}} Constant returns to scale.
#' \item{\code{"variable"}} Variable returns to scale (default).
#' }
#'
#' @param direction
#' Direction of the vector to project on the frontier. Only applied if \code{measure = "ddf"}.
#' \itemize{
#' \item{\code{"mean"}}  Projection vector given by the average value of inputs and outputs of all DMUs. Applied if \code{direction = NULL}.
#' \item{\code{"briec"}} Projection vector given by the value of inputs and outputs of the evaluated DMU.
#' }
#'
#' @param weights
#' Weights for the additive model. Only applied if \code{measure = "wam"}.
#' \itemize{
#' \item{\code{"ONE"}} Weighted Additive Model proposed by \insertCite{charnes1985;textual}{aces}.
#' \item{\code{"MIP"}} Measure of Inefficiency Proportions proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{"NOR"}} Normalized Weighted Additive Model proposed by \insertCite{lovell1995;textual}{aces}.
#' \item{\code{"RAM"}} Range Adjusted Measure proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{"BAM"}} Bounded Adjusted Measure proposed by \insertCite{cooper2011;textual}{aces}.
#' }
#'
#' @param digits
#' Number of decimal places to which efficiency scores are rounded.
#'
#' @importFrom dplyr summarise %>% mutate_if
#' @importFrom stats median quantile sd
#'
#' @references
#'
#' \insertRef{espana2024rf}{aces} \cr \cr
#' \insertRef{banker1984}{aces} \cr \cr
#' \insertRef{chambers1998}{aces} \cr \cr
#' \insertRef{fare1978}{aces} \cr \cr
#' \insertRef{charnes1985}{aces} \cr \cr
#' \insertRef{cooper1999}{aces} \cr \cr
#' \insertRef{lovell1995}{aces} \cr \cr
#' \insertRef{cooper2011}{aces}
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
    method = "rf_aces",
    measure = "rad_out",
    returns = "variable",
    direction = NULL,
    weights = NULL,
    digits = 3
    ) {

  # error handling
  display_errors_scores (
    data = eval_data,
    x = x,
    y = y,
    object = object,
    method = method,
    measure = measure,
    returns = returns,
    direction = direction,
    weights = weights,
    digits = digits
  )

  # number of inputs
  nX <- length(x)

  # number of outputs
  nY <- length(y)

  # =================== #
  # Data for technology #
  # =================== #

  # matrix of inputs
  tech_xmat <- as.matrix(object[["technology"]][[method]][["xmat"]])

  # matrix of outputs
  tech_ymat <- as.matrix(object[["technology"]][[method]][["ymat"]])

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
      convexity = TRUE,
      returns = returns
    )

  } else if (measure == "rad_inp") {

    scores <- rad_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns
    )

  } else if (measure == "ddf") {

    if (is.null(direction)) direction <- "mean"

    scores <- ddf (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      direction = direction,
      convexity = TRUE,
      returns = returns
    )

  } else if (measure == "rsl_out") {

    scores <- rsl_out (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns
    )

  } else if (measure == "rsl_inp") {

    scores <- rsl_inp (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns
    )

  } else if (measure == "wam") {

    scores <- wam (
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      weights = weights,
      convexity = TRUE,
      returns = returns
    )

  } else {

    y_hat_point <- predict (
      object = object,
      newdata = eval_data,
      x = x,
      method = method
    )

    # ratio y_hat / y_obs
    ratios <- as.data.frame(y_hat_point / eval_data[, y])

    # score: minimum ratio by row
    scores <- apply(ratios, 1, min)

  }

  # model name
  model <- ifelse(returns == "constant", "aces_crt", "aces_vrt")
  model <- paste(model, measure, sep = "_")

  # transform scores as a data.frame format
  scores <- as.data.frame(scores)

  # column names of the new data.frame
  colnames(scores) <- model

  # row names of the new data.frame
  rownames(scores) <- row.names(eval_data)

  return(round(scores, digits))

}
