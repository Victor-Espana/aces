#' @title Compute Output-Oriented Radial Scores
#'
#' @description
#' Solves the output-oriented radial envelopment model for each evaluated DMU.
#' All outputs are expanded by a common factor while inputs remain fixed; a DMU
#' on the frontier has an expansion factor of one.
#'
#' @param tech_xmat
#' Matrix of inputs for the reference technology, with one row per reference DMU.
#'
#' @param tech_ymat
#' Matrix of outputs for the reference technology, with one row per reference DMU.
#'
#' @param eval_xmat
#' Matrix of inputs for the evaluated DMUs.
#'
#' @param eval_ymat
#' Matrix of outputs for the evaluated DMUs.
#'
#' @param convexity
#' If \code{TRUE}, use a convex DEA technology. If \code{FALSE}, use binary
#' intensity variables for an FDH technology.
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param type
#' Value to return: \code{"objective"} for the objective value or
#' \code{"variables"} for the efficiency variable \eqn{\phi}.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds set.rhs set.mat get.objective get.variables
#'
#' @return
#' A one-column matrix with the requested objective value or radial expansion
#' factor for each evaluated DMU.

rad_out <- function(
  tech_xmat,
  tech_ymat,
  eval_xmat,
  eval_ymat,
  convexity,
  returns,
  type = "objective"
) {
  # FDH-VRS scores: closed-form solution by enumeration over dominating units.
  if (!convexity && returns == "variable") {
    return(
      rad_out_fdh(
        tech_xmat = tech_xmat,
        tech_ymat = tech_ymat,
        eval_xmat = eval_xmat,
        eval_ymat = eval_ymat
      )
    )
  }

  # number of DMUs in the technology
  tech_dmu <- nrow(tech_xmat)

  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs and outputs
  nX <- ncol(tech_xmat)
  nY <- ncol(tech_ymat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  # objective function: max phi
  objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
  objVal[1] <- 1

  # ======================== #
  # Build LP model structure #
  # ======================== #

  lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
  lp.control(lps, sense = "max")
  set.objfn(lps, objVal)

  # input constraints: 0*phi + sum(lambda_j * x_ji) <= x_di
  for (xi in 1:nX) {
    add.constraint(lps, xt = c(0, tech_xmat[, xi]), "<=", rhs = 0)
  }

  # output constraints: -y_dr*phi + sum(lambda_j * y_jr) >= 0
  for (yi in 1:nY) {
    add.constraint(lps, xt = c(-1, tech_ymat[, yi]), ">=", rhs = 0)
  }

  # technology constraint
  if (returns == "variable") {
    add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)

    if (!convexity) {
      set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
    }
  }

  # ============================= #
  # Solve for each evaluated DMU  #
  # ============================= #

  for (d in 1:eval_dmu) {

    # update RHS of input constraints
    for (xi in 1:nX) {
      set.rhs(lps, b = eval_xmat[d, xi], constraints = xi)
    }

    # update coefficient of phi in output constraints
    for (yi in 1:nY) {
      set.mat(lps, nX + yi, 1, -eval_ymat[d, yi])
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

#' @title Compute Output-Oriented FDH Scores
#'
#' @description
#' Computes output-oriented radial scores under a variable-returns-to-scale Free
#' Disposal Hull (FDH) technology:
#' \deqn{\phi_d = \max_{j : x_j \le x_d} \min_r (y_{jr} / y_{dr})}
#' Each DMU is compared with observed reference units that use no more of any
#' input; convex combinations of reference units are not allowed.
#'
#' @param tech_xmat
#' Matrix of inputs for the reference technology, with one row per reference DMU.
#'
#' @param tech_ymat
#' Matrix of outputs for the reference technology, with one row per reference DMU.
#'
#' @param eval_xmat
#' Matrix of inputs for the evaluated DMUs.
#'
#' @param eval_ymat
#' Matrix of outputs for the evaluated DMUs.
#'
#' @return
#' A one-column matrix with one efficiency score per evaluated DMU.

rad_out_fdh <- function(
  tech_xmat,
  tech_ymat,
  eval_xmat,
  eval_ymat
) {
  # number of DMUs to be evaluated
  eval_dmu <- nrow(eval_xmat)

  # number of inputs
  nX <- ncol(tech_xmat)

  # transposed technology inputs (for vectorized dominance check)
  t_tech_xmat <- t(tech_xmat)

  # initialize vector of scores
  scores <- matrix(nrow = eval_dmu, ncol = 1)

  for (d in 1:eval_dmu) {
    # dominating units: x_j <= x_d component-wise
    dominates <- colSums(t_tech_xmat <= eval_xmat[d, ]) == nX

    # infeasible if no unit dominates the evaluated DMU
    if (!any(dominates)) {
      scores[d, 1] <- NA
      next
    }

    # output ratios y_jr / y_dr for dominating units
    ratios <- tech_ymat[dominates, , drop = FALSE] /
      matrix(
        eval_ymat[d, ],
        nrow = sum(dominates),
        ncol = ncol(tech_ymat),
        byrow = TRUE
      )

    # phi_d = max_j min_r (y_jr / y_dr)
    scores[d, 1] <- max(apply(ratios, 1, min))
  }

  return(scores)
}

#' @title Compute Input-Oriented Radial Scores
#'
#' @description
#' Solves the input-oriented radial envelopment model for each evaluated DMU.
#' All inputs are contracted by a common factor while outputs remain fixed; a DMU
#' on the frontier has a contraction factor of one.
#'
#' @param tech_xmat
#' Matrix of inputs for the reference technology, with one row per reference DMU.
#'
#' @param tech_ymat
#' Matrix of outputs for the reference technology, with one row per reference DMU.
#'
#' @param eval_xmat
#' Matrix of inputs for the evaluated DMUs.
#'
#' @param eval_ymat
#' Matrix of outputs for the evaluated DMUs.
#'
#' @param convexity
#' If \code{TRUE}, use a convex DEA technology. If \code{FALSE}, use binary
#' intensity variables for an FDH technology.
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param type
#' A \code{character} string specifying the return value:
#' \code{"objective"} for the optimal objective value, or
#' \code{"variables"} for the optimal value of the efficiency variable \eqn{\theta}.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds set.rhs set.mat get.objective get.variables
#'
#' @return
#' A one-column matrix with the requested objective value or radial contraction
#' factor for each evaluated DMU.

rad_inp <- function(
  tech_xmat,
  tech_ymat,
  eval_xmat,
  eval_ymat,
  convexity,
  returns,
  type = "objective"
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

  # objective function: min theta
  objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
  objVal[1] <- 1

  # ======================== #
  # Build LP model structure #
  # ======================== #

  lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
  lp.control(lps, sense = "min")
  set.objfn(lps, objVal)

  # input constraints: -x_di*theta + sum(lambda_j * x_ji) <= 0
  for (xi in 1:nX) {
    add.constraint(lps, xt = c(-1, tech_xmat[, xi]), "<=", rhs = 0)
  }

  # output constraints: sum(lambda_j * y_jr) >= y_dr
  for (yi in 1:nY) {
    add.constraint(lps, xt = c(0, tech_ymat[, yi]), ">=", rhs = 0)
  }

  # technology constraint
  if (returns == "variable") {
    add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)

    if (!convexity) {
      set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
    }
  }

  # ============================= #
  # Solve for each evaluated DMU  #
  # ============================= #

  for (d in 1:eval_dmu) {

    # update coefficient of theta in input constraints
    for (xi in 1:nX) {
      set.mat(lps, xi, 1, -eval_xmat[d, xi])
    }

    # update RHS of output constraints
    for (yi in 1:nY) {
      set.rhs(lps, b = eval_ymat[d, yi], constraints = nX + yi)
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


#' @title Compute Directional Distance Scores
#'
#' @description
#' Solves the directional distance envelopment model for each evaluated DMU. The
#' score is the greatest feasible movement that contracts inputs and expands
#' outputs in the proportions specified by \code{direction}.
#'
#' @param tech_xmat
#' Matrix of inputs for the reference technology, with one row per reference DMU.
#'
#' @param tech_ymat
#' Matrix of outputs for the reference technology, with one row per reference DMU.
#'
#' @param eval_xmat
#' Matrix of inputs for the evaluated DMUs.
#'
#' @param eval_ymat
#' Matrix of outputs for the evaluated DMUs.
#'
#' @param direction
#' A \code{matrix} or \code{data.frame} with \code{N} rows and
#' \code{nX + nY} columns specifying the direction vector for each DMU.
#' Columns \code{1:nX} correspond to input directions (\eqn{g_x}) and
#' columns \code{(nX+1):(nX+nY)} to output directions (\eqn{g_y}).
#'
#' @param convexity
#' If \code{TRUE}, use a convex DEA technology. If \code{FALSE}, use binary
#' intensity variables for an FDH technology.
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param type
#' A \code{character} string specifying the return value:
#' \code{"objective"} for the optimal objective value, or
#' \code{"variables"} for the optimal value of \eqn{\beta}.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds set.rhs set.mat get.objective get.variables
#'
#' @return
#' A one-column matrix with the requested objective value or directional distance
#' for each evaluated DMU.

ddf <- function(
  tech_xmat,
  tech_ymat,
  eval_xmat,
  eval_ymat,
  direction,
  convexity,
  returns,
  type = "objective"
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

  # objective function: max beta
  objVal <- matrix(ncol = 1 + tech_dmu, nrow = 1)
  objVal[1] <- 1

  # =============================== #
  # Build LP model structure #
  # =============================== #

  lps <- make.lp(nrow = 0, ncol = 1 + tech_dmu)
  lp.control(lps, sense = "max")
  set.objfn(lps, objVal)

  # input constraints: g_xi*beta + sum(lambda_j * x_ji) <= x_di
  for (xi in 1:nX) {
    add.constraint(lps, xt = c(1, tech_xmat[, xi]), "<=", rhs = 0)
  }

  # output constraints: -g_yr*beta + sum(lambda_j * y_jr) >= y_dr
  for (yi in 1:nY) {
    add.constraint(lps, xt = c(-1, tech_ymat[, yi]), ">=", rhs = 0)
  }

  # technology constraint
  if (returns == "variable") {
    add.constraint(lprec = lps, xt = c(0, rep(1, tech_dmu)), type = "=", rhs = 1)

    if (!convexity) {
      set.type(lps, columns = 1:tech_dmu + 1, type = c("binary"))
    }
  }

  # lower bound for beta: -Inf (allow negative efficiency)
  set.bounds(lps, lower = c(-Inf, rep(0, tech_dmu)))

  # ============================== #
  # Solve for each evaluated DMU   #
  # ============================== #

  for (d in 1:eval_dmu) {

    # update input constraints: coefficient of beta and RHS
    for (xi in 1:nX) {
      set.mat(lps, xi, 1, G_x[d, xi])
      set.rhs(lps, b = eval_xmat[d, xi], constraints = xi)
    }

    # update output constraints: coefficient of beta and RHS
    for (yi in 1:nY) {
      set.mat(lps, nX + yi, 1, -G_y[d, yi])
      set.rhs(lps, b = eval_ymat[d, yi], constraints = nX + yi)
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

#' @title Compute Output-Oriented Russell Scores
#'
#' @description
#' Solves the output-oriented Russell envelopment model for each evaluated DMU.
#' Unlike a radial model, it allows a different expansion factor for each output
#' and can therefore identify non-proportional output inefficiency.
#'
#' @param tech_xmat
#' Matrix of inputs for the reference technology, with one row per reference DMU.
#'
#' @param tech_ymat
#' Matrix of outputs for the reference technology, with one row per reference DMU.
#'
#' @param eval_xmat
#' Matrix of inputs for the evaluated DMUs.
#'
#' @param eval_ymat
#' Matrix of outputs for the evaluated DMUs.
#'
#' @param convexity
#' If \code{TRUE}, use a convex DEA technology. If \code{FALSE}, use binary
#' intensity variables for an FDH technology.
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param type
#' A \code{character} string specifying the return value:
#' \code{"objective"} for the optimal objective value, or
#' \code{"variables"} for the optimal values of \eqn{\phi_r} (one per output).
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds set.rhs set.mat get.objective get.variables
#'
#' @return
#' A matrix with one row per evaluated DMU. It contains one objective value when
#' \code{type = "objective"}, or one expansion factor per output when
#' \code{type = "variables"}.

rsl_out <- function(
  tech_xmat,
  tech_ymat,
  eval_xmat,
  eval_ymat,
  convexity,
  returns,
  type = "objective"
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

  # objective function: max (1/nY) * sum(phi_r)
  objVal <- matrix(ncol = nY + tech_dmu, nrow = 1)
  objVal[1:nY] <- 1 / nY

  # ======================== #
  # Build LP model structure #
  # ======================== #

  lps <- make.lp(nrow = 0, ncol = tech_dmu + nY)
  lp.control(lps, sense = "max")
  set.objfn(lps, objVal)

  # input constraints: sum(lambda_j * x_ji) <= x_di
  for (xi in 1:nX) {
    add.constraint(lps, xt = c(rep(0, nY), tech_xmat[, xi]), "<=", rhs = 0)
  }

  # output constraints: -y_dr*phi_r + sum(lambda_j * y_jr) >= 0
  for (yi in 1:nY) {
    phi <- rep(0, nY)
    phi[yi] <- -1
    add.constraint(lps, xt = c(phi, tech_ymat[, yi]), ">=", rhs = 0)
  }

  # lower bounds: phi >= 1
  set.bounds(lps, lower = rep(1, nY), columns = 1:nY)

  # technology constraint
  if (returns == "variable") {
    add.constraint(lprec = lps, xt = c(rep(0, nY), rep(1, tech_dmu)), type = "=", rhs = 1)

    if (!convexity) {
      set.type(lps, columns = 1:tech_dmu + nY, type = c("binary"))
    }
  }

  # ============================= #
  # Solve for each evaluated DMU  #
  # ============================= #

  for (d in 1:eval_dmu) {

    # update RHS of input constraints
    for (xi in 1:nX) {
      set.rhs(lps, b = eval_xmat[d, xi], constraints = xi)
    }

    # update coefficient of phi_r in output constraints
    for (yi in 1:nY) {
      set.mat(lps, nX + yi, yi, -eval_ymat[d, yi])
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

#' @title Compute Input-Oriented Russell Scores
#'
#' @description
#' Solves the input-oriented Russell envelopment model for each evaluated DMU.
#' Unlike a radial model, it allows a different contraction factor for each input
#' and can therefore identify non-proportional input inefficiency.
#'
#' @param tech_xmat
#' Matrix of inputs for the reference technology, with one row per reference DMU.
#'
#' @param tech_ymat
#' Matrix of outputs for the reference technology, with one row per reference DMU.
#'
#' @param eval_xmat
#' Matrix of inputs for the evaluated DMUs.
#'
#' @param eval_ymat
#' Matrix of outputs for the evaluated DMUs.
#'
#' @param convexity
#' If \code{TRUE}, use a convex DEA technology. If \code{FALSE}, use binary
#' intensity variables for an FDH technology.
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param type
#' A \code{character} string specifying the return value:
#' \code{"objective"} for the optimal objective value, or
#' \code{"variables"} for the optimal values of \eqn{\theta_i} (one per input).
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds set.rhs set.mat get.objective get.variables
#'
#' @return
#' A matrix with one row per evaluated DMU. It contains one objective value when
#' \code{type = "objective"}, or one contraction factor per input when
#' \code{type = "variables"}.

rsl_inp <- function(
  tech_xmat,
  tech_ymat,
  eval_xmat,
  eval_ymat,
  convexity,
  returns,
  type = "objective"
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
    scores <- matrix(nrow = eval_dmu, ncol = nX)
  }

  # objective function: min (1/nX) * sum(theta_i)
  objVal <- matrix(ncol = nX + tech_dmu, nrow = 1)
  objVal[1:nX] <- 1 / nX

  # ======================== #
  # Build LP model structure #
  # ======================== #

  lps <- make.lp(nrow = 0, ncol = nX + tech_dmu)
  lp.control(lps, sense = "min")
  set.objfn(lps, objVal)

  # input constraints: -x_di*theta_i + sum(lambda_j * x_ji) <= 0
  for (xi in 1:nX) {
    theta <- rep(0, nX)
    theta[xi] <- -1
    add.constraint(lps, xt = c(theta, tech_xmat[, xi]), "<=", rhs = 0)
  }

  # output constraints: sum(lambda_j * y_jr) >= y_dr
  for (yi in 1:nY) {
    add.constraint(lps, xt = c(rep(0, nX), tech_ymat[, yi]), ">=", rhs = 0)
  }

  # upper bounds: theta <= 1
  set.bounds(lps, upper = rep(1, nX), columns = 1:nX)

  # technology constraint
  if (returns == "variable") {
    add.constraint(lprec = lps, xt = c(rep(0, nX), rep(1, tech_dmu)), type = "=", rhs = 1)

    if (!convexity) {
      set.type(lps, columns = 1:tech_dmu + nX, type = c("binary"))
    }
  }

  # ============================= #
  # Solve for each evaluated DMU  #
  # ============================= #

  for (d in 1:eval_dmu) {

    # update coefficient of theta_i in input constraints
    for (xi in 1:nX) {
      set.mat(lps, xi, xi, -eval_xmat[d, xi])
    }

    # update RHS of output constraints
    for (yi in 1:nY) {
      set.rhs(lps, b = eval_ymat[d, yi], constraints = nX + yi)
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

#' @title Compute Weighted Additive Scores
#'
#' @description
#' Solves a weighted additive envelopment model for each evaluated DMU. The model
#' measures input excesses and output shortfalls directly through slacks, with a
#' weighting rule that makes variables with different scales comparable.
#'
#' @param tech_xmat
#' Matrix of inputs for the reference technology, with one row per reference DMU.
#'
#' @param tech_ymat
#' Matrix of outputs for the reference technology, with one row per reference DMU.
#'
#' @param eval_xmat
#' Matrix of inputs for the evaluated DMUs.
#'
#' @param eval_ymat
#' Matrix of outputs for the evaluated DMUs.
#'
#' @param weights Weighting rule for the additive model:
#' \itemize{
#' \item{\code{"wam_mip"}} Measure of Inefficiency Proportions.
#' \item{\code{"wam_nor"}} Normalized Weighted Additive Model.
#' \item{\code{"wam_ram"}} Range Adjusted Measure.
#' \item{\code{"wam_bam"}} Bounded Adjusted Measure.
#' }
#'
#' @param convexity
#' If \code{TRUE}, use a convex DEA technology. If \code{FALSE}, use binary
#' intensity variables for an FDH technology.
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param type
#' A \code{character} string specifying the return value:
#' \code{"objective"} for the optimal objective value, or
#' \code{"variables"} for the optimal slack values.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds set.rhs set.mat get.objective get.variables
#'
#' @return
#' A matrix with one row per evaluated DMU. It contains one objective value when
#' \code{type = "objective"}, or the input and output slacks when
#' \code{type = "variables"}.

wam <- function(
  tech_xmat,
  tech_ymat,
  eval_xmat,
  eval_ymat,
  weights,
  convexity,
  returns,
  type = "objective"
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

  # precompute weights that do not change across DMUs
  objVal <- matrix(ncol = nX + nY + tech_dmu, nrow = 1)

  if (weights == "wam_nor") {
    objVal[1:(nX + nY)] <- c(1 / apply(eval_xmat, 2, sd), 1 / apply(eval_ymat, 2, sd))
  } else if (weights == "wam_ram") {
    xranges <- apply(eval_xmat, 2, max) - apply(eval_xmat, 2, min)
    yranges <- apply(eval_ymat, 2, max) - apply(eval_ymat, 2, min)
    objVal[1:(nX + nY)] <- c(1 / ((nX + nY) * xranges), 1 / ((nX + nY) * yranges))
  }

  # ======================== #
  # Build LP model structure #
  # ======================== #

  lps <- make.lp(nrow = 0, ncol = nX + nY + tech_dmu)
  lp.control(lps, sense = "max")
  set.objfn(lps, objVal)

  # input constraints: s_xi + sum(lambda_j * x_ji) = x_di
  for (xi in 1:nX) {
    x_slack <- rep(0, nX)
    x_slack[xi] <- 1
    slacks <- c(x_slack, rep(0, nY))

    add.constraint(lps, xt = c(slacks, tech_xmat[, xi]), "=", rhs = 0)
  }

  # output constraints: -s_yr + sum(lambda_j * y_jr) = y_dr
  for (yi in 1:nY) {
    y_slack <- rep(0, nY)
    y_slack[yi] <- -1
    slacks <- c(rep(0, nX), y_slack)

    add.constraint(lps, xt = c(slacks, tech_ymat[, yi]), "=", rhs = 0)
  }

  # technology constraint
  if (returns == "variable") {
    add.constraint(lprec = lps, xt = c(rep(0, nX + nY), rep(1, tech_dmu)), type = "=", rhs = 1)

    if (!convexity) {
      set.type(lps, columns = 1:tech_dmu + (nX + nY), type = c("binary"))
    }
  }

  # ============================= #
  # Solve for each evaluated DMU  #
  # ============================= #

  for (d in 1:eval_dmu) {

    # update objective function for DMU-specific weights
    if (weights == "wam_mip") {
      objVal[1:(nX + nY)] <- c(1 / eval_xmat[d, ], 1 / eval_ymat[d, ])
      set.objfn(lps, objVal)

    } else if (weights == "wam_bam") {
      p1 <- eval_xmat[d, ] - apply(eval_xmat, 2, min)
      p2 <- apply(eval_ymat, 2, max) - eval_ymat[d, ]
      objVal[1:(nX + nY)] <- c(1 / ((nX + nY) * p1), 1 / ((nX + nY) * p2))
      set.objfn(lps, objVal)
    }

    # update RHS of input constraints
    for (xi in 1:nX) {
      set.rhs(lps, b = eval_xmat[d, xi], constraints = xi)
    }

    # update RHS of output constraints
    for (yi in 1:nY) {
      set.rhs(lps, b = eval_ymat[d, yi], constraints = nX + yi)
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

#' @title Compute Efficiency Scores
#'
#' @description
#' Evaluates each decision-making unit (DMU) against the production technology
#' constructed by ACES or RF-ACES. The stored input-output reference points are
#' enveloped under the selected returns-to-scale assumption. The function
#' supports proportional, directional, non-radial, and slack-based measures, so
#' the score can reflect different notions of technical inefficiency. If
#' \code{method = NULL}, the standard model for the supplied object is used.
#'
#' @param eval_data
#' A \code{data.frame} or \code{matrix} containing the DMUs to evaluate.
#'
#' @param x
#' Column indexes of input variables in \code{eval_data}.
#'
#' @param y
#' Column indexes of output variables in \code{eval_data}.
#'
#' @param relevant
#' If \code{TRUE}, compute the efficiency result with only the original inputs
#' represented in at least one selected basis function. Inputs used through an
#' interaction are also retained.
#'
#' @param object
#' An \code{aces} or \code{rf_aces} object.
#'
#' @param method
#' Fitted model used to define the technology. When \code{NULL}, use \code{"aces"}
#' for an \code{aces} object and \code{"rf_aces"} for an \code{rf_aces} object:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_cubic"}}: Cubic Smooth Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_quintic"}}: Quintic Smooth Adaptive Constrained Enveloping Splines.
#' \item{\code{"rf_aces"}}: Random Forest Adaptive Constrained Enveloping Splines.
#' \item{\code{"rf_aces_cubic"}}: Cubic-smoothed RF-ACES.
#' \item{\code{"rf_aces_quintic"}}: Quintic-smoothed RF-ACES.
#' }
#'
#' @param measure
#' Efficiency measure to compute:
#' \itemize{
#' \item{\code{"rad_out"}: Proportional expansion of all outputs while inputs
#' remain fixed.}
#' \item{\code{"rad_inp"}: Proportional contraction of all inputs while outputs
#' remain fixed.}
#' \item{\code{"ddf"}: Simultaneous input contraction and output expansion along
#' the supplied direction.}
#' \item{\code{"rsl_out"}: Separate proportional adjustments for each output.}
#' \item{\code{"rsl_inp"}: Separate proportional adjustments for each input.}
#' \item{\code{"wam_mip"}: Slack-based Measure of Inefficiency Proportions.}
#' \item{\code{"wam_nor"}: Slack-based normalized weighted additive measure.}
#' \item{\code{"wam_ram"}: Slack-based range-adjusted measure.}
#' \item{\code{"wam_bam"}: Slack-based bounded adjusted measure.}
#' \item{\code{"rf_aces_rad_out"}: Output-oriented radial scores based directly
#' on RF-ACES predictions. Available only for \code{rf_aces} objects.}
#' }
#'
#' @param returns
#' Returns-to-scale assumption used to construct the reference technology.
#' Choose \code{"constant"} for a conical technology or \code{"variable"} for a
#' convex technology with a convexity constraint.
#'
#' @param direction
#' Direction vectors used when \code{measure = "ddf"}. Supply a matrix or data
#' frame with one row per evaluated DMU. Input directions must come first,
#' followed by output directions. The direction determines the relative input
#' contractions and output expansions represented by the score.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{espana2024rf}{aces} \cr \cr
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
#' A data frame with one row and one efficiency score for each evaluated DMU.
#' Row names are copied from \code{eval_data}; the column name identifies the
#' returns-to-scale assumption and efficiency measure.
#'
#' @details
#' The \code{relevant} argument applies a structural rule, not an importance
#' threshold. For an ACES object, the retained inputs come from the selected
#' basis functions of \code{method}. For an RF-ACES object, the function uses
#' the union of the inputs selected by any learner. The reduced input set is
#' used only in the envelopment problem; the model is not refitted. In
#' particular, \code{aces_varimp()} scores do not determine this set. The
#' argument has no effect on \code{measure = "rf_aces_rad_out"}, which is based
#' directly on prediction-to-observation ratios.
#'
#' @export

get_scores <- function(
  eval_data,
  x,
  y,
  relevant = FALSE,
  object,
  method = NULL,
  measure = "rad_out",
  returns = "variable",
  direction = NULL
  ) {

  # default method based on object class
  if (is.null(method)) {
    method <- if (inherits(object, "rf_aces")) "rf_aces" else "aces"
  }

  # handle errors:
  display_errors_scores(
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
    if (inherits(object, "rf_aces")) {
      # aggregate knots across all trees
      all_xi <- unique(unlist(lapply(object[["forest"]], function(tree) {
        tree[["methods"]][[method]][["knots"]]$xi
      })))
    } else {
      all_xi <- unique(object[["methods"]][[method]][["knots"]]$xi)
    }

    # variable degree
    xi_degree <- object[["control"]][["xi_degree"]]
    colnames(xi_degree) <- names(object[["control"]][["kn_grid"]])

    # check participating variables
    participating_vars <- intersect(xi_degree[1, ], all_xi)

    # names of participant variables
    participating_vars_names <- colnames(xi_degree)[participating_vars]

    # extract the original variables from interaction names such as "x1:x2"
    split_vars <- unique(unlist(strsplit(
      participating_vars_names,
      ":",
      fixed = TRUE
    )))

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

  # ============ #
  # Compute scores
  # ============ #

  if (measure == "rf_aces_rad_out") {
    # RF-ACES specific: predict with the forest and compute ratio
    y_hat_point <- rf_aces_predict(
      object = object,
      eval_data = eval_data,
      x = x,
      method = method
    )

    # ratio y_hat / y_obs
    ratios <- as.data.frame(y_hat_point / eval_data[, y])

    # score: minimum ratio by row
    scores <- apply(ratios, 1, min)

  } else if (measure == "rad_out") {
    scores <- rad_out(
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )
  } else if (measure == "rad_inp") {
    scores <- rad_inp(
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )
  } else if (measure == "ddf") {
    scores <- ddf(
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
    scores <- rsl_out(
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )
  } else if (measure == "rsl_inp") {
    scores <- rsl_inp(
      tech_xmat = tech_xmat,
      tech_ymat = tech_ymat,
      eval_xmat = eval_xmat,
      eval_ymat = eval_ymat,
      convexity = TRUE,
      returns = returns,
      type = "objective"
    )
  } else if (grepl("wam", measure)) {
    scores <- wam(
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
