#' @title Estimate Coefficients for Adaptive Constrained Enveloping Splines (ACES)
#'
#' @description
#'
#' This function estimates a vector of coefficients for the Adaptive Constrained Enveloping Splines (ACES) model. These coefficients enforce the desired shape properties, such as monotonicity and concavity to fit the data within the ACES framework.
#'
#' @param B
#' A \code{matrix} of linear basis functions derived from input variables.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_scores
#' A \code{matrix} containing DEA-VRS efficiency scores, calculated using an output-oriented radial model. For models with multiple outputs, each column corresponds to the scores for one specific output.
#'
#' @param fdh_scores
#' A \code{matrix} containing FDH efficiency scores, calculated using an output-oriented radial model. For models with multiple outputs, each column corresponds to the scores for one specific output.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity.
#'
#' @return
#'
#' A \code{vector} containing the estimated coefficients for the ACES model.

estimate_coefficients <- function(
  B,
  y_obs,
  dea_scores,
  fdh_scores,
  it_list,
  Bp_list,
  shape
) {
  coefs <- estimate_coefficients_high_speed(
    B = B,
    y_obs = y_obs,
    dea_scores = dea_scores,
    fdh_scores = fdh_scores,
    it_list = it_list,
    Bp_list = Bp_list,
    shape = shape
  )

  return(coefs)
}

#' @title Estimate Coefficients in Adaptive Constrained Enveloping Splines: Envelopment Version
#'
#' @description
#'
#' This function estimates a vector of coefficients for the Adaptive Constrained Enveloping Splines (ACES) model in its envelopment version. It solves a mathematical programming problem to determine the coefficients that best fit the observed data while respecting certain shape constraints such as monotonicity and/or concavity. The model aims to envelop the data points, particularly when estimating production functions or efficiency frontiers.
#'
#' @param B
#' A \code{matrix} of linear basis functions derived from input variables.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_scores
#' A \code{matrix} containing DEA-VRS efficiency scores, calculated using an output-oriented radial model. For models with multiple outputs, each column corresponds to the scores for one specific output.
#'
#' @param fdh_scores
#' A \code{matrix} containing FDH efficiency scores, calculated using an output-oriented radial model. For models with multiple outputs, each column corresponds to the scores for one specific output.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity.
#'
#' @return
#'
#' A \code{vector} containing the estimated coefficients for the ACES model.

estimate_coefficients_envelopment <- function(
  B,
  y_obs,
  dea_scores,
  fdh_scores,
  it_list,
  Bp_list,
  shape
) {
  # monotonicity
  mono <- shape[["mono"]]

  # concavity
  conc <- shape[["conc"]]

  # dea weights
  deaw <- as.matrix(1 / dea_scores)

  # sample size
  N <- nrow(B)

  # number of basis functions (inputs)
  p <- ncol(B)

  # number of outputs
  nY <- ncol(y_obs)

  # number of inputs
  nX <- length(Bp_list)

  # structure of coefficients
  coefs <- matrix(NA, nrow = p, ncol = nY)

  for (out in 1:nY) {
    # envelopment indicator
    env_ind <- fdh_scores[, out] < 1 + 1e-5
    num_env <- sum(env_ind)

    # select an output
    y_ind <- y_obs[env_ind, out]

    # ==================================================== #
    # vars: c(coef_0, coef_1, ..., coef_P, e_1, ... , e_n) #
    # ==================================================== #

    objVal <- c(rep(0, p), deaw[env_ind, out])

    # ==================================================== #
    # A: envelopment + concavity + monotonicity            #
    # ==================================================== #

    # envelopment: y_hat - e = y
    EMat <- cbind(B[env_ind, ], diag(rep(-1, num_env), num_env))

    # matrix of coefficients
    Amat <- rbind(EMat)

    # generate non-decreasing monotonicity matrix
    if (mono) {
      MMat <- monotonocity_matrix(
        it_list = it_list,
        Bp_list = Bp_list,
        N = num_env
      )

      # add MMat to Amat
      Amat <- rbind(Amat, MMat)
    }

    # generate concavity matrix
    if (conc) {
      CMat <- concavity_matrix(
        Bp_list = Bp_list,
        N = num_env
      )

      # add CMat to Amat
      Amat <- rbind(Amat, CMat)
    }

    # ==================================================== #
    # b: envelopment + concavity + monotonicity            #
    # ==================================================== #

    bvec <- c(y_ind, rep(0, nrow(Amat) - num_env))

    # ==================================================== #
    # directions of inequalities                           #
    # ==================================================== #

    dirs <- c(rep("==", num_env), rep(">=", nrow(Amat) - num_env))

    # ==================================================== #
    # lower and upper bounds                               #
    # ==================================================== #

    bnds <- list(lower = list(ind = 1:p, val = rep(-Inf, p)))

    # ==================================================== #
    # solution of the optimization problem                 #
    # ==================================================== #

    sols <- Rglpk_solve_LP(
      obj = objVal,
      mat = Amat,
      dir = dirs,
      rhs = bvec,
      bounds = bnds
    )

    if (sols[["status"]] == 0) {
      coefs[, out] <- sols$solution[1:p]
    } else {
      coefs[, out] <- rep(0, p)
    }
  }

  return(coefs)
}

#' @title Estimate Coefficients in ACES: High-Speed Envelopment Version
#'
#' @description
#'
#' Mathematically equivalent reformulation of \code{estimate_coefficients_envelopment}.
#' The N error variables are eliminated algebraically: since e = B coefs - y with
#' e >= 0, minimizing sum(w * e) is equivalent to minimizing (w'B) coefs subject to
#' B coefs >= y plus the shape constraints. The resulting LP has only p variables.
#' It is solved through its DUAL (p rows, N + k columns), which is much faster
#' because the simplex bases are p x p instead of (N + k) x (N + k). The optimal
#' primal coefficients are recovered from the dual values of the p equality
#' constraints. A self-check validates feasibility and strong duality; on any
#' failure the reduced primal is solved with GLPK as a fallback.
#'
#' Note: when the LP has multiple optimal solutions, the returned vertex may
#' differ from the legacy implementation, with identical objective value and
#' identical fitted values.
#'
#' @param B
#' A \code{matrix} of linear basis functions derived from input variables.
#'
#' @param y_obs
#' A \code{matrix} of the observed output data.
#'
#' @param dea_scores
#' A \code{matrix} containing DEA-VRS efficiency scores.
#'
#' @param fdh_scores
#' A \code{matrix} containing FDH efficiency scores.
#'
#' @param it_list
#' A \code{list} containing the set of intervals by input.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions by input.
#'
#' @param shape
#' A \code{list} indicating whether to impose monotonicity and/or concavity.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint get.objective
#'
#' @return
#'
#' A \code{matrix} containing the estimated coefficients for the ACES model.

estimate_coefficients_high_speed <- function(
  B,
  y_obs,
  dea_scores,
  fdh_scores,
  it_list,
  Bp_list,
  shape
) {
  # monotonicity
  mono <- shape[["mono"]]

  # concavity
  conc <- shape[["conc"]]

  # dea weights
  deaw <- as.matrix(1 / dea_scores)

  # number of basis functions (columns of B)
  p <- ncol(B)

  # number of outputs
  nY <- ncol(y_obs)

  # structure of coefficients
  coefs <- matrix(NA, nrow = p, ncol = nY)

  # shape constraint matrix (k x p), without error-variable padding
  SMat <- NULL

  if (mono) {
    SMat <- rbind(SMat, monotonocity_matrix_hs(it_list = it_list, Bp_list = Bp_list))
  }

  if (conc) {
    SMat <- rbind(SMat, concavity_matrix_hs(Bp_list = Bp_list))
  }

  k <- if (is.null(SMat)) 0L else nrow(SMat)

  for (out in 1:nY) {
    # envelopment indicator
    env_ind <- fdh_scores[, out] < 1 + 1e-5
    num_env <- sum(env_ind)

    # data restricted to enveloped observations
    y_env <- y_obs[env_ind, out]
    w_env <- deaw[env_ind, out]
    B_env <- B[env_ind, , drop = FALSE]

    # ==================================================== #
    # reduced primal:                                      #
    #   min  (w'B) beta                                    #
    #   s.t. [B_env; SMat] beta >= [y_env; 0], beta free   #
    # ==================================================== #

    cvec <- as.vector(crossprod(B_env, w_env))
    Amat <- if (k > 0) rbind(B_env, SMat) else B_env
    bvec <- c(y_env, rep(0, k))

    # ==================================================== #
    # dual:                                                #
    #   max  b'lambda                                      #
    #   s.t. A'lambda = c, lambda >= 0                     #
    # (p equality rows; num_env + k non-negative columns)  #
    # ==================================================== #

    beta <- tryCatch({
      lps <- make.lp(nrow = 0, ncol = num_env + k)
      invisible(lp.control(lps, sense = "max"))
      set.objfn(lps, bvec)

      # r-th dual constraint = r-th column of Amat
      for (r in 1:p) {
        add.constraint(lps, xt = Amat[, r], type = "=", rhs = cvec[r])
      }

      status <- solve(lps)

      if (status != 0) {
        NULL
      } else {
        # primal coefficients = dual values of the p equality constraints
        duals <- lpSolveAPI::get.sensitivity.rhs(lps)[["duals"]][1:p]
        obj_dual <- get.objective(lps)

        # self-check (covers lp_solve dual sign conventions):
        # the candidate must be primal-feasible and satisfy strong duality
        sel <- NULL

        for (cand in list(duals, -duals)) {
          feas_gap <- min(Amat %*% cand - bvec)
          dual_gap <- abs(sum(cvec * cand) - obj_dual)

          tol_feas <- -1e-6 * (1 + max(abs(bvec)))
          tol_dual <- 1e-6 * (1 + abs(obj_dual))

          if (feas_gap >= tol_feas && dual_gap <= tol_dual) {
            sel <- cand
            break
          }
        }

        sel
      }
    }, error = function(e) NULL)

    # fallback: reduced primal solved with GLPK
    if (is.null(beta)) {
      beta <- estimate_coefficients_reduced(
        Amat = Amat,
        bvec = bvec,
        cvec = cvec,
        p = p
      )
    }

    coefs[, out] <- beta
  }

  return(coefs)
}

#' @title Solve the Reduced Primal LP of ACES with GLPK
#'
#' @description
#' Fallback solver for \code{estimate_coefficients_high_speed}: solves the
#' reduced primal (p variables) directly with GLPK.
#'
#' @param Amat
#' Constraint \code{matrix}: envelopment rows plus shape rows.
#'
#' @param bvec
#' Right-hand side \code{vector}.
#'
#' @param cvec
#' Objective \code{vector} (w'B).
#'
#' @param p
#' Number of coefficients.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return
#' A \code{vector} of estimated coefficients.

estimate_coefficients_reduced <- function(
  Amat,
  bvec,
  cvec,
  p
) {
  sols <- Rglpk_solve_LP(
    obj = cvec,
    mat = Amat,
    dir = rep(">=", nrow(Amat)),
    rhs = bvec,
    bounds = list(lower = list(ind = 1:p, val = rep(-Inf, p)))
  )

  if (sols[["status"]] == 0) {
    return(sols$solution[1:p])
  } else {
    return(rep(0, p))
  }
}
