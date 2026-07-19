#' @title Estimate ACES Coefficients
#'
#' @description
#'
#' Solves the reduced linear program for the coefficients of a linear ACES model.
#'
#' Error variables are removed analytically so that only the basis coefficients
#' need to be optimized. The function first solves the smaller dual problem and
#' recovers the primal coefficients from its solution. If the recovered solution
#' does not pass the feasibility checks, the reduced primal problem is solved
#' directly with GLPK.
#'
#' @param B
#' Matrix of evaluated linear basis functions.
#'
#' @param y_obs
#' Matrix of observed outputs.
#'
#' @param dea_scores
#' A matrix of output-oriented DEA-VRS scores, with one column per output.
#'
#' @param fdh_scores
#' A matrix of output-oriented FDH scores, with one column per output.
#'
#' @param it_list
#' Intervals and active basis functions, grouped by input.
#'
#' @param Bp_list
#' Basis functions grouped by input.
#'
#' @param shape
#' A list with logical elements \code{mono} and \code{conc}.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint get.objective
#'
#' @return
#'
#' A matrix of estimated coefficients, with one column per output.

estimate_coefficients <- function(
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
    SMat <- rbind(SMat, monotonicity_matrix(it_list = it_list, Bp_list = Bp_list))
  }

  if (conc) {
    SMat <- rbind(SMat, concavity_matrix(Bp_list = Bp_list))
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

#' @title Solve the Reduced ACES Primal Problem
#'
#' @description
#' Solves the reduced primal problem directly with GLPK. This is the fallback
#' used by \code{estimate_coefficients}.
#'
#' @param Amat
#' Constraint matrix containing the envelopment and shape rows.
#'
#' @param bvec
#' Right-hand-side vector.
#'
#' @param cvec
#' Objective vector.
#'
#' @param p
#' Number of coefficients to return.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#'
#' @return
#' A numeric vector of estimated coefficients.

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
