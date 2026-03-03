#' @title Super-Efficiency Ranking for ACES models
#'
#' @description
#' This function calculates super-efficiency scores by performing a Leave-One-Out (LOO) approach. For each DMU, the technology is re-defined by re-estimating the frontier coefficients excluding that specific unit. The excluded unit is then evaluated  against this peer-frontier using standard DEA measures. This approach allows for breaking ties among efficient units and provides a rigorous ranking.
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
#' \item{\code{rad_out}} The output-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{rad_inp}} The input-oriented radial measure proposed by \insertCite{banker1984;textual}{aces}.
#' \item{\code{ddf}}     The directional distance function proposed by \insertCite{chambers1998;textual}{aces}.
#' \item{\code{rsl_out}} The output-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{rsl_inp}} The input-oriented Russell measure proposed by \insertCite{fare1978;textual}{aces}.
#' \item{\code{wam}}     A weighted additive model.
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
#' Direction of the vector to project on the frontier. Only applied if \code{measure = "ddf"}. It can be a \code{matrix} of custom orientation coefficients or one of the following predefined strings:
#' \itemize{
#'   \item{\code{"mean"}}: Projection vector given by the average value of inputs and outputs of all DMUs. Applied if \code{direction = NULL}.
#'   \item{\code{"briec"}}: Projection vector given by the observed values of the evaluated DMU.
#' }
#'
#' @param weights
#' Weights for the additive model. Only applied if \code{measure = "wam"}.
#' \itemize{
#' \item{\code{"ONE"}} Weighted Additive Model proposed by \insertCite{charnes1985;textual}{aces}. Applied if \code{weights = NULL}.
#' \item{\code{"MIP"}} Measure of Inefficiency Proportions proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{"NOR"}} Normalized Weighted Additive Model proposed by \insertCite{lovell1995;textual}{aces}.
#' \item{\code{"RAM"}} Range Adjusted Measure proposed by \insertCite{cooper1999;textual}{aces}.
#' \item{\code{"BAM"}} Bounded Adjusted Measure proposed by \insertCite{cooper2011;textual}{aces}.
#' }
#'
#' @param digits
#' Number of decimal places to which efficiency scores are rounded.
#'
#' @references
#'
#' \insertRef{espana2024}{aces} \cr \cr
#' \insertRef{banker1984}{aces} \cr \cr
#' \insertRef{chambers1998}{aces} \cr \cr
#' \insertRef{fare1978}{aces} \cr \cr
#' \insertRef{charnes1985}{aces} \cr \cr
#' \insertRef{cooper1999}{aces} \cr \cr
#' \insertRef{lovell1995}{aces} \cr \cr
#' \insertRef{cooper2011}{aces}
#'
#' @return
#' A data frame with DMU names and their corresponding super-efficiency scores.
#'
#' @export
ranking_aces <- function (
    eval_data,
    x,
    y,
    relevant = FALSE,
    object,
    method = "aces",
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

  # extract model data & structure
  model_obj <- object[["methods"]][[method]]

  if (method == "aces") {

    # list submodels
    submodels <- object[["methods"]][["aces"]][["aces_submodels"]]

    # find best submodel based on gcv
    best_gcv <- object[["methods"]][["aces"]][["GCV"]]

    # identify the optimal submodel to get its specific Bp_list
    best_idx  <- which(sapply(submodels, function(m) abs(m$GCV - best_gcv) < 1e-9))[1]

    # bp_list
    Bp_list   <- submodels[[best_idx]][["Bp_list"]]

    # build it_list
    it_list   <- set_intervals(Bp_list)

  } else {

    # for smoothed models, structural info comes from the knots data frame
    knots_df  <- model_obj[["knots"]]

  }

}
