% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_aces_algorithm.R
\name{rf_aces_algorithm}
\alias{rf_aces_algorithm}
\title{Algorithm of Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).}
\usage{
rf_aces_algorithm(
  data,
  inps,
  outs,
  nets,
  y_type,
  model_type,
  error_type,
  degree,
  metric,
  shape,
  nterms,
  nvars,
  err_red,
  hd_cost,
  minspan,
  endspan,
  kn_grid,
  wc,
  wq
)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} containing the variables in the model.}

\item{inps}{Column indexes of input variables in \code{data}.}

\item{outs}{Column indexes of output variables in \code{data}.}

\item{nets}{Column indexes of netput (outputs evaluated as inputs) variables in \code{data}. These variables are treated as inputs during prediction computation and as outputs when computing efficiency scores.}

\item{y_type}{A \code{character} string that determines the prediction approach for \code{y}.}

\item{model_type}{A \code{character} string specifying the nature of the production frontier that the function estimates.}

\item{error_type}{A \code{character} string specifying the error structure that the function will use when fitting the model.}

\item{degree}{Maximum degree of interaction between variables.}

\item{metric}{A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance.}

\item{shape}{A \code{list} indicating whether to impose monotonicity and/or concavity and/or passing through the origin.}

\item{nterms}{Maximum number of terms created before pruning.}

\item{nvars}{An \code{integer} indicating the number of variables randomly chosen at each split.}

\item{err_red}{Minimum reduced error rate for the addition of a new pair of 1-degree basis functions.}

\item{hd_cost}{Minimum percentage of improvement over the best 1 degree basis function to incorporate a higher degree basis function.}

\item{minspan}{Minimum number of observations between two adjacent knots.}

\item{endspan}{Minimum number of observations before the first and after the final knot.}

\item{kn_grid}{Grid design for knots placement in ACES.}

\item{wc}{Hyperparameter for side knot distances in the cubic smoothing procedure.}

\item{wq}{Hyperparameter for the side knot distances in the quintic smoothing procedure.}
}
\value{
An \code{aces} object.
}
\description{
This function implements the Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) algorithm, which estimates a production frontier satisfying classical production theory axioms like monotonicity and concavity. It offers both stochastic and envelopment versions. These estimations are based on the adaptation of the Multivariate Adaptive Regression Splines (MARS) technique developed by \insertCite{friedman1991;textual}{aces} and the Random Forest methodology \insertCite{breiman2001}{aces}. For details, see \insertCite{espana2024;textual}{aces}.
}
\references{
\insertRef{friedman1991}{aces} \cr \cr
\insertRef{breiman2001}{aces}
}
