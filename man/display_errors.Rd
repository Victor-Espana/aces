% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/display_errors.R
\name{display_errors}
\alias{display_errors}
\title{Error Messaging in Adaptive Constrained Enveloping Splines (ACES) functions.}
\usage{
display_errors(
  caller,
  data,
  x,
  y,
  y_type,
  model_type,
  error_type,
  nvars,
  degree,
  metric,
  nterms,
  err_red,
  hd_cost,
  minspan,
  endspan,
  kn_grid,
  d,
  wc,
  wq,
  object,
  measure,
  convexity,
  returns,
  direction,
  digits
)
}
\arguments{
\item{caller}{A \code{character} string specifying the function that calls \code{display_errors}.}

\item{data}{A \code{data.frame} or \code{matrix} containing the variables in the model.}

\item{x}{Column indexes of input variables in \code{data}.}

\item{y}{Column indexes of output variables in \code{data}.}

\item{y_type}{A \code{character} string that determines the prediction approach for \code{y}. It can be either: \code{"ind"} or \code{"all"}.}

\item{model_type}{A \code{character} string specifying the nature of the production frontier that the function estimates. It can be either: \code{"env"} or \code{"sto"}.}

\item{error_type}{A \code{character} string specifying the error structure that the function will use when fitting the model. It can be either: \code{"add"} or \code{"mul"}.}

\item{nvars}{An \code{integer} indicating the number of variables randomly chosen at each split in RF-ACES.}

\item{degree}{A \code{list} with input indexes for interaction of variables, or a \code{numeric} specifying the maximum degree of interaction.}

\item{metric}{A \code{character} string specifying the lack-of-fit criterion to evaluate the model performance. Options are: \code{"mae"}, \code{"mape"}, \code{"mse"}, \code{"rmse"}, \code{"nrmse1"} or \code{"nrmse2"}.}

\item{nterms}{A positive \code{integer} specifying the maximum number of terms created before pruning.}

\item{err_red}{A \code{numeric} value specifying the minimum reduced error rate for the addition of a new pair of 1-degree basis functions.}

\item{hd_cost}{A \code{numeric} value specifying the minimum percentage of improvement over the best 1-degree basis function to incorporate a higher degree basis function.}

\item{minspan}{A \code{numeric} value specifying the minimum number of observations between two adjacent knots. Options are: \code{"-2"}, \code{"-1"} or \code{"m"}.}

\item{endspan}{A \code{numeric} value specifying the minimum number of observations before the first and after the final knot. Options are: \code{"-2"}, \code{"-1"} or \code{"m"}.}

\item{kn_grid}{Grid of knots to perform ACES. It can be: \code{-1} or a \code{list}.}

\item{d}{A positive \code{numeric} value specifying the Generalized Cross Validation (GCV) penalty per knot.}

\item{wc}{A numeric value used for the cubic smoothing procedure.}

\item{wq}{A numeric value used for the quintic smoothing procedure.}

\item{object}{An \code{aces} object.}

\item{measure}{Mathematical programming model to compute the efficiency scores.}

\item{convexity}{A \code{logical} value indicating if a convex technology is assumed.}

\item{returns}{Type of returns to scale for computing the efficiency scores.}

\item{direction}{Vector direction for DMU projection in Directional Distance Function when computing the efficiency scores.}

\item{digits}{Number of digits to round efficiency scores.}
}
\value{
This function return error messages if hyperparameters are incorrectly specified.
}
\description{
This function displays error messages if hyperparameters are incorrectly specified in ACES functions.
}
