% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efficiency_scores.R
\name{rsl_out}
\alias{rsl_out}
\title{The output-oriented Russell model in the envelopment format}
\usage{
rsl_out(tech_xmat, tech_ymat, eval_xmat, eval_ymat, convexity, returns)
}
\arguments{
\item{tech_xmat}{A \code{matrix} containing the observed input variables used to define the technology.}

\item{tech_ymat}{A \code{matrix} containing the observed output variables used to define the technology.}

\item{eval_xmat}{A \code{matrix} containing the input data of the Decision-Making Units (DMUs) to be evaluated.}

\item{eval_ymat}{A \code{matrix} containing the output data of the Decision-Making Units (DMUs) to be evaluated.}

\item{convexity}{A \code{logical} value indicating if a convex technology is assumed.}

\item{returns}{Type of returns to scale.}
}
\value{
A \code{vector} of \code{"numeric"} scores computed through the output-oriented Russell model.
}
\description{
This function computes efficiency scores through the output-oriented Russell model in the envelopment format.
}
