% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ACES.R
\name{set_knots_grid}
\alias{set_knots_grid}
\title{Set the Grid of Knots}
\usage{
set_knots_grid(data, nX, inps, kn_grid)
}
\arguments{
\item{data}{A \code{matrix} containing the variables in the model.}

\item{nX}{Number of inputs (with interactions and netputs).}

\item{inps}{Number of original inputs (without interactions).}

\item{kn_grid}{Grid of knots to perform ACES.}
}
\value{
A \code{list} with the available vector of knots for each variable.
}
\description{
This function sets the grid of knots to perform ACES.
}
