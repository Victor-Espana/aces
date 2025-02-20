% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ACES.R
\name{set_knots_grid}
\alias{set_knots_grid}
\title{Define the Grid of Knots}
\usage{
set_knots_grid(data, n_input_1, n_input_2, kn_grid, quick_aces, dea_scores)
}
\arguments{
\item{data}{A \code{matrix} containing the variables in the model.}

\item{n_input_1}{Number of inputs (excluding interactions) and contextual variables.}

\item{n_input_2}{Number of inputs (including interactions) and contextual variables.}

\item{kn_grid}{A \code{list} providing a custom grid of knots for each variable. If not supplied, the function automatically generates a grid of knots for each variable based on the data.}

\item{quick_aces}{A \code{logical} indicating whether to use the fast version of ACES.}

\item{dea_scores}{A \code{matrix} containing DEA-VRS efficiency scores, calculated using an output-oriented radial model.}
}
\value{
A \code{list} where each element contains the generated or provided vector of knots for the corresponding variable in the model.
}
\description{
This function defines the grid of knots to perform Adaptive Constrained Enveloping Splines (ACES).
}
