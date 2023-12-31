% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_coefficients.R
\name{concavity_matrix}
\alias{concavity_matrix}
\title{Build Constraint Matrix for Enforcing Concavity}
\usage{
concavity_matrix(Bp_list, N)
}
\arguments{
\item{Bp_list}{A \code{list} containing the set of basis functions by input.}

\item{N}{Sample size.}
}
\value{
A \code{matrix} to enforce concavity constraint on the estimated function.
}
\description{
This function constructs the constraint matrix (often denoted as matrix "A") to enforce concavity.
}
