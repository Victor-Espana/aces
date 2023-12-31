% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backward_algorithm.R
\name{mapping_basis}
\alias{mapping_basis}
\title{Mapping between Basis Function ID and Matrix Column Index by Variable}
\usage{
mapping_basis(Bp_list, nX)
}
\arguments{
\item{Bp_list}{A \code{list} containing the set of basis functions by input.}

\item{nX}{Number of inputs.}
}
\value{
A \code{matrix} that establishes the connection between the basis function ID and its column index in the basis function matrix by variable.
}
\description{
This function creates a mapping between the ID of a basis function and its corresponding column index in the basis function matrix (B matrix) by variable.
}
