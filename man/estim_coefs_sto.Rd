% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_coefficients_linear.R
\name{estim_coefs_sto}
\alias{estim_coefs_sto}
\title{Estimate Coefficients in Adaptive Constrained Enveloping Splines Fitting: Stochastic Version}
\usage{
estim_coefs_sto(B, y_obs, it_list, Bp_list, mono, conc, ptto)
}
\arguments{
\item{B}{A \code{matrix} of linear basis functions.}

\item{y_obs}{A \code{matrix} of the observed output data.}

\item{it_list}{A \code{list} containing the set of intervals by input.}

\item{Bp_list}{A \code{list} containing the set of basis functions by input.}

\item{mono}{A \code{logical} value indicating whether to enforce the constraint of non-decreasing monotonicity in the estimator.}

\item{conc}{A \code{logical} value indicating whether to enforce the constraint of concavity in the estimator.}

\item{ptto}{A \code{logical} value indicating whether the estimator should satisfy f(0) = 0.}
}
\value{
A \code{vector} of estimated coefficients.
}
\description{
This function solves a Quadratic Programming problem to obtain a vector of coefficients that impose the required shaped on the Adaptive Constrained Enveloping Splines estimator in the stochastic version
}
