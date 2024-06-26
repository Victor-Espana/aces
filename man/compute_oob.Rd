% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_aces_algorithm.R
\name{compute_oob}
\alias{compute_oob}
\title{Compute Out-of-Bag Error in Random Forest Adaptive Constrained Enveloping Splines (RF-ACES).}
\usage{
compute_oob(data, y, oob_idxs, oob_pred, model)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} containing the variables in the model.}

\item{y}{Column indexes of output variables in \code{data}.}

\item{oob_pred}{A \code{list} with the Out-Of-Bag predictions for each model.}

\item{model}{Index indicating the model of the ensamble.}

\item{x}{Column indexes of input variables in \code{data}.}

\item{object}{An \code{aces} that has been trained using bagging as in Random Forest.}

\item{oob_data}{A \code{data.frame} containing the out-of-bag samples. These are the sample that were not used during the training of the ACES model and will be used to evaluate its performance.}

\item{inb_idxs}{A \code{vector} of integers representing the indices of the samples that are in the bag.}

\item{method}{Model for prediction:
\itemize{
\item{\code{"aces_forward"}}: Random Forest Forward Adaptive Constrained Enveloping Splines model.
\item{\code{"aces_cubic"}}: Random Forest Cubic Smoothed Adaptive Constrained Enveloping Splines model.
\item{\code{"aces_quintic"}}: Random Forest Quintic Smoothed Adaptive Constrained Enveloping Splines model.
}}
}
\value{
This function returns the Out-of-Bag error metric, which provides an estimate of the model's prediction error on unseen data.
}
\description{
This function calculates the Out-Of-Bag (OOB) error for a given Random Forest Adaptive Constrained Enveloping Splines (RF-ACES) with a certain number of base learners.
}
