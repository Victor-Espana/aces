#' @title Generate the Set of Unidimensional Intervals
#'
#' @description
#' This function determines the boundaries for the intervals and identifies what basis functions are active in each interval.
#'
#' @param Bp_list
#' A \code{list} containing the set of basis functions for each input.
#'
#' @return
#' A \code{list} with the boundaries of the intervals and the IDs of the active basis functions within each interval.

set_intervals <- function (
    Bp_list
    ) {

  # number of inputs
  nX <- length(Bp_list)

  # set of intervals by variable.
  it_list <- vector("list", nX)

  for (v in 1:nX) {

    # get vector of knots in ascending order
    t_vector <- c(0, Inf)

    for (side in c("paired", "right", "left")) {
      if (!is.null(Bp_list[[v]][[side]])) {
        t_vector <- c(t_vector, sapply(Bp_list[[v]][[side]], "[[", "t"))
      }
    }

    t_vector <- sort(t_vector)

    for (j in 1:(length(t_vector) - 1)) {

      # lower bound in the interval
      lb <- t_vector[[j]]

      # upper bound in the interval
      ub <- t_vector[[j + 1]]

      # add lower and upper bounds
      it_list[[v]][[j]] <- list("Lb" = lb, "Ub" = ub)

      # basis function index
      Bp_xi <- c()

      # status of the basis function
      status <- c()

      for (side in c("paired", "right", "left")) {

        if (!is.null(Bp_list[[v]][[side]])) {

          for (l in 1:length(Bp_list[[v]][[side]])) {

            # knot
            t_value <- Bp_list[[v]][[side]][[l]][["t"]]

            # index
            Bp_xi_value <- Bp_list[[v]][[side]][[l]][["Bp_xi"]]

            if (side == "paired") {

              if (t_value <= lb) {
                Bp_xi <- c(Bp_xi, Bp_xi_value[[1]])
                status  <- c(status, "right")

              } else if (t_value >= ub) {
                Bp_xi <- c(Bp_xi, Bp_xi_value[[2]])
                status  <- c(status, "left")

              }

            } else if (side == "right") {

              if (t_value <= lb) {
                Bp_xi <- c(Bp_xi, Bp_xi_value)
                status  <- c(status, "right")

              }

            } else {

              if (t_value >= ub) {
                Bp_xi <- c(Bp_xi, Bp_xi_value)
                status  <- c(status, "left")

              }

            }
          }
        }
      }

      if (is.null(Bp_xi)) {
        it_list[[v]][[j]] <- NULL
      } else {
        it_list[[v]][[j]][["Bp"]] <- data.frame(Bp_xi, status)
      }

    }
  }

  return(it_list)

}
