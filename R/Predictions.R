#' @title Model Prediction for Multivariate Adaptive Frontier Splines.
#'
#' @description This function predicts the expected output by a \code{MAFS} object.
#'
#' @param object A \code{MAFS} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param class Model for prediction. \code{1} for MAFS before pruning, \code{2} for MAFS, \code{3} for Smooth MAFS and \code{4} for Smooth Forward MAFS.
#'
#' @return \code{data.frame} with the predicted values.
#'
#' @export
predict.MAFS <- function(object, newdata, x, class = 2) {

  train_names <- object[["data"]][["input_names"]]
   test_names <- names(newdata)[x]

  if (!identical(sort(train_names), sort(test_names))) {
    stop("Different variable names in training and test sets.")
  }

  # Select variables and reorder as in training data
  newdata <- newdata[, x, drop = FALSE][train_names]
        N <- nrow(newdata)
        y <- object[["data"]][["y"]]

  # B1
  B <- matrix(rep(1, N), nrow = N)

  if (class == 1) {
    model <- object[["Forward.MAFS"]]
    # Always in pairs
    knots <- model[["knots"]]

    for (i in 1:nrow(knots)) {

      xi <- knots[i, 1]
       t <- knots[i, 2]

      # Basis functions for the new data
      hinge1 <- ifelse(newdata[, xi] > t, newdata[, xi] - t, 0)
      hinge2 <- ifelse(newdata[, xi] < t, t - newdata[, xi], 0)

      # Update B
      B <- cbind(B, hinge1, hinge2)
    }

  } else if (class == 2) {
    model <- object[["MAFS"]]
    knots <- model[["knots"]]

    if (nrow(knots) == 0) {
      B <- B

    } else {
      for (i in 1:nrow(knots)) {
        xi <- knots[i, 1]
         t <- knots[i, 2]

        # Basis functions for the new data

        if (knots[i, 3] == "R") {
          hinge1 <- ifelse(newdata[, xi] > t, newdata[, xi] - t, 0)

          # Update B
          B <- cbind(B, hinge1)

        } else {
          hinge2 <- ifelse(newdata[, xi] < t, t - newdata[, xi], 0)

          # Update B
          B <- cbind(B, hinge2)
        }
      }
    }

  } else if (class %in% c(3, 4)) {

    if (class == 3) {
      model <- object[["Smooth.MAFS"]]

    } else {
      model <- object[["Smooth.MAFS.Forward"]]

    }

    knots <- model[["knots"]]

    if (nrow(knots) == 0) {
      # Only possible in Smooth MAFS (from backward)
      B <- B

    } else {

      for (st in c("paired", "not paired")) {
        for (i in 1:length(model[["cubic_knots"]])) {

          # This variables does not split the data
          if (is.null(model[["cubic_knots"]][[i]])) next

          for (j in 1:length(model[["cubic_knots"]][[i]])) {

            if (model[["cubic_knots"]][[i]][[j]][["status"]] != st) next

            # Cubic Knot
            Ct <- model[["cubic_knots"]][[i]][[j]][["t"]]

            # Knot
            t <- Ct[2]

            # Sides of that knot
            if (class == 3) {
              side <- knots[knots[, "t"] == t, "side"]

            } else {
              side <- c("L", "R")

            }

            # Basis functions for the new data
            B <- CreateCubicBF(newdata, i, Ct, B, side)
          }
        }
      }
    }
  }

  predictions <- as.data.frame(B %*% model[["alpha"]])
  names(predictions) <- paste(object[["data"]][["output_names"]], "_pred", sep = "")

  return(predictions)
}
