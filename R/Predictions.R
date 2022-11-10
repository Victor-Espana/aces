#' @title Model Prediction for Additive Adaptive Frontier Splines.
#'
#' @description This function predicts the expected output by a \code{AAFS} object.
#'
#' @param object A \code{AAFS} object.
#' @param newdata \code{data.frame}. Set of input variables to predict on.
#' @param x Inputs index.
#' @param class Model for prediction. \code{1} for AAFS before pruning, \code{2} for AAFS, \code{3} for CSAAFS (cubic version) and \code{4} for QSAAFS (quintic version) and .
#'
#' @return \code{data.frame} with the predicted values.
#'
#' @export
predict.AAFS <- function(object, newdata, x, class = 2) {

  train_names <- object[["data"]][["input_names"]]
   test_names <- colnames(newdata)[x]

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
    model <- object[["Forward.AAFS"]]
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
    model <- object[["AAFS"]]
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

  } else if (class == 4){

    model <- object[["QSAAFS"]]
    knots <- model[["Qknots"]]

    if (is.null(knots)) {
      # Only possible in Smooth AAFS (from backward)
      B <- B

    } else {
      for (st in c("paired", "unpaired")) {
        # Variable
        for (i in 1:length(knots)) {
          # This variables does not expand the data
          if (is.null(knots[[i]])) next
          for (j in 1:length(knots[[i]])) {

            if (knots[[i]][[j]][["status"]] != st) next

            # Quintic Knot
            Qt <- knots[[i]][[j]][["t"]]

            t0 <- Qt[1] # t-
            t1 <- Qt[2] # t
            t2 <- Qt[3] # t+

            d  <- t2 - t0 # t+ - t-
            d1 <- t2 - t1 # t+ - t
            d2 <- t1 - t0 # t  - t-

            alpha1 <- (6 * d1 - 4 * d2) / d ^ 3
            alpha2 <- (4 * d1 - 6 * d2) / d ^ 3

            beta1 <- (7 * d2 - 8 * d1) / d ^ 4
            beta2 <- (7 * d1 - 8 * d2) / d ^ 4

            gamma1 <- gamma2 <- (3 * d1 - 3 * d2) / d ^ 5

            # Both or right
            if (st == "paired") {

              term1 <- alpha1 * (newdata[, i] - t0) ^ 3
              term2 <- beta1  * (newdata[, i] - t0) ^ 4
              term3 <- gamma1 * (newdata[, i] - t0) ^ 5

              # Q1
              Q1 <- ifelse(data[, i] <= t0,
                           0,
                           (ifelse((newdata[, i] > t0) & (newdata[, i] < t2),
                                   term1 + term2 + term3,
                                   newdata[, i] - t1)))

              B <- cbind(B, Q1)
            }

            # Both or left
            if (st == "paired" | st == "unpaired") {

              term1 <- alpha2 * (newdata[, i] - t2) ^ 3
              term2 <- beta2  * (newdata[, i] - t2) ^ 4
              term3 <- gamma2 * (newdata[, i] - t2) ^ 5

              # Q2
              Q2 <- ifelse(newdata[, i] <= t0,
                           t1 - newdata[, i],
                           (ifelse((newdata[, i] > t0) & (newdata[, i] < t2),
                                   term1 + term2 + term3,
                                   0)))
              B <- cbind(B, Q2)
            }
          }
        }
      }
    }
  } else {

    model <- object[["CSAAFS"]]
    knots <- model[["Cknots"]]

    if (is.null(knots)) {
      # Only possible in Smooth AAFS (from backward)
      B <- B

    } else {
      for (st in c("paired", "unpaired")) {
        # Variable
        for (i in 1:length(knots)) {
          # This variables does not expand the data
          if (is.null(knots[[i]])) next
          for (j in 1:length(knots[[i]])) {

            if (knots[[i]][[j]][["status"]] != st) next

            # Cubic Knot
            Ct <- knots[[i]][[j]][["t"]]

            t0 <- Ct[1] # t-
            t1 <- Ct[2] # t
            t2 <- Ct[3] # t+

            d <- t1 - t0 # t  - t-
            e <- t2 - t1 # t+ - t

            p1 <- (2 * e - d) / (e + d) ^ 2
            r1 <- (d - e) / (e + d) ^ 3

            p2 <- (2 * d - e) / (- e - d) ^ 2
            r2 <- (e - d) / (- e - d) ^ 3

            # Both or right
            if (st == "paired") {

              term1 <- p1 * (newdata[, i] - t0) ^ 2
              term2 <- r1 * (newdata[, i] - t0) ^ 3

              # C1
              C1 <- ifelse(newdata[, i] <= t0,
                           0,
                           (ifelse((newdata[, i] > t0) & (newdata[, i] < t2),
                                   term1 + term2,
                                   newdata[, i] - t1)))


              B <- cbind(B, C1)
            }

            # Both or left
            if (st == "paired" | st == "unpaired") {

              term1 <- p2 * (newdata[, i] - t2) ^ 2
              term2 <- r2 * (newdata[, i] - t2) ^ 3

              # C2
              C2 <- ifelse(newdata[, i] <= t0,
                           t1 - newdata[, i],
                           (ifelse((newdata[, i] > t0) & (newdata[, i] < t2),
                                   term1 + term2,
                                   0)))

              B <- cbind(B, C2)
            }
          }
        }
      }
    }
  }

  # Predictions
  y_hat <- matrix(NA, nrow = N, ncol = length(y))

  for (out in 1:length(y)) {
    y_hat[, out] <- pmax(0, B %*% model[["coefs"]][, out])
  }

  predictions <- as.data.frame(y_hat)
  names(predictions) <- paste(object[["data"]][["output_names"]], "_pred", sep = "")

  return(predictions)
}
