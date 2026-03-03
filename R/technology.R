#' @title Generate a Technology Set
#'
#' @description
#' This function generates a technology set based on input data and an ACES model output.
#' The technology set represents the feasible combinations of inputs and outputs.
#'
#' @param var_names
#' A \code{string} with variables names for technology set.
#'
#' @param tech_xmat
#' A \code{matrix} representing the input data.
#'
#' @param tech_ymat1
#' A \code{matrix} representing the output data.
#'
#' @param tech_ymat2
#' A \code{matrix} representing the output data generated from an ACES model.
#'
#' @param table_scores
#' A \code{matrix} containing radial output scores using all outputs and using just each individual output.
#'
#' @return
#' A \code{matrix} representing the technology set.

generate_technology <- function (
    var_names,
    tech_xmat,
    tech_ymat1,
    tech_ymat2,
    table_scores
    ) {

  # define technology
  tech <- list()

  # matrix of inputs for the technology
  tech[["xmat"]] <- as.matrix(tech_xmat)
  colnames(tech[["xmat"]]) <- var_names[1:ncol(tech_xmat)]

  # matrix of outputs for the technology
  tech_ymat <- as.data.frame (
    matrix (
      NA,
      nrow = nrow(table_scores),
      ncol = ncol(table_scores) - 1
      )
    )

  colnames(tech_ymat) <- colnames(tech_ymat1)

  for (i in 1:nrow(table_scores)) {

    for (j in 2:ncol(table_scores)) {

      if (abs(table_scores[i, 1] - table_scores[i, j]) <= 0.05) {

        tech_ymat[i, j - 1] <- tech_ymat2[i, j - 1]

      } else {

        tech_ymat[i, j - 1] <- tech_ymat1[i, j - 1]

      }
    }
  }

  tech[["ymat"]] <- as.matrix(tech_ymat)
  colnames(tech[["ymat"]]) <- var_names[(ncol(tech_xmat) + 1):length(var_names)]

  return(tech)

}

#' @title Update Technology Set
#'
#' @description
#' This function constructs a new production technology by modifying the output matrix based on a refinement procedure that replaces overestimated outputs. Additionally, it allows incorporating the origin into the technology set.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param tech_xmat
#' A \code{matrix} representing the input data (optional). Required if \code{psi} is provided.
#'
#' @param tech_ymat
#' A \code{matrix} representing the output data (optional). Required if \code{psi} is provided.
#'
#' @param psi
#' A \code{numeric} threshold controlling the refinement of predicted outputs (optional). If the predicted value for a given output exceeds the observed value by more than \code{psi}, the observed value is used instead.
#'
#' @param pto
#' A \code{logical} indicating if (0,0) should be added to the technology set. Default is \code{FALSE}.
#'
#' @return
#' An \code{aces} object with updated technologies.
#'
#' @export

change_technology <- function (
    object,
    tech_xmat = NULL,
    tech_ymat = NULL,
    psi = NULL,
    pto = FALSE
    ) {

  # ==================================================== #
  # 0. DATA SCALING                                      #
  # ==================================================== #

  # scaling parameters
  scaling <- object[["control"]][["scale"]]

  # copy of the original variables
  new_x <- if(!is.null(tech_xmat)) as.matrix(tech_xmat) else NULL
  new_y <- if(!is.null(tech_ymat)) as.matrix(tech_ymat) else NULL

  # if the model was scaled, we need to scale new data too
  if (!is.null(scaling) && scaling$is_scaled) {

    if (!is.null(new_x)) {
      new_x <- sweep(new_x, 2, scaling$mean_x, "/")
    }

    if (!is.null(new_y)) {
      new_y <- sweep(new_y, 2, scaling$mean_y, "/")
    }

  }

  # ==================================================== #
  # 1. REFINEMENT LOGIC (PSI)                            #
  # ==================================================== #

  # Execute only if data and threshold are provided
  if (!is.null(psi) && !is.null(new_x) && !is.null(new_y)) {

    # sample size
    N <- nrow(new_x)

    # number of outputs
    nY <- ncol(new_y)

    # initialize table of scores
    table_scores <- matrix (
      ncol = nY + 1,
      nrow = N,
      dimnames = list(NULL, c("y_all", paste("y", 1:nY, sep = "")))
    ) %>% as.data.frame()

    # calculate global scores
    table_scores[, 1] <- rad_out (
      tech_xmat = as.matrix(new_x),
      tech_ymat = as.matrix(new_y),
      eval_xmat = as.matrix(new_x),
      eval_ymat = as.matrix(new_y),
      convexity = TRUE,
      returns = "variable",
      type = "objective"
    )[, 1]

    # calculate individual scores per output
    for (out in 1:nY) {
      table_scores[, 1 + out] <- rad_out (
        tech_xmat = as.matrix(new_x),
        tech_ymat = as.matrix(new_y[, out]),
        eval_xmat = as.matrix(new_x),
        eval_ymat = as.matrix(new_y[, out]),
        convexity = TRUE,
        returns = "variable",
        type = "objective"
      )[, 1]

    }

    # update technology for each method
    for (m in names(object[["technology"]])) {

      # predictions from the model
      Bmat <- object[["methods"]][[m]][["Bmatx"]]
      coef <- object[["methods"]][[m]][["coefs"]]

      tech_ymat2 <- Bmat %*% coef
      colnames(tech_ymat2) <- colnames(tech_ymat)

      # initialize refined output matrix with observed values
      tech_out <- as.matrix(new_y)

      # refinement procedure
      for (i in 1:nrow(table_scores)) {
        for (j in 1:nY) {
          # compare global score (col 1) vs individual score (col j+1)
          if (abs(table_scores[i, 1] - table_scores[i, j + 1]) <= psi) {
            # keep predicted value
            tech_out[i, j] <- tech_ymat2[i, j]
          }
          # else: keep observed value (already in tech_out)
        }
      }

      # update the object
      object[["technology"]][[m]][["xmat"]] <- as.matrix(new_x)
      object[["technology"]][[m]][["ymat"]] <- tech_out

    }
  }

  object[["control"]][["psi"]] <- psi

  # ==================================================== #
  # 2. PASS THROUGH THE ORIGIN (PTO) LOGIC               #
  # ==================================================== #

  # Execute sequentially after refinement (if any)
  if (pto) {
    for (m in names(object[["technology"]])) {

      current_x <- object[["technology"]][[m]][["xmat"]]
      current_y <- object[["technology"]][[m]][["ymat"]]

      # create zero rows
      zero_x <- rep(0, ncol(current_x))
      zero_y <- rep(0, ncol(current_y))

      # append (0,0) to the technology set
      object[["technology"]][[m]][["xmat"]] <- rbind(current_x, unname(zero_x))
      object[["technology"]][[m]][["ymat"]] <- rbind(current_y, unname(zero_y))
    }
  }

  return(object)

}

#' @title Get Technology Set
#'
#' @description
#' Retrieves the matrix constituting the estimated technology for a specific method from an ACES object.
#'
#' @param object
#' An \code{aces} object.
#' @param method
#' Model for prediction:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_cubic"}}: Cubic Smoothed Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_quintic"}}: Quintic Smoothed Adaptive Constrained Enveloping Splines.
#' }
#'
#' @return A matrix containing the inputs (xmat) and the estimated output (ymat).
#'
#' @export
#'
get_technology <- function (
    object,
    method = "aces"
    ) {

  # 1. check if the requested method exists in the object structure
  if (is.null(object[["technology"]]) || is.null(object[["technology"]][[method]])) {
    stop (
      sprintf (
        "Method '%s' not found in the 'technology' slot of the provided object.",
        method
        )
      )
  }

  # 2. access the specific technology list
  tech_data <- object[["technology"]][[method]]

  # 3. extract xmat and ymat
  x_mat <- tech_data[["xmat"]]
  y_mat <- tech_data[["ymat"]]

  # 4. restore original scale if needed
  scaling <- object[["control"]][["scale"]]

  if (!is.null(scaling) && scaling$is_scaled) {
    # multiply by scaling factors (Means) to recover original magnitude
    x_mat <- sweep(x_mat, 2, scaling$mean_x, "*")
    y_mat <- sweep(y_mat, 2, scaling$mean_y, "*")
  }

  # 5. combine results into a single matrix
  tech <- cbind(x_mat, y_mat)

  return(tech)

}
