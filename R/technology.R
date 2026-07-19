#' @title Build an Estimated Technology
#'
#' @description
#' Builds the reference sample used to define an ACES production technology.
#' For each DMU and output, the estimated maximum output is retained when the
#' joint and output-specific radial scores agree within the refinement threshold.
#' Otherwise, the observed output is kept to avoid using a marginal prediction
#' that may distort the jointly feasible output vector.
#'
#' With one output, the joint and output-specific scores are identical, so the
#' reference sample uses the estimated maximum attainable outputs. With multiple outputs,
#' it can contain both predicted and observed components. The original inputs and
#' these refined output vectors define the DEA-type technology used for efficiency
#' measurement.
#'
#' @param var_names
#' Character vector of input and output names.
#'
#' @param tech_xmat
#' Matrix of inputs that define the technology.
#'
#' @param tech_ymat1
#' Matrix of observed outputs, treated as feasible production levels.
#'
#' @param tech_ymat2
#' Matrix of maximum attainable outputs estimated by ACES.
#'
#' @param table_scores
#' Matrix of output-oriented radial scores computed jointly and separately for
#' each output. A predicted component is retained when the absolute difference
#' between its joint and output-specific scores is at most \code{0.05}; otherwise,
#' the observed component is used.
#'
#' @return
#' A list with \code{xmat}, the original input vectors, and \code{ymat}, the
#' refined output vectors. Together they form the reference sample for the
#' estimated production technology.

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

#' @title Update an Estimated Technology
#'
#' @description
#' Rebuilds the reference samples stored in an ACES object from supplied data and
#' a new output-refinement threshold. This changes the points used to construct
#' the production technologies without refitting the spline models. The origin
#' can also be added as an explicit reference point.
#'
#' @param object
#' An \code{aces} object.
#'
#' @param tech_xmat
#' Optional matrix of inputs used to rebuild the technologies.
#'
#' @param tech_ymat
#' Optional matrix of observed outputs used to rebuild the technologies.
#'
#' @param psi
#' Optional output-refinement threshold. When a joint and a single-output score
#' differ by more than this value, the observed output is retained; otherwise,
#' the estimated maximum output is used. Larger values therefore retain more
#' predicted components.
#'
#' @param pto
#' If \code{TRUE}, add the origin to each technology.
#'
#' @return
#' The supplied \code{aces} object with rebuilt technology components. Fitted
#' models and other controls are unchanged.
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

#' @title Extract an Estimated Technology
#'
#' @description
#' Extracts the input and output points that define one estimated technology.
#' These are the reference points used by \code{get_scores()} and
#' \code{get_targets()}, rather than every feasible point in the production set
#' or fitted values for a new data set.
#'
#' @param object
#' An \code{aces} object.
#' @param method
#' Fitted model whose technology is returned:
#' \itemize{
#' \item{\code{"aces_forward"}}: Forward Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces"}}: Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_cubic"}}: Cubic Smoothed Adaptive Constrained Enveloping Splines.
#' \item{\code{"aces_quintic"}}: Quintic Smoothed Adaptive Constrained Enveloping Splines.
#' }
#'
#' @return A matrix containing the technology inputs followed by its outputs, on
#' the original data scale.
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
