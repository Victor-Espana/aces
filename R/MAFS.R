#' @title Multivariate Adaptive Frontier Splines
#'
#' @description This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity, which is based upon the adaptation of the machine learning technique known as Multivariate Adaptive Regression Splines (MARS).
#'
#' @name MAFS
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param nterms Maximum number of terms created by the forward algorithm (before pruning).
#' @param Kp Maximum degree of interaction allowed. Default is \code{1}.
#' @param d Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}. If it is set to \code{-1}, \code{GCV = RSS / n}.
#' @param err_red Minimum reduced error rate for the addition of two new basis functions. Default is \code{0.01}.
#' @param minspan Minimum number of observations between knots. When \code{minspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param endspan Minimum number of observations before the first and after the final knot. When \code{endspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param linpreds \code{logical}. If \code{TRUE}, predictors can enter linearly.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#'
#' @return A \code{MAFS} object.
#'
#' @export
MAFS <- function(data, x, y, nterms, Kp = 1, d = 2, err_red = 0.01,
                  minspan = 0, endspan = 0, linpreds = FALSE,
                  na.rm = TRUE) {

  # Data in data[x, y] format.
  data <- preProcess(data = data, x = x, y = y, na.rm = na.rm)

  # Samples in data
  N <- nrow(data)

  # Number of inputs
  nX <- length(x)

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  # nterms
  # Max value between hyperparameter selected by the user and
  # C(M) < N during backward
  nterms <- min(nterms, floor((2 * N + d) / (2 + d)))

  # ================= #
  # FORWARD ALGORITHM #
  # ================= #

  # basis function
    #     id: index
    # status: intercept / paired / not paired
    #   side: E (entire) / R (right) / L (left)
    #     Bp: basis function
    #     xi: variable for splitting
    #      t: knot for splitting
    #      R: mean squared error between true data and predicted data (B %*% alpha)
    #  alpha: regression coefficients
  bf <- list(
    "id" = 1,
    "status" = "intercept",
    "side" = "E",
    "Bp" = rep(1, N),
    "xi" = c(-1),
    "t" = c(-1),
    "R" = mse(data[, y], matrix(rep(1, N)) %*% max(data[, y])),
    "alpha" = max(data[, y])
  )

  # Set of knots. It saves indexes of data used as knots.
  knots_list <- vector("list", nX)

  # Set of basis functions and B Matrix
  ForwardModel <- list(BF = list(bf),
                       B = matrix(rep(1, N)))

  # Endspan
  if (endspan == 0){
    Le <- ceiling(3 - log2(0.05 / nX))
  } else {
    Le <- endspan
  }

  # Error of first basis function
  err <- bf[["R"]]

  while(length(ForwardModel[["BF"]]) + 2 < nterms) {

    # Divide the node
    B_BF_knots_err <- AddBF(data, x, y, ForwardModel, knots_list,
                            Kp, minspan, Le, linpreds, err_min = err)

     ForwardModel[["B"]] <- B_BF_knots_err[[1]]
    ForwardModel[["BF"]] <- B_BF_knots_err[[2]]
              knots_list <- B_BF_knots_err[[3]]
                 New_err <- B_BF_knots_err[[4]]

    # New minimun error
    if (New_err < err * (1 - err_red)) {
      err <- New_err

    } else {
      break
    }
  }

  # ==
  # Forward MAFS
  # ==

     xi <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["xi"]]))
      t <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["t"]]))
  alpha <- ForwardModel[["BF"]][[length(ForwardModel[["BF"]])]][["alpha"]]

  knotsForward <- unique(cbind(xi, t))

  MAFS_Forward = list(
       BF = ForwardModel[["BF"]],
        B = ForwardModel[["B"]],
    knots = knotsForward,
    alpha = alpha)

  # ================== #
  # BACKWARD ALGORITHM #
  # ================== #

  SubModels_MAFS <- PruningMAFS(data, y, MAFS_Forward, d)

  # Lack-of-fit at each model
  LOFs <- sapply(SubModels_MAFS, function(x) x[["LOF"]])

  # Model with minimum error
  BackwardModel <- SubModels_MAFS[[which.min(LOFs)]]

  # ==
  # MAFS Model
  # ==

  knotsBackward <- do.call(rbind.data.frame, BackwardModel[["t"]])

  MAFS_Model <- list(
    MAFS_SubModels = SubModels_MAFS,
                 B = BackwardModel[["B"]],
             knots = knotsBackward,
             alpha = BackwardModel[["alpha"]],
               LOF = BackwardModel[["LOF"]])

  # ======================== #
  # TRUNCATED CUBIC FUNCTION #
  # ======================== #

  # ==
  # Smooth MAFS (backward)
  # ==

  if (nrow(knotsBackward) == 0) {
    Smooth_Backward_MAFS <- MAFS_Model
    warning("Smoothing not available because there are no knots. MAFS model is used instead.")

  } else {
    Smooth_Backward_MAFS <- SmoothMAFS(data, nX, knotsBackward, y)
  }

  # ==
  # Smooth MAFS (forward)
  # ==

  Smooth_Forward_MAFS <- SmoothForwardMAFS(data, nX, knotsForward, data[, y])

  # MAFS object
  MAFS <- MAFS_object(data, x, y, rownames(data), nterms, Kp,
                      d, err_red, minspan, endspan, na.rm,
                      MAFS_Forward, MAFS_Model,
                      Smooth_Backward_MAFS,
                      Smooth_Forward_MAFS)

  return(MAFS)
}

#' @title Create a MAFS object
#'
#' @description This function saves information about the Multivariate Adapative Frontier Splines model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param rownames \code{string}. Data rownames.
#' @param nterms Maximum number of terms created by the forward algorithm (before pruning).
#' @param Kp Maximum degree of interaction allowed. Default is \code{1}.
#' @param d Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}. If set to \code{-1}, \code{GCV = RSS / n}.
#' @param err_red Minimun reduced error rate for the addition of two new basis functions. Default is \code{0.01}.
#' @param minspan Minimum number of observations between knots. When \code{minspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param endspan Minimum number of observations before the first and after the final knot. When \code{endspan = 0} (default), it is calculated as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' @param Depth_MAFS_Model \code{list} containing the Multivariate Adaptive Frontier Splines model before pruning.
#' @param MAFS_Model \code{list} containing the Multivariate Adaptive Frontier Splines model.
#' @param Smooth_Backward_MAFS \code{list} containing the Smooth Multivariate Adaptive Frontier Splines model.
#' @param Smooth_Forward_MAFS \code{list} containing the Smooth Multivariate Adaptive Frontier Splines model (from forward).
#'
#' @return A \code{MAFS} object.
MAFS_object <- function(data, x, y, rownames, nterms, Kp, d, err_red, minspan, endspan, na.rm,
                        Depth_MAFS_Model, MAFS_Model, Smooth_Backward_MAFS, Smooth_Forward_MAFS) {

  MAFS_object <- list("data" = list(df = data,
                                    x = x,
                                    y = y,
                                    input_names = names(data)[x],
                                    output_names = names(data)[y],
                                    row_names = rownames),
                       "control" = list(nterms = nterms,
                                        Kp = Kp,
                                        d = d,
                                        err_red = err_red,
                                        minspan = minspan,
                                        endspan = endspan,
                                        na.rm = na.rm),
                       "Forward.MAFS" = Depth_MAFS_Model,
                       "MAFS" = MAFS_Model,
                       "Smooth.MAFS" = Smooth_Backward_MAFS,
                       "Smooth.MAFS.Forward" = Smooth_Forward_MAFS)

  class(MAFS_object) <- "MAFS"

  return(MAFS_object)

}
