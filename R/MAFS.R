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
#' @param err.red Minimum reduced error rate for the addition of two new basis functions. Default is \code{0.01}.
#' @param minspan Minimum number of observations between knots.
#' \itemize{
#' \item{\code{minspan = -2}} Minspan computed as in Zhang (1994).
#' \item{\code{minspan = -1}} Minspan computed as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' \item{\code{minspan =  m}} \code{m} observations between knots.
#' }
#' @param endspan Minimum number of observations before the first and after the final knot.
#' \itemize{
#' \item{\code{endspan = -2}} Endspan computed as in Zhang (1994).
#' \item{\code{endspan = -1}} Endspan computed as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' \item{\code{endspan = +m}} \code{m} observations between knots.
#' }
#' @param knotsGrid Grid of knots to perform MAFS:
#' \itemize{
#' \item{\code{knotsGrid = -1}} The original approach based on the observed data is used.
#' \item{\code{knotsGrid =  0}} A grid of equidistant knots is created from the scattering of the data.
#' \item{\code{knotsGrid = +p}} \code{p} percentiles are computed to perform the grid of knots.
#' }
#' @param linpreds \code{logical}. If \code{TRUE}, predictors can enter linearly.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#'
#' @return A \code{MAFS} object.
#'
#' @export
MAFS <- function(data, x, y, nterms, Kp = 1, d = 2,
                 err.red = 0.01, minspan = -2,
                 endspan = -1, knotsGrid = -1,
                 linpreds = FALSE, na.rm = TRUE) {

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
    # status: intercept / paired / unpaired
    #   side: E (entire) / R (right) / L (left)
    #     Bp: basis function
    #     xi: variable for splitting
    #      t: knot for splitting
    #      R: mean squared error between true data and predicted data (B %*% alpha)
    #  alpha: regression coefficients
  bf <- list(
    "id"     = 1,
    "status" = "intercept",
    "side"   = "E",
    "Bp"     = rep(1, N),
    "xi"     = c(-1),
    "t"      = c(-1),
    "R"      = mse(data[, y, drop = F], matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, max)),
    "alpha"  = apply(data[, y, drop = FALSE], 2, max)
  )

  # Set of knots. It saves indexes of data used as knots.
  knots.list <- vector("list", nX)

  # Set of basis functions and B Matrix
  ForwardModel <- list(BF = list(bf), B = matrix(rep(1, N)))

  # Error of first basis function
  err <- bf[["R"]]

  # Minimum span and end span
  L.Le <- computeSpan(data, minspan, endspan, nX)

  # Grid of knots
  if (knotsGrid == 0) {
    knotsGrid <- vector("list", nX)

    for(var in 1:nX) {
      median.diff      <- median(diff(data[order(data[, var]), var]))
      knotsGrid[[var]] <- seq(min(data[, var]), max(data[, var]), median.diff)
    }

  } else if (knotsGrid > 0) {
    gr        <- knotsGrid
    knotsGrid <- vector("list", nX)

    for(var in 1:nX) {
      knotsGrid[[var]] <- quantile(data[, var], probs = seq(0, 1, 1 / (gr - 1)))
      names(knotsGrid[[var]]) <- NULL
    }

  } else {
    knotsGrid <- NULL
  }

  while(length(ForwardModel[["BF"]]) + 2 < nterms) {

    # Divide the node
    B.BF.knots.err <- AddBF(data, x, y,
                            ForwardModel,
                            knots.list, Kp,
                            L.Le[[1]], L.Le[[2]],
                            knotsGrid, linpreds,
                            err.min = err)

    ForwardModel[["B"]]  <- B.BF.knots.err[[1]]
    ForwardModel[["BF"]] <- B.BF.knots.err[[2]]
    knots.list           <- B.BF.knots.err[[3]]
    new.err              <- B.BF.knots.err[[4]]

    # New minimun error
    if (new.err < err * (1 - err.red)) {
      err <- new.err

    } else {
      break
    }
  }

  # ==
  # Forward MAFS
  # ==

  xi    <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["xi"]]))
  t     <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["t"]]))
  alpha <- ForwardModel[["BF"]][[length(ForwardModel[["BF"]])]][["alpha"]]

  knotsForward <- unique(cbind(xi, t))

  MAFS_Forward = list(
    BF    = ForwardModel[["BF"]],
    B     = ForwardModel[["B"]],
    knots = knotsForward,
    alpha = alpha
    )

  # ================== #
  # BACKWARD ALGORITHM #
  # ================== #

  SubModels_MAFS <- PruningMAFS(data, y, MAFS_Forward, d)

  # Lack-of-fit at each model
  LOFs <- sapply(SubModels_MAFS, function(x) x[["LOF"]])

  # Model with minimum error
  BackwardModel <- SubModels_MAFS[[which.min(LOFs)]]
  knotsBackward <- do.call(rbind.data.frame, BackwardModel[["t"]])

  # Interaction of variables
  if (Kp > 1 & any(knotsBackward$status == "unpaired")) {
    BackwardModel <- InteractionModel(data, x, BackwardModel, knotsBackward, err.red, Kp, d)
  }

  # ==
  # MAFS Model
  # ==

  MAFS_Model <- list(
    MAFS_SubModels = SubModels_MAFS,
    B              = BackwardModel[["B"]],
    knots          = knotsBackward,
    alpha          = BackwardModel[["alpha"]],
    LOF            = BackwardModel[["LOF"]]
    )

  # =================== #
  # SMOOTHING PROCEDURE #
  # =================== #

  # ==
  # Smooth MAFS
  # ==

  if (nrow(knotsBackward) == 0) {
    Smooth_MAFS <- MAFS_Model
    warning("Smoothing not available because there are no knots. MAFS model is used instead.")

  } else {
    QSmooth_MAFS <- QSmoothMAFS(data, nX, knotsBackward, y)
    CSmooth_MAFS <- CSmoothMAFS(data, nX, knotsBackward, y)
  }

  # MAFS object
  MAFS <- MAFS_object(data, x, y,
                      nterms, Kp, d,
                      err.red, na.rm,
                      MAFS_Forward,
                      MAFS_Model,
                      QSmooth_MAFS,
                      CSmooth_MAFS)

  return(MAFS)
}

#' @title Compute Span
#'
#' @description This function computes the minimum span and the end span.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param minspan Minimum number of observations between knots.
#' \itemize{
#' \item{\code{minspan = -2}} Minspan computed as in Zhang (1994).
#' \item{\code{minspan = -1}} Minspan computed as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' \item{\code{minspan =  m}} \code{m} observations between knots.
#' }
#' @param endspan Minimum number of observations before the first and after the final knot.
#' \itemize{
#' \item{\code{endspan = -2}} Endspan computed as in Zhang (1994).
#' \item{\code{endspan = -1}} Endspan computed as in Friedman's MARS paper section 3.8 with alpha = 0.05.
#' \item{\code{endspan =  m}} \code{m} observations between knots.
#' }
#'
#' @param nX Number of inputs.
#'
#' @return \code{list} with the minimum span and the end span.
computeSpan <- function(data, minspan, endspan, nX) {

  N <- nrow(data)

  # Minimum span (L)
  if (minspan == - 2) {
    # Zhang approach
    L <- c()

    for (var in 1:nX) {
      max3 <- data[order(data[, var], decreasing = TRUE), var][1:3]
      min3 <- data[order(data[, var]), var][1:3]

      m1 <- - (max3[1] - min3[1]) / (2.5 * (N - 1)) * log2(- (1 / N) * log(0.95))
      m2 <- (1 / N) * sum(max3 - min3)

      Lvar <- ceiling(max(m1, m2))

      L <- c(L, Lvar)
    }

  } else if (minspan == - 1) {
    # Friedman approach (this value must be computed later)
    L <- -1

  } else {
    L <- minspan
  }

  # End span (Le)
  if (endspan == - 2) {
    # Zhang approach
    Le <- c()

    for (var in 1:nX) {
      max3 <- data[order(data[, var], decreasing = TRUE), var][1:3]
      min3 <- data[order(data[, var]), var][1:3]

      m1 <- - (max3[1] - min3[1]) / (2.5 * (N - 1)) * log2(- (1 / N) * log(0.95))
      m2 <- (1 / N) * sum(max3 - min3)

      Levar <- ceiling(max(m1, m2))

      Le <- c(Le, Levar)
    }

  } else if (endspan == - 1) {
    # Friedman approach
    Le <- ceiling(3 - log2(0.05 / nX))

  } else {
    Le <- endspan
  }

  return(list(L, Le))

}

#' @title Create a MAFS object
#'
#' @description This function saves information about the Multivariate Adapative Frontier Splines model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param nterms Maximum number of terms created by the forward algorithm (before pruning).
#' @param Kp Maximum degree of interaction allowed. Default is \code{1}.
#' @param d Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}. If set to \code{-1}, \code{GCV = RSS / n}.
#' @param err.red Minimun reduced error rate for the addition of two new basis functions. Default is \code{0.01}.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' @param Forward.MAFS \code{list} containing the Forward Multivariate Adaptive Frontier Splines.
#' @param MAFS \code{list} containing the Multivariate Adaptive Frontier Splines model.
#' @param QSmooth.MAFS \code{list} containing the Smooth Multivariate Adaptive Frontier Splines model (quintic version).
#' @param CSmooth.MAFS \code{list} containing the Smooth Multivariate Adaptive Frontier Splines model (cubic version).
#'
#' @importFrom dplyr %>% filter select
#'
#' @return A \code{MAFS} object.
MAFS_object <- function(data, x, y, nterms, Kp, d, err.red, na.rm,
                        Forward.MAFS, MAFS, QSmooth.MAFS, CSmooth.MAFS) {

  MAFS_object <- list("data" = list(df = data,
                                    x  = x,
                                    y  = y,
                                    input_names  = colnames(data)[x],
                                    output_names = colnames(data)[y],
                                    row_names    = rownames(data)),
                      "control" = list(nterms = nterms,
                                        Kp = Kp,
                                        d  = d,
                                        err.red = err.red,
                                        na.rm = na.rm),
                      "Forward.MAFS" = Forward.MAFS,
                      "MAFS"         = MAFS,
                      "QSmooth.MAFS" = QSmooth.MAFS,
                      "CSmooth.MAFS" = CSmooth.MAFS)

  class(MAFS_object) <- "MAFS"

  # Data for Efficiency Models
  knots <- unique(MAFS[["knots"]][, c("xi", "t")])
  knots.list <- vector("list", length(x))

  for (v in 1:length(x)) {
    knots.list[[v]] <- knots %>% filter(xi == v) %>% select(t) %>% t()
    knots.list[[v]] <- sort(c(knots.list[[v]], min(data[, v]), max(data[, v])))
  }

  aknots <- expand.grid(knots.list)
  names(aknots) <- colnames(data)[x]

  yknots <- predict(MAFS_object, aknots, x, 2)

  MAFS_object[["Efficiency"]] <- list("aknots" = aknots,
                                      "yknots" = yknots)

  return(MAFS_object)
}
