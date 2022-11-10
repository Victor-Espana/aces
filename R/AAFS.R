#' @title Additive Adaptive Frontier Splines: 1 output
#'
#' @description This function estimates a production frontier satisfying some classical production theory axioms, such as monotonicity and concavity, which is based upon the adaptation of the machine learning technique known as Multivariate Adaptive Regression Splines (MARS) developed by \insertCite{friedman1991;textual}{aafs}.
#'
#' @name AAFS
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param nterms For the Forward algorithm. Maximum number of terms created before pruning.
#' @param Kp For the Forward algorithm. Maximum degree of interaction allowed. Only \code{1} is available.
#' @param err.red For the Forward algorithm. Minimum reduced error rate for the addition of a new pair of basis functions. Default is \code{0.01}.
#' @param minspan For the Forward algorithm. Minimum number of observations between two adjacent knots.
#' \itemize{
#' \item{\code{minspan = -2}} Computed as in \insertCite{zhang1994;textual}{aafs}.
#' \item{\code{minspan = -1}} Computed as in \insertCite{friedman1991;textual}{aafs}.
#' \item{\code{minspan = +m}}
#' }
#' @param endspan For the Forward algorithm. Minimum number of observations before the first and after the final knot.
#' \itemize{
#' \item{\code{endspan = -2}} Computed as in \insertCite{zhang1994;textual}{aafs}.
#' \item{\code{endspan = -1}} Computed as in \insertCite{friedman1991;textual}{aafs}.
#' \item{\code{endspan = +m}}
#' }
#' @param knotsGrid For the Forward algorithm. If \code{knotsGrid == -1} (default), the original approach of \insertCite{friedman1991;textual}{aafs} based on the observed data is used. Alternatively, a \code{list} with the grid of knots to perform AAFS can be entered. In this case, each element of the \code{list} contains a vector with the knots of one variable (e.g., the second element of the list contains the knot values of the second variable).
#' @param d For the Backward algorithm. Generalized Cross Validation (GCV) penalty per knot. Default is \code{2}.
#' @param w1 For the cubic smoothing procedure \insertCite{friedman1991}{aafs}. The distance between the central knot and a side knot is at most \code{w} times the distance between that central knot and the other side knot for each truncated cubic basis function. It must be set between 1 and 2. If a \code{vector} is entered, the \code{w} value that most reduce the lack-of-fit is selected.
#' @param w2 For the quintic smoothing procedure \insertCite{chen1999}{aafs}. The distance between the central knot and a side knot is at most \code{w} times the distance between that central knot and the other side knot for each truncated quintic basis function. It must be set between 8/7 and 1.5. If a \code{vector} is entered, the \code{w} value that most reduce the lack-of-fit is selected.
#' @param weights Weights for outputs.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#'
#' @references
#' \insertRef{zhang1994}{aafs} \cr
#' \cr
#' \insertRef{friedman1991}{aafs} \cr
#' \cr
#' \insertRef{chen1999}{aafs}
#'
#' @importFrom Rdpack reprompt
#'
#' @return An \code{AAFS} object.
#'
#' @export
AAFS <- function(data, x, y, nterms, Kp = 1, err.red = 0.01,
                 minspan = -2, endspan = -2, knotsGrid = -1,
                 d = 2, w1 = seq(1, 2, length.out = 5),
                 w2 = seq(8 / 7, 1.5, length.out = 5),
                 na.rm = TRUE) {

  # Hyperparameters errors
  if (Kp != 1) stop("Kp can only be 1.")
  if (!all(w1 >= 1 & w1 <= 2)) stop("w1 must be between 1 and 2")
  if (!all(w2 >= 8 / 7 & w2 <= 1.5)) stop("w2 must be between 8/7 and 1.5")

  # Data in data[x, y] format.
  data <- preProcess(data = data, x = x, y = y, na.rm = na.rm)

  # Samples in data
  N <- nrow(data)

  # Number of inputs
  nX <- length(x)

  # Reorder index 'x' and 'y' in data
  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  # nterms: max value between nterms and C(M) < N during backward
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
    #      R: mean squared error between true data and predicted data (B %*% coefs)
    #  coefs: regression coefficients
  bf <- list(
    "id"     = 1,
    "status" = "intercept",
    "side"   = "E",
    "Bp"     = rep(1, N),
    "xi"     = c(-1),
    "t"      = c(-1),
    "R"      = me(data[, y, drop = F], matrix(rep(1, N)) %*% apply(data[, y, drop = F], 2, max)),
    "coefs"  = apply(data[, y, drop = FALSE], 2, max)
  )

  # Set of knots. It saves indexes of data used as knots.
  knots.list <- vector("list", nX)

  # Set of basis functions and B Matrix
  ForwardModel <- list(BF = list(bf), B = matrix(rep(1, N)))

  # Error of the first basis function
  err <- bf[["R"]]
  max.err.y <- which.max(err) # maximum error id
  try.y <- rep(TRUE, length(y))

  # Minimum span (minspan) and end span (endspan)
  L.Le <- computeSpan(data, minspan, endspan, nX)

  if (!is.list(knotsGrid)) {
    knotsGrid <- lapply(1:nX, function(i) data[, i])
  }

  while(length(ForwardModel[["BF"]]) + 2 < nterms) {

    # Add new BFs
    B.BF.knots.err <- AddBF(data, x, y,
                            ForwardModel,
                            knots.list, Kp,
                            L.Le[[1]], L.Le[[2]],
                            knotsGrid,
                            err.min = err,
                            err.idx = max.err.y)

    new.err <- B.BF.knots.err[[4]]

    # New minimum error
    if (mean(unlist(new.err)) < mean(unlist(err)) * (1 - err.red)) {

      # If the reduction of the error is enough...
      # Update B
      ForwardModel[["B"]]  <- B.BF.knots.err[[1]]
      # Update basis functions
      ForwardModel[["BF"]] <- B.BF.knots.err[[2]]
      # Update the knot list
      knots.list <- B.BF.knots.err[[3]]
      # Update error
      err <- new.err
      max.err.y <- which.max(err) # maximum error id
      try.y <- rep(TRUE, length(y))

    } else if (any(try.y)){
      try.y[max.err.y] <- FALSE

      if (!any(try.y)) {
        break
      } else {
        max.err.y <- c(1:length(y))[try.y][1]
      }
    }
  }

  # ==
  # Forward AAFS
  # ==

  xi <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["xi"]]))
   t <- unlist(sapply(ForwardModel[["BF"]], function(x) if(all(x[["xi"]] != -1)) x[["t"]]))
  coefs <- ForwardModel[["BF"]][[length(ForwardModel[["BF"]])]][["coefs"]]

  knotsForward <- unique(cbind(xi, t))

  AAFS.Forward = list(
       BF = ForwardModel[["BF"]],
        B = ForwardModel[["B"]],
    knots = knotsForward,
    coefs = coefs
    )

  # ================== #
  # BACKWARD ALGORITHM #
  # ================== #
  AAFS.Submodels <- PruningAAFS(data, y, AAFS.Forward, d)

  # Lack-of-fit at each model
  GCVs <- sapply(AAFS.Submodels, function(x) x[["GCV"]])

  # Model with minimum error
  AAFS.Backward <- AAFS.Submodels[[which.min(GCVs)]]
  knotsBackward <- do.call(rbind.data.frame, AAFS.Backward[["t"]])

  # ==
  # AAFS Model
  # ==

  AAFS <- list(
    AAFS_SubModels = AAFS.Submodels,
                 B = AAFS.Backward[["B"]],
             knots = knotsBackward,
             coefs = AAFS.Backward[["coefs"]],
               GCV = AAFS.Backward[["GCV"]]
    )

  # =================== #
  # SMOOTHING PROCEDURE #
  # =================== #

  # ==
  # Smooth AAFS
  # ==

  if (nrow(knotsBackward) == 0) {
    CSmooth.AAFS <- QSmooth.AAFS <- AAFS.Model
    warning("Smoothing not available because there are no knots. AAFS model is used instead.")

  } else {
    # List for x.space
    x.space <- vector("list", nX)

    for (xi in 1:nX) {
      # Knots for the xi variable
      kt_xi <- sort(unique(knotsBackward[knotsBackward[, 1] == xi, 2]))

      if (length(kt_xi) == 0) next

      # Add the initial and the end observation. They cannot be used as knots.
      x.space[[xi]] <- c(min(data[, xi]), kt_xi, max(data[, xi]))

      # Calculate Midpoints
      x.space[[xi]] <- sort(c(x.space[[xi]], x.space[[xi]][- length(x.space[[xi]])] + diff(x.space[[xi]]) / 2))
    }

    CSmooth.AAFS <- CSAAFS(data, nX, knotsBackward, x.space, y, w1)
    QSmooth.AAFS <- QSAAFS(data, nX, knotsBackward, x.space, y, w2)
  }

  # AAFS object
  AAFS <- AAFS_object(data, x, y, nterms, Kp, err.red,
                      minspan, endspan, knotsGrid, d,
                      na.rm, AAFS.Forward, AAFS,
                      QSmooth.AAFS, CSmooth.AAFS)

  return(AAFS)
}

#' @title Compute Span
#'
#' @description This function computes the minimum span and the end span.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param minspan Minimum number of observations between two adjacent knots.
#' \itemize{
#' \item{\code{minspan = -2}} Computed as in \insertCite{zhang1994;textual}{aafs}.
#' \item{\code{minspan = -1}} Computed as in \insertCite{friedman1991;textual}{aafs}.
#' \item{\code{minspan = +m}}
#' }
#' @param endspan Minimum number of observations before the first and after the final knot.
#' \itemize{
#' \item{\code{endspan = -2}} Computed as in \insertCite{zhang1994;textual}{aafs}.
#' \item{\code{endspan = -1}} Computed as in \insertCite{friedman1991;textual}{aafs}.
#' \item{\code{endspan = +m}}
#' }
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

#' @title Create an AAFS object
#'
#' @description This function saves information about the Additive Adaptive Frontier Splines model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param nterms Maximum number of terms created by the forward algorithm (before pruning).
#' @param Kp Maximum degree of interaction allowed. Only \code{1} is available.
#' @param err.red Minimum reduced error rate for the addition of two new basis functions.
#' @param minspan Minimum number of observations between two adjacent knots.
#' @param endspan Minimum number of observations before the first and after the final knot.
#' @param knotsGrid Grid of knots to perform AAFS.
#' @param d Generalized Cross Validation (GCV) penalty per knot.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' @param Forward.AAFS \code{list} containing the Forward Additive Adaptive Frontier Splines.
#' @param AAFS \code{list} containing the Additive Adaptive Frontier Splines model.
#' @param QSmooth.AAFS \code{list} containing the Smooth Additive Adaptive Frontier Splines model (quintic version).
#' @param CSmooth.AAFS \code{list} containing the Smooth Additive Adaptive Frontier Splines model (cubic version).
#'
#' @importFrom dplyr %>% filter select
#'
#' @return A \code{AAFS} object.
AAFS_object <- function(data, x, y, nterms, Kp, err.red,
                        minspan, endspan, knotsGrid, d,
                        na.rm, Forward.AAFS, AAFS, QSmooth.AAFS,
                        CSmooth.AAFS) {

  AAFS.object <- list("data" = list(df = data,
                                    x  = x,
                                    y  = y,
                                    input_names  = colnames(data)[x],
                                    output_names = colnames(data)[y],
                                    row_names    = rownames(data)),
                      "control" = list(nterms = nterms,
                                       Kp = Kp,
                                       err.red = err.red,
                                       minspan = minspan,
                                       endspan = endspan,
                                       knotsGrid = knotsGrid,
                                       d  = d,
                                       na.rm = na.rm),
                      "Forward.AAFS" = Forward.AAFS,
                      "AAFS"         = AAFS,
                      "CSAAFS"       = CSmooth.AAFS,
                      "QSAAFS"       = QSmooth.AAFS)

  class(AAFS.object) <- "AAFS"

  return(AAFS.object)
}
