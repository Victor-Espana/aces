#' @title The output-oriented Banker, Charnes and Cooper Programming Model for Additive Adaptive Frontier Splines
#'
#' @description This function computes efficiency scores through a BCC model with output orientation for AAFS.
#'
#' @param dmu Number of DMUs in the sample.
#' @param xmat \code{data.frame} or \code{matrix} containing the inputs in the model.
#' @param ymat \code{data.frame} or \code{matrix} containing the outputs in the model.
#' @param aknot Knots coordinates.
#' @param yknot Predictions for knots.
#' @param nX Number of inputs.
#' @param nY Number of outputs.
#' @param nk Number of knots.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A numerical vector with the efficiency scores computed.
AAFS_BCC.OUT <- function(dmu, xmat, ymat, aknot, yknot, nX, nY, nk) {

  scores <- matrix(nrow = dmu, ncol = 1)

  for(d in 1:dmu){

    objVal <- matrix(ncol = nk + 1, nrow = 1)
    objVal[1] <- 1

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = nk + 1)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX) {
      add.constraint(lps, xt = c(0, aknot[, xi]), "<=",  rhs = xmat[d, xi])
    }

    for(yi in 1:nY) {
      add.constraint(lps, xt = c(- ymat[d, yi], yknot[, yi]), ">=", rhs = 0)
    }

    # Constrain 2.3 - phi = 1
    add.constraint(lprec = lps, xt = c(0, rep(1, nk)), type = "=", rhs = 1)

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The input-oriented Banker, Charnes and Cooper Programming Model for Additive Adaptive Frontier Splines
#'
#' @description This function computes efficiency scores through a BCC model with input orientation for AAFS.
#'
#' @param dmu Number of DMUs in the sample.
#' @param xmat \code{data.frame} or \code{matrix} containing the inputs in the model.
#' @param ymat \code{data.frame} or \code{matrix} containing the outputs in the model.
#' @param aknot Knots coordinates.
#' @param yknot Predictions for knots.
#' @param nX Number of inputs.
#' @param nY Number of outputs.
#' @param nk Number of knots.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A numerical vector with the efficiency scores computed.
AAFS_BCC.INP <- function(dmu, xmat, ymat, aknot, yknot, nX, nY, nk) {

  scores <- matrix(nrow = dmu, ncol = 1)

  for(d in 1:dmu){

    objVal <- matrix(ncol = nk + 1, nrow = 1)
    objVal[1] <- 1

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = nk + 1)
    lp.control(lps, sense = 'min')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX) {
      add.constraint(lps, xt = c(- xmat[d, xi], aknot[, xi]), "<=",  rhs = 0)
    }

    for(yi in 1:nY) {
      add.constraint(lps, xt = c(0, yknot[, yi]), ">=", rhs = ymat[d, yi])
    }

    # Constrain 2.3 - lambda = 1
    add.constraint(lprec = lps, xt = c(0, rep(1, nk)), type = "=", rhs = 1)

    # Constrain 2.4
    set.type(lps, columns = 1:nk + 1, type = c("binary"))

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The Directional Distance Function Programming Model for Additive Adaptive Frontier Splines
#'
#' @description This function computes efficiency scores through a DDF model for AAFS.
#'
#' @param dmu Number of DMUs in the sample.
#' @param xmat \code{data.frame} or \code{matrix} containing the inputs in the model.
#' @param ymat \code{data.frame} or \code{matrix} containing the outputs in the model.
#' @param aknot Knots coordinates.
#' @param yknot Predictions for knots.
#' @param nX Number of inputs.
#' @param nY Number of outputs.
#' @param nk Number of knots.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A numerical vector with the efficiency scores computed.
AAFS_DDF <- function(dmu, xmat, ymat, aknot, yknot, nX, nY, nk) {

  scores <- matrix(nrow = dmu, ncol = 1)

  for(d in 1:dmu){

    objVal <- matrix(ncol = 1 + nk, nrow = 1) # beta + lambdas
    objVal[1] <- 1 # beta

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = nk + 1)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX) { # beta[g-], a, <=, x
      add.constraint(lps, xt = c(xmat[d, xi], aknot[, xi]), "<=",  rhs = xmat[d, xi])
    }

    for(yi in 1:nY) { # - y, d(a), >=, beta[g+]
      add.constraint(lps, xt = c(- ymat[d, yi], yknot[, yi]), ">=", rhs = ymat[d, yi])
    }

    # Constrain 2.3 - lambda = 1
    add.constraint(lprec = lps, xt = c(0, rep(1, nk)), type = "=", rhs = 1)

    # Constrain 2.4
    set.type(lps, columns = 1:nk + 1, type = c("binary"))

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The input-oriented Russell Model for Additive Adaptive Frontier Splines
#'
#' @description This function computes efficiency scores through a Russell model with input orientation for AAFS.
#'
#' @param dmu Number of DMUs in the sample.
#' @param xmat \code{data.frame} or \code{matrix} containing the inputs in the model.
#' @param ymat \code{data.frame} or \code{matrix} containing the outputs in the model.
#' @param aknot Knots coordinates.
#' @param yknot Predictions for knots.
#' @param nX Number of inputs.
#' @param nY Number of outputs.
#' @param nk Number of knots.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A numerical vector with efficiency scores.
AAFS_RSL.INP <- function(dmu, xmat, ymat, aknot, yknot, nX, nY, nk) {

  scores <- matrix(nrow = dmu, ncol = 1)

  for(d in 1:dmu){

    objVal <- matrix(ncol = nk + nX, nrow = 1)
    objVal[1:nX] <- 1 / nX

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = nk + nX)
    lp.control(lps, sense = 'min')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX) {
      vec                     <- c()
      vec[xi]                 <- - xmat[d, xi]
      vec[(1:nX)[- xi]]       <- 0
      vec[(nX + 1):(nX + nk)] <- aknot[, xi]

      add.constraint(lps, xt = vec, "<=",  rhs = 0)
    }

    for(yi in 1:nY) {
      add.constraint(lps, xt = c(rep(0, nX), yknot[, yi]), ">=", rhs = ymat[d, yi])
    }

    # Constrain 2.3 - lambda = 1
    add.constraint(lprec = lps, xt = c(rep(0, nX), rep(1, nk)), type = "=", rhs = 1)

    # Constrain 2.4
    set.type(lps, columns = 1:nk + nX, type = c("binary"))

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The output-oriented Russell Model for Additive Adaptive Frontier Splines
#'
#' @description This function computes efficiency scores through a Russell model with output orientation for AAFS.
#'
#' @param dmu Number of DMUs in the sample.
#' @param xmat \code{data.frame} or \code{matrix} containing the inputs in the model.
#' @param ymat \code{data.frame} or \code{matrix} containing the outputs in the model.
#' @param aknot Knots coordinates.
#' @param yknot Predictions for knots.
#' @param nX Number of inputs.
#' @param nY Number of outputs.
#' @param nk Number of knots.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A numerical vector with efficiency scores.
AAFS_RSL.OUT <- function(dmu, xmat, ymat, aknot, yknot, nX, nY, nk) {

  scores <- matrix(nrow = dmu, ncol = 1)

  for(d in 1:dmu){

    objVal <- matrix(ncol = nk + nY, nrow = 1)
    objVal[1:nY] <- 1 / nY

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = nk + nY)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX) {
      add.constraint(lps, xt = c(rep(0, nY), aknot[, xi]), "<=",  rhs = xmat[d, xi])
    }

    for(yi in 1:nY) {
      vec                     <- c()
      vec[yi]                 <- - ymat[d, yi]
      vec[(1:nY)[- yi]]       <- 0
      vec[(nY + 1):(nY + nk)] <- yknot[, yi]

      add.constraint(lps, xt = vec, ">=", rhs = 0)
    }

    # Constrain 2.3 - lambda = 1
    add.constraint(lprec = lps, xt = c(rep(0, nY), rep(1, nk)), type = "=", rhs = 1)

    # Constrain 2.4
    set.type(lps, columns = 1:nk + nY, type = c("binary"))

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title The Weighted Additive Model for Additive Adaptive Frontier Splines
#'
#' @description This function computes efficiency scores through a WAM model for AAFS.
#'
#' @param dmu Number of DMUs in the sample.
#' @param xmat \code{data.frame} or \code{matrix} containing the inputs in the model.
#' @param ymat \code{data.frame} or \code{matrix} containing the outputs in the model.
#' @param aknot Knots coordinates.
#' @param yknot Predictions for knots.
#' @param nX Number of inputs.
#' @param nY Number of outputs.
#' @param nk Number of knots.
#' @param weights Character. \code{"MIP"} for Measure of Inefficiency Proportion or \code{"RAM"} for Range Adjusted Measure of Inefficiency.
#'
#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint set.type set.bounds get.objective
#'
#' @return A numerical vector with efficiency scores.
AAFS_WAM <- function(dmu, xmat, ymat, aknot, yknot, nX, nY, nk, weights) {

  # Range for RAM measures
  if (weights == "RAM") {
    InpRanges <- apply(xmat, 2, max) - apply(xmat, 2, min)
    OutRanges <- apply(ymat, 2, max) - apply(ymat, 2, min)

    ranges    <- c(InpRanges, OutRanges) / (nX + nY)
  }

  scores <- matrix(nrow = dmu, ncol = 1)

  for(d in 1:dmu){

    objVal <- matrix(ncol = nX + nY + nk, nrow = 1)

    if (weights == "MIP") {
      objVal[1:(nX + nY)] <- c(1 / xmat[d, ], 1 / ymat[d, ])

    } else if (weights == "RAM"){
      objVal[1:(nX + nY)] <- ranges

    }

    # structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = nX + nY + nk)
    lp.control(lps, sense = 'max')
    set.objfn(lps, objVal)

    # constrain 2.1 and 2.2
    for(xi in 1:nX) {
      vec                               <- c()
      vec[xi]                           <- 1
      vec[(1:nX)[- xi]]                 <- 0
      vec[(nX + 1):(nX + nY)]           <- 0
      vec[(nX + nY + 1):(nY + nX + nk)] <- aknot[, xi]

      add.constraint(lps, xt = vec, "=",  rhs = xmat[d, xi])
    }

    for(yi in 1:nY) {
      vec                               <- c()
      vec[1:nX]                         <- 0
      vec[nX + yi]                      <- - 1
      vec[((nX + 1):(nX + nY))[- yi]]   <- 0
      vec[(nX + nY + 1):(nY + nX + nk)] <- yknot[, yi]

      add.constraint(lps, xt = vec, "=", rhs = ymat[d, yi])
    }

    # Constrain 2.3 - lambda = 1
    add.constraint(lprec = lps, xt = c(rep(0, nY + nX), rep(1, nk)), type = "=", rhs = 1)

    # Constrain 2.4
    set.type(lps, columns = 1:nk + (nX + nY), type = c("binary"))

    solve(lps)
    scores[d, ] <- get.objective(lps)
  }

  return(scores)
}

#' @title Efficiency Scores computed through a Additive Adaptive Frontier Splines model.
#'
#' @description This function computes the efficiency scores for each DMU through a Additive Adaptive Frontier Splines model.
#'
#' @param data \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column input indexes in \code{data}.
#' @param y Column output indexes in \code{data}.
#' @param object An \code{AAFS} object.
#' @param Gdepth Number of partitions into which the data is divided (in each dimension of \code{x}) to create the data sample to perform DEA.
#' @param ModelPred Model to predict the virtual data.
#' @param scores_model Mathematical programming model to calculate scores.
#' \itemize{
#' \item{\code{BCC.OUT}} BCC model. Output-oriented. Efficiency level at 1.
#' \item{\code{BCC.INP}} BCC model. Input-oriented. Efficiency level at 1.
#' \item{\code{DDF}}     Directional Distance Function. Efficiency level at 0.
#' \item{\code{RSL.OUT}} Russell model. Output-oriented. Efficiency level at 1.
#' \item{\code{RSL.INP}} Russell model. Input-oriented. Efficiency level at 1.
#' \item{\code{WAM.MIP}} Weighted Additive Model. Measure of Inefficiency Proportions. Efficiency level at 0.
#' \item{\code{WAM.RAM}} Weighted Additive Model. Range Adjusted Measure of Inefficiency. Efficiency level at 0.
#' }
#' @param digits Decimal units for scores.
#' @param DEA \code{logical}. If \code{TRUE}, DEA scores are also computed with the programming model selected in \code{scores_model}.
#' @param print.table \code{logical}. If \code{TRUE}, a summary descriptive table of the efficiency scores is displayed.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#'
#' @importFrom dplyr summarise %>% mutate_if
#' @importFrom stats median quantile sd
#'
#' @export
#'
#' @return A \code{data.frame} with the efficiency scores computed through a Additive Adaptive Frontier Splines model.
AAFS.scores <- function(data, x, y, object, Gdepth, ModelPred,
                        scores_model, digits = 3, DEA = TRUE,
                        print.table = FALSE, na.rm = TRUE) {

  # Possible errors
  if (!is(object, "AAFS")) {
    stop(paste(deparse(substitute(object)), "must be an AAFS object."))

  } else if (digits < 0) {
    stop(paste('digits =', digits, 'must be greater than 0.'))

  } else if (!scores_model %in% c("BCC.OUT", "BCC.INP", "DDF", "RSL.OUT", "RSL.INP",
                                  "WAM.MIP","WAM.RAM")) {
    stop(paste(scores_model, "is not available. Please, check help(\"AAFS.scores\")"))
  }

  data <- preProcess(data, x, y, na.rm = na.rm)

  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  train_names <- c(object[["data"]][["input_names"]], object[["data"]][["output_names"]])

  # Possible errors
  if (!identical(sort(train_names), sort(colnames(data)))) {
    stop("Different variable names in training and data.")
  }

  dmu  <- nrow(data)
  xmat <- as.matrix(data[, x])
  ymat <- as.matrix(data[, y])
  nX   <- length(x)
  nY   <- length(y)

  # Data for Efficiency Models: from forward model
  knots.list <- vector("list", length(x))

  for (v in 1:length(x)) {
    knots.list[[v]] <- seq(min(data[, v]), max(data[, v]), length.out = Gdepth)
  }

  aknots <- expand.grid(knots.list)
  names(aknots) <- colnames(data)[x]

  # Model prediction
  yknots <- predict(object, aknots, x, ModelPred)

  nk <- nrow(aknots)

  if (scores_model == "BCC.OUT") {
    scores <- AAFS_BCC.OUT(dmu, xmat, ymat, aknots, yknots, nX, nY, nk)
    model  <- "AAFS_BCC.OUT"

    if (DEA == TRUE) {
      DEA.scores <- AAFS_BCC.OUT(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu)
      DEA.model  <- "DEA_BCC.OUT"
    }

  } else if (scores_model == "BCC.INP"){
    scores <- AAFS_BCC.INP(dmu, xmat, ymat, aknots, yknots, nX, nY, nk)
    model  <- "AAFS_BCC.INP"

    if (DEA == TRUE) {
      DEA.scores <- AAFS_BCC.INP(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu)
      DEA.model  <- "DEA_BCC.INP"
    }

  } else if (scores_model == "DDF"){
    scores <- AAFS_DDF(dmu, xmat, ymat, aknots, yknots, nX, nY, nk)
    model  <- "AAFS_DDF"

    if (DEA == TRUE) {
      DEA.scores <- AAFS_DDF(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu)
      DEA.model  <- "DEA_DDF"
    }

  } else if (scores_model == "RSL.OUT"){
    scores <- AAFS_RSL.OUT(dmu, xmat, ymat, aknots, yknots, nX, nY, nk)
    model  <- "AAFS_RSL.OUT"

    if (DEA == TRUE) {
      DEA.scores <- AAFS_RSL.OUT(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu)
      DEA.model  <- "DEA_RSL.OUT"
    }

  } else if (scores_model == "RSL.INP"){
    scores <- AAFS_RSL.INP(dmu, xmat, ymat, aknots, yknots, nX, nY, nk)
    model  <- "AAFS_RSL.INP"

    if (DEA == TRUE) {
      DEA.scores <- AAFS_RSL.INP(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu)
      DEA.model  <- "DEA_RSL.INP"
    }

  } else if (scores_model == "WAM.MIP"){
    scores <- AAFS_WAM(dmu, xmat, ymat, aknots, yknots, nX, nY, nk, "MIP")
    model  <- "AAFS_WAM.MIP"

    if (DEA == TRUE) {
      DEA.scores <- AAFS_WAM(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu, "MIP")
      DEA.model  <- "DEA_WAM.MIP"
    }

  } else if (scores_model == "WAM.RAM") {
    scores <- AAFS_WAM(dmu, xmat, ymat, aknots, yknots, nX, nY, nk, "RAM")
    model  <- "FAAFS_WAM.RAM"

    if (FDH == TRUE){
      DEA.scores <- AAFS_WAM(dmu, xmat, ymat, xmat, ymat, nX, nY, dmu, "RAM")
      DEA.model <- "DEA_WAM.RAM"
    }
  }

  scores <- as.data.frame(scores)
  names(scores) <- model
  rownames(scores) <- row.names(data)

  descriptive <- scores %>%
    summarise("Model"     = "AAFS",
              "Mean"      = mean(scores[, 1]),
              "Std. Dev." = sd(scores[, 1]),
              "Min"       = min(scores[, 1]),
              "Q1"        = quantile(scores[, 1])[[2]],
              "Median"    = median(scores[, 1]),
              "Q3"        = quantile(scores[, 1])[[3]],
              "Max"       = max(scores[, 1])) %>%
    mutate_if(is.numeric, round, digits)

  if (DEA == TRUE){

    DEA.scores          <- as.data.frame(DEA.scores)
    names(DEA.scores)   <- DEA.model
    rownames(DEA.model) <- row.names(data)

    descriptive[2, ] <- DEA.scores %>%
      summarise("Model"     = "DEA",
                "Mean"      = mean(DEA.scores[, 1]),
                "Std. Dev." = sd(DEA.scores[, 1]),
                "Min"       = min(DEA.scores[, 1]),
                "Q1"        = quantile(DEA.scores[, 1])[[2]],
                "Median"    = median(DEA.scores[, 1]),
                "Q3"        = quantile(DEA.scores[, 1])[[3]],
                "Max"       = max(DEA.scores[, 1])) %>%
      mutate_if(is.numeric, round, digits)

    scores.df <- cbind(data, round(scores, digits), round(DEA.scores, digits))

    if (print.table == TRUE) {
      print(descriptive, row.names = FALSE)
      cat("\n")
    }

    return(scores.df[, c(ncol(scores.df) - 1, ncol(scores.df))])

  } else {

    scores.df <- cbind(data, round(scores, digits))

    if (print.table == TRUE) {
      print(descriptive, row.names = FALSE)
      cat("\n")
    }

    return(round(scores.df[, ncol(scores.df)], digits))
  }
}

