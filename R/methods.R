#' @export
print.aces <- function(
    x,
    method = c("aces"),
    digits = 2,
    ...
    ) {

  model <- x[["methods"]][[method]]
  if (is.null(model)) {
    cat(sprintf("No modelel '%s' stored in this aces object.\n", method))
    return(invisible(x))
  }

  coefs <- as.matrix(model[["coefs"]])
  knots <- model[["knots"]]

  # ---- Output names (columns of the coefficient matrix)
  ynames <- x[["data"]][["ynames"]]
  if (is.null(colnames(coefs))) {
    colnames(coefs) <- if (!is.null(ynames)) ynames else paste0("y", seq_len(ncol(coefs)))
  }

  # ---- Expanded input names (originals + interactions) from kn_grid
  xnames_ext <- names(x[["control"]][["kn_grid"]])

  # Helper: pretty-print a single hinge with max() using expanded names
  hinge_str <- function(xi, t, side) {

    vnm <- if (!is.null(xnames_ext) && xi >= 1 && xi <= length(xnames_ext)) {
      xnames_ext[xi]
    } else {
      paste0("x", xi)
    }

    if (side == "R") {
      sprintf("max(0, %s - %.*f)", vnm, digits, t)
    } else if (side == "L") {
      sprintf("max(0, %.*f - %s)", digits, t, vnm)
    } else {
      "1"
    }

  }

  # Column names for printed coefficients
  coef_names <- paste0("Coef_", colnames(coefs))

  # Initialize output table with fixed columns to avoid name-mismatch in rbind
  out <- data.frame(
    BF = character(0),
    BaseFunction = character(0),
    matrix(numeric(0), ncol = length(coef_names)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(out) <- c("BF","BaseFunction", coef_names)

  # ---- Intercept row (single row, all outputs as columns)
  inter_vals <- round(coefs[1, , drop = FALSE], digits)
  inter_row  <- data.frame(
    BF = "BF0 (Intercept)",
    BaseFunction = "1",
    as.list(drop(inter_vals)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(inter_row)[-(1:2)] <- coef_names
  out <- rbind(out, inter_row)

  # ---- One row per BF
  for (i in seq_len(nrow(knots))) {

    krow <- knots[i, ]

    expr <- hinge_str (
      xi   = as.integer(krow$xi),
      t    = as.numeric(krow$t),
      side = as.character(krow$side)
    )

    coef_vals <- round(coefs[i + 1, , drop = FALSE], digits)

    term_row <- data.frame (
      BF = sprintf("BF%d", i),
      BaseFunction = expr,
      as.list(drop(coef_vals)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    names(term_row)[-(1:2)] <- coef_names

    out <- rbind(out, term_row)

    }

  print(out, row.names = FALSE, right = FALSE)
  invisible(x)
}
