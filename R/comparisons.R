#' Build the result object shared by every hypothesis test
#'
#' @param test Short machine name of the test.
#' @param method Human-readable method name.
#' @param comparisons Data frame, one row per treatment concentration.
#' @param reference Manual section.
#' @param ... Further named elements specific to the test.
#' @return An object of class `tox_comparison`.
#' @noRd
new_comparison <- function(test, method, comparisons, reference, ...) {
  endpoints <- derive_noec_loec(
    comparisons$conc,
    comparisons$significant
  )

  structure(
    c(
      list(
        test = test,
        method = method,
        comparisons = comparisons,
        noec = endpoints$noec,
        loec = endpoints$loec,
        monotone = endpoints$monotone,
        reference = reference
      ),
      list(...)
    ),
    class = "tox_comparison"
  )
}

#' Derive the no- and lowest-observed-effect concentrations
#'
#' @param conc Treatment concentrations, control excluded.
#' @param significant Logical, one per concentration.
#' @return A list with `noec`, `loec` and `monotone`.
#' @noRd
derive_noec_loec <- function(conc, significant) {
  ordering <- order(conc)
  conc <- conc[ordering]
  significant <- significant[ordering]

  if (!any(significant)) {
    # Nothing was detected, so the highest concentration tested is the NOEC and
    # the LOEC lies above the tested range.
    return(list(noec = max(conc), loec = NA_real_, monotone = TRUE))
  }

  loec <- min(conc[significant])
  below <- conc[conc < loec]

  # Section 9.6.5.1 warns about a non-significant concentration lying between
  # two significant ones but prescribes no rule for it. The pattern is reported
  # as observed rather than forced monotone, and `monotone` records whether it
  # needs comment.
  list(
    noec = if (length(below)) max(below) else NA_real_,
    loec = loec,
    monotone = all(significant[conc >= loec])
  )
}

#' Pooled analysis of variance quantities plus per-concentration means
#'
#' @param x A `tox_data` object.
#' @return A list with `sw`, `df`, `conc`, `n`, `mean` and `control_index`.
#' @noRd
comparison_parts <- function(x) {
  y <- analysis_response(x)
  g <- factor(x$replicates$conc)
  fit <- stats::lm(y ~ g)
  conc <- as.numeric(levels(g))

  list(
    sw = stats::sigma(fit),
    df = stats::df.residual(fit),
    conc = conc,
    n = as.integer(table(g)),
    mean = as.numeric(tapply(y, g, mean)),
    control_index = which(conc == x$control),
    values = split(y, g)
  )
}

#' Sign that makes a toxic effect positive
#'
#' The EPA tests are one-sided, looking for the response to move away from the
#' control in the direction toxicity would move it.
#'
#' @param direction Either "decreasing" or "increasing".
#' @return 1 or -1.
#' @noRd
effect_sign <- function(direction) {
  if (direction == "decreasing") 1 else -1
}

#' @export
print.tox_comparison <- function(x, ...) {
  cat(x$method, "\n", sep = "")
  cat("  ", x$reference, "\n\n", sep = "")

  out <- x$comparisons
  numeric_columns <- vapply(out, is.numeric, logical(1))
  out[numeric_columns] <- lapply(out[numeric_columns], signif, digits = 4)
  print(out, row.names = FALSE)

  cat("\n")
  cat(
    "  NOEC ",
    format_endpoint(x$noec, x$comparisons$conc, "noec"),
    "\n",
    sep = ""
  )
  cat(
    "  LOEC ",
    format_endpoint(x$loec, x$comparisons$conc, "loec"),
    "\n",
    sep = ""
  )
  if (!x$monotone) {
    cat(
      "\n  Note: a concentration below the LOEC is significant while one\n",
      "  above it is not. Section 9.6.5.1 advises that such results be used\n",
      "  with caution.\n",
      sep = ""
    )
  }
  invisible(x)
}

#' Format an endpoint that may fall outside the tested range
#'
#' @param value The endpoint, possibly NA.
#' @param conc The tested concentrations.
#' @param which Either "noec" or "loec".
#' @return A string.
#' @noRd
format_endpoint <- function(value, conc, which) {
  if (!is.na(value)) {
    return(format(value))
  }
  if (which == "noec") {
    paste0("< ", format(min(conc)))
  } else {
    paste0("> ", format(max(conc)))
  }
}
