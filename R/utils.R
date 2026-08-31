#' Coerce to a toxcalc_data object
#'
#' Lets every user-facing function accept either a prepared [toxcalc_data()]
#' object or a bare data frame plus the arguments needed to prepare one.
#'
#' @param x A `toxcalc_data` object or a data frame.
#' @param ... Passed to [toxcalc_data()] when `x` is a data frame.
#' @return A `toxcalc_data` object.
#' @noRd
as_toxcalc_data <- function(x, ...) {
  if (inherits(x, "toxcalc_data")) {
    chk::chk_unused(...)
    return(x)
  }
  toxcalc_data(x, ...)
}

#' The response the hypothesis-testing flowchart operates on
#'
#' A measurement for a continuous endpoint, and the replicate proportion for a
#' quantal one. The EPA flowchart works on replicate proportions rather than
#' counts, which is why this distinction is made in one place rather than at
#' every call site.
#'
#' @param x A `toxcalc_data` object.
#' @return A numeric vector, one element per replicate.
#' @noRd
analysis_response <- function(x) {
  if (x$type == "quantal") x$replicates$proportion else x$replicates$response
}

#' Pooled within-group centred residuals
#'
#' Every observation minus the mean of its own concentration. This is what the
#' EPA manuals test for normality -- not the raw values, and not each
#' concentration separately (Appendix B, section 2.3).
#'
#' @param x A `toxcalc_data` object.
#' @return A numeric vector, one element per replicate.
#' @noRd
centred_residuals <- function(x) {
  y <- analysis_response(x)
  y - stats::ave(y, factor(x$replicates$conc), FUN = mean)
}

#' Build a result object shared by the assumption and hypothesis tests
#'
#' @param ... Named elements of the result.
#' @return An object of class `toxcalc_htest`.
#' @noRd
new_htest <- function(...) {
  structure(list(...), class = "toxcalc_htest")
}

#' @export
print.toxcalc_htest <- function(x, ...) {
  cat(x$method, "\n", sep = "")
  cat("  ", x$reference, "\n\n", sep = "")
  cat("  ", x$statistic_name, " = ", format(x$statistic, digits = 5), sep = "")
  if (!is.null(x$p_value) && !is.na(x$p_value)) {
    cat(", p = ", format(x$p_value, digits = 4), sep = "")
  }
  if (!is.null(x$critical) && !is.na(x$critical)) {
    cat(", critical value = ", format(x$critical, digits = 5), sep = "")
  }
  cat("\n")
  cat("  alpha = ", x$alpha, "\n", sep = "")
  cat("  ", x$conclusion, "\n", sep = "")
  invisible(x)
}
