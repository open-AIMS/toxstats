#' Minimum significant difference
#'
#' The smallest difference from the control that the test could have detected
#' as significant. Reporting it alongside a no-observed-effect concentration is
#' what distinguishes a test that found no effect from a test that could not
#' have found one.
#'
#' @details
#' The minimum significant difference is
#'
#' \deqn{MSD = d\, S_w \sqrt{1/n_1 + 1/n_i}}
#'
#' where \eqn{d} is the one-sided critical value of the multiple-comparison
#' procedure, \eqn{S_w} is the square root of the within mean square from a
#' one-way analysis of variance, \eqn{n_1} is the number of control replicates
#' and \eqn{n_i} the number at concentration \eqn{i} (Appendix C, section 1.10).
#'
#' The critical value for Dunnett's procedure is computed from the multivariate
#' t distribution by numerical integration. Because all `k` comparisons share
#' the control, the correlation matrix has the one-factor form
#' \eqn{\rho_{ij} = \lambda_i \lambda_j} with
#' \eqn{\lambda_i = \sqrt{n_i / (n_1 + n_i)}}, and conditioning on that single
#' factor reduces the problem to a two-dimensional integral.
#'
#' The manuals instead reproduce a table credited to Miller (1981). Computation
#' is used here for three reasons: it agrees with the table, returning 2.3561
#' against a printed 2.36 on the Appendix C example; it extends to unbalanced
#' designs the table does not cover; and it avoids reproducing a third-party
#' table whose permission does not transfer.
#'
#' Integration is used in preference to `mvtnorm::qmvt()`, which is a randomised
#' quasi-Monte Carlo method and returns a slightly different answer on every
#' call. The same data must give the same critical value every time.
#'
#' When the design is unbalanced the difference detectable is not the same at
#' every concentration, so a value is returned for each.
#'
#' @inheritParams toxcalc_params
#' @param test The multiple-comparison procedure whose critical value is used.
#'   `"dunnett"` is the EPA method for a balanced design; `"bonferroni"` and
#'   `"sidak"` correspond to the t test with Bonferroni's or Dunn-Sidak's
#'   adjustment, used when replication is unequal.
#'
#' @return An object of class `toxcalc_msd`, a list with elements `msd` (named
#'   by concentration), `critical`, `sw`, `df`, `n_control`, `n`, `test`,
#'   `alpha` and `balanced`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix C, section 1.10.
#'
#' @seealso [pmsd()] to express this as a percentage of the control mean.
#'
#' @examples
#' # Appendix C worked example; the manual prints MSD = 0.162
#' msd(fathead_c1, response = "weight")
#'
#' @export
msd <- function(
  x,
  ...,
  alpha = 0.05,
  test = c("dunnett", "bonferroni", "sidak")
) {
  test <- match.arg(test)
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))

  x <- as_toxcalc_data(x, ...)
  parts <- anova_parts(x)

  n_control <- parts$n[parts$is_control]
  n_treat <- parts$n[!parts$is_control]
  k <- length(n_treat)
  if (!k) {
    chk::abort_chk("At least one non-control concentration is needed.")
  }

  critical <- critical_value(
    test = test,
    alpha = alpha,
    df = parts$df,
    n_control = n_control,
    n_treat = n_treat
  )

  values <- critical * parts$sw * sqrt(1 / n_control + 1 / n_treat)
  names(values) <- as.character(parts$conc[!parts$is_control])

  structure(
    list(
      msd = values,
      critical = critical,
      sw = parts$sw,
      df = parts$df,
      n_control = n_control,
      n = n_treat,
      test = test,
      alpha = alpha,
      balanced = length(unique(parts$n)) == 1L
    ),
    class = "toxcalc_msd"
  )
}

#' Percent minimum significant difference
#'
#' The minimum significant difference expressed as a percentage of the control
#' mean, which is how the EPA method-variability criteria are stated. A value
#' outside the bounds published for a test method indicates a test that was
#' either unusually insensitive or unusually sensitive.
#'
#' @details
#' `PMSD = 100 * MSD / control mean`.
#'
#' EPA publishes lower and upper bounds for three sublethal endpoints, shipped
#' here as [epa_pmsd_bounds]. Pass `bounds` to compare against one of them.
#'
#' The manuals require the percent minimum significant difference to be
#' computed parametrically **even when the flowchart selected a non-parametric
#' test**, so this function does not consult the flowchart.
#'
#' @inheritParams msd
#' @param bounds Optional. Either a row name of [epa_pmsd_bounds] such as
#'   `"fathead_growth"`, or a length-two numeric vector giving the lower and
#'   upper bounds directly. `NULL` reports the value without comparison.
#'
#' @return An object of class `toxcalc_pmsd`, a list with elements `pmsd`
#'   (named by concentration), `control_mean`, `bounds`, `status` and the
#'   underlying `msd` object.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, section 10.2.8 and Table 6.
#'
#' @examples
#' # Appendix C worked example; the manual reports about 24 per cent
#' pmsd(fathead_c1, response = "weight", bounds = "fathead_growth")
#'
#' @export
pmsd <- function(
  x,
  ...,
  alpha = 0.05,
  test = c("dunnett", "bonferroni", "sidak"),
  bounds = NULL
) {
  test <- match.arg(test)
  x <- as_toxcalc_data(x, ...)
  fit <- msd(x, alpha = alpha, test = test)

  control_mean <- x$pooled$mean[x$pooled$conc == x$control]
  values <- 100 * fit$msd / control_mean

  bounds <- resolve_bounds(bounds)
  status <- if (is.null(bounds)) {
    rep(NA_character_, length(values))
  } else {
    ifelse(
      values < bounds[1],
      "below_lower",
      ifelse(values > bounds[2], "above_upper", "within")
    )
  }
  names(status) <- names(values)

  structure(
    list(
      pmsd = values,
      control_mean = control_mean,
      bounds = bounds,
      status = status,
      msd = fit
    ),
    class = "toxcalc_pmsd"
  )
}

#' Within mean square, degrees of freedom and replication
#'
#' @param x A `toxcalc_data` object.
#' @return A list with `sw`, `df`, `conc`, `n` and `is_control`.
#' @noRd
anova_parts <- function(x) {
  y <- analysis_response(x)
  g <- factor(x$replicates$conc)
  fit <- stats::lm(y ~ g)
  conc <- as.numeric(levels(g))

  list(
    sw = stats::sigma(fit),
    df = stats::df.residual(fit),
    conc = conc,
    n = as.integer(table(g)),
    is_control = conc == x$control
  )
}

#' Resolve the bounds argument to a length-two numeric vector or NULL
#'
#' @param bounds A string naming a row of epa_pmsd_bounds, a numeric pair, or
#'   NULL.
#' @return A length-two numeric vector, or NULL.
#' @noRd
resolve_bounds <- function(bounds) {
  if (is.null(bounds)) {
    return(NULL)
  }
  if (is.character(bounds)) {
    chk::chk_string(bounds)
    chk::chk_subset(bounds, toxcalc::epa_pmsd_bounds$id)
    row <- toxcalc::epa_pmsd_bounds[
      toxcalc::epa_pmsd_bounds$id == bounds,
      ,
      drop = FALSE
    ]
    return(c(row$lower, row$upper))
  }
  chk::chk_numeric(bounds)
  if (length(bounds) != 2L) {
    chk::abort_chk("`bounds` must be a length two numeric vector.")
  }
  if (bounds[1] >= bounds[2]) {
    chk::abort_chk("`bounds` must be increasing: lower then upper.")
  }
  bounds
}

#' @export
print.toxcalc_msd <- function(x, ...) {
  cat("Minimum significant difference\n")
  cat("  EPA-821-R-02-013 Appendix C, section 1.10\n\n")
  cat(
    "  ",
    x$test,
    " critical value = ",
    format(x$critical, digits = 5),
    " (one-sided, alpha = ",
    x$alpha,
    ", df = ",
    x$df,
    ")\n",
    sep = ""
  )
  cat(
    "  within mean square root Sw = ",
    format(x$sw, digits = 4),
    "\n\n",
    sep = ""
  )
  print(x$msd)
  if (!x$balanced) {
    cat("\n  Design is unbalanced, so the MSD differs by concentration.\n")
  }
  invisible(x)
}

#' @export
print.toxcalc_pmsd <- function(x, ...) {
  cat("Percent minimum significant difference\n")
  cat("  EPA-821-R-02-013 section 10.2.8\n\n")
  cat("  control mean = ", format(x$control_mean, digits = 4), "\n\n", sep = "")
  print(round(x$pmsd, 1))
  if (!is.null(x$bounds)) {
    cat(
      "\n  EPA bounds: ",
      x$bounds[1],
      " to ",
      x$bounds[2],
      " per cent\n",
      sep = ""
    )
    print(x$status)
  }
  invisible(x)
}
