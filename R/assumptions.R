#' Test the normality assumption
#'
#' Tests whether the pooled within-group centred residuals are normally
#' distributed, the first branch point of the EPA hypothesis-testing flowchart.
#' The parametric tests downstream -- Dunnett's procedure and the t test with
#' Bonferroni's adjustment -- assume normality; failing this test sends the
#' analysis to the non-parametric branch.
#'
#' @details
#' The manuals are specific about what is tested. Each observation is centred
#' by subtracting the mean of its own concentration, and the resulting
#' residuals are pooled across all concentrations and tested together
#' (Appendix B, section 2.3). Testing the raw values, or testing each
#' concentration separately, gives a different answer.
#'
#' Two statistics are provided, and `"auto"` picks between them as the manuals
#' direct: Shapiro-Wilk for 50 or fewer observations, and the Kolmogorov "D"
#' statistic above that.
#'
#' The Shapiro-Wilk statistic is computed by [stats::shapiro.test()], which
#' uses Royston's AS R94 algorithm. The manuals instead tabulate Conover's
#' coefficients and quantiles, a hand-calculation aid from before this was
#' available in software. The two agree closely: on the manual's own Appendix B
#' worked example, Royston returns `W = 0.9601` where the manual prints
#' `0.959`, and the two agree exactly once the manual's rounding of the
#' concentration means is applied. Royston's algorithm is used here because it
#' reproduces the published result, is not restricted to the tabulated sample
#' sizes, and returns a p-value rather than a single fixed cutoff.
#'
#' The Kolmogorov statistic is Stephens' (1974) modified form,
#' `D* = D (sqrt(n) - 0.01 + 0.85 / sqrt(n))`, compared against the critical
#' values in the manuals' Table B.11. It has no p-value, only a decision at the
#' tabulated alpha levels of 0.010, 0.025, 0.050, 0.100 and 0.150.
#'
#' @inheritParams toxcalc_params
#' @param method Which statistic to use. `"auto"` (the default) follows the
#'   manuals: Shapiro-Wilk for 50 or fewer observations, Kolmogorov above that.
#'
#' @return An object of class `toxcalc_htest`, a list with elements `method`,
#'   `statistic`, `statistic_name`, `p_value`, `critical`, `alpha`, `n`,
#'   `normal`, `conclusion` and `reference`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix B, sections 2.1 and 2.10.
#'
#' Royston P (1995) A remark on Algorithm AS 181: the W test for normality.
#' *Applied Statistics* 44:547-551.
#'
#' Stephens MA (1974) EDF statistics for goodness of fit and some comparisons.
#' *Journal of the American Statistical Association* 69:730-737.
#'
#' @examples
#' # Appendix B worked example; the manual prints W = 0.959 and concludes
#' # the data are normally distributed.
#' epa_normality(fathead_b1, response = "weight")
#'
#' @export
epa_normality <- function(
  x,
  ...,
  method = c("auto", "shapiro_wilk", "kolmogorov"),
  alpha_assumption = 0.01
) {
  method <- match.arg(method)
  chk::chk_number(alpha_assumption)
  chk::chk_range(alpha_assumption, c(0, 1))

  x <- as_toxcalc_data(x, ...)
  residuals <- centred_residuals(x)
  n <- length(residuals)

  if (method == "auto") {
    method <- if (n <= 50) "shapiro_wilk" else "kolmogorov"
  }

  if (method == "shapiro_wilk") {
    shapiro_wilk_result(residuals, alpha_assumption)
  } else {
    kolmogorov_result(residuals, alpha_assumption)
  }
}

#' @noRd
shapiro_wilk_result <- function(residuals, alpha) {
  if (length(residuals) < 3) {
    chk::abort_chk(
      "Shapiro-Wilk needs at least 3 observations; ",
      length(residuals),
      " were supplied."
    )
  }
  fit <- stats::shapiro.test(residuals)
  normal <- fit$p.value > alpha

  new_htest(
    method = "Shapiro-Wilk test for normality (Royston AS R94)",
    statistic = unname(fit$statistic),
    statistic_name = "W",
    p_value = fit$p.value,
    critical = NA_real_,
    alpha = alpha,
    n = length(residuals),
    normal = normal,
    conclusion = if (normal) {
      "The residuals are consistent with a normal distribution."
    } else {
      "The residuals are not normally distributed."
    },
    reference = "EPA-821-R-02-013 Appendix B, section 2.1"
  )
}

#' @noRd
kolmogorov_result <- function(residuals, alpha) {
  # EPA Table B.11. These five values are Stephens' (1974) critical points for
  # the modified statistic and are the only alpha levels the manuals tabulate,
  # so an alpha outside them cannot be honoured.
  critical_values <- c(
    "0.01" = 1.035,
    "0.025" = 0.955,
    "0.05" = 0.895,
    "0.1" = 0.819,
    "0.15" = 0.775
  )
  key <- as.character(alpha)
  if (!key %in% names(critical_values)) {
    chk::abort_chk(
      "The Kolmogorov test is tabulated only at alpha 0.01, 0.025, 0.05, ",
      "0.1 and 0.15 (Table B.11); `alpha_assumption` was ",
      alpha,
      "."
    )
  }

  n <- length(residuals)
  ordered <- sort(residuals)
  z <- (ordered - mean(residuals)) / stats::sd(residuals)
  p <- stats::pnorm(z)
  d_plus <- max(seq_len(n) / n - p)
  d_minus <- max(p - (seq_len(n) - 1) / n)
  d <- max(d_plus, d_minus)
  d_star <- d * (sqrt(n) - 0.01 + 0.85 / sqrt(n))

  critical <- unname(critical_values[key])
  normal <- d_star <= critical

  new_htest(
    method = "Kolmogorov \"D\" test for normality (Stephens 1974)",
    statistic = d_star,
    statistic_name = "D*",
    p_value = NA_real_,
    critical = critical,
    alpha = alpha,
    n = n,
    normal = normal,
    conclusion = if (normal) {
      "The residuals are consistent with a normal distribution."
    } else {
      "The residuals are not normally distributed."
    },
    reference = "EPA-821-R-02-013 Appendix B, section 2.10"
  )
}

#' Test the homogeneity of variance assumption
#'
#' Tests whether the variance of the response is the same at every
#' concentration, the second branch point of the EPA hypothesis-testing
#' flowchart. Dunnett's procedure pools the variance across concentrations, so
#' it is only valid when they are equal; failing this test sends the analysis
#' to the non-parametric branch.
#'
#' @details
#' Bartlett's test is what the manuals specify (Appendix B, section 3), and is
#' the default here. It is computed by [stats::bartlett.test()], which uses the
#' same formula the manuals set out, including the `C` correction factor.
#'
#' Two robust alternatives are offered because Bartlett's test is sensitive to
#' departures from normality, which is awkward given it is applied immediately
#' after a normality test that the data may only just have passed. Levene's
#' test is an analysis of variance on absolute deviations from the group
#' median, and the Fligner-Killeen test is its rank-based counterpart. Neither
#' is an EPA method, and using one is recorded as a departure.
#'
#' Note that the manual's own printed Bartlett statistic for the Appendix B
#' example, `B = 7.691`, is not reproducible from the raw data, which give
#' `6.836`. The difference arises because the manual computes the statistic
#' from its printed table of variances rounded to four decimal places, one of
#' which (`0.0020` for a variance of `0.002055`) is itself rounded incorrectly.
#' Both values lead to the same conclusion. This package computes from the
#' unrounded data.
#'
#' @inheritParams toxcalc_params
#' @param method Which test to use. `"bartlett"` is the EPA method and the
#'   default.
#'
#' @return An object of class `toxcalc_htest`, a list with elements `method`,
#'   `statistic`, `statistic_name`, `p_value`, `critical`, `alpha`, `df`,
#'   `homogeneous`, `conclusion` and `reference`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix B, section 3.
#'
#' Snedecor GW, Cochran WG (1980) *Statistical Methods*, 7th edition. Iowa
#' State University Press.
#'
#' @examples
#' # Appendix B worked example; the variances are not significantly different.
#' epa_variance(fathead_b1, response = "weight")
#'
#' @export
epa_variance <- function(
  x,
  ...,
  method = c("bartlett", "levene", "fligner"),
  alpha_assumption = 0.01
) {
  method <- match.arg(method)
  chk::chk_number(alpha_assumption)
  chk::chk_range(alpha_assumption, c(0, 1))

  x <- as_toxcalc_data(x, ...)
  y <- analysis_response(x)
  g <- factor(x$replicates$conc)

  if (nlevels(g) < 2) {
    chk::abort_chk(
      "At least two concentrations are needed to compare variances."
    )
  }
  if (any(table(g) < 2)) {
    chk::abort_chk(
      "Every concentration needs at least two replicates to have a variance."
    )
  }

  fit <- switch(
    method,
    bartlett = stats::bartlett.test(y, g),
    fligner = stats::fligner.test(y, g),
    levene = levene_test(y, g)
  )

  homogeneous <- fit$p.value > alpha_assumption
  label <- switch(
    method,
    bartlett = "Bartlett's test for homogeneity of variance",
    fligner = "Fligner-Killeen test for homogeneity of variance (not an EPA method)",
    levene = "Levene's test for homogeneity of variance (not an EPA method)"
  )

  new_htest(
    method = label,
    statistic = unname(fit$statistic),
    statistic_name = if (method == "levene") "F" else "B",
    p_value = unname(fit$p.value),
    critical = if (method == "levene") {
      stats::qf(1 - alpha_assumption, fit$parameter[1], fit$parameter[2])
    } else {
      stats::qchisq(1 - alpha_assumption, unname(fit$parameter))
    },
    alpha = alpha_assumption,
    df = unname(fit$parameter),
    homogeneous = homogeneous,
    conclusion = if (homogeneous) {
      "The variances are not significantly different."
    } else {
      "The variances are heterogeneous."
    },
    reference = if (method == "bartlett") {
      "EPA-821-R-02-013 Appendix B, section 3"
    } else {
      "Not an EPA method; a robust alternative to Bartlett's test"
    }
  )
}

#' Levene's test on absolute deviations from the group median
#'
#' @param y Numeric response.
#' @param g Grouping factor.
#' @return A list shaped like an htest, with `statistic`, `parameter` (a length
#'   two vector of degrees of freedom) and `p.value`.
#' @noRd
levene_test <- function(y, g) {
  # Deviations are taken from the median rather than the mean, which is the
  # Brown-Forsythe form and is the more robust of the two in common use.
  deviations <- abs(y - stats::ave(y, g, FUN = stats::median))
  fit <- stats::anova(stats::lm(deviations ~ g))
  list(
    statistic = c(F = fit[["F value"]][1]),
    parameter = c(df1 = fit[["Df"]][1], df2 = fit[["Df"]][2]),
    p.value = fit[["Pr(>F)"]][1]
  )
}
