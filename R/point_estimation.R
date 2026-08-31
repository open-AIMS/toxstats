#' Prepare quantal data for point estimation
#'
#' Pools counts within each concentration, smooths the proportions to be
#' monotone non-decreasing, and adjusts for control response with Abbott's
#' formula. All four EPA point-estimation methods begin here.
#'
#' @param x A `tox_data` object of type "quantal".
#' @return A list with `conc`, `x` (log10 concentration), `n`, `observed`,
#'   `smoothed`, `adjusted` and `control_smoothed`, the first six excluding the
#'   control.
#' @noRd
quantal_prep <- function(x) {
  if (x$type != "quantal") {
    chk::abort_chk(
      "Point estimation needs quantal data; call `tox_data()` with ",
      "`type = \"quantal\"` and an `n_exposed` column."
    )
  }

  pooled <- x$pooled[order(x$pooled$conc), , drop = FALSE]
  control <- which(pooled$conc == x$control)
  if (control != 1L) {
    chk::abort_chk(
      "The control must be the lowest concentration for point estimation, ",
      "because the smoothing step orders the response by concentration."
    )
  }
  if (any(pooled$conc[-control] <= 0)) {
    chk::abort_chk(
      "Concentrations must be positive apart from the control, because the ",
      "methods work on the log10 of concentration."
    )
  }

  smoothed <- smooth_monotone(pooled$proportion, direction = "increasing")
  adjusted <- abbott(smoothed, p_control = smoothed[control])
  treatments <- seq_len(nrow(pooled))[-control]

  list(
    conc = pooled$conc[treatments],
    x = log10(pooled$conc[treatments]),
    n = pooled$n_exposed[treatments],
    observed = pooled$proportion[treatments],
    smoothed = smoothed[treatments],
    adjusted = adjusted[treatments],
    control_observed = pooled$proportion[control],
    control_smoothed = smoothed[control]
  )
}

#' Build the result object shared by the point-estimation methods
#'
#' @noRd
new_estimate <- function(method, estimates, working, reference, ...) {
  structure(
    c(
      list(
        method = method,
        estimates = estimates,
        working = working,
        reference = reference
      ),
      list(...)
    ),
    class = "tox_estimate"
  )
}

#' The multiplier EPA applies to a standard error at a given confidence level
#'
#' Section 11.2.3.3 step 6 writes the interval as `m +/- 2.0 sqrt(V(m))`. Two
#' is the manual's rounding of 1.96, and using it is what reproduces the
#' printed intervals, so it is used at the 95 per cent level throughout and the
#' exact normal quantile is used elsewhere.
#'
#' @noRd
epa_multiplier <- function(ci_level) {
  if (isTRUE(all.equal(ci_level, 0.95))) {
    2
  } else {
    stats::qnorm(1 - (1 - ci_level) / 2)
  }
}

#' Estimate the LC50 by the graphical method
#'
#' Interpolates between the two concentrations that bracket 50 per cent
#' response, on a logarithmic concentration scale. The EPA flowchart selects
#' this method when there are no partial responses, so no model can be fitted.
#'
#' @details
#' The proportions are smoothed and adjusted for control response first, then
#' the estimate is read off the straight line joining the two bracketing points
#' on semi-logarithmic axes, which is linear interpolation in `log10(conc)`.
#'
#' No confidence interval is available. The manual gives none, and with no
#' partial response there is no information from which to construct one.
#'
#' @inheritParams toxstats_params
#'
#' @return An object of class `tox_estimate`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-012, section 11.2.2.
#'
#' @examples
#' # Acute Table 20, all-or-nothing column; the manual reads 35 per cent off
#' # the plot, and interpolation gives 35.4.
#' graphical_lc50(
#'   acute_table20,
#'   response = "graphical", n_exposed = "exposed", type = "quantal"
#' )
#'
#' @export
graphical_lc50 <- function(x, ...) {
  x <- as_tox_data(x, ...)
  parts <- quantal_prep(x)

  if (!any(parts$adjusted < 0.5) || !any(parts$adjusted >= 0.5)) {
    chk::abort_chk(
      "The smoothed, adjusted response proportions must bracket 0.5; they ",
      "run from ",
      signif(min(parts$adjusted), 3),
      " to ",
      signif(max(parts$adjusted), 3),
      "."
    )
  }

  below <- max(which(parts$adjusted < 0.5))
  estimate <- 10^(parts$x[below] +
    (0.5 - parts$adjusted[below]) *
      (parts$x[below + 1] - parts$x[below]) /
      (parts$adjusted[below + 1] - parts$adjusted[below]))

  new_estimate(
    method = "Graphical method",
    estimates = data.frame(
      endpoint = "LC50",
      p = 50,
      estimate = estimate,
      lower = NA_real_,
      upper = NA_real_,
      ci_method = NA_character_,
      stringsAsFactors = FALSE
    ),
    working = working_frame(parts),
    reference = "EPA-821-R-02-012 section 11.2.2"
  )
}

#' Estimate the LC50 by the Spearman-Karber method
#'
#' Estimates the mean of the distribution of the log10 tolerance, which is the
#' median when that distribution is symmetric. The EPA flowchart selects this
#' method when partial responses occur but the probit model does not fit, and
#' the response runs from zero at the lowest concentration to complete at the
#' highest.
#'
#' @details
#' After smoothing and Abbott adjustment,
#'
#' \deqn{m = \sum_{i=1}^{k-1} (p^a_{i+1} - p^a_i)\,(X_i + X_{i+1})/2}
#'
#' with \eqn{X_i} the log10 concentration, and
#'
#' \deqn{V(m) = \sum_{i=2}^{k-1} \frac{p^a_i (1 - p^a_i)(X_{i+1} -
#'   X_{i-1})^2}{4(n_i - 1)}}
#'
#' The interval is `m ± 2 sqrt(V(m))` on the log scale, back-transformed. The
#' manual writes the multiplier as 2.0 rather than 1.96; that rounding is kept
#' because it is what reproduces the printed interval.
#'
#' The method requires the adjusted proportion to be zero at the lowest
#' concentration and one at the highest. When it is not, the flowchart directs
#' the analysis to [trimmed_spearman_karber()].
#'
#' @inheritParams toxstats_params
#'
#' @inherit graphical_lc50 return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-012, section 11.2.3.
#'
#' Finney DJ (1978) *Statistical Method in Biological Assay*, 3rd edition.
#' Charles Griffin, London.
#'
#' @examples
#' # Acute Table 20; the manual gives m = 1.656527, V(m) = 0.0010977 and an
#' # LC50 of 45.3 per cent with limits 38.9 and 52.8.
#' spearman_karber(
#'   acute_table20,
#'   response = "spearman_karber", n_exposed = "exposed", type = "quantal"
#' )
#'
#' @export
spearman_karber <- function(x, ..., ci_level = 0.95) {
  chk::chk_number(ci_level)
  chk::chk_range(ci_level, c(0, 1))

  x <- as_tox_data(x, ...)
  parts <- quantal_prep(x)
  k <- length(parts$adjusted)

  if (k < 3) {
    chk::abort_chk(
      "At least three concentrations besides the control are needed."
    )
  }
  if (parts$adjusted[1] != 0 || parts$adjusted[k] != 1) {
    chk::abort_chk(
      "The Spearman-Karber method requires a smoothed, adjusted proportion ",
      "of 0 at the lowest concentration and 1 at the highest; they are ",
      signif(parts$adjusted[1], 3),
      " and ",
      signif(parts$adjusted[k], 3),
      ". Use `trimmed_spearman_karber()` instead."
    )
  }

  m <- sum(diff(parts$adjusted) * (parts$x[-k] + parts$x[-1]) / 2)
  interior <- 2:(k - 1)
  variance <- sum(
    parts$adjusted[interior] *
      (1 - parts$adjusted[interior]) *
      (parts$x[interior + 1] - parts$x[interior - 1])^2 /
      (4 * (parts$n[interior] - 1))
  )
  half <- epa_multiplier(ci_level) * sqrt(variance)

  new_estimate(
    method = "Spearman-Karber method",
    estimates = data.frame(
      endpoint = "LC50",
      p = 50,
      estimate = 10^m,
      lower = 10^(m - half),
      upper = 10^(m + half),
      ci_method = "normal on the log10 scale",
      stringsAsFactors = FALSE
    ),
    working = working_frame(parts),
    reference = "EPA-821-R-02-012 section 11.2.3",
    m = m,
    variance = variance,
    ci_level = ci_level
  )
}

#' Estimate the LC50 by the trimmed Spearman-Karber method
#'
#' Estimates a trimmed mean of the log10 tolerance distribution, so that the
#' estimate can be formed when the response does not run all the way from zero
#' to complete. The EPA flowchart selects this method when the probit model
#' does not fit and the untrimmed Spearman-Karber requirements are not met.
#'
#' @details
#' ## The trim
#'
#' Section 11.2.4.3 step 4 defines it as
#'
#' \deqn{\mathrm{trim} = \max(p^a_1,\; 1 - p^a_k)}
#'
#' where the proportions are the smoothed, Abbott-adjusted ones. The ordering
#' matters: smoothing first, then the Abbott adjustment, then the trim.
#' Computing it from the raw proportions gives a different answer. This
#' definition reproduces the 20.51 per cent printed in the acute manual and the
#' 20.41 per cent printed in Appendix K of the chronic manual.
#'
#' `trim = NULL`, the default, is the automatic trim ToxCalc advertised.
#'
#' ## The estimate
#'
#' The concentrations at which the adjusted response crosses the trim and its
#' complement are found by interpolation on the log10 scale, the interior
#' points are rescaled to run from zero to one, and the Spearman-Karber formula
#' is applied to that set. On the acute Table 20 data this gives 77.1105
#' against the 77.11 the manual's program prints.
#'
#' ## The interval
#'
#' **This is the one place where the package does not reproduce the manual.**
#' The manual delegates the interval to a program whose source is not
#' available, and it does not state the formula. Hamilton et al. (1977) gave a
#' variance expression that Hamilton et al. (1978) then corrected, and the
#' correction is not retrievable from any public source found.
#'
#' The interval here is therefore obtained by the delta method applied to the
#' estimator actually implemented, which is exact for that estimator rather
#' than adopted on faith. On the acute Table 20 data it gives 69.5 to 85.6
#' against the 69.74 to 85.26 the manual prints, a difference in the third
#' significant figure. The point estimate is unaffected.
#'
#' @inheritParams toxstats_params
#' @param trim The proportion to trim from each tail, or `NULL` for the
#'   automatic trim defined above.
#'
#' @inherit graphical_lc50 return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-012, section 11.2.4.
#'
#' Hamilton MA, Russo RC, Thurston RV (1977) Trimmed Spearman-Karber method for
#' estimating median lethal concentrations in toxicity bioassays.
#' *Environmental Science & Technology* 11:714-719, with the correction at
#' 12:417 (1978).
#'
#' @examples
#' # Acute Table 20; the manual's program prints a trim of 20.51 per cent and
#' # an LC50 of 77.11.
#' trimmed_spearman_karber(
#'   acute_table20,
#'   response = "trimmed", n_exposed = "exposed", type = "quantal"
#' )
#'
#' @export
trimmed_spearman_karber <- function(x, ..., trim = NULL, ci_level = 0.95) {
  chk::chk_number(ci_level)
  chk::chk_range(ci_level, c(0, 1))
  chk::chk_null_or(trim, vld = chk::vld_number)

  x <- as_tox_data(x, ...)
  parts <- quantal_prep(x)

  automatic <- is.null(trim)
  if (automatic) {
    trim <- auto_trim(parts$adjusted)
  }
  chk::chk_range(trim, c(0, 0.5))
  if (trim >= 0.5) {
    chk::abort_chk(
      "The trim must be below 0.5; it is ",
      signif(trim, 4),
      ", which means the response never brackets 50 per cent."
    )
  }

  m <- tsk_estimate(parts$x, parts$adjusted, trim)

  # The delta method on the estimator as implemented. The smoothing and the
  # Abbott adjustment are not differentiable everywhere, so the derivative is
  # taken numerically with respect to the observed counts.
  gradient <- vapply(
    seq_along(parts$observed),
    function(i) {
      h <- 1e-6
      # A proportion of exactly 0 or 1 is common in these designs, and
      # perturbing past the boundary is not a proportion at all. The step is
      # clipped and the divisor taken from the step actually used, so the
      # difference becomes one-sided at a boundary rather than wrong.
      upper <- min(parts$observed[i] + h, 1)
      lower <- max(parts$observed[i] - h, 0)
      if (upper == lower) {
        return(0)
      }
      at <- function(value) {
        perturbed <- parts$observed
        perturbed[i] <- value
        recompute_tsk(
          perturbed,
          parts$control_observed,
          parts$x,
          trim,
          automatic
        )
      }
      (at(upper) - at(lower)) / (upper - lower)
    },
    numeric(1)
  )
  variance <- sum(
    gradient^2 * parts$observed * (1 - parts$observed) / parts$n
  )
  half <- epa_multiplier(ci_level) * sqrt(variance)

  new_estimate(
    method = paste0(
      "Trimmed Spearman-Karber method (",
      if (automatic) "automatic " else "",
      "trim ",
      signif(100 * trim, 4),
      " per cent)"
    ),
    estimates = data.frame(
      endpoint = "LC50",
      p = 50,
      estimate = 10^m,
      lower = 10^(m - half),
      upper = 10^(m + half),
      ci_method = "delta method on the log10 scale",
      stringsAsFactors = FALSE
    ),
    working = working_frame(parts),
    reference = "EPA-821-R-02-012 section 11.2.4",
    trim = trim,
    automatic_trim = automatic,
    m = m,
    variance = variance,
    ci_level = ci_level
  )
}

#' The automatic trim of section 11.2.4.3 step 4
#'
#' @param adjusted Smoothed, Abbott-adjusted proportions, control excluded.
#' @return A single number.
#' @noRd
auto_trim <- function(adjusted) {
  max(adjusted[1], 1 - adjusted[length(adjusted)])
}

#' The trimmed Spearman-Karber point estimate on the log10 scale
#'
#' @param x Log10 concentrations.
#' @param adjusted Smoothed, Abbott-adjusted proportions.
#' @param trim The trim proportion.
#' @return The estimate of log10(LC50).
#' @noRd
tsk_estimate <- function(x, adjusted, trim) {
  crossing <- function(target) {
    below <- which(adjusted <= target)
    if (!length(below)) {
      return(x[1])
    }
    j <- max(below)
    if (j == length(adjusted) || adjusted[j] == target) {
      return(x[j])
    }
    x[j] +
      (target - adjusted[j]) *
        (x[j + 1] - x[j]) /
        (adjusted[j + 1] - adjusted[j])
  }

  lower <- crossing(trim)
  upper <- crossing(1 - trim)
  interior <- adjusted > trim & adjusted < 1 - trim

  points_x <- c(lower, x[interior], upper)
  points_p <- c(trim, adjusted[interior], 1 - trim)
  rescaled <- (points_p - trim) / (1 - 2 * trim)

  sum(
    diff(rescaled) *
      (points_x[-length(points_x)] + points_x[-1]) /
      2
  )
}

#' Recompute the trimmed estimate from perturbed proportions
#'
#' Used by the delta method, so it must repeat every step the estimator takes,
#' including the smoothing and the Abbott adjustment.
#'
#' @noRd
recompute_tsk <- function(observed, control, x, trim, automatic) {
  smoothed <- smooth_monotone(c(control, observed), direction = "increasing")
  adjusted <- abbott(smoothed, p_control = smoothed[1])[-1]
  if (automatic) {
    trim <- auto_trim(adjusted)
  }
  if (trim >= 0.5) {
    return(NA_real_)
  }
  tsk_estimate(x, adjusted, trim)
}

#' Estimate lethal concentrations by the probit method
#'
#' Fits a probit regression of response on log10 concentration and reads the
#' requested lethal concentrations off it. The EPA flowchart selects this
#' method when at least two partial responses occur and the fit is adequate.
#'
#' @details
#' The model is `glm(family = binomial(link = "probit"))` on `log10(conc)`.
#' Where the control response is not zero the proportions are Abbott-adjusted
#' first, as the manual's own output does.
#'
#' Adequacy of fit is judged by the **Pearson** chi-square statistic for
#' heterogeneity, compared with the chi-square distribution on `k - 2` degrees
#' of freedom. On the manual's worked example this gives 3.076 against the
#' 3.076 printed, confirming that the statistic is Pearson's and not the
#' deviance, which is 3.859 for the same fit.
#'
#' Confidence limits are Fieller's. That reproduces the manual's printed limits
#' exactly on the worked example: 18.787 to 27.846 for the LC50 and 4.147 to
#' 10.959 for the LC1. The delta method does not, giving 19.04 to 27.47 and
#' 5.17 to 12.14, so the choice is settled by evidence rather than convention.
#'
#' When the heterogeneity chi-square is significant the manual's practice is to
#' inflate the variance by `chi-square / df`. The worked example does not
#' exercise this, its chi-square being well below the tabular value, so the
#' behaviour is offered under `heterogeneity` and defaults to applying it.
#'
#' @inheritParams toxstats_params
#' @param heterogeneity Should the variance be inflated by `chi-square / df`
#'   when the heterogeneity chi-square is significant at `alpha`? A flag.
#'
#' @inherit graphical_lc50 return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-012, section 11.2.5.
#'
#' Finney DJ (1978) *Statistical Method in Biological Assay*, 3rd edition.
#'
#' @examples
#' # Acute Table 20; the manual prints LC50 22.872 (18.787, 27.846) and
#' # LC1 7.924 (4.147, 10.959), with a chi-square of 3.076.
#' probit_lc(
#'   acute_table20,
#'   response = "probit", n_exposed = "exposed", type = "quantal",
#'   p = c(1, 50)
#' )
#'
#' @export
probit_lc <- function(
  x,
  ...,
  p = c(1, 50),
  ci_level = 0.95,
  alpha = 0.05,
  heterogeneity = TRUE
) {
  chk::chk_numeric(p)
  chk::chk_range(p, c(0, 100))
  chk::chk_number(ci_level)
  chk::chk_range(ci_level, c(0, 1))
  chk::chk_flag(heterogeneity)

  x <- as_tox_data(x, ...)
  parts <- quantal_prep(x)

  # Abbott adjustment of the observed proportions, as the manual's own output
  # reports it. With no control response this leaves the data unchanged.
  adjusted <- abbott(parts$observed, p_control = parts$control_observed)
  # Fitted as proportions with the group sizes as weights rather than as a
  # two-column count matrix, because the Abbott adjustment makes the implied
  # counts non-integer and the count form then warns on every call.
  separation <- FALSE
  fit <- withCallingHandlers(
    stats::glm(
      adjusted ~ parts$x,
      family = stats::binomial(link = "probit"),
      weights = parts$n
    ),
    warning = function(w) {
      message <- conditionMessage(w)
      # Two warnings are expected here and say nothing the user can act on.
      # The EPA designs routinely include a concentration with no response and
      # one with complete response, which produces the separation warning on
      # almost every worked example; it is recorded on the result instead. And
      # the Abbott adjustment makes the implied number of responses
      # non-integer whenever there is any control response, which is exactly
      # what the adjustment is for.
      if (grepl("fitted probabilities numerically 0 or 1", message)) {
        separation <<- TRUE
        invokeRestart("muffleWarning")
      }
      if (grepl("non-integer", message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  coefficients <- stats::coef(fit)
  covariance <- stats::vcov(fit)
  df <- stats::df.residual(fit)
  chisq <- sum(stats::residuals(fit, type = "pearson")^2)

  inflate <- heterogeneity && df > 0 && chisq > stats::qchisq(1 - alpha, df)
  factor <- if (inflate) chisq / df else 1

  estimates <- do.call(
    rbind,
    lapply(p, function(pct) {
      limits <- fieller_limits(
        coefficients,
        covariance * factor,
        pct / 100,
        ci_level
      )
      data.frame(
        endpoint = paste0("LC", pct),
        p = pct,
        estimate = 10^limits$estimate,
        lower = 10^limits$lower,
        upper = 10^limits$upper,
        ci_method = "Fieller on the log10 scale",
        stringsAsFactors = FALSE
      )
    })
  )

  new_estimate(
    method = "Probit method",
    estimates = estimates,
    working = working_frame(parts),
    reference = "EPA-821-R-02-012 section 11.2.5",
    fit = fit,
    chisq = chisq,
    df = df,
    chisq_critical = if (df > 0) stats::qchisq(1 - alpha, df) else NA_real_,
    heterogeneity_applied = inflate,
    separation = separation,
    ci_level = ci_level
  )
}

#' Fieller's confidence limits for a lethal concentration
#'
#' @param coefficients The fitted intercept and slope.
#' @param covariance Their covariance matrix, already inflated if required.
#' @param p The response proportion.
#' @param ci_level Confidence level.
#' @return A list with `estimate`, `lower` and `upper`, on the log10 scale.
#' @noRd
fieller_limits <- function(coefficients, covariance, p, ci_level) {
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  slope <- coefficients[2]
  estimate <- unname((stats::qnorm(p) - coefficients[1]) / slope)

  v11 <- covariance[1, 1]
  v12 <- covariance[1, 2]
  v22 <- covariance[2, 2]
  g <- unname(z^2 * v22 / slope^2)

  # g >= 1 means the slope is not significantly different from zero, and the
  # interval is unbounded rather than merely wide.
  if (g >= 1) {
    return(list(estimate = estimate, lower = -Inf, upper = Inf))
  }

  centre <- estimate + g * (estimate + v12 / v22) / (1 - g)
  half <- unname(
    (z / (slope * (1 - g))) *
      sqrt(
        v11 +
          2 * estimate * v12 +
          estimate^2 * v22 -
          g * (v11 - v12^2 / v22)
      )
  )

  list(estimate = estimate, lower = centre - half, upper = centre + half)
}

#' The audit table of what each step did to the response
#'
#' @noRd
working_frame <- function(parts) {
  data.frame(
    conc = parts$conc,
    log10_conc = parts$x,
    n = parts$n,
    observed = parts$observed,
    smoothed = parts$smoothed,
    adjusted = parts$adjusted,
    stringsAsFactors = FALSE
  )
}

#' @export
print.tox_estimate <- function(x, ...) {
  cat(x$method, "\n", sep = "")
  cat("  ", x$reference, "\n\n", sep = "")

  out <- x$estimates
  numeric_columns <- vapply(out, is.numeric, logical(1))
  out[numeric_columns] <- lapply(out[numeric_columns], signif, digits = 5)
  print(out, row.names = FALSE)

  if (!is.null(x$chisq)) {
    cat(
      "\n  Chi-square for heterogeneity: ",
      signif(x$chisq, 4),
      " on ",
      x$df,
      " df (critical ",
      signif(x$chisq_critical, 4),
      ")\n",
      sep = ""
    )
    if (x$heterogeneity_applied) {
      cat("  Variance inflated by chi-square / df.\n")
    }
  }
  invisible(x)
}
