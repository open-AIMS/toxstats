#' Estimate an inhibition concentration by linear interpolation
#'
#' Estimates the concentration causing a given percentage reduction in a
#' sublethal response such as growth or reproduction, by interpolating between
#' the two concentrations whose smoothed mean responses bracket it. This is the
#' EPA linear interpolation method, implemented in the manual by the ICPIN
#' program.
#'
#' @details
#' ## The estimate
#'
#' The concentration means are smoothed to be monotone non-increasing, then
#'
#' \deqn{IC_p = C_J + [M_1 (1 - p/100) - M_J]\frac{C_{J+1} - C_J}{M_{J+1} -
#'   M_J}}
#'
#' where \eqn{C_J} and \eqn{C_{J+1}} bracket the target response and \eqn{M_1}
#' is the smoothed control mean. On the Appendix M worked example this gives
#' 8.5716 for the IC25 against the 8.5715 the manual's program prints, and
#' 10.893 for the IC50 against 10.89.
#'
#' Smoothing is by pool-adjacent-violators with **equal weight on each
#' concentration mean**, not weighted by the number of replicates. That is what
#' reproduces the manual's smoothed means of 28.75 across the control and the
#' three lowest concentrations. With unequal replication the manual does not say
#' what it does, and this choice is a documented assumption.
#'
#' Where the target response falls outside the range of the smoothed means the
#' estimate is reported as an inequality: greater than the highest
#' concentration, or less than the lowest.
#'
#' ## The interval
#'
#' Ordinary interval methods do not apply to an interpolated estimate, so the
#' manual uses a bootstrap. Replicate values are resampled **with replacement
#' within each concentration**, the control included, and the means are
#' recomputed, **re-smoothed** and re-interpolated on every iteration.
#'
#' The re-smoothing is the detail most easily missed. The manual's step list
#' places smoothing before resampling, which reads as though it happens once.
#' Omitting it from the loop does not merely widen the interval; it biases the
#' estimate badly. On the Appendix M data the resampled estimates then average
#' about 10.4 against a point estimate of 8.57, with a standard deviation of
#' 0.49 rather than 0.14, because unsmoothed resampled means are not monotone
#' and the interpolation selects the wrong bracketing pair.
#'
#' Limits are the empirical order statistics the manual describes, "the second
#' smallest and second largest" of 80 resamples, generalised to
#' `floor((nboot + 1) * (1 - ci_level) / 2)` in from each end.
#'
#' **The interval is not exactly reproducible from the manual.** Its printed
#' limits of 8.3112 and 9.0418 come from 80 resamples drawn with the seed
#' -641671986 by a Turbo Pascal generator, so no other program can return them.
#' The estimate, the smoothed means and the bootstrap standard deviation are
#' reproducible, and are what the test suite checks.
#'
#' The default of 200 resamples departs from the manual's 80. The manual itself
#' warns that "confidence limits based on the empirical quantiles of a
#' bootstrap distribution of 80 samples may be unstable", and permits up to
#' 1000.
#'
#' ## Expanded limits
#'
#' The manual reports a second, wider interval when there are fewer than seven
#' replicates per concentration. Its definition appears only in the ICPIN
#' version 2.0 program documentation, which is not publicly retrievable, so it
#' is **not** implemented. A warning is raised when the design would have
#' triggered it.
#'
#' @inheritParams toxcalc_params
#' @param p The percentage reductions to estimate, on a 1 to 99 scale.
#'
#' @return An object of class `toxcalc_estimate`, with an additional `boot`
#'   element holding the resampled estimates, the count that could not be
#'   computed, and the seed used.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix M.
#'
#' Norberg-King TJ (1993) *A linear interpolation method for sublethal
#' toxicity: the inhibition concentration (ICp) approach*, version 2.0.
#' EPA/600/M-91/037.
#'
#' Efron B (1982) *The Jackknife, the Bootstrap, and Other Resampling Plans*.
#' Society for Industrial and Applied Mathematics, Philadelphia.
#'
#' @examples
#' # Appendix M worked example; the manual gives IC25 = 8.5715 and
#' # IC50 = 10.89.
#' icp(ceriodaphnia_m1, response = "young", p = c(25, 50), seed = 42)
#'
#' @export
icp <- function(
  x,
  ...,
  p = 25,
  nboot = 200,
  seed = NULL,
  ci_level = 0.95
) {
  chk::chk_numeric(p)
  chk::chk_range(p, c(1, 99))
  chk::chk_whole_number(nboot)
  chk::chk_gte(nboot, 20)
  chk::chk_number(ci_level)
  chk::chk_range(ci_level, c(0, 1))
  chk::chk_null_or(seed, vld = chk::vld_whole_number)

  x <- as_toxcalc_data(x, ...)
  response <- analysis_response(x)
  group <- factor(x$replicates$conc)
  conc <- as.numeric(levels(group))

  if (which(conc == x$control) != 1L) {
    chk::abort_chk(
      "The control must be the lowest concentration, because the smoothing ",
      "step orders the response by concentration."
    )
  }
  if (nrow(x$pooled) < 3) {
    chk::abort_chk(
      "At least three concentrations including the control are needed to ",
      "interpolate."
    )
  }
  if (any(x$pooled$n_rep < 7)) {
    warning(
      "At least one concentration has fewer than seven replicates. The ",
      "manual reports an additional expanded interval in that case, whose ",
      "definition appears only in the ICPIN program documentation and is not ",
      "implemented here. Only the ordinary interval is reported.",
      call. = FALSE
    )
  }

  values <- split(response, group)
  means <- vapply(values, mean, numeric(1))
  smoothed <- smooth_monotone(means, direction = "decreasing")

  if (!is.null(seed)) {
    set.seed(seed)
  }
  resampled <- vapply(
    seq_len(nboot),
    function(i) {
      drawn <- lapply(values, function(v) sample(v, length(v), replace = TRUE))
      resmoothed <- smooth_monotone(
        vapply(drawn, mean, numeric(1)),
        direction = "decreasing"
      )
      vapply(
        p,
        function(pct) interpolate_icp(conc, resmoothed, pct)$estimate,
        numeric(1)
      )
    },
    numeric(length(p))
  )
  resampled <- matrix(resampled, nrow = length(p))

  index <- max(1L, floor((nboot + 1) * (1 - ci_level) / 2))
  estimates <- do.call(
    rbind,
    lapply(seq_along(p), function(i) {
      point <- interpolate_icp(conc, smoothed, p[i])
      draws <- sort(resampled[i, ][is.finite(resampled[i, ])])
      enough <- length(draws) >= 2 * index + 1
      data.frame(
        endpoint = paste0("IC", p[i]),
        p = p[i],
        estimate = point$estimate,
        lower = if (enough) draws[index] else NA_real_,
        upper = if (enough) draws[length(draws) - index + 1L] else NA_real_,
        ci_method = "bootstrap order statistics",
        bound = point$bound,
        boot_mean = if (length(draws)) mean(draws) else NA_real_,
        boot_sd = if (length(draws) > 1) stats::sd(draws) else NA_real_,
        n_undefined = nboot - length(draws),
        stringsAsFactors = FALSE
      )
    })
  )

  undefined <- max(estimates$n_undefined) / nboot
  if (undefined > 0.05) {
    warning(
      round(100 * undefined),
      " per cent of bootstrap resamples did not bracket the target response, ",
      "so the interval rests on fewer draws than requested.",
      call. = FALSE
    )
  }

  new_estimate(
    method = "Linear interpolation method",
    estimates = estimates,
    working = data.frame(
      conc = conc,
      n_rep = as.integer(table(group)),
      mean = as.numeric(means),
      sd = vapply(values, stats::sd, numeric(1)),
      smoothed = as.numeric(smoothed),
      stringsAsFactors = FALSE
    ),
    reference = "EPA-821-R-02-013 Appendix M",
    boot = list(replicates = resampled, nboot = nboot, seed = seed),
    ci_level = ci_level
  )
}

#' Interpolate one inhibition concentration from smoothed means
#'
#' @param conc Concentrations in increasing order, control first.
#' @param smoothed Smoothed, monotone non-increasing means.
#' @param p The percentage reduction.
#' @return A list with `estimate` and `bound`, the latter naming the side the
#'   estimate falls outside when it cannot be interpolated.
#' @noRd
interpolate_icp <- function(conc, smoothed, p) {
  target <- smoothed[1] * (1 - p / 100)

  # Section 4.3: if the response never falls to the target the estimate lies
  # above the highest concentration tested, and if the first concentration
  # already responds beyond it the estimate lies below the lowest.
  if (all(smoothed > target)) {
    return(list(estimate = Inf, bound = "above the highest concentration"))
  }
  if (smoothed[2] <= target && length(smoothed) > 1) {
    below <- 1L
  } else {
    below <- max(which(smoothed > target))
  }
  if (below >= length(smoothed)) {
    return(list(estimate = Inf, bound = "above the highest concentration"))
  }

  slope <- smoothed[below + 1L] - smoothed[below]
  if (slope == 0) {
    return(list(estimate = NA_real_, bound = "indeterminate"))
  }

  list(
    estimate = conc[below] +
      (target - smoothed[below]) * (conc[below + 1L] - conc[below]) / slope,
    bound = NA_character_
  )
}
