#' Williams' test
#'
#' Compares each concentration with the control while assuming the underlying
#' concentration-response relationship is monotone. Using that assumption makes
#' the test more powerful than Dunnett's procedure when it holds, and
#' misleading when it does not.
#'
#' @details
#' ## This is not an EPA method
#'
#' Williams' test is **not** on the EPA flowchart. Section 9.4.1.2 of the
#' chronic manual mentions it only as an alternative that "requires additional
#' assumptions", and neither manual gives a worked example or a table of
#' critical values for it. It is provided here because ToxCalc offered it, and
#' every result is labelled as an extension so it cannot be mistaken for the
#' analysis the manual prescribes.
#'
#' ## The method
#'
#' The concentration means are replaced by their isotonic estimates under the
#' restriction that the response is monotone in concentration. **The control is
#' not part of that restriction**; it is estimated freely. The statistic for
#' concentration `i` is then
#'
#' \deqn{\bar{t}_i = \frac{\bar{Y}_0 - M_i}{s\sqrt{1/n_i + 1/n_0}}}
#'
#' with `s` the square root of the within mean square from the full analysis of
#' variance, and `M_i` the isotonic estimate at concentration `i` formed from
#' concentrations 1 to `i` only.
#'
#' The test proceeds **downwards** from the highest concentration. Once a
#' concentration is not significant, every lower one is declared not
#' significant without being tested. That is what makes the test a trend
#' procedure rather than a set of pairwise comparisons, and it is recorded in
#' the `tested` column of the result.
#'
#' ## Critical values are simulated
#'
#' Williams tabulated his critical values, and those tables are reproduced from
#' *Biometrics*, so they are not transcribed here (see the package vignette).
#' They are obtained instead by simulating the null distribution of
#' \eqn{\bar{t}_i} for the design in hand: the concentration means and the
#' pooled variance are drawn under the hypothesis of no effect, the isotonic
#' estimate is formed, and the upper `alpha` quantile is taken.
#'
#' Simulation also removes the two things that go wrong with the tables. The
#' critical value depends on `i`, the position of the concentration being
#' tested, not on the total number of concentrations, which is the commonest
#' implementation error. And the tables are tabulated only at particular
#' degrees of freedom and require interpolation in `1/nu`, which simulation
#' avoids entirely.
#'
#' `seed` defaults to a fixed value so that the same design always returns the
#' same critical value. The Monte Carlo standard error is reported alongside
#' it; raise `nsim` to reduce it.
#'
#' ## Validation
#'
#' With no EPA worked example and no retrievable table, this implementation is
#' checked against mathematical identities rather than published output. With a
#' single concentration there is no order restriction, so the critical value
#' must equal `qt(1 - alpha, nu)`, and the simulation reproduces it. Beyond
#' one concentration the critical value must fall below Dunnett's for the same
#' design, because the order restriction is additional information, and it
#' must increase with `i`. All three are asserted in the test suite.
#'
#' @inheritParams toxstats_params
#' @param nsim Number of simulations used for each critical value.
#' @param seed Integer seed for the simulation. The random number stream is
#'   restored afterwards, so calling this function does not disturb any other
#'   analysis.
#'
#' @inherit dunnett return
#'
#' @references
#' Williams DA (1971) A test for differences between treatment means when
#' several dose levels are compared with a zero dose control. *Biometrics*
#' 27:103-117.
#'
#' Williams DA (1972) The comparison of several dose levels with a zero dose
#' control. *Biometrics* 28:519-531.
#'
#' @examples
#' # Not an EPA method; compare with dunnett() on the same data.
#' williams(fathead_c1, response = "weight", nsim = 2000)
#'
#' @export
williams <- function(x, ..., alpha = 0.05, nsim = 20000, seed = 1L) {
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))
  chk::chk_whole_number(nsim)
  chk::chk_gte(nsim, 1000)
  chk::chk_whole_number(seed)

  x <- as_tox_data(x, ...)
  parts <- comparison_parts(x)
  control <- parts$control_index
  treatments <- setdiff(seq_along(parts$conc), control)
  k <- length(treatments)
  if (!k) {
    chk::abort_chk("At least one non-control concentration is needed.")
  }

  # Work on a scale where a toxic effect lowers the response, so the order
  # restriction is always "non-increasing" whatever the endpoint direction.
  sign <- effect_sign(x$direction)
  oriented <- sign * parts$mean
  control_mean <- oriented[control]
  dose_means <- oriented[treatments]
  dose_n <- parts$n[treatments]

  monotone <- is_monotone(dose_means, direction = "decreasing")
  if (!monotone$monotone) {
    warning(
      "Williams' test assumes a monotone concentration-response relationship. ",
      "The observed means depart from it at position(s) ",
      paste(monotone$violations, collapse = ", "),
      ". The isotonic step will absorb the departure, but the assumption ",
      "should be considered before the result is used.",
      call. = FALSE
    )
  }

  seed_state <- preserve_seed(seed)
  on.exit(seed_state(), add = TRUE)

  # One simulation serves every position. Deriving all k critical values from
  # the same draws is not only k times faster, it makes them coherent: drawn
  # independently they can come out very slightly non-monotone in i, which is
  # Monte Carlo noise but looks like an error.
  simulated <- williams_null(
    n_control = parts$n[control],
    n_doses = dose_n,
    df = parts$df,
    nsim = nsim
  )

  isotonic <- vapply(
    seq_len(k),
    function(i) isotonic_last(dose_means[seq_len(i)], dose_n[seq_len(i)]),
    numeric(1)
  )
  statistic <- (control_mean - isotonic) /
    (parts$sw * sqrt(1 / dose_n + 1 / parts$n[control]))
  critical <- vapply(
    seq_len(k),
    function(i) unname(stats::quantile(simulated[, i], 1 - alpha)),
    numeric(1)
  )
  monte_carlo_se <- vapply(
    seq_len(k),
    function(i) monte_carlo_error(simulated[, i], alpha),
    numeric(1)
  )
  # The simulated null distribution gives a p-value directly, which is more
  # informative than reporting only the side of the critical value.
  p_value <- vapply(
    seq_len(k),
    function(i) mean(simulated[, i] >= statistic[i]),
    numeric(1)
  )

  # Step down from the highest concentration. The first that is not
  # significant stops the procedure, and everything below it is declared not
  # significant without being tested.
  significant <- logical(k)
  tested <- logical(k)
  for (i in rev(seq_len(k))) {
    tested[i] <- TRUE
    if (statistic[i] > critical[i]) {
      significant[i] <- TRUE
    } else {
      break
    }
  }

  comparisons <- data.frame(
    conc = parts$conc[treatments],
    n = dose_n,
    mean = parts$mean[treatments],
    isotonic = sign * isotonic,
    statistic = statistic,
    critical = critical,
    mc_se = monte_carlo_se,
    p_value = p_value,
    tested = tested,
    significant = significant,
    stringsAsFactors = FALSE
  )
  rownames(comparisons) <- NULL

  new_comparison(
    test = "williams",
    method = paste0(
      "Williams' test (not an EPA method; ",
      format(nsim, big.mark = ","),
      " simulations)"
    ),
    comparisons = comparisons,
    reference = paste0(
      "Williams (1971) Biometrics 27:103-117; ",
      "not on the EPA-821-R-02-013 Figure 2 flowchart"
    ),
    alpha = alpha,
    sw = parts$sw,
    df = parts$df,
    nsim = nsim,
    seed = seed,
    monotone_observed = monotone$monotone
  )
}

#' The isotonic estimate at the highest of a run of concentrations
#'
#' For a non-increasing isotonic fit the estimate at the last position is the
#' smallest weighted average of any trailing run, which is the max-min formula
#' reduced to its final index.
#'
#' @param y Oriented concentration means, lowest concentration first.
#' @param w Weights, the number of replicates at each concentration.
#' @return A single number.
#' @noRd
isotonic_last <- function(y, w) {
  numerator <- rev(cumsum(rev(w * y)))
  denominator <- rev(cumsum(rev(w)))
  min(numerator / denominator)
}

#' Simulate the null distribution of Williams' statistic
#'
#' Under the hypothesis of no effect the concentration means and the pooled
#' variance are independent, so there is no need to simulate individual
#' observations: the means are drawn directly and the variance from its
#' chi-squared distribution.
#'
#' @param n_control Replicates in the control.
#' @param n_doses Replicates at concentrations 1 to i.
#' @param df Residual degrees of freedom from the full analysis of variance.
#' @param nsim Number of simulations.
#' @return A matrix with `nsim` rows and one column per position `i`, holding
#'   the simulated statistic for the test at concentration `i`.
#' @noRd
williams_null <- function(n_control, n_doses, df, nsim) {
  k <- length(n_doses)
  scale <- sqrt(stats::rchisq(nsim, df) / df)
  control <- stats::rnorm(nsim, sd = 1 / sqrt(n_control))

  means <- vapply(
    n_doses,
    function(n) stats::rnorm(nsim, sd = 1 / sqrt(n)),
    numeric(nsim)
  )
  means <- matrix(means, nrow = nsim)

  out <- matrix(NA_real_, nrow = nsim, ncol = k)
  for (i in seq_len(k)) {
    # The isotonic estimate at position i, computed for every simulation at
    # once: accumulate trailing weighted sums over concentrations 1 to i and
    # keep the running minimum.
    isotonic <- rep(Inf, nsim)
    numerator <- numeric(nsim)
    denominator <- 0
    for (u in rev(seq_len(i))) {
      numerator <- numerator + n_doses[u] * means[, u]
      denominator <- denominator + n_doses[u]
      isotonic <- pmin(isotonic, numerator / denominator)
    }
    out[, i] <- (control - isotonic) /
      (scale * sqrt(1 / n_doses[i] + 1 / n_control))
  }
  out
}

#' Standard error of a simulated quantile
#'
#' @param values The simulated statistics.
#' @param alpha The upper tail probability.
#' @return A single number.
#' @noRd
monte_carlo_error <- function(values, alpha) {
  n <- length(values)
  quantile <- stats::quantile(values, 1 - alpha)
  density <- stats::density(values, from = quantile, to = quantile, n = 1)$y
  if (!is.finite(density) || density <= 0) {
    return(NA_real_)
  }
  sqrt(alpha * (1 - alpha) / n) / density
}

#' Set a seed and return a function that restores the previous stream
#'
#' Simulating a critical value must not disturb a user's own random number
#' stream, which it would otherwise do silently.
#'
#' @param seed The seed to set.
#' @return A function that restores the stream when called.
#' @noRd
preserve_seed <- function(seed) {
  existed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  previous <- if (existed) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  set.seed(seed)

  function() {
    if (existed) {
      assign(".Random.seed", previous, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    invisible(NULL)
  }
}
