# Tests for the point-estimation methods in R/point_estimation.R
#
# Reference values come from Table 20 and sections 11.2.2 to 11.2.5 of
# EPA-821-R-02-012, and from Appendix K of EPA-821-R-02-013.

# helpers -------------------------------------------------------------------

table20 <- function(column) {
  tox_data(
    acute_table20,
    response = column,
    n_exposed = "exposed",
    type = "quantal"
  )
}

# graphical -----------------------------------------------------------------

test_that("graphical_lc50 reproduces the section 11.2.2 example", {
  fit <- graphical_lc50(table20("graphical"))

  # The manual smooths to 0.0125 across the control and the lower three
  # concentrations, adjusts to 0, 0, 0, 1, 1, and reads 35 per cent off the
  # plot. Interpolating rather than reading gives 35.36.
  expect_equal(round(fit$working$smoothed, 4), c(0.0125, 0.0125, 0.0125, 1, 1))
  expect_equal(fit$working$adjusted, c(0, 0, 0, 1, 1))
  expect_equal(round(fit$estimates$estimate, 3), 35.355)
  expect_true(is.na(fit$estimates$lower))
})

test_that("the graphical interpolation is linear in log10 concentration", {
  # 25 per cent and 50 per cent bracket the estimate, and the geometric mean
  # of the two is what a semi-logarithmic plot gives.
  fit <- graphical_lc50(table20("graphical"))
  expect_equal(fit$estimates$estimate, sqrt(25 * 50))
})

test_that("graphical_lc50 requires the response to bracket 0.5", {
  low <- acute_table20
  low$graphical <- c(0, 0, 0, 0, 1, 2)
  expect_error(
    graphical_lc50(
      low,
      response = "graphical",
      n_exposed = "exposed",
      type = "quantal"
    ),
    regexp = "must bracket 0.5"
  )
})

# spearman-karber -----------------------------------------------------------

test_that("spearman_karber reproduces the section 11.2.3 example", {
  fit <- spearman_karber(table20("spearman_karber"))

  # The manual prints smoothed 0.025 across the control and lower three
  # concentrations, adjusted 0, 0, 0, 0.641, 1.
  expect_equal(round(fit$working$smoothed, 4), c(0.025, 0.025, 0.025, 0.65, 1))
  expect_equal(round(fit$working$adjusted, 3), c(0, 0, 0, 0.641, 1))

  # m = 1.656527, V(m) = 0.0010977, LC50 45.3 with limits 38.9 and 52.8.
  expect_equal(round(fit$m, 5), 1.65652)
  expect_equal(signif(fit$variance, 4), 0.001098)
  expect_equal(round(fit$estimates$estimate, 1), 45.3)
  expect_equal(round(fit$estimates$lower, 1), 38.9)
  expect_equal(round(fit$estimates$upper, 1), 52.8)
})

test_that("the Spearman-Karber interval uses the manual's multiplier of 2", {
  # Section 11.2.3.3 step 6 writes m +/- 2.0 sqrt(V(m)). That is the manual's
  # rounding of 1.96, and keeping it is what reproduces the printed limits.
  fit <- spearman_karber(table20("spearman_karber"))
  expect_equal(
    log10(fit$estimates$upper) - fit$m,
    2 * sqrt(fit$variance)
  )
  expect_equal(toxstats:::epa_multiplier(0.95), 2)
  expect_equal(toxstats:::epa_multiplier(0.99), qnorm(0.995))
})

test_that("spearman_karber refuses data that do not meet its requirements", {
  # The trimmed column does not reach a complete response, which is exactly
  # why the manual analyses it with the trimmed method.
  expect_error(
    spearman_karber(table20("trimmed")),
    regexp = "Use `trimmed_spearman_karber\\(\\)` instead"
  )
})

# trimmed spearman-karber ---------------------------------------------------

test_that("the automatic trim reproduces both published trims", {
  # Acute Table 20 prints 20.51 per cent.
  acute <- trimmed_spearman_karber(table20("trimmed"))
  expect_equal(round(100 * acute$trim, 2), 20.51)
  expect_true(acute$automatic_trim)

  # Chronic Appendix K prints 20.41 per cent: 40 organisms per concentration
  # with 2, 0, 2, 0, 0 and 32 responding.
  appendix_k <- data.frame(
    conc = c(0, 6.25, 12.5, 25, 50, 100),
    dead = c(2, 0, 2, 0, 0, 32),
    exposed = 40
  )
  chronic <- trimmed_spearman_karber(
    appendix_k,
    response = "dead",
    n_exposed = "exposed",
    type = "quantal"
  )
  expect_equal(round(100 * chronic$trim, 2), 20.41)
})

test_that("the trim is computed after smoothing and Abbott adjustment", {
  # Order matters. Taking the trim from the raw proportions gives a different
  # answer, and it is the commonest way to get this method wrong.
  fit <- trimmed_spearman_karber(table20("trimmed"))

  expect_equal(
    fit$trim,
    max(fit$working$adjusted[1], 1 - fit$working$adjusted[5])
  )
  raw <- acute_table20$trimmed[-1] / 20
  expect_false(isTRUE(all.equal(fit$trim, max(raw[1], 1 - raw[5]))))
})

test_that("trimmed_spearman_karber reproduces the printed LC50", {
  fit <- trimmed_spearman_karber(table20("trimmed"))
  # The manual's program prints 77.11.
  expect_equal(round(fit$estimates$estimate, 2), 77.11)
})

test_that("the trimmed interval is close to, but not identical with, the manual", {
  # Documented departure: the manual delegates the interval to a program whose
  # variance formula it does not state, and Hamilton's corrected expression is
  # not retrievable. The delta method is used instead. The manual prints
  # 69.74 to 85.26.
  fit <- trimmed_spearman_karber(table20("trimmed"))

  expect_equal(fit$estimates$ci_method, "delta method on the log10 scale")
  expect_equal(round(fit$estimates$lower), 70)
  expect_equal(round(fit$estimates$upper), 85)
  expect_lt(abs(fit$estimates$lower - 69.74), 0.5)
  expect_lt(abs(fit$estimates$upper - 85.26), 0.5)
})

test_that("an explicit trim overrides the automatic one", {
  fit <- trimmed_spearman_karber(table20("trimmed"), trim = 0.25)
  expect_equal(fit$trim, 0.25)
  expect_false(fit$automatic_trim)
  expect_match(fit$method, "trim 25 per cent")
})

test_that("a trim of one half or more is refused", {
  expect_error(
    trimmed_spearman_karber(table20("trimmed"), trim = 0.5),
    regexp = "must be below 0.5"
  )
})

test_that("the trimmed method agrees with the untrimmed one at zero trim", {
  # With a full response range the automatic trim is zero, and the estimator
  # must then reduce to plain Spearman-Karber.
  sk <- spearman_karber(table20("spearman_karber"))
  tsk <- trimmed_spearman_karber(table20("spearman_karber"))

  expect_equal(tsk$trim, 0)
  expect_equal(tsk$m, sk$m)
})

# probit --------------------------------------------------------------------

test_that("probit_lc reproduces the section 11.2.5 example", {
  fit <- probit_lc(table20("probit"), p = c(1, 50))

  # The manual prints LC1 7.924 (4.147, 10.959) and LC50 22.872
  # (18.787, 27.846).
  expect_equal(round(fit$estimates$estimate, 3), c(7.924, 22.872))
  expect_equal(round(fit$estimates$lower, 3), c(4.147, 18.787))
  expect_equal(round(fit$estimates$upper, 3), c(10.959, 27.846))
})

test_that("the heterogeneity statistic is Pearson's, not the deviance", {
  # The manual prints 3.076 against a tabular 7.815. Pearson gives 3.076; the
  # deviance for the same fit is 3.859, so the two are distinguishable and the
  # manual is using Pearson.
  fit <- probit_lc(table20("probit"))

  expect_equal(round(fit$chisq, 3), 3.076)
  expect_equal(fit$df, 3L)
  expect_equal(round(fit$chisq_critical, 3), 7.815)
  expect_equal(round(deviance(fit$fit), 3), 3.859)
  expect_false(fit$heterogeneity_applied)
})

test_that("the probit limits are Fieller's, not the delta method", {
  # The delta method gives 19.04 to 27.47 for the LC50, which does not match
  # the manual. This test pins the choice.
  fit <- probit_lc(table20("probit"), p = 50)
  expect_equal(fit$estimates$ci_method, "Fieller on the log10 scale")
  expect_equal(round(fit$estimates$lower, 3), 18.787)
})

test_that("the variance is inflated only when heterogeneity is significant", {
  adequate <- probit_lc(table20("probit"))
  expect_false(adequate$heterogeneity_applied)

  # Turning the check off must not change an adequate fit.
  off <- probit_lc(table20("probit"), heterogeneity = FALSE)
  expect_equal(adequate$estimates$lower, off$estimates$lower)
})

test_that("probit_lc suppresses only the two expected glm warnings", {
  # Every EPA design has a concentration with no response and one with a
  # complete response, so the separation warning would fire on almost every
  # call. It is recorded rather than printed.
  expect_no_warning(probit_lc(table20("probit")))
  expect_type(probit_lc(table20("probit"))$separation, "logical")
})

# shared behaviour ----------------------------------------------------------

test_that("the point-estimation methods need quantal data", {
  expect_error(
    graphical_lc50(fathead_c1, response = "weight"),
    regexp = "needs quantal data"
  )
})

test_that("the working table records what each step did", {
  fit <- spearman_karber(table20("spearman_karber"))
  expect_equal(
    names(fit$working),
    c("conc", "log10_conc", "n", "observed", "smoothed", "adjusted")
  )
  expect_equal(fit$working$conc, c(6.25, 12.5, 25, 50, 100))
  expect_equal(fit$working$log10_conc, log10(c(6.25, 12.5, 25, 50, 100)))
})

test_that("the estimates print with their reference", {
  expect_output(print(spearman_karber(table20("spearman_karber"))), "11.2.3")
  expect_output(print(probit_lc(table20("probit"))), "Chi-square")
})
