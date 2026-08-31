# Tests for the assumption tests in R/assumptions.R
#
# Reference values come from the worked examples in EPA-821-R-02-013.

# epa_normality: Shapiro-Wilk ------------------------------------------------

test_that("epa_normality reproduces the Appendix B Shapiro-Wilk example", {
  fit <- epa_normality(fathead_b1, response = "weight")

  # The manual prints W = 0.959. Royston's algorithm gives 0.9601 from the
  # unrounded data, and 0.9594 once the manual's rounding of the concentration
  # means is applied, so agreement is to three decimal places.
  expect_equal(round(fit$statistic, 3), 0.960)
  expect_equal(fit$statistic_name, "W")
  expect_true(fit$normal)
  expect_equal(fit$n, 20L)
  expect_match(fit$method, "Shapiro-Wilk")
})

test_that("epa_normality tests the pooled within-group centred residuals", {
  # Not the raw values: the raw weights span five concentration means and are
  # not expected to be normal even when the residuals are.
  x <- tox_data(fathead_b1, response = "weight")
  residuals <- fathead_b1$weight -
    ave(fathead_b1$weight, fathead_b1$conc, FUN = mean)

  expect_equal(
    epa_normality(x)$statistic,
    unname(shapiro.test(residuals)$statistic)
  )
  expect_false(isTRUE(all.equal(
    epa_normality(x)$statistic,
    unname(shapiro.test(fathead_b1$weight)$statistic)
  )))
})

test_that("epa_normality decides against the assumption alpha", {
  x <- tox_data(fathead_b1, response = "weight")
  # p is about 0.55 here, so the conclusion is stable across sensible alphas.
  expect_true(epa_normality(x, alpha_assumption = 0.01)$normal)
  expect_true(epa_normality(x, alpha_assumption = 0.05)$normal)
  expect_equal(epa_normality(x, alpha_assumption = 0.05)$alpha, 0.05)
})

# epa_normality: Kolmogorov --------------------------------------------------

test_that("epa_normality reproduces the Appendix B Kolmogorov example", {
  # 60 observations, so "auto" must select the Kolmogorov statistic.
  fit <- epa_normality(ceriodaphnia_b7, response = "young")

  expect_match(fit$method, "Kolmogorov")
  expect_equal(fit$statistic_name, "D*")
  # The manual prints D* = 0.4684 and a critical value of 1.035.
  # The correct value is 0.4572. The manual prints 0.4684 because it rounds
  # each z to two decimal places to look the probability up in a printed
  # normal table; see the test below.
  expect_equal(round(fit$statistic, 4), 0.4572)
  expect_equal(fit$critical, 1.035)
  expect_true(fit$normal)
  expect_equal(fit$n, 60L)
})

test_that("the Kolmogorov critical values match Table B.11", {
  x <- tox_data(ceriodaphnia_b7, response = "young")
  expected <- c(
    "0.01" = 1.035,
    "0.025" = 0.955,
    "0.05" = 0.895,
    "0.1" = 0.819,
    "0.15" = 0.775
  )
  for (a in names(expected)) {
    fit <- epa_normality(
      x,
      method = "kolmogorov",
      alpha_assumption = as.numeric(a)
    )
    expect_equal(fit$critical, unname(expected[a]))
  }
})

test_that("the Kolmogorov test rejects an untabulated alpha", {
  # Table B.11 gives five alpha levels and no interpolation rule, so an
  # untabulated level cannot be honoured.
  expect_error(
    epa_normality(
      ceriodaphnia_b7,
      response = "young",
      method = "kolmogorov",
      alpha_assumption = 0.02
    ),
    regexp = "tabulated only at alpha"
  )
})

test_that("method selection follows the 50 observation rule", {
  small <- tox_data(fathead_b1, response = "weight")
  large <- tox_data(ceriodaphnia_b7, response = "young")

  expect_match(epa_normality(small)$method, "Shapiro-Wilk")
  expect_match(epa_normality(large)$method, "Kolmogorov")
  # An explicit method overrides the rule.
  expect_match(
    epa_normality(large, method = "shapiro_wilk")$method,
    "Shapiro-Wilk"
  )
})

# epa_variance ---------------------------------------------------------------

test_that("epa_variance reproduces Bartlett on the Appendix B example", {
  fit <- epa_variance(fathead_b1, response = "weight")

  # The manual prints B = 7.691, computed from its own table of variances
  # rounded to four decimal places, one of which is itself misrounded. From
  # the unrounded data the statistic is 6.836. Both are far below the 13.277
  # critical value, so the conclusion is unchanged.
  expect_equal(round(fit$statistic, 3), 6.836)
  expect_equal(fit$df, 4L)
  expect_equal(round(fit$critical, 3), 13.277)
  expect_true(fit$homogeneous)
  expect_match(fit$method, "Bartlett")
})

test_that("the manual's Bartlett value follows from its rounded variances", {
  # Documents why the package value differs, so the discrepancy cannot be
  # mistaken for an implementation error later.
  printed <- c(0.0018, 0.0020, 0.0001, 0.0059, 0.0037)
  nu <- rep(3, 5)
  s2bar <- sum(nu * printed) / sum(nu)
  # The manual also rounds C to 1.133 before dividing; using the unrounded
  # 1.133333 gives 7.689 rather than the printed 7.691.
  cc <- 1.133
  b <- (sum(nu) * log(s2bar) - sum(nu * log(printed))) / cc

  expect_equal(round(b, 3), 7.691)
})

test_that("epa_variance offers robust alternatives, flagged as non-EPA", {
  x <- tox_data(fathead_b1, response = "weight")

  levene <- epa_variance(x, method = "levene")
  fligner <- epa_variance(x, method = "fligner")

  expect_match(levene$method, "not an EPA method")
  expect_match(fligner$method, "not an EPA method")
  expect_true(levene$homogeneous)
  expect_true(fligner$homogeneous)
  expect_equal(levene$statistic_name, "F")
})

test_that("epa_variance rejects designs it cannot test", {
  single <- data.frame(conc = c(0, 0, 1, 1), y = c(1, 2, 3, 4))
  expect_error(
    epa_variance(single[single$conc == 0, ], response = "y"),
    regexp = "two concentrations"
  )

  one_rep <- data.frame(conc = c(0, 1), y = c(1, 2))
  expect_error(
    epa_variance(one_rep, response = "y"),
    regexp = "at least two replicates"
  )
})

# printing -------------------------------------------------------------------

test_that("assumption results print their reference", {
  expect_output(
    print(epa_normality(fathead_b1, response = "weight")),
    "Appendix B"
  )
  expect_output(
    print(epa_variance(fathead_b1, response = "weight")),
    "Appendix B"
  )
})

# reconciling with the manual's printed statistics ---------------------------

test_that("the manual's Kolmogorov value follows from rounding z to 2 dp", {
  # The manual looks each standardised residual up in a printed normal table,
  # which means rounding z to two decimal places first. Doing the same
  # reproduces its D = 0.0597 and D* = 0.4684; working at full precision gives
  # 0.0583 and 0.4572. Both are far below the 1.035 critical value.
  x <- tox_data(ceriodaphnia_b7, response = "young")
  residuals <- ceriodaphnia_b7$young -
    ave(ceriodaphnia_b7$young, ceriodaphnia_b7$conc, FUN = function(z) {
      round(mean(z), 1)
    })

  n <- length(residuals)
  z <- round(sort(residuals) / sd(residuals), 2)
  p <- pnorm(z)
  d <- max(max(seq_len(n) / n - p), max(p - (seq_len(n) - 1) / n))

  expect_equal(round(d, 4), 0.0597)
  expect_equal(round(d * (sqrt(n) - 0.01 + 0.85 / sqrt(n)), 3), 0.468)
})
