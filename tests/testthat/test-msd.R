# Tests for the minimum significant difference in R/msd.R
#
# Reference values come from Appendix C of EPA-821-R-02-013.

# msd ------------------------------------------------------------------------

test_that("msd reproduces the Appendix C worked example", {
  fit <- msd(fathead_c1, response = "weight")

  # The manual prints Sw = 0.097 on 15 degrees of freedom, a Dunnett critical
  # value of 2.36, and MSD = 0.162.
  expect_equal(round(fit$sw, 3), 0.097)
  expect_equal(fit$df, 15L)
  expect_equal(round(fit$critical, 2), 2.36)
  expect_true(all(round(fit$msd, 3) == 0.162))
  expect_true(fit$balanced)
})

test_that("the computed Dunnett critical value matches the EPA table", {
  # EPA Table C.5 gives 2.36 for a one-sided test at alpha = 0.05 with 15
  # degrees of freedom and four concentrations excluding the control. The
  # exact multivariate t value is 2.3561, which is the same to the two decimal
  # places the table prints.
  fit <- msd(fathead_c1, response = "weight")
  expect_equal(round(fit$critical, 4), 2.3561)
})

test_that("msd names its values by concentration", {
  fit <- msd(fathead_c1, response = "weight")
  expect_equal(names(fit$msd), c("32", "64", "128", "256"))
  expect_length(fit$msd, 4)
})

test_that("msd handles an unbalanced design per concentration", {
  unbalanced <- fathead_c1[-1, ]
  fit <- msd(unbalanced, response = "weight")

  expect_false(fit$balanced)
  expect_equal(fit$n_control, 3L)
  # The control has fewer replicates, so every difference is harder to detect
  # than in the balanced case.
  expect_true(all(fit$msd > msd(fathead_c1, response = "weight")$msd))
})

test_that("the Bonferroni and Sidak critical values bracket sensibly", {
  x <- toxcalc_data(fathead_c1, response = "weight")
  dunnett <- msd(x, test = "dunnett")$critical
  bonferroni <- msd(x, test = "bonferroni")$critical
  sidak <- msd(x, test = "sidak")$critical

  # Dunnett accounts for the correlation between comparisons sharing a
  # control, so it is the least conservative of the three. Sidak is exact for
  # independent comparisons and so sits just inside Bonferroni.
  expect_lt(dunnett, sidak)
  expect_lt(sidak, bonferroni)
})

test_that("a single treatment concentration falls back to the t quantile", {
  two_level <- fathead_c1[fathead_c1$conc %in% c(0, 32), ]
  fit <- msd(two_level, response = "weight")
  expect_equal(fit$critical, qt(0.95, fit$df))
})

test_that("msd validates its arguments", {
  expect_error(
    msd(fathead_c1, response = "weight", alpha = 1.5),
    regexp = "between 0 and 1"
  )
  expect_error(
    msd(fathead_c1, response = "weight", test = "tukey"),
    regexp = "should be one of"
  )
})

# pmsd -----------------------------------------------------------------------

test_that("pmsd reproduces the Appendix C percentage", {
  fit <- pmsd(fathead_c1, response = "weight")

  # The manual describes the MSD as about 24 per cent of the control mean.
  expect_equal(round(fit$control_mean, 3), 0.677)
  expect_true(all(round(fit$pmsd, 1) == 23.9))
})

test_that("pmsd compares against the EPA Table 6 bounds", {
  fit <- pmsd(fathead_c1, response = "weight", bounds = "fathead_growth")

  expect_equal(fit$bounds, c(12, 30))
  expect_true(all(fit$status == "within"))
})

test_that("pmsd flags values outside the bounds in both directions", {
  x <- toxcalc_data(fathead_c1, response = "weight")

  expect_true(all(pmsd(x, bounds = c(30, 40))$status == "below_lower"))
  expect_true(all(pmsd(x, bounds = c(5, 10))$status == "above_upper"))
})

test_that("pmsd reports no status when no bounds are given", {
  fit <- pmsd(fathead_c1, response = "weight")
  expect_null(fit$bounds)
  expect_true(all(is.na(fit$status)))
})

test_that("pmsd validates the bounds argument", {
  expect_error(
    pmsd(fathead_c1, response = "weight", bounds = "trout_growth"),
    regexp = "trout_growth"
  )
  expect_error(
    pmsd(fathead_c1, response = "weight", bounds = c(30, 12)),
    regexp = "increasing"
  )
  expect_error(
    pmsd(fathead_c1, response = "weight", bounds = c(1, 2, 3)),
    regexp = "length two"
  )
})

# the shipped bounds table ---------------------------------------------------

test_that("epa_pmsd_bounds matches EPA Table 6", {
  expect_equal(nrow(epa_pmsd_bounds), 3L)
  expect_equal(
    epa_pmsd_bounds$lower,
    c(12, 13, 9.1)
  )
  expect_equal(
    epa_pmsd_bounds$upper,
    c(30, 47, 29)
  )
  expect_equal(epa_pmsd_bounds$method, c("1000.0", "1002.0", "1003.0"))
})

# printing -------------------------------------------------------------------

test_that("msd and pmsd print their reference", {
  expect_output(print(msd(fathead_c1, response = "weight")), "Appendix C")
  expect_output(
    print(pmsd(fathead_c1, response = "weight", bounds = "fathead_growth")),
    "EPA bounds"
  )
})

# determinism ----------------------------------------------------------------

test_that("the Dunnett critical value is identical on repeated calls", {
  # mvtnorm::qmvt is randomised and returned values between 2.3552 and 2.3572
  # for this design. A regulatory analysis must be reproducible, which is why
  # the critical value is obtained by numerical integration instead.
  values <- vapply(
    1:5,
    function(i) msd(fathead_c1, response = "weight")$critical,
    numeric(1)
  )
  expect_equal(length(unique(values)), 1L)
})

test_that("the computed critical value agrees with mvtnorm within its error", {
  skip_if_not_installed("mvtnorm")
  fit <- msd(fathead_c1, response = "weight")

  rho <- matrix(0.5, 4, 4)
  diag(rho) <- 1
  # Averaged over several draws, because a single one carries exactly the
  # Monte Carlo error this comparison is meant to bound: individual values
  # between 2.3552 and 2.3585 have been observed for this design. The
  # tolerance accommodates that spread rather than asserting agreement to
  # more places than the reference method can deliver.
  reference <- mean(vapply(
    1:10,
    function(i) {
      abs(mvtnorm::qmvt(0.95, tail = "lower", df = 15, corr = rho)$quantile)
    },
    numeric(1)
  ))

  expect_equal(fit$critical, reference, tolerance = 5e-3)
})
