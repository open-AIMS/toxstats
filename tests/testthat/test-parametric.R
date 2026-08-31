# Tests for the parametric comparisons in R/parametric.R
#
# Reference values come from Appendices C and D of EPA-821-R-02-013.

# fixtures ------------------------------------------------------------------

# Appendix D uses the Appendix C data with the third replicate of the
# 256 ug/L concentration presumed lost, which is row 19.
appendix_d <- function() fathead_c1[-19, ]

# dunnett -------------------------------------------------------------------

test_that("dunnett reproduces the Appendix C worked example", {
  fit <- dunnett(fathead_c1, response = "weight")

  # The manual prints t = 1.487, 0.248, 1.633, 3.251 against a critical value
  # of 2.36, and concludes NOEC 128, LOEC 256.
  expect_equal(
    round(fit$comparisons$statistic, 3),
    c(1.486, 0.248, 1.635, 3.248)
  )
  expect_equal(round(fit$critical, 2), 2.36)
  expect_equal(round(fit$sw, 3), 0.097)
  expect_equal(fit$df, 15L)
  expect_equal(fit$comparisons$significant, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(fit$noec, 128)
  expect_equal(fit$loec, 256)
  expect_true(fit$monotone)
})

test_that("dunnett reports adjusted p-values consistent with the decision", {
  fit <- dunnett(fathead_c1, response = "weight")
  # An adjusted p-value at or below alpha must agree with the comparison of
  # the statistic against the familywise critical value.
  expect_equal(fit$comparisons$p_value <= 0.05, fit$comparisons$significant)
  # The comparison whose statistic equals the critical value has p = alpha.
  at_critical <- 1 -
    toxstats:::dunnett_prob(fit$critical, rep(sqrt(0.5), 4), 15)
  expect_equal(at_critical, 0.05, tolerance = 1e-6)
})

test_that("dunnett refuses an unbalanced design", {
  # The manual sends unequal replication to the t test with Bonferroni's
  # adjustment, so returning a Dunnett result here would be wrong.
  expect_error(
    dunnett(appendix_d(), response = "weight"),
    regexp = "same number of replicates"
  )
})

# bonferroni_t --------------------------------------------------------------

test_that("bonferroni_t reproduces the Appendix D worked example", {
  fit <- bonferroni_t(appendix_d(), response = "weight")

  # The manual prints t = 1.623, 0.220, 1.782, 4.022 against a critical value
  # of 2.510, and concludes NOEC 128, LOEC 256.
  #
  # The 64 ug/L value is a slip in the manual: from its own means and pooled
  # standard deviation the statistic is 0.270, and the three others agree to
  # within rounding. Both are far below the critical value.
  expect_equal(
    round(fit$comparisons$statistic, 3),
    c(1.622, 0.270, 1.785, 4.028)
  )
  expect_equal(round(fit$critical, 3), 2.510)
  expect_equal(fit$df, 14L)
  expect_equal(fit$noec, 128)
  expect_equal(fit$loec, 256)
})

test_that("the Appendix D mean square is misprinted", {
  # Table D.3 gives SSW = 0.111 on 14 degrees of freedom and then a mean
  # square of 0.0029. 0.111 / 14 is 0.0079, which is what the manual's own
  # printed t values were computed from.
  fit <- bonferroni_t(appendix_d(), response = "weight")
  expect_equal(round(fit$sw^2, 4), 0.0079)
  expect_equal(round(fit$sw^2 * fit$df, 3), 0.111)
})

test_that("the Bonferroni critical value reproduces EPA Table D.5", {
  # Table D.5 is exactly qt(1 - alpha / k, df); its row for infinite degrees
  # of freedom at K = 4 prints 2.242 against a computed 2.2414.
  expect_equal(round(qt(1 - 0.05 / 4, 14), 3), 2.510)
  expect_equal(round(qnorm(1 - 0.05 / 4), 3), 2.241)
})

test_that("bonferroni_t accepts a balanced design as well", {
  # Nothing prevents it; it is simply more conservative than Dunnett's.
  balanced <- bonferroni_t(fathead_c1, response = "weight")
  expect_gt(
    balanced$critical,
    dunnett(fathead_c1, response = "weight")$critical
  )
  expect_equal(balanced$noec, 128)
})

# dunn_sidak_t --------------------------------------------------------------

test_that("dunn_sidak_t sits between Dunnett and Bonferroni", {
  x <- tox_data(fathead_c1, response = "weight")
  expect_lt(dunnett(x)$critical, dunn_sidak_t(x)$critical)
  expect_lt(dunn_sidak_t(x)$critical, bonferroni_t(x)$critical)
})

test_that("dunn_sidak_t is flagged as off the EPA flowchart", {
  fit <- dunn_sidak_t(fathead_c1, response = "weight")
  expect_match(fit$reference, "Not on the EPA flowchart")
})

# welch_t -------------------------------------------------------------------

test_that("welch_t gives each comparison its own degrees of freedom", {
  fit <- welch_t(fathead_c1, response = "weight")

  expect_equal(nrow(fit$comparisons), 4L)
  expect_true(all(fit$comparisons$df < 15))
  expect_match(fit$method, "not an EPA multi-concentration method")
})

test_that("welch_t agrees with stats::t.test on a single comparison", {
  two_level <- fathead_c1[fathead_c1$conc %in% c(0, 256), ]
  fit <- welch_t(two_level, response = "weight", adjust = "none")
  reference <- t.test(
    two_level$weight[two_level$conc == 0],
    two_level$weight[two_level$conc == 256],
    alternative = "greater"
  )
  expect_equal(fit$comparisons$statistic, unname(reference$statistic))
  expect_equal(fit$comparisons$p_value, reference$p.value)
})

# direction -----------------------------------------------------------------

test_that("an increasing endpoint reverses the tail", {
  # The same data read as an endpoint that toxicity increases must give the
  # mirror image of the decreasing case.
  decreasing <- dunnett(fathead_c1, response = "weight")
  increasing <- dunnett(
    fathead_c1,
    response = "weight",
    direction = "increasing"
  )
  expect_equal(
    increasing$comparisons$statistic,
    -decreasing$comparisons$statistic
  )
  expect_true(all(!increasing$comparisons$significant))
})

# NOEC and LOEC derivation --------------------------------------------------

test_that("endpoints fall outside the range when nothing or everything is significant", {
  none <- toxstats:::derive_noec_loec(c(1, 2, 4), c(FALSE, FALSE, FALSE))
  expect_equal(none$noec, 4)
  expect_true(is.na(none$loec))

  all_of_them <- toxstats:::derive_noec_loec(c(1, 2, 4), c(TRUE, TRUE, TRUE))
  expect_true(is.na(all_of_them$noec))
  expect_equal(all_of_them$loec, 1)
})

test_that("a non-monotone significance pattern is reported, not smoothed", {
  # Section 9.6.5.1 gives no rule for this, so the LOEC is the lowest
  # significant concentration and the pattern is flagged.
  broken <- toxstats:::derive_noec_loec(c(1, 2, 4), c(FALSE, TRUE, FALSE))
  expect_equal(broken$loec, 2)
  expect_equal(broken$noec, 1)
  expect_false(broken$monotone)
})

# printing ------------------------------------------------------------------

test_that("comparisons print their endpoints and reference", {
  expect_output(print(dunnett(fathead_c1, response = "weight")), "Appendix C")
  expect_output(print(dunnett(fathead_c1, response = "weight")), "NOEC 128")
  expect_output(print(dunnett(fathead_c1, response = "weight")), "LOEC 256")
})

test_that("an endpoint outside the tested range prints as an inequality", {
  x <- tox_data(fathead_c1, response = "weight", direction = "increasing")
  expect_output(print(dunnett(x)), "LOEC > 256")
})
