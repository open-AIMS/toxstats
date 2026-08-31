# Tests for Williams' test in R/williams.R
#
# Neither EPA manual gives a worked example or a table of critical values for
# this test, and Williams' own tables are not retrievable, so this
# implementation is validated against mathematical identities rather than
# published output.

# the identity that pins the simulation ---------------------------------------

test_that("one concentration gives the ordinary t critical value", {
  # With a single concentration there is no order restriction, so the isotonic
  # estimate is the mean itself and the statistic is an ordinary one-sided t.
  # The simulation must therefore reproduce qt(1 - alpha, nu). This is the
  # check that the simulation machinery is right; everything else rests on it.
  set.seed(1)
  for (df in c(10, 15, 30)) {
    simulated <- toxstats:::williams_null(
      n_control = 4,
      n_doses = 4,
      df = df,
      nsim = 200000
    )
    expect_equal(
      unname(quantile(simulated[, 1], 0.95)),
      qt(0.95, df),
      tolerance = 0.02
    )
  }
})

test_that("the critical value rises with position and stays below Dunnett's", {
  # Adding concentrations makes the maximum of the isotonic statistics larger,
  # so the critical value must increase with i. And the order restriction is
  # information Dunnett's procedure does not use, so Williams' value must be
  # the smaller of the two for the same design.
  fit <- suppressWarnings(
    williams(fathead_c1, response = "weight", nsim = 40000)
  )
  reference <- dunnett(fathead_c1, response = "weight")$critical

  expect_equal(fit$comparisons$critical[1], qt(0.95, fit$df), tolerance = 0.02)
  expect_true(all(diff(fit$comparisons$critical) > -0.01))
  expect_true(all(fit$comparisons$critical < reference))
})

test_that("the same design always gives the same critical value", {
  once <- suppressWarnings(williams(
    fathead_c1,
    response = "weight",
    nsim = 5000
  ))
  twice <- suppressWarnings(williams(
    fathead_c1,
    response = "weight",
    nsim = 5000
  ))
  expect_equal(once$comparisons$critical, twice$comparisons$critical)
})

test_that("simulating does not disturb the caller's random stream", {
  # A user's own analysis must not shift because a critical value was
  # simulated inside a function they called.
  set.seed(123)
  expected <- rnorm(3)

  set.seed(123)
  suppressWarnings(williams(fathead_c1, response = "weight", nsim = 2000))
  expect_equal(rnorm(3), expected)
})

# the isotonic step -----------------------------------------------------------

test_that("the isotonic estimate is the smallest trailing weighted average", {
  # The max-min formula reduced to its final index.
  y <- c(10, 8, 9, 4)
  w <- rep(1, 4)
  expect_equal(
    toxstats:::isotonic_last(y, w),
    min(mean(y), mean(y[2:4]), mean(y[3:4]), y[4])
  )
  # An already monotone sequence is left alone at the last position.
  expect_equal(toxstats:::isotonic_last(c(10, 8, 4), rep(1, 3)), 4)
})

test_that("the isotonic step respects weights", {
  y <- c(10, 4)
  expect_equal(toxstats:::isotonic_last(y, c(1, 1)), 4)
  # min over trailing runs: the whole run gives (1*4 + 3*10)/4 = 8.5, the
  # trailing single value gives 10, so 8.5 is the estimate.
  expect_equal(toxstats:::isotonic_last(c(4, 10), c(1, 3)), 8.5)
})

test_that("the isotonic means absorb a departure from monotonicity", {
  # The Appendix C data rise from 0.575 at 32 to 0.660 at 64, so the isotonic
  # estimate at position 2 is the average of the two.
  fit <- suppressWarnings(
    williams(fathead_c1, response = "weight", nsim = 2000)
  )
  expect_equal(
    round(fit$comparisons$isotonic[2], 4),
    round(mean(c(0.575, 0.660)), 4),
    tolerance = 1e-3
  )
  # Where the means already decrease the isotonic estimate is the mean itself.
  expect_equal(fit$comparisons$isotonic[4], fit$comparisons$mean[4])
})

# the step-down procedure -----------------------------------------------------

test_that("testing stops at the first non-significant concentration", {
  fit <- suppressWarnings(
    williams(fathead_c1, response = "weight", nsim = 20000)
  )

  # Only 256 is significant, so 128 is tested, fails, and the two below it are
  # never tested at all. That is what makes this a trend procedure.
  expect_equal(fit$comparisons$significant, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(fit$comparisons$tested, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(fit$noec, 128)
  expect_equal(fit$loec, 256)
})

test_that("Williams agrees with Dunnett on the Appendix C data", {
  # It is more powerful, so its NOEC can be lower but never higher.
  williams_fit <- suppressWarnings(
    williams(fathead_c1, response = "weight", nsim = 20000)
  )
  dunnett_fit <- dunnett(fathead_c1, response = "weight")
  expect_lte(williams_fit$noec, dunnett_fit$noec)
})

test_that("a p-value is reported from the simulated null", {
  fit <- suppressWarnings(
    williams(fathead_c1, response = "weight", nsim = 20000)
  )
  # A p-value at or below alpha must agree with the critical value comparison
  # for every concentration that was actually tested.
  tested <- fit$comparisons$tested
  expect_equal(
    fit$comparisons$p_value[tested] <= 0.05,
    fit$comparisons$significant[tested]
  )
  expect_true(all(fit$comparisons$p_value >= 0))
  expect_true(all(fit$comparisons$p_value <= 1))
})

# the monotonicity assumption -------------------------------------------------

test_that("a non-monotone pattern of means warns", {
  # The Appendix C means rise from 32 to 64, which is a departure at
  # position 2.
  expect_warning(
    williams(fathead_c1, response = "weight", nsim = 2000),
    regexp = "monotone concentration-response"
  )
})

test_that("a monotone pattern does not warn", {
  monotone <- data.frame(
    conc = rep(c(0, 1, 2, 4), each = 4),
    y = c(10, 11, 9, 10, 8, 9, 7, 8, 6, 5, 7, 6, 3, 2, 4, 3)
  )
  expect_no_warning(williams(monotone, response = "y", nsim = 2000))
})

# labelling -------------------------------------------------------------------

test_that("the result is labelled as a non-EPA extension", {
  fit <- suppressWarnings(
    williams(fathead_c1, response = "weight", nsim = 2000)
  )
  expect_match(fit$method, "not an EPA method")
  expect_match(fit$reference, "not on the EPA-821-R-02-013 Figure 2 flowchart")
  expect_match(fit$reference, "Williams \\(1971\\)")
})

test_that("the flowchart never selects Williams, but it can be forced", {
  # It is not a terminal of the chart, so reaching it always counts as an
  # override and always warns.
  expect_false("williams" %in% unname(toxstats:::flowchart_terminals()))

  # Two warnings fire here: the override and the monotonicity note. Both are
  # collected rather than matched one at a time, so neither masks the other.
  raised <- character()
  fit <- withCallingHandlers(
    tox_test(fathead_c1, response = "weight", test = "williams"),
    warning = function(w) {
      raised <<- c(raised, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("flowchart selected", raised)))
  expect_true(any(grepl("monotone", raised)))
  expect_true(fit$overridden)
  expect_equal(fit$comparison$test, "williams")
  expect_equal(fit$selected, "dunnett")
})

# validation ------------------------------------------------------------------

test_that("williams validates its arguments", {
  expect_error(
    williams(fathead_c1, response = "weight", alpha = 2),
    regexp = "between 0 and 1"
  )
  expect_error(
    williams(fathead_c1, response = "weight", nsim = 10),
    regexp = "greater than or equal to 1000"
  )
})

test_that("williams reports its simulation settings", {
  fit <- suppressWarnings(
    williams(fathead_c1, response = "weight", nsim = 3000, seed = 7)
  )
  expect_equal(fit$nsim, 3000)
  expect_equal(fit$seed, 7)
  expect_false(fit$monotone_observed)
  # The Monte Carlo standard error is reported so the precision is visible.
  expect_true(all(fit$comparisons$mc_se > 0))
  expect_true(all(fit$comparisons$mc_se < 0.1))
})
