# Tests for the linear interpolation method in R/icp.R
#
# Reference values come from Appendix M of EPA-821-R-02-013 and from the ICPIN
# program output reproduced there as Figure 9.

# the estimate ---------------------------------------------------------------

test_that("icp reproduces the Appendix M smoothed means", {
  fit <- icp(ceriodaphnia_m1, response = "young", nboot = 40, seed = 1)

  # Table M.2 prints 28.75 across the control and the three lowest
  # concentrations, then 9.4 and 0.
  expect_equal(fit$working$smoothed, c(28.75, 28.75, 28.75, 28.75, 9.4, 0))
  expect_equal(fit$working$mean, c(22.4, 26.3, 34.6, 31.7, 9.4, 0))
  # The ICPIN output prints these standard deviations.
  expect_equal(
    round(fit$working$sd, 3),
    c(6.931, 8.001, 4.835, 2.946, 3.893, 0)
  )
})

test_that("icp reproduces the Appendix M point estimates", {
  fit <- icp(
    ceriodaphnia_m1,
    response = "young",
    p = c(25, 50),
    nboot = 40,
    seed = 1
  )

  # The manual's program prints 8.5715 for the IC25 and 10.89 for the IC50.
  expect_equal(round(fit$estimates$estimate, 4), c(8.5715, 10.8931))
  expect_equal(fit$estimates$endpoint, c("IC25", "IC50"))
})

test_that("the interpolation follows the manual's own arithmetic", {
  # Section 6.3: a 25 per cent reduction gives a target of 21.56 young, which
  # falls between the smoothed 28.75 at 6.25 per cent and 9.4 at 12.5 per cent.
  smoothed <- c(28.75, 28.75, 28.75, 28.75, 9.4, 0)
  conc <- c(0, 1.56, 3.12, 6.25, 12.5, 25)
  target <- smoothed[1] * 0.75

  expect_equal(round(target, 2), 21.56)
  expect_equal(
    toxstats:::interpolate_icp(conc, smoothed, 25)$estimate,
    6.25 + (target - 28.75) * (12.5 - 6.25) / (9.4 - 28.75)
  )
})

test_that("smoothing uses equal weights on the concentration means", {
  # 28.75 is the unweighted mean of 22.4, 26.3, 34.6 and 31.7. Weighting by
  # replication would give the same here because the design is balanced, so
  # the check is made directly on the smoother.
  expect_equal(mean(c(22.4, 26.3, 34.6, 31.7)), 28.75)
  expect_equal(
    smooth_monotone(c(22.4, 26.3, 34.6, 31.7, 9.4, 0)),
    c(28.75, 28.75, 28.75, 28.75, 9.4, 0)
  )
})

# the bootstrap --------------------------------------------------------------

test_that("a seed makes the interval reproducible", {
  once <- icp(ceriodaphnia_m1, response = "young", nboot = 60, seed = 99)
  twice <- icp(ceriodaphnia_m1, response = "young", nboot = 60, seed = 99)

  expect_equal(once$estimates$lower, twice$estimates$lower)
  expect_equal(once$estimates$upper, twice$estimates$upper)
  expect_equal(once$boot$replicates, twice$boot$replicates)
})

test_that("the interval brackets the estimate and is close to the manual", {
  # The manual's limits of 8.3112 and 9.0418 came from 80 resamples drawn with
  # a seed no other program can reproduce, so they cannot be matched exactly.
  # What can be checked is that the interval is of the right size and position.
  fit <- icp(ceriodaphnia_m1, response = "young", nboot = 1000, seed = 3)

  expect_lt(fit$estimates$lower, fit$estimates$estimate)
  expect_gt(fit$estimates$upper, fit$estimates$estimate)
  expect_lt(abs(fit$estimates$lower - 8.3112), 0.2)
  expect_lt(abs(fit$estimates$upper - 9.0418), 0.3)
})

test_that("omitting the re-smoothing biases the bootstrap badly", {
  # The manual's step list places smoothing before resampling, which reads as
  # though it happens once. Leaving it out of the loop is not a small error:
  # unsmoothed resampled means are not monotone, so the interpolation selects
  # the wrong bracketing pair and the resampled estimates drift away from the
  # point estimate entirely. This is asserted rather than described, because
  # it is the detail of the method most easily got wrong.
  values <- split(ceriodaphnia_m1$young, ceriodaphnia_m1$conc)
  conc <- as.numeric(names(values))

  set.seed(11)
  drawn <- replicate(
    500,
    lapply(values, function(v) sample(v, length(v), replace = TRUE)),
    simplify = FALSE
  )
  resmoothed <- vapply(
    drawn,
    function(d) {
      means <- smooth_monotone(
        vapply(d, mean, numeric(1)),
        direction = "decreasing"
      )
      toxstats:::interpolate_icp(conc, means, 25)$estimate
    },
    numeric(1)
  )
  unsmoothed <- vapply(
    drawn,
    function(d) {
      toxstats:::interpolate_icp(conc, vapply(d, mean, numeric(1)), 25)$estimate
    },
    numeric(1)
  )

  # Re-smoothing keeps the bootstrap distribution centred on the estimate.
  expect_lt(abs(mean(resmoothed) - 8.5715), 0.1)
  # Omitting it shifts the centre by about two units and roughly triples the
  # spread.
  expect_gt(mean(unsmoothed), 10)
  expect_gt(sd(unsmoothed), 3 * sd(resmoothed))
})

test_that("the order statistic rule follows the manual's description", {
  # "The second smallest and second largest" of 80 resamples, which is
  # floor((80 + 1) * 0.025) = 2 in from each end.
  fit <- icp(ceriodaphnia_m1, response = "young", nboot = 80, seed = 5)
  draws <- sort(fit$boot$replicates[1, ])

  expect_equal(fit$estimates$lower, draws[2])
  expect_equal(fit$estimates$upper, draws[79])
})

test_that("a wider confidence level takes fewer order statistics in", {
  narrow <- icp(ceriodaphnia_m1, response = "young", nboot = 200, seed = 5)
  wide <- icp(
    ceriodaphnia_m1,
    response = "young",
    nboot = 200,
    seed = 5,
    ci_level = 0.99
  )

  expect_lt(wide$estimates$lower, narrow$estimates$lower)
  expect_gt(wide$estimates$upper, narrow$estimates$upper)
})

test_that("the bootstrap summary is reported", {
  fit <- icp(ceriodaphnia_m1, response = "young", nboot = 100, seed = 2)

  expect_equal(fit$boot$nboot, 100)
  expect_equal(fit$boot$seed, 2)
  expect_equal(dim(fit$boot$replicates), c(1L, 100L))
  expect_equal(fit$estimates$n_undefined, 0)
  expect_gt(fit$estimates$boot_sd, 0)
})

# edge cases -----------------------------------------------------------------

test_that("a response that never reaches the target is reported as a bound", {
  # Removing the two concentrations that show any inhibition leaves nothing
  # below 75 per cent of the control.
  result <- toxstats:::interpolate_icp(
    c(0, 1.56, 3.12, 6.25),
    rep(28.75, 4),
    25
  )
  expect_equal(result$estimate, Inf)
  expect_match(result$bound, "above the highest")
})

test_that("fewer than seven replicates warns about the expanded interval", {
  # The manual reports a second, wider interval in that case, whose definition
  # is not publicly retrievable and is therefore not implemented.
  thin <- ceriodaphnia_m1[ceriodaphnia_m1$replicate %in% as.character(1:5), ]
  expect_warning(
    icp(thin, response = "young", nboot = 40, seed = 1),
    regexp = "expanded interval"
  )
})

test_that("icp validates its arguments", {
  expect_error(
    icp(ceriodaphnia_m1, response = "young", p = 0),
    regexp = "between 1 and 99"
  )
  expect_error(
    icp(ceriodaphnia_m1, response = "young", nboot = 5),
    regexp = "greater than or equal to 20"
  )
  expect_error(
    icp(ceriodaphnia_m1, response = "young", ci_level = 2),
    regexp = "between 0 and 1"
  )
})

test_that("icp needs at least three concentrations", {
  two <- ceriodaphnia_m1[ceriodaphnia_m1$conc %in% c(0, 25), ]
  expect_error(
    icp(two, response = "young", nboot = 40),
    regexp = "three concentrations"
  )
})

# printing -------------------------------------------------------------------

test_that("icp prints its reference", {
  fit <- icp(ceriodaphnia_m1, response = "young", nboot = 40, seed = 1)
  expect_output(print(fit), "Appendix M")
  expect_output(print(fit), "IC25")
})
