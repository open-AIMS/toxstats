# Tests for monotone smoothing in R/smooth.R

# smooth_monotone -----------------------------------------------------------

test_that("smooth_monotone reproduces the EPA Appendix K smoothing", {
  # EPA-821-R-02-013 Appendix K: mortality proportions 0.05, 0.00, 0.05, 0.00,
  # 0.00 collapse to a single smoothed value of 0.02. This is the case that
  # repeated pairwise averaging gets wrong, because five values pool at once.
  observed <- c(0.05, 0.00, 0.05, 0.00, 0.00)
  expect_equal(
    smooth_monotone(observed, direction = "increasing"),
    rep(0.02, 5)
  )
})

test_that("smooth_monotone pools a run of three at once", {
  # Pairwise averaging would give 6, 6, 3 here; the correct pooled mean is 5.
  expect_equal(
    smooth_monotone(c(9, 3, 3), direction = "increasing"),
    rep(5, 3)
  )
})

test_that("smooth_monotone leaves an already monotone sequence unchanged", {
  increasing <- c(0, 0.1, 0.4, 0.9, 1)
  decreasing <- c(10, 8, 8, 3, 1)
  expect_equal(
    smooth_monotone(increasing, direction = "increasing"),
    increasing
  )
  expect_equal(
    smooth_monotone(decreasing, direction = "decreasing"),
    decreasing
  )
})

test_that("smooth_monotone enforces the requested direction", {
  y <- c(10, 12, 5)
  smoothed <- smooth_monotone(y, direction = "decreasing")
  expect_true(all(diff(smoothed) <= 0))
  expect_equal(smoothed, c(11, 11, 5))

  smoothed_up <- smooth_monotone(c(5, 2, 9), direction = "increasing")
  expect_true(all(diff(smoothed_up) >= 0))
  expect_equal(smoothed_up, c(3.5, 3.5, 9))
})

test_that("smooth_monotone preserves the weighted mean", {
  y <- c(4, 9, 2, 7, 1)
  expect_equal(mean(smooth_monotone(y, direction = "decreasing")), mean(y))

  w <- c(1, 2, 3, 2, 1)
  smoothed <- smooth_monotone(y, w = w, direction = "decreasing")
  expect_equal(
    stats::weighted.mean(smoothed, w),
    stats::weighted.mean(y, w)
  )
})

test_that("smooth_monotone respects weights", {
  # Equal weights pool 10 and 4 to 7; weighting the second value three times
  # as heavily pulls the pooled value down to 5.5.
  y <- c(10, 4)
  expect_equal(smooth_monotone(y, direction = "increasing"), c(7, 7))
  expect_equal(
    smooth_monotone(y, w = c(1, 3), direction = "increasing"),
    c(5.5, 5.5)
  )
})

test_that("smooth_monotone handles degenerate input", {
  expect_equal(smooth_monotone(5, direction = "decreasing"), 5)
  expect_equal(
    smooth_monotone(numeric(0), direction = "decreasing"),
    numeric(0)
  )
})

test_that("smooth_monotone rejects invalid input", {
  expect_error(smooth_monotone(c(1, NA)), regexp = "missing")
  expect_error(smooth_monotone("a"), regexp = "numeric")
  expect_error(smooth_monotone(c(1, 2), w = 1), regexp = "same length")
  expect_error(smooth_monotone(c(1, 2), w = c(1, 0)), regexp = "positive")
  expect_error(
    smooth_monotone(c(1, 2), direction = "sideways"),
    regexp = "should be one of"
  )
})

# is_monotone ---------------------------------------------------------------

test_that("is_monotone detects direction and reports violations", {
  flat <- is_monotone(c(10, 8, 8, 3), direction = "decreasing")
  expect_true(flat$monotone)
  expect_length(flat$violations, 0)

  broken <- is_monotone(c(10, 8, 9, 4), direction = "decreasing")
  expect_false(broken$monotone)
  expect_equal(broken$violations, 3L)
})

test_that("is_monotone treats ties as monotone in both directions", {
  expect_true(is_monotone(c(1, 1, 1), direction = "decreasing")$monotone)
  expect_true(is_monotone(c(1, 1, 1), direction = "increasing")$monotone)
})

test_that("is_monotone handles degenerate input", {
  expect_true(is_monotone(5)$monotone)
  expect_true(is_monotone(numeric(0))$monotone)
})
