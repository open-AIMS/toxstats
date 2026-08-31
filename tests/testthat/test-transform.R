# Tests for transformations in R/transform.R
#
# Reference values are taken from the worked examples printed in the EPA WET
# method manuals and are quoted to the precision the manuals print.

# arcsine_sqrt --------------------------------------------------------------

test_that("arcsine_sqrt reproduces the EPA Appendix B 4.2 worked values", {
  # EPA-821-R-02-013 Appendix B, section 4.2 prints 0.8861 for RP = 0.60.
  expect_equal(round(arcsine_sqrt(0.60), 4), 0.8861)

  # The same section prints 0.1120 for RP = 0 and 1.4588 for RP = 1, both with
  # 20 organisms per replicate. These depend on the endpoint adjustments; a
  # bare asin(sqrt(p)) gives 0 and pi/2 instead.
  expect_equal(round(arcsine_sqrt(0, n = 20), 4), 0.1120)
  expect_equal(round(arcsine_sqrt(1, n = 20), 4), 1.4588)
})

test_that("arcsine_sqrt applies the endpoint adjustments and nothing else", {
  expect_equal(arcsine_sqrt(0, n = 20), asin(sqrt(1 / 80)))
  expect_equal(arcsine_sqrt(1, n = 20), asin(sqrt(1 - 1 / 80)))
  # An interior proportion is untouched by the adjustments.
  expect_equal(arcsine_sqrt(0.5, n = 20), asin(sqrt(0.5)))
})

test_that("arcsine_sqrt is vectorised and accepts n per element", {
  expect_equal(
    arcsine_sqrt(c(0, 0.5, 1), n = c(10, 10, 40)),
    c(asin(sqrt(1 / 40)), asin(sqrt(0.5)), asin(sqrt(1 - 1 / 160)))
  )
  expect_length(arcsine_sqrt(c(0.1, 0.2, 0.3)), 3)
})

test_that("arcsine_sqrt requires n only when an endpoint is present", {
  # No 0 or 1, so n is genuinely unnecessary.
  expect_no_error(arcsine_sqrt(c(0.25, 0.75)))
  expect_error(arcsine_sqrt(c(0, 0.5)), regexp = "must be supplied")
  expect_error(arcsine_sqrt(c(0.5, 1)), regexp = "must be supplied")
})

test_that("arcsine_sqrt rejects invalid input", {
  expect_error(arcsine_sqrt(c(0.5, 1.5)), regexp = "between 0 and 1")
  expect_error(arcsine_sqrt(c(0.5, -0.1)), regexp = "between 0 and 1")
  expect_error(arcsine_sqrt(c(0.5, NA)), regexp = "missing")
  expect_error(arcsine_sqrt("a"), regexp = "numeric")
  expect_error(arcsine_sqrt(0, n = c(10, 20)), regexp = "length 1")
})

# inv_arcsine_sqrt ----------------------------------------------------------

test_that("inv_arcsine_sqrt inverts arcsine_sqrt for interior proportions", {
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_equal(inv_arcsine_sqrt(arcsine_sqrt(p)), p)
})

test_that("inv_arcsine_sqrt returns the adjusted value at the endpoints", {
  # Documented behaviour: the adjustment is not undone, because the adjusted
  # value is the one the analysis used.
  expect_equal(inv_arcsine_sqrt(arcsine_sqrt(0, n = 20)), 1 / 80)
  expect_equal(inv_arcsine_sqrt(arcsine_sqrt(1, n = 20)), 1 - 1 / 80)
})

test_that("inv_arcsine_sqrt rejects values outside the transformed range", {
  expect_error(inv_arcsine_sqrt(pi), regexp = "must be between")
  expect_error(inv_arcsine_sqrt(-0.1), regexp = "must be between")
})

# abbott --------------------------------------------------------------------

test_that("abbott adjusts for control response", {
  # Appendix K control proportion after smoothing is 0.02; the top
  # concentration smooths to 0.80, giving (0.80 - 0.02) / 0.98.
  expect_equal(abbott(0.80, p_control = 0.02), 0.78 / 0.98)
  expect_equal(abbott(0.02, p_control = 0.02), 0)
})

test_that("abbott leaves a zero control response unchanged", {
  p <- c(0, 0.25, 1)
  expect_equal(abbott(p, p_control = 0), p)
})

test_that("abbott clamps below the control by default", {
  expect_equal(abbott(0.01, p_control = 0.10), 0)
  expect_lt(abbott(0.01, p_control = 0.10, clamp = FALSE), 0)
})

test_that("abbott rejects a fully responding control", {
  # Undefined rather than infinite: such a control carries no information.
  expect_error(abbott(0.5, p_control = 1), regexp = "less than 1")
})

test_that("abbott rejects invalid input", {
  expect_error(abbott(1.5, p_control = 0.1), regexp = "between 0 and 1")
  expect_error(abbott(0.5, p_control = -0.1), regexp = "between 0 and 1")
  expect_error(abbott(0.5, p_control = c(0.1, 0.2)), regexp = "number")
  expect_error(abbott(0.5, p_control = 0.1, clamp = "yes"), regexp = "flag")
})
