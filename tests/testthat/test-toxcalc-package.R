# Tests for package-level metadata in R/toxcalc-package.R

# scaffold ------------------------------------------------------------------

# Placeholder while the package has no exported functions. It exists so that
# `R CMD check` runs a test suite rather than failing with "No test files
# found", and is replaced by real tests as each phase lands.
test_that("package metadata is available", {
  expect_true(nzchar(utils::packageDescription("toxcalc")$Version))
})
