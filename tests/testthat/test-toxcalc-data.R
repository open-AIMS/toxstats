# Tests for the input contract in R/toxcalc_data.R

# fixtures ------------------------------------------------------------------
#
# These are SYNTHETIC, used only to exercise the input contract. They are not
# EPA data and carry no published expected values. The real worked-example
# fixtures are extracted from the manuals in Phase 2, where they are needed to
# validate the statistical methods.

growth_data <- function() {
  data.frame(
    conc = rep(c(0, 32, 64, 128, 256), each = 4),
    growth = c(
      0.711,
      0.662,
      0.054,
      0.785,
      0.517,
      0.501,
      0.723,
      0.760,
      0.602,
      0.669,
      0.694,
      0.706,
      0.348,
      0.400,
      0.041,
      0.512,
      0.216,
      0.277,
      0.328,
      0.347
    )
  )
}

survival_data <- function() {
  data.frame(
    conc = rep(c(0, 25, 50, 100), each = 2),
    dead = c(0, 1, 1, 0, 4, 5, 10, 10),
    exposed = rep(10, 8)
  )
}

# continuous ----------------------------------------------------------------

test_that("toxcalc_data accepts a continuous endpoint", {
  x <- toxcalc_data(growth_data(), response = "growth")

  expect_s3_class(x, "toxcalc_data")
  expect_equal(x$type, "continuous")
  expect_equal(x$direction, "decreasing")
  expect_equal(x$control, 0)
  expect_equal(nrow(x$replicates), 20)
  expect_equal(nrow(x$pooled), 5)
  expect_false("proportion" %in% names(x$replicates))
})

test_that("pooled summaries match a direct calculation", {
  raw <- growth_data()
  x <- toxcalc_data(raw, response = "growth")
  control <- raw$growth[raw$conc == 0]

  expect_equal(x$pooled$conc, c(0, 32, 64, 128, 256))
  expect_equal(x$pooled$n_rep, rep(4L, 5))
  expect_equal(x$pooled$mean[1], mean(control))
  expect_equal(x$pooled$sd[1], stats::sd(control))
  expect_equal(x$pooled$var[1], stats::var(control))
})

# quantal -------------------------------------------------------------------

test_that("toxcalc_data accepts a quantal endpoint", {
  x <- toxcalc_data(
    survival_data(),
    response = "dead",
    n_exposed = "exposed",
    type = "quantal"
  )

  expect_equal(x$type, "quantal")
  expect_equal(x$replicates$proportion, c(0, 0.1, 0.1, 0, 0.4, 0.5, 1, 1))
  # Pooled counts are what the point-estimation methods use.
  expect_equal(x$pooled$n_exposed, rep(20, 4))
  expect_equal(x$pooled$n_affected, c(1, 1, 9, 20))
  expect_equal(x$pooled$proportion, c(0.05, 0.05, 0.45, 1))
  # The mean is the mean of replicate proportions, which the hypothesis
  # flowchart uses; it need not equal the pooled proportion when replicates
  # differ in size.
  expect_equal(x$pooled$mean, c(0.05, 0.05, 0.45, 1))
})

test_that("quantal data requires n_exposed", {
  expect_error(
    toxcalc_data(survival_data(), response = "dead", type = "quantal"),
    regexp = "n_exposed"
  )
})

test_that("quantal data rejects impossible counts", {
  bad <- survival_data()
  bad$dead[1] <- 11
  expect_error(
    toxcalc_data(
      bad,
      response = "dead",
      n_exposed = "exposed",
      type = "quantal"
    ),
    regexp = "must not exceed"
  )

  negative <- survival_data()
  negative$dead[1] <- -1
  expect_error(
    toxcalc_data(
      negative,
      response = "dead",
      n_exposed = "exposed",
      type = "quantal"
    ),
    regexp = "non-negative"
  )
})

# missing values ------------------------------------------------------------

test_that("rows with a missing response are dropped and counted", {
  raw <- growth_data()
  raw$growth[c(2, 5)] <- NA
  x <- toxcalc_data(raw, response = "growth")

  expect_equal(x$n_dropped, 2L)
  expect_equal(nrow(x$replicates), 18)
  # The affected concentrations now have three replicates, which is what the
  # flowchart's balance check must see.
  expect_equal(x$pooled$n_rep, c(3L, 3L, 4L, 4L, 4L))
})

test_that("a wholly missing response is an error", {
  raw <- growth_data()
  raw$growth <- NA_real_
  expect_error(
    toxcalc_data(raw, response = "growth"),
    regexp = "no rows with a non-missing response"
  )
})

# validation ----------------------------------------------------------------

test_that("the control must be a concentration present in the data", {
  expect_error(
    toxcalc_data(growth_data(), response = "growth", control = 5),
    regexp = "must be one of the concentrations"
  )
})

test_that("column names must exist in the data", {
  expect_error(
    toxcalc_data(growth_data(), response = "biomass"),
    regexp = "biomass"
  )
  expect_error(
    toxcalc_data(growth_data(), conc = "dose", response = "growth"),
    regexp = "dose"
  )
})

test_that("negative concentrations and empty data are rejected", {
  negative <- growth_data()
  negative$conc[1] <- -1
  expect_error(
    toxcalc_data(negative, response = "growth"),
    regexp = "greater than or equal to"
  )
  expect_error(
    toxcalc_data(growth_data()[0, ], response = "growth"),
    regexp = "at least one row"
  )
})

test_that("direction must be one of the two recognised values", {
  expect_error(
    toxcalc_data(growth_data(), response = "growth", direction = "flat"),
    regexp = "should be one of"
  )
  x <- toxcalc_data(
    growth_data(),
    response = "growth",
    direction = "increasing"
  )
  expect_equal(x$direction, "increasing")
})

# replicate labels ----------------------------------------------------------

test_that("replicate labels are generated when not supplied", {
  x <- toxcalc_data(growth_data(), response = "growth")
  expect_equal(x$replicates$replicate[1:4], as.character(1:4))
})

test_that("supplied replicate labels are retained", {
  raw <- growth_data()
  raw$rep <- rep(LETTERS[1:4], times = 5)
  x <- toxcalc_data(raw, response = "growth", replicate = "rep")
  expect_equal(sort(unique(x$replicates$replicate)), LETTERS[1:4])
})

# methods -------------------------------------------------------------------

test_that("methods return the documented pieces", {
  x <- toxcalc_data(growth_data(), response = "growth")

  expect_identical(summary(x), x$pooled)
  expect_identical(as.data.frame(x), x$replicates)
  expect_output(print(x), "toxcalc_data")
  expect_output(print(x), "balanced")
  # print returns its input invisibly
  capture.output(returned <- print(x))
  expect_identical(returned, x)
})

test_that("print reports an unbalanced design as such", {
  raw <- growth_data()
  raw$growth[1] <- NA
  expect_output(print(toxcalc_data(raw, response = "growth")), "unbalanced")
})
