# Tests for the flowchart driver in R/tox_test.R and the engine in R/flowchart.R

# fixtures ------------------------------------------------------------------

# A synthetic quantal design, used only to exercise the transform branch of the
# flowchart. It is not EPA data and carries no published expected values.
quantal_design <- function() {
  data.frame(
    conc = rep(c(0, 10, 20, 40), each = 4),
    affected = c(
      0,
      1,
      0,
      1,
      1,
      2,
      1,
      0,
      3,
      4,
      3,
      5,
      9,
      10,
      8,
      9
    ),
    exposed = rep(10, 16)
  )
}

# end to end: the parametric branch ------------------------------------------

test_that("the flowchart selects Dunnett on the Appendix C data", {
  fit <- tox_test(fathead_c1, response = "weight")

  expect_s3_class(fit, "tox_test")
  expect_equal(fit$selected, "dunnett")
  expect_false(fit$overridden)
  expect_equal(fit$transform, "none")
  # The manual's own conclusion for these data.
  expect_equal(fit$noec, 128)
  expect_equal(fit$loec, 256)
})

test_that("the Appendix C decision trail records every branch point", {
  trail <- decisions(tox_test(fathead_c1, response = "weight"))

  expect_equal(
    trail$node,
    c("transform", "normality", "variance", "replication_p", "test")
  )
  expect_equal(trail$step, 1:5)
  expect_equal(trail$answer[trail$node == "normality"], "yes")
  expect_equal(trail$answer[trail$node == "variance"], "yes")
  expect_equal(trail$answer[trail$node == "replication_p"], "yes")
  # Every row cites the manual, which is the point of the trail.
  expect_true(all(nzchar(trail$reference)))
  expect_true(all(grepl("EPA-821", trail$reference)))
})

test_that("the two significance levels are applied where they belong", {
  fit <- tox_test(fathead_c1, response = "weight")
  trail <- decisions(fit)

  expect_equal(fit$alpha, 0.05)
  expect_equal(fit$alpha_assumption, 0.01)
  # The assumption rows quote 0.01, not 0.05.
  expect_match(trail$criterion[trail$node == "normality"], "0.01")
  expect_match(trail$criterion[trail$node == "variance"], "0.01")
})

# end to end: the non-parametric branch --------------------------------------

test_that("the flowchart selects Steel on the Appendix E data", {
  # The 50 per cent concentration is dropped because Appendix G found survival
  # there significantly reduced.
  fit <- tox_test(ceriodaphnia_e1, response = "young", exclude = 50)

  expect_equal(fit$selected, "steel")
  expect_equal(fit$noec, 3)
  expect_equal(fit$loec, 6)
})

test_that("Steel is reached because normality fails, not by assumption", {
  fit <- tox_test(ceriodaphnia_e1, response = "young", exclude = 50)
  trail <- decisions(fit)

  expect_equal(trail$answer[trail$node == "normality"], "no")
  expect_lt(trail$p_value[trail$node == "normality"], 0.01)
  # The variance test is never reached: the normality arm goes straight to the
  # replication questions.
  expect_false("variance" %in% trail$node)
  expect_true("min_replicates" %in% trail$node)
})

# exclusions -----------------------------------------------------------------

test_that("caller-supplied exclusions are recorded and applied", {
  fit <- tox_test(ceriodaphnia_e1, response = "young", exclude = 50)

  expect_equal(fit$excluded$conc, 50)
  expect_match(fit$excluded$reason, "9.5.2")
  expect_false(50 %in% fit$working$pooled$conc)
  # The original data is retained untouched, because point estimation uses it.
  expect_true(50 %in% fit$data$pooled$conc)
})

test_that("a completely responding concentration is excluded automatically", {
  raw <- quantal_design()
  raw$affected[raw$conc == 40] <- 10
  fit <- tox_test(
    raw,
    response = "affected",
    n_exposed = "exposed",
    type = "quantal"
  )

  expect_equal(fit$excluded$conc, 40)
  expect_match(fit$excluded$reason, "complete response")
})

test_that("the control cannot be excluded", {
  expect_error(
    tox_test(fathead_c1, response = "weight", exclude = 0),
    regexp = "control concentration cannot be excluded"
  )
})

test_that("an unknown exclusion is rejected", {
  expect_error(
    tox_test(fathead_c1, response = "weight", exclude = 999),
    regexp = "999"
  )
})

# the transform branch -------------------------------------------------------

test_that("quantal data is arc sine transformed before the assumption tests", {
  fit <- tox_test(
    quantal_design(),
    response = "affected",
    n_exposed = "exposed",
    type = "quantal"
  )

  expect_equal(fit$transform, "arcsine_sqrt")
  expect_equal(decisions(fit)$answer[1], "yes")
  expect_match(decisions(fit)$outcome[1], "arc sine")
  # The working response is the transformed one, on the radian scale.
  expect_equal(fit$working$type, "continuous")
  expect_true(all(fit$working$replicates$response <= pi / 2))
  expect_equal(
    fit$working$replicates$response,
    arcsine_sqrt(fit$working$replicates$proportion, n = 10)
  )
})

test_that("continuous data is not transformed", {
  fit <- tox_test(fathead_c1, response = "weight")
  expect_equal(fit$transform, "none")
  expect_equal(decisions(fit)$answer[1], "no")
  expect_equal(fit$working$replicates$response, fit$data$replicates$response)
})

# forcing a test -------------------------------------------------------------

test_that("forcing a different test warns and is recorded", {
  expect_warning(
    fit <- tox_test(fathead_c1, response = "weight", test = "steel"),
    regexp = "flowchart selected"
  )

  expect_true(fit$overridden)
  expect_equal(fit$selected, "dunnett")
  expect_equal(fit$comparison$test, "steel")
  expect_true("override" %in% decisions(fit)$node)
  expect_match(
    decisions(fit)$outcome[decisions(fit)$node == "override"],
    "user forced steel"
  )
})

test_that("forcing the test the flowchart chose is not an override", {
  expect_no_warning(
    fit <- tox_test(fathead_c1, response = "weight", test = "dunnett")
  )
  expect_false(fit$overridden)
  expect_false("override" %in% decisions(fit)$node)
})

test_that("an unknown test name is rejected", {
  expect_error(
    tox_test(fathead_c1, response = "weight", test = "tukey"),
    regexp = "tukey"
  )
})

# the PMSD lower-bound override ----------------------------------------------

test_that("no override is applied without bounds", {
  fit <- tox_test(fathead_c1, response = "weight")
  expect_false("pmsd_override" %in% decisions(fit)$node)
  expect_equal(fit$loec, 256)
})

test_that("a significant concentration below the lower bound is reversed", {
  # Section 10.2.8.2.5. The 256 ug/L concentration differs from the control by
  # about 33 per cent, so a lower bound above that must reverse it.
  plain <- tox_test(fathead_c1, response = "weight")
  expect_equal(plain$loec, 256)

  overridden <- tox_test(
    fathead_c1,
    response = "weight",
    pmsd_bounds = c(40, 60)
  )
  expect_true("pmsd_override" %in% decisions(overridden)$node)
  expect_true(is.na(overridden$loec))
  expect_equal(overridden$noec, 256)
  expect_match(
    decisions(overridden)$reference[
      decisions(overridden)$node == "pmsd_override"
    ],
    "10.2.8.2.5"
  )
})

test_that("bounds that do not bite leave the result alone", {
  fit <- tox_test(
    fathead_c1,
    response = "weight",
    pmsd_bounds = "fathead_growth"
  )
  expect_false("pmsd_override" %in% decisions(fit)$node)
  expect_equal(fit$loec, 256)
  expect_equal(unname(fit$pmsd$status[1]), "within")
})

# the no-valid-test terminal -------------------------------------------------

test_that("a design with no EPA-sanctioned test is refused explicitly", {
  # Non-normal residuals with three replicates: section 9.4.5.2 rules out the
  # non-parametric branch, and the parametric branch is unavailable, so the
  # flowchart has no terminal.
  set.seed(42)
  skewed <- data.frame(
    conc = rep(c(0, 1, 2, 4), each = 3),
    y = c(
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      10,
      1000
    )
  )
  expect_error(
    tox_test(skewed, response = "y"),
    regexp = "no\\s+EPA-sanctioned test is available"
  )
})

# sensitivity ----------------------------------------------------------------

test_that("MSD and PMSD are reported for the Appendix C example", {
  fit <- tox_test(
    fathead_c1,
    response = "weight",
    pmsd_bounds = "fathead_growth"
  )

  expect_equal(round(unname(fit$pmsd$msd$msd[1]), 3), 0.162)
  expect_equal(round(unname(fit$pmsd$pmsd[1]), 1), 23.9)
})

test_that("MSD is computed even when a non-parametric test was selected", {
  # Section 10.2.8 asks for it parametrically regardless of the branch taken.
  fit <- tox_test(ceriodaphnia_e1, response = "young", exclude = 50)
  expect_equal(fit$selected, "steel")
  expect_false(is.null(fit$pmsd))
  expect_gt(unname(fit$pmsd$msd$msd[1]), 0)
})

# methods --------------------------------------------------------------------

test_that("as.data.frame returns one tidy row per endpoint", {
  out <- as.data.frame(
    tox_test(fathead_c1, response = "weight", pmsd_bounds = "fathead_growth")
  )

  expect_equal(out$endpoint, c("NOEC", "LOEC", "MSD", "PMSD"))
  expect_equal(out$value[out$endpoint == "NOEC"], 128)
  expect_equal(out$value[out$endpoint == "LOEC"], 256)
  expect_equal(round(out$value[out$endpoint == "MSD"], 3), 0.162)
  expect_equal(round(out$value[out$endpoint == "PMSD"], 1), 23.9)
  expect_true(all(out$method == "dunnett"))
  expect_true(all(out$selected_by_flowchart == "dunnett"))
  expect_false(any(out$overridden))
})

test_that("print is brief and summary carries the trail", {
  fit <- tox_test(fathead_c1, response = "weight")

  expect_output(print(fit), "EPA WET hypothesis test")
  expect_output(print(fit), "NOEC 128")
  expect_output(print(fit), "summary\\(\\) for the decision trail")

  expect_output(summary(fit), "Appendix B, section 2.1")
  expect_output(summary(fit), "Is replication equal")
})

test_that("summary reports an override prominently", {
  suppressWarnings(
    fit <- tox_test(fathead_c1, response = "weight", test = "steel")
  )
  expect_output(print(fit), "OVERRIDDEN")
})

test_that("the Appendix C summary is stable", {
  expect_snapshot(
    summary(tox_test(
      fathead_c1,
      response = "weight",
      pmsd_bounds = "fathead_growth"
    ))
  )
})

test_that("the Appendix E summary is stable", {
  expect_snapshot(
    summary(tox_test(ceriodaphnia_e1, response = "young", exclude = 50))
  )
})

# the engine -----------------------------------------------------------------

test_that("the chart is well formed", {
  chart <- toxstats:::flowchart_hypothesis_multi()
  terminals <- names(toxstats:::flowchart_terminals())

  # Every destination is either another node or a known terminal.
  destinations <- unique(c(chart$yes, chart$no))
  known <- c(chart$node, terminals, "@no_valid_test")
  expect_true(all(destinations %in% known))
  # Every node is reachable from the first.
  expect_true(all(chart$node[-1] %in% destinations))
})

test_that("walk_flowchart rejects an unknown node", {
  chart <- toxstats:::flowchart_hypothesis_multi()
  chart$yes[chart$node == "transform"] <- "nowhere"
  chart$no[chart$node == "transform"] <- "nowhere"
  expect_error(
    toxstats:::walk_flowchart(
      chart,
      list(
        working = tox_data(fathead_c1, response = "weight"),
        alpha_assumption = 0.01,
        assumptions = list()
      )
    ),
    regexp = "Unknown flowchart node"
  )
})

# running both branches ------------------------------------------------------

test_that("the default runs the hypothesis branch only", {
  fit <- tox_test(fathead_c1, response = "weight")

  expect_equal(fit$branch, "hypothesis")
  expect_null(fit$point)
  expect_equal(as.data.frame(fit)$endpoint, c("NOEC", "LOEC", "MSD", "PMSD"))
})

test_that("a continuous endpoint takes the interpolation branch", {
  fit <- tox_test(
    ceriodaphnia_m1,
    response = "young",
    branch = "both",
    nboot = 60,
    seed = 42
  )

  expect_s3_class(fit$point, "tox_estimate")
  expect_match(fit$point$method, "Linear interpolation")
  # The Appendix M estimate, reached through the driver rather than directly.
  expect_equal(round(fit$point$estimates$estimate, 4), 8.5715)
})

test_that("a quantal endpoint takes the LC50 branch", {
  raw <- quantal_design()
  fit <- tox_test(
    raw,
    response = "affected",
    n_exposed = "exposed",
    type = "quantal",
    branch = "both"
  )

  expect_s3_class(fit$point, "tox_lc50")
  expect_true(fit$point$selected %in% unname(toxstats:::lc50_terminals()))
  expect_equal(nrow(fit$point$estimate$estimates), 2L)
})

test_that("a design with one replicate per concentration is refused clearly", {
  # acute_table20 is pooled, one row per concentration, so there is no
  # within-concentration variation and no hypothesis test is possible. The
  # error should say so rather than surfacing one from shapiro.test.
  expect_error(
    tox_test(
      acute_table20,
      response = "probit",
      n_exposed = "exposed",
      type = "quantal"
    ),
    regexp = "no variation, so normality cannot be tested"
  )
  # The point-estimation methods handle it, which is what the message says.
  expect_no_error(
    lc50(
      acute_table20,
      response = "probit",
      n_exposed = "exposed",
      type = "quantal"
    )
  )
})

test_that("point estimation uses all the data, exclusions and all", {
  # Section 9.5.2 drops a completely responding concentration from the NOEC
  # but keeps it for point estimation, so the two branches see different data.
  raw <- quantal_design()
  raw$affected[raw$conc == 40] <- 10
  fit <- tox_test(
    raw,
    response = "affected",
    n_exposed = "exposed",
    type = "quantal",
    branch = "both"
  )

  expect_equal(fit$excluded$conc, 40)
  expect_false(40 %in% fit$working$pooled$conc)
  # The point branch still sees it.
  expect_true(40 %in% fit$point$data$pooled$conc)
})

test_that("both branches appear in the tidy output together", {
  fit <- tox_test(
    ceriodaphnia_m1,
    response = "young",
    branch = "both",
    nboot = 60,
    seed = 42
  )
  out <- as.data.frame(fit)

  expect_equal(out$endpoint, c("NOEC", "LOEC", "MSD", "PMSD", "IC25"))
  # Only the point estimate carries an interval.
  expect_true(all(is.na(out$lower[out$endpoint != "IC25"])))
  expect_false(is.na(out$lower[out$endpoint == "IC25"]))
  expect_equal(out$method[out$endpoint == "IC25"], "icp")
})

test_that("p selects which point estimates are reported", {
  fit <- tox_test(
    ceriodaphnia_m1,
    response = "young",
    branch = "both",
    p = c(10, 25, 50),
    nboot = 60,
    seed = 42
  )
  expect_equal(fit$point$estimates$endpoint, c("IC10", "IC25", "IC50"))
})

test_that("print shows the point estimate when there is one", {
  fit <- tox_test(
    ceriodaphnia_m1,
    response = "young",
    branch = "both",
    nboot = 60,
    seed = 42
  )
  expect_output(print(fit), "IC25")
  expect_output(print(fit), "NOEC")
})
