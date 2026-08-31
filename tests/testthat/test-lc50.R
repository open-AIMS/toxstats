# Tests for the point-estimation flowchart in R/lc50.R
#
# The four columns of EPA Table 20 were each constructed to exercise one
# terminal of the Figure 6 chart, so together they test the whole chart.

table20 <- function(column) {
  tox_data(
    acute_table20,
    response = column,
    n_exposed = "exposed",
    type = "quantal"
  )
}

# routing --------------------------------------------------------------------

test_that("each Table 20 column routes to the method it was built for", {
  expect_equal(lc50(table20("graphical"))$selected, "graphical_lc50")
  expect_equal(lc50(table20("spearman_karber"))$selected, "spearman_karber")
  expect_equal(lc50(table20("trimmed"))$selected, "trimmed_spearman_karber")
  expect_equal(lc50(table20("probit"))$selected, "probit_lc")
})

test_that("all-or-nothing data reach the graphical method", {
  trail <- decisions(lc50(table20("graphical")))

  expect_equal(trail$answer[trail$node == "partial_two"], "no")
  expect_equal(trail$answer[trail$node == "partial_one"], "no")
  expect_equal(trail$statistic[trail$node == "partial_one"], 0)
  # The probit branch is never reached.
  expect_false("probit_fits" %in% trail$node)
})

test_that("partial responses are counted after smoothing and adjustment", {
  # The Spearman-Karber column shows two partial responses in the raw data but
  # only one after smoothing and Abbott adjustment. Counting the raw ones
  # sends the analysis to the probit branch, which is not what the manual
  # does, so this is the check that pins the choice.
  trail <- decisions(lc50(table20("spearman_karber")))

  expect_equal(trail$statistic[trail$node == "partial_two"], 1)
  raw <- acute_table20$spearman_karber[-1] / 20
  expect_equal(sum(raw > 0 & raw < 1), 2)
})

test_that("the Spearman-Karber requirements are checked before that method", {
  trail <- decisions(lc50(table20("spearman_karber")))

  expect_equal(trail$answer[trail$node == "partial_one"], "yes")
  expect_equal(trail$answer[trail$node == "sk_requirements"], "yes")
  expect_match(
    trail$outcome[trail$node == "sk_requirements"],
    "runs 0 to 1"
  )
})

test_that("a failing probit fit falls through to the trimmed method", {
  trail <- decisions(lc50(table20("trimmed")))

  expect_equal(trail$answer[trail$node == "partial_two"], "yes")
  expect_equal(trail$answer[trail$node == "probit_fits"], "no")
  expect_match(
    trail$outcome[trail$node == "probit_fits"],
    "does not fit"
  )
  # And the untrimmed requirements are not met either.
  expect_equal(trail$answer[trail$node == "sk_requirements"], "no")
})

test_that("an adequate probit fit is taken", {
  trail <- decisions(lc50(table20("probit")))

  expect_equal(trail$answer[trail$node == "probit_fits"], "yes")
  expect_equal(round(trail$statistic[trail$node == "probit_fits"], 3), 3.076)
  expect_false("sk_requirements" %in% trail$node)
})

# results --------------------------------------------------------------------

test_that("the driver returns the same estimate as the method called directly", {
  driven <- lc50(table20("probit"), p = c(1, 50))
  direct <- probit_lc(table20("probit"), p = c(1, 50))
  expect_equal(driven$estimate$estimates, direct$estimates)
})

test_that("every Table 20 column reproduces its published estimate", {
  expect_equal(
    round(lc50(table20("graphical"))$estimate$estimates$estimate, 1),
    35.4
  )
  expect_equal(
    round(lc50(table20("spearman_karber"))$estimate$estimates$estimate, 1),
    45.3
  )
  expect_equal(
    round(lc50(table20("trimmed"))$estimate$estimates$estimate, 2),
    77.11
  )
  expect_equal(
    round(lc50(table20("probit"))$estimate$estimates$estimate, 3),
    c(7.924, 22.872)
  )
})

# overriding -----------------------------------------------------------------

test_that("forcing a different method warns and is recorded", {
  expect_warning(
    fit <- lc50(table20("probit"), method = "spearman_karber"),
    regexp = "flowchart selected"
  )
  expect_true(fit$overridden)
  expect_equal(fit$selected, "probit_lc")
  expect_true("override" %in% decisions(fit)$node)
})

test_that("forcing the selected method is not an override", {
  expect_no_warning(fit <- lc50(table20("probit"), method = "probit_lc"))
  expect_false(fit$overridden)
})

test_that("an unknown method is rejected", {
  expect_error(lc50(table20("probit"), method = "logit"), regexp = "logit")
})

# the chart ------------------------------------------------------------------

test_that("the Figure 6 chart is well formed", {
  chart <- toxstats:::flowchart_lc50()
  destinations <- unique(c(chart$yes, chart$no))
  known <- c(chart$node, names(toxstats:::lc50_terminals()))

  expect_true(all(destinations %in% known))
  expect_true(all(chart$node[-1] %in% destinations))
  expect_setequal(
    unname(toxstats:::lc50_terminals()),
    c(
      "probit_lc",
      "spearman_karber",
      "trimmed_spearman_karber",
      "graphical_lc50"
    )
  )
})

# methods --------------------------------------------------------------------

test_that("as.data.frame carries the selection alongside the estimate", {
  out <- as.data.frame(lc50(table20("probit")))

  expect_equal(nrow(out), 2L)
  expect_equal(out$selected_by_flowchart, rep("probit_lc", 2))
  expect_false(any(out$overridden))
  expect_equal(out$conc_units, rep("%", 2))
})

test_that("print is brief and summary carries the trail", {
  fit <- lc50(table20("probit"))

  expect_output(print(fit), "EPA WET point estimate")
  expect_output(print(fit), "summary\\(\\) for the decision trail")
  expect_output(summary(fit), "Figure 6")
  expect_output(summary(fit), "Is the probit model appropriate")
})

test_that("the Table 20 summaries are stable", {
  expect_snapshot({
    summary(lc50(table20("spearman_karber")))
    summary(lc50(table20("trimmed")))
    summary(lc50(table20("probit")))
  })
})
