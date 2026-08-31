# Tests for the non-parametric comparisons in R/nonparametric.R
#
# Reference values come from Appendices E, F and G of EPA-821-R-02-013.

# fixtures ------------------------------------------------------------------

# The 50 per cent concentration is excluded from the reproduction analysis
# because survival there was significantly reduced; see the Appendix G example.
appendix_e <- function() ceriodaphnia_e1[ceriodaphnia_e1$conc < 50, ]

# Appendix F presumes two males occurred, one in the control and one at
# 12 per cent, giving unequal replication.
appendix_f <- function() {
  x <- appendix_e()
  x[
    !(x$conc == 0 & x$replicate == "1") &
      !(x$conc == 12 & x$replicate == "5"),
  ]
}

# steel ---------------------------------------------------------------------

test_that("steel reproduces the Appendix E rank sums", {
  fit <- steel(appendix_e(), response = "young")

  # The manual prints 84, 64, 76 and 55. Its own Table E.3 ranks for the
  # 6 per cent concentration sum to 63.5, so the printed 64 is a slip.
  expect_equal(fit$comparisons$rank_sum, c(84, 63.5, 76, 55))
  expect_equal(fit$comparisons$conc, c(3, 6, 12, 25))
})

test_that("steel reproduces the Appendix E conclusion", {
  fit <- steel(appendix_e(), response = "young")

  # The manual's critical rank sum is 76 and a rank sum at or below it is
  # significant, so 6, 12 and 25 per cent are significant and 3 per cent is
  # not, giving NOEC 3 and LOEC 6.
  expect_equal(fit$comparisons$significant, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(fit$noec, 3)
  expect_equal(fit$loec, 6)
  expect_true(fit$monotone)
})

test_that("the 12 per cent comparison is marginal, as the manual has it", {
  # Its rank sum equals the tabulated critical value exactly, and the adjusted
  # p-value is just below 0.05. This is the case most likely to diverge if the
  # multiplicity adjustment were wrong, so it is asserted explicitly.
  fit <- steel(appendix_e(), response = "young")
  p <- fit$comparisons$p_value[fit$comparisons$conc == 12]
  expect_lt(p, 0.05)
  expect_gt(p, 0.04)
})

test_that("steel refuses an unbalanced design", {
  expect_error(
    steel(appendix_f(), response = "young"),
    regexp = "same number of replicates"
  )
})

test_that("steel refuses fewer than four replicates", {
  # Section 9.4.5.2 requires at least four for the non-parametric branch.
  thin <- appendix_e()[appendix_e()$replicate %in% c("1", "2", "3"), ]
  expect_error(
    steel(thin, response = "young"),
    regexp = "at least four replicates"
  )
})

# wilcoxon_rank_sum ---------------------------------------------------------

test_that("wilcoxon_rank_sum reproduces the Appendix F worked example", {
  fit <- wilcoxon_rank_sum(appendix_f(), response = "young")

  # The manual prints rank sums of 79, 57, 58 and 55, and concludes that
  # 6, 12 and 25 per cent differ from the control, giving NOEC 3 and LOEC 6.
  expect_equal(fit$comparisons$rank_sum, c(79, 57, 58, 55))
  expect_equal(fit$comparisons$n, c(10L, 10L, 9L, 10L))
  expect_equal(fit$comparisons$significant, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(fit$noec, 3)
  expect_equal(fit$loec, 6)
})

test_that("wilcoxon_rank_sum applies the Bonferroni adjustment", {
  fit <- wilcoxon_rank_sum(appendix_f(), response = "young")
  expect_equal(fit$per_comparison_alpha, 0.05 / 4)
  # The 3 per cent comparison is significant at 0.05 but not after adjustment,
  # which is what keeps the NOEC at 3.
  p_three <- fit$comparisons$p_value[fit$comparisons$conc == 3]
  expect_lt(p_three, 0.05)
  expect_gt(p_three, 0.05 / 4)
})

test_that("wilcoxon_rank_sum does not warn about ties", {
  # Reproduction counts tie constantly, so the exact p-value is unavailable
  # and stats::wilcox.test warns. The warning is muffled deliberately.
  expect_no_warning(wilcoxon_rank_sum(appendix_f(), response = "young"))
})

test_that("wilcoxon_rank_sum accepts a balanced design", {
  fit <- wilcoxon_rank_sum(appendix_e(), response = "young")
  expect_equal(fit$comparisons$rank_sum, c(84, 63.5, 76, 55))
  expect_equal(fit$noec, 3)
})

# fisher_exact --------------------------------------------------------------

test_that("fisher_exact reproduces the Appendix G worked example", {
  fit <- fisher_exact(
    ceriodaphnia_g2,
    response = "dead",
    n_exposed = "exposed",
    type = "quantal"
  )

  # The manual concludes that only 25 per cent differs from the control,
  # giving NOEC 12 and LOEC 25.
  expect_equal(fit$comparisons$conc, c(1, 3, 6, 12, 25))
  expect_equal(
    fit$comparisons$significant,
    c(FALSE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_equal(fit$noec, 12)
  expect_equal(fit$loec, 25)
})

test_that("fisher_exact applies no multiplicity adjustment by default", {
  # Appendix G section 1 asks for a pairwise error rate, because the test is
  # itself conservative.
  fit <- fisher_exact(
    ceriodaphnia_g2,
    response = "dead",
    n_exposed = "exposed",
    type = "quantal"
  )
  expect_equal(fit$adjust, "none")
  expect_equal(fit$per_comparison_alpha, 0.05)

  adjusted <- fisher_exact(
    ceriodaphnia_g2,
    response = "dead",
    n_exposed = "exposed",
    type = "quantal",
    adjust = "bonferroni"
  )
  expect_equal(adjusted$per_comparison_alpha, 0.05 / 5)
  # The 25 per cent result survives the adjustment comfortably.
  expect_equal(adjusted$loec, 25)
})

test_that("fisher_exact requires quantal data", {
  expect_error(
    fisher_exact(fathead_c1, response = "weight"),
    regexp = "needs quantal data"
  )
})

test_that("fisher_exact uses the pooled counts, control included", {
  fit <- fisher_exact(
    ceriodaphnia_g2,
    response = "dead",
    n_exposed = "exposed",
    type = "quantal"
  )
  expect_equal(fit$comparisons$n_exposed, rep(10, 5))
  expect_equal(fit$comparisons$n_affected, c(0, 0, 0, 0, 10))
  expect_equal(fit$comparisons$proportion, c(0, 0, 0, 0, 1))
})

# printing ------------------------------------------------------------------

test_that("non-parametric results print their reference", {
  expect_output(print(steel(appendix_e(), response = "young")), "Appendix E")
  expect_output(
    print(wilcoxon_rank_sum(appendix_f(), response = "young")),
    "Appendix F"
  )
  expect_output(
    print(fisher_exact(
      ceriodaphnia_g2,
      response = "dead",
      n_exposed = "exposed",
      type = "quantal"
    )),
    "Appendix G"
  )
})
