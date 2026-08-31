#' Steel's Many-One Rank Test
#'
#' Compares each concentration with the control by ranking the two together,
#' holding the overall error rate at `alpha` across all the comparisons. This
#' is the test the EPA flowchart selects when the residuals are not normal, or
#' the variances are not homogeneous, and replication is equal.
#'
#' @details
#' For each control and concentration pair the observations are combined and
#' ranked, and the ranks belonging to the concentration are summed. A small
#' rank sum indicates a response below the control.
#'
#' The manual compares that rank sum with a tabulated critical value
#' (Table E.5). This package instead obtains a multiplicity-adjusted p-value
#' from [kSamples::Steel.test()], which evaluates the joint null distribution
#' of all `k` comparisons directly. That avoids reproducing a third-party
#' table, extends beyond the tabulated designs, and reproduces the manual's
#' conclusion on its own worked example, including the borderline comparison
#' whose rank sum equals the tabulated critical value exactly.
#'
#' Steel's test requires equal replication. The manual directs unequal
#' replication to [wilcoxon_rank_sum()], and this function signals an error
#' rather than returning a value the manual would not accept. It also requires
#' at least four replicates per concentration (section 9.4.5.2).
#'
#' @inheritParams toxstats_params
#' @param method How the p-value is obtained: `"asymptotic"` (the default),
#'   `"simulated"`, or `"exact"`. Passed to [kSamples::Steel.test()].
#' @param nsim Number of simulations when `method = "simulated"`.
#'
#' @inherit dunnett return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix E.
#'
#' Steel RGD (1959) A multiple comparison rank sum test: treatments versus
#' control. *Biometrics* 15:560-572.
#'
#' @examples
#' # Appendix E worked example. The 50 per cent concentration is excluded from
#' # the reproduction analysis because survival there was significantly
#' # reduced; see section 9.5.2 and the Appendix G example.
#' repro <- ceriodaphnia_e1[ceriodaphnia_e1$conc < 50, ]
#' steel(repro, response = "young")
#'
#' @export
steel <- function(
  x,
  ...,
  alpha = 0.05,
  method = c("asymptotic", "simulated", "exact"),
  nsim = 10000
) {
  method <- match.arg(method)
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))

  x <- as_tox_data(x, ...)
  parts <- comparison_parts(x)
  control <- parts$control_index
  treatments <- setdiff(seq_along(parts$conc), control)

  if (length(unique(parts$n)) != 1L) {
    chk::abort_chk(
      "Steel's Many-One Rank Test requires the same number of replicates at ",
      "every concentration; replication here is ",
      paste(parts$n, collapse = ", "),
      ". The EPA flowchart uses `wilcoxon_rank_sum()` when replication is ",
      "unequal."
    )
  }
  check_nonparametric_replication(parts$n)

  # kSamples expects the control first, which is also the order its per
  # comparison output follows.
  samples <- c(parts$values[control], parts$values[treatments])
  fit <- kSamples::Steel.test(
    samples,
    method = method,
    alternative = if (effect_sign(x$direction) > 0) "less" else "greater",
    Nsim = nsim
  )

  adjusted <- if (method == "asymptotic") {
    fit$pval.asympt.adj
  } else {
    fit$pval.adj
  }

  comparisons <- data.frame(
    conc = parts$conc[treatments],
    n = parts$n[treatments],
    mean = parts$mean[treatments],
    rank_sum = vapply(
      parts$values[treatments],
      function(v) rank_sum(parts$values[[control]], v),
      numeric(1)
    ),
    statistic = as.numeric(fit$Wstand),
    p_value = as.numeric(adjusted),
    stringsAsFactors = FALSE
  )
  comparisons$significant <- comparisons$p_value <= alpha
  rownames(comparisons) <- NULL

  new_comparison(
    test = "steel",
    method = paste0("Steel's Many-One Rank Test (", method, " p-values)"),
    comparisons = comparisons,
    reference = "EPA-821-R-02-013 Appendix E",
    alpha = alpha,
    n_ties = fit$n.ties
  )
}

#' Wilcoxon Rank Sum Test with Bonferroni's adjustment
#'
#' Compares each concentration with the control by ranking the two together,
#' dividing the significance level by the number of comparisons. This is the
#' test the EPA flowchart selects when the assumptions fail and replication is
#' unequal, which rules out Steel's test.
#'
#' @details
#' The manual compares each rank sum with a tabulated critical value
#' (Table F.5). This package computes a p-value with [stats::wilcox.test()] and
#' compares it with `alpha / k`, which avoids reproducing the table and
#' reproduces the manual's conclusion on its own worked example.
#'
#' Where the data contain ties an exact p-value is unavailable and the normal
#' approximation with a continuity correction is used, which
#' [stats::wilcox.test()] signals with a warning. Ties are common in
#' reproduction counts, so the warning is suppressed here and the number of
#' tied comparisons is recorded in the result instead.
#'
#' Bonferroni's adjustment bounds the overall error rate rather than fixing it,
#' so this test is more conservative than Steel's and correspondingly less
#' powerful. It requires at least four replicates per concentration
#' (section 9.4.5.2).
#'
#' @inheritParams toxstats_params
#'
#' @inherit dunnett return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix F.
#'
#' @examples
#' # Appendix F uses the Appendix E data with two males removed, one from the
#' # control and one from the 12 per cent concentration.
#' repro <- ceriodaphnia_e1[ceriodaphnia_e1$conc < 50, ]
#' unequal <- repro[!(repro$conc == 0 & repro$replicate == "1") &
#'   !(repro$conc == 12 & repro$replicate == "5"), ]
#' wilcoxon_rank_sum(unequal, response = "young")
#'
#' @export
wilcoxon_rank_sum <- function(x, ..., alpha = 0.05) {
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))

  x <- as_tox_data(x, ...)
  parts <- comparison_parts(x)
  control <- parts$control_index
  treatments <- setdiff(seq_along(parts$conc), control)
  check_nonparametric_replication(parts$n)

  sign <- effect_sign(x$direction)
  per_comparison <- alpha / length(treatments)

  fits <- lapply(parts$values[treatments], function(v) {
    withCallingHandlers(
      stats::wilcox.test(
        v,
        parts$values[[control]],
        alternative = if (sign > 0) "less" else "greater"
      ),
      warning = function(w) invokeRestart("muffleWarning")
    )
  })

  comparisons <- data.frame(
    conc = parts$conc[treatments],
    n = parts$n[treatments],
    mean = parts$mean[treatments],
    rank_sum = vapply(
      parts$values[treatments],
      function(v) rank_sum(parts$values[[control]], v),
      numeric(1)
    ),
    statistic = vapply(fits, function(f) unname(f$statistic), numeric(1)),
    p_value = vapply(fits, function(f) f$p.value, numeric(1)),
    stringsAsFactors = FALSE
  )
  comparisons$significant <- comparisons$p_value <= per_comparison
  rownames(comparisons) <- NULL

  new_comparison(
    test = "wilcoxon",
    method = "Wilcoxon Rank Sum Test with Bonferroni's adjustment",
    comparisons = comparisons,
    reference = "EPA-821-R-02-013 Appendix F",
    alpha = alpha,
    per_comparison_alpha = per_comparison
  )
}

#' Fisher's Exact Test
#'
#' Compares the proportion responding at each concentration with the control
#' proportion, using the hypergeometric distribution. The EPA manuals apply it
#' to survival, whose outcome is counted rather than measured.
#'
#' @details
#' Requires quantal data, so `tox_data(type = "quantal")` with `n_exposed`
#' supplied.
#'
#' **No multiplicity adjustment is applied by default.** Section 1 of Appendix
#' G states that because Fisher's Exact Test is itself conservative, "a
#' pair-wise comparison error rate of 0.05 is suggested rather than an
#' experiment-wise error rate". Set `adjust` to depart from that.
#'
#' The manual works the test by looking up a table of significance levels
#' (Table G.5). This package uses [stats::fisher.test()], which computes the
#' same probability directly.
#'
#' @inheritParams toxstats_params
#' @param adjust How to adjust for multiplicity. `"none"` is the EPA default
#'   for this test.
#'
#' @inherit dunnett return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix G.
#'
#' @examples
#' # Appendix G worked example; the manual reports NOEC 12 and LOEC 25
#' fisher_exact(
#'   ceriodaphnia_g2,
#'   response = "dead",
#'   n_exposed = "exposed",
#'   type = "quantal"
#' )
#'
#' @export
fisher_exact <- function(
  x,
  ...,
  alpha = 0.05,
  adjust = c("none", "bonferroni")
) {
  adjust <- match.arg(adjust)
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))

  x <- as_tox_data(x, ...)
  if (x$type != "quantal") {
    chk::abort_chk(
      "Fisher's Exact Test needs quantal data; call `tox_data()` with ",
      "`type = \"quantal\"` and an `n_exposed` column."
    )
  }

  pooled <- x$pooled
  control <- which(pooled$conc == x$control)
  treatments <- setdiff(seq_len(nrow(pooled)), control)
  per_comparison <- if (adjust == "bonferroni") {
    alpha / length(treatments)
  } else {
    alpha
  }

  sign <- effect_sign(x$direction)
  p_value <- vapply(
    treatments,
    function(i) {
      table_2x2 <- matrix(
        c(
          pooled$n_affected[control],
          pooled$n_exposed[control] - pooled$n_affected[control],
          pooled$n_affected[i],
          pooled$n_exposed[i] - pooled$n_affected[i]
        ),
        nrow = 2,
        byrow = TRUE
      )
      # `response` counts organisms affected, so a toxic effect at a
      # decreasing-response endpoint means more affected in the treatment row,
      # which is the "less" tail of the control-first table.
      stats::fisher.test(
        table_2x2,
        alternative = if (sign > 0) "less" else "greater"
      )$p.value
    },
    numeric(1)
  )

  comparisons <- data.frame(
    conc = pooled$conc[treatments],
    n_exposed = pooled$n_exposed[treatments],
    n_affected = pooled$n_affected[treatments],
    proportion = pooled$proportion[treatments],
    p_value = p_value,
    significant = p_value <= per_comparison,
    stringsAsFactors = FALSE
  )
  rownames(comparisons) <- NULL

  new_comparison(
    test = "fisher",
    method = paste0(
      "Fisher's Exact Test",
      if (adjust == "bonferroni") ", Bonferroni adjusted" else ""
    ),
    comparisons = comparisons,
    reference = "EPA-821-R-02-013 Appendix G",
    alpha = alpha,
    adjust = adjust,
    per_comparison_alpha = per_comparison
  )
}

#' Sum of the treatment ranks when control and treatment are ranked together
#'
#' @param control Control values.
#' @param treatment Treatment values.
#' @return A single number.
#' @noRd
rank_sum <- function(control, treatment) {
  combined <- rank(c(control, treatment))
  sum(combined[-seq_along(control)])
}

#' The non-parametric branch needs at least four replicates
#'
#' @param n Integer vector of replicate counts.
#' @return Invisibly TRUE, or an error.
#' @noRd
check_nonparametric_replication <- function(n) {
  if (any(n < 4)) {
    chk::abort_chk(
      "The non-parametric tests require at least four replicates at every ",
      "concentration (EPA-821-R-02-013 section 9.4.5.2); replication here ",
      "is ",
      paste(n, collapse = ", "),
      "."
    )
  }
  invisible(TRUE)
}
