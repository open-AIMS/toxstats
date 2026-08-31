#' Dunnett's procedure
#'
#' Compares each concentration mean with the control mean, holding the overall
#' error rate at `alpha` across all the comparisons. This is the parametric
#' test the EPA flowchart selects when the residuals are normal, the variances
#' are homogeneous, and every concentration has the same number of replicates.
#'
#' @details
#' The statistic for concentration `i` is
#'
#' \deqn{t_i = \frac{\bar{Y}_1 - \bar{Y}_i}{S_w \sqrt{1/n_1 + 1/n_i}}}
#'
#' with `Sw` the square root of the within mean square from a one-way analysis
#' of variance. It is compared with a one-sided critical value computed as
#' described in [msd()]; the manual reads the same value from a table credited
#' to Miller (1981).
#'
#' Dunnett's procedure requires equal replication. When replication is unequal
#' the manual directs the analysis to [bonferroni_t()] instead, and this
#' function signals an error rather than silently returning a value the manual
#' would not accept.
#'
#' @inheritParams toxstats_params
#'
#' @return An object of class `tox_comparison`, a list with elements
#'   `test`, `method`, `comparisons` (a data frame with one row per treatment
#'   concentration), `noec`, `loec`, `monotone`, `critical`, `sw`, `df`,
#'   `alpha` and `reference`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix C.
#'
#' Dunnett CW (1955) A multiple comparison procedure for comparing several
#' treatments with a control. *Journal of the American Statistical
#' Association* 50:1096-1121.
#'
#' @examples
#' # Appendix C worked example; the manual reports NOEC 128 and LOEC 256
#' dunnett(fathead_c1, response = "weight")
#'
#' @export
dunnett <- function(x, ..., alpha = 0.05) {
  x <- as_tox_data(x, ...)
  parts <- comparison_parts(x)

  if (length(unique(parts$n)) != 1L) {
    chk::abort_chk(
      "Dunnett's procedure requires the same number of replicates at every ",
      "concentration; replication here is ",
      paste(parts$n, collapse = ", "),
      ". The EPA flowchart uses `bonferroni_t()` when replication is unequal."
    )
  }

  many_one_t(x, parts, alpha, test = "dunnett")
}

#' The t test with Bonferroni's adjustment
#'
#' Compares each concentration mean with the control mean, dividing the
#' significance level by the number of comparisons. This is the parametric test
#' the EPA flowchart selects when the assumptions are met but replication is
#' unequal, which rules out Dunnett's procedure.
#'
#' @details
#' The statistic is the same as in [dunnett()]; only the critical value
#' differs, being `qt(1 - alpha / k, df)` for `k` treatment concentrations.
#' This reproduces the manual's Table D.5 exactly: 2.510 for `alpha = 0.05`
#' with 14 degrees of freedom and four concentrations, and 2.241 in the row
#' for infinite degrees of freedom against a printed 2.242.
#'
#' Bonferroni's adjustment sets an upper bound of `alpha` on the overall error
#' rate rather than fixing it there, so it is more conservative than Dunnett's
#' procedure and correspondingly less powerful.
#'
#' @inheritParams toxstats_params
#'
#' @inherit dunnett return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix D.
#'
#' @examples
#' # Appendix D uses the Appendix C data with one replicate lost
#' lost <- fathead_c1[-19, ]
#' bonferroni_t(lost, response = "weight")
#'
#' @export
bonferroni_t <- function(x, ..., alpha = 0.05) {
  x <- as_tox_data(x, ...)
  many_one_t(x, comparison_parts(x), alpha, test = "bonferroni")
}

#' The t test with Dunn-Sidak's adjustment
#'
#' As [bonferroni_t()], but using Dunn-Sidak's adjustment, which is exact when
#' the comparisons are independent and so slightly less conservative than
#' Bonferroni's.
#'
#' @details
#' The critical value is `qt(1 - (1 - (1 - alpha)^(1 / k)), df)`. ToxCalc
#' offered this test; the EPA flowchart does not select it, so using it is a
#' documented departure from the manual.
#'
#' @inheritParams toxstats_params
#'
#' @inherit dunnett return
#'
#' @references
#' Sidak Z (1967) Rectangular confidence regions for the means of multivariate
#' normal distributions. *Journal of the American Statistical Association*
#' 62:626-633.
#'
#' @examples
#' dunn_sidak_t(fathead_c1, response = "weight")
#'
#' @export
dunn_sidak_t <- function(x, ..., alpha = 0.05) {
  x <- as_tox_data(x, ...)
  many_one_t(x, comparison_parts(x), alpha, test = "sidak")
}

#' Welch's t test against the control
#'
#' Compares each concentration mean with the control mean without pooling the
#' variance, so it remains valid when the variances differ. Each comparison
#' carries its own degrees of freedom, from the Welch-Satterthwaite
#' approximation.
#'
#' @details
#' This is not on the EPA multi-concentration flowchart, which sends
#' heterogeneous variances to the non-parametric branch. It appears in the
#' single-concentration chart of the acute manual (Figure 12) and was offered
#' by ToxCalc, so it is provided here and flagged as a departure.
#'
#' @inheritParams toxstats_params
#' @param adjust How to adjust for multiplicity across the `k` comparisons.
#'
#' @inherit dunnett return
#'
#' @references
#' US EPA (2002) EPA-821-R-02-012, Figure 12.
#'
#' @examples
#' welch_t(fathead_c1, response = "weight")
#'
#' @export
welch_t <- function(
  x,
  ...,
  alpha = 0.05,
  adjust = c("bonferroni", "sidak", "none")
) {
  adjust <- match.arg(adjust)
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))

  x <- as_tox_data(x, ...)
  parts <- comparison_parts(x)
  control <- parts$control_index
  treatments <- setdiff(seq_along(parts$conc), control)
  sign <- effect_sign(x$direction)

  per_comparison <- switch(
    adjust,
    bonferroni = alpha / length(treatments),
    sidak = 1 - (1 - alpha)^(1 / length(treatments)),
    none = alpha
  )

  fits <- lapply(treatments, function(i) {
    stats::t.test(
      parts$values[[control]],
      parts$values[[i]],
      alternative = if (sign > 0) "greater" else "less",
      var.equal = FALSE
    )
  })

  comparisons <- data.frame(
    conc = parts$conc[treatments],
    n = parts$n[treatments],
    mean = parts$mean[treatments],
    statistic = vapply(fits, function(f) unname(f$statistic), numeric(1)),
    df = vapply(fits, function(f) unname(f$parameter), numeric(1)),
    p_value = vapply(fits, function(f) f$p.value, numeric(1)),
    stringsAsFactors = FALSE
  )
  comparisons$critical <- stats::qt(1 - per_comparison, comparisons$df)
  comparisons$significant <- comparisons$p_value <= per_comparison

  new_comparison(
    test = "welch_t",
    method = paste0(
      "Welch's t test against the control, ",
      adjust,
      " adjusted (not an EPA multi-concentration method)"
    ),
    comparisons = comparisons,
    reference = "EPA-821-R-02-012 Figure 12; not on the Figure 2 flowchart",
    alpha = alpha,
    adjust = adjust,
    per_comparison_alpha = per_comparison
  )
}

#' Many-to-one comparison against a pooled variance estimate
#'
#' Shared by Dunnett's procedure and the two adjusted t tests, which differ
#' only in the critical value they compare the statistic with.
#'
#' @param x A `tox_data` object.
#' @param parts The output of `comparison_parts()`.
#' @param alpha Familywise significance level.
#' @param test One of "dunnett", "bonferroni", "sidak".
#' @return A `tox_comparison`.
#' @noRd
many_one_t <- function(x, parts, alpha, test) {
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))

  control <- parts$control_index
  treatments <- setdiff(seq_along(parts$conc), control)
  if (!length(treatments)) {
    chk::abort_chk("At least one non-control concentration is needed.")
  }
  sign <- effect_sign(x$direction)

  critical <- critical_value(
    test = test,
    alpha = alpha,
    df = parts$df,
    n_control = parts$n[control],
    n_treat = parts$n[treatments]
  )

  statistic <- sign *
    (parts$mean[control] - parts$mean[treatments]) /
    (parts$sw * sqrt(1 / parts$n[control] + 1 / parts$n[treatments]))

  # The reported p-value is adjusted for multiplicity, so that comparing it
  # with alpha gives the same answer as comparing the statistic with the
  # critical value. Reporting the unadjusted p-value beside a familywise
  # critical value would invite the two to be read together and disagree.
  raw <- stats::pt(statistic, parts$df, lower.tail = FALSE)
  k <- length(treatments)
  adjusted <- switch(
    test,
    dunnett = vapply(
      statistic,
      function(t) {
        1 -
          dunnett_prob(
            t,
            sqrt(
              parts$n[treatments] / (parts$n[control] + parts$n[treatments])
            ),
            parts$df
          )
      },
      numeric(1)
    ),
    bonferroni = pmin(1, k * raw),
    sidak = 1 - (1 - raw)^k
  )

  comparisons <- data.frame(
    conc = parts$conc[treatments],
    n = parts$n[treatments],
    mean = parts$mean[treatments],
    statistic = statistic,
    critical = critical,
    p_value = adjusted,
    significant = statistic > critical,
    stringsAsFactors = FALSE
  )

  labels <- c(
    dunnett = "Dunnett's procedure",
    bonferroni = "t test with Bonferroni's adjustment",
    sidak = "t test with Dunn-Sidak's adjustment"
  )
  references <- c(
    dunnett = "EPA-821-R-02-013 Appendix C",
    bonferroni = "EPA-821-R-02-013 Appendix D",
    sidak = "Not on the EPA flowchart; offered by ToxCalc"
  )

  new_comparison(
    test = test,
    method = labels[[test]],
    comparisons = comparisons,
    reference = references[[test]],
    critical = critical,
    sw = parts$sw,
    df = parts$df,
    alpha = alpha
  )
}
