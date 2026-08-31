#' The EPA multi-concentration hypothesis-testing flowchart
#'
#' The chart is held as data rather than as nested `if` statements, so that each
#' branch point can be tested on its own, the traversal can be recorded without
#' the recording being tangled into the logic, and a second chart can be added
#' without touching the engine.
#'
#' Chronic Figure 2 and acute Figure 13 are the same chart. Nodes whose name
#' begins with `@` are terminals naming the test to run.
#'
#' @return A data frame with columns `node`, `question`, `rule`, `yes`, `no`
#'   and `reference`.
#' @noRd
flowchart_hypothesis_multi <- function() {
  data.frame(
    node = c(
      "transform",
      "normality",
      "variance",
      "min_replicates",
      "replication_p",
      "replication_np"
    ),
    question = c(
      "Is the response a proportion requiring transformation?",
      "Are the pooled within-group residuals normally distributed?",
      "Are the variances homogeneous across concentrations?",
      "Are there at least four replicates at every concentration?",
      "Is replication equal across all concentrations?",
      "Is replication equal across all concentrations?"
    ),
    rule = c(
      "transform",
      "normality",
      "variance",
      "min_replicates",
      "balanced",
      "balanced"
    ),
    # The transform node leads to the normality test whichever way it is
    # answered; it records what was done to the response, not a branch.
    yes = c(
      "normality",
      "variance",
      "replication_p",
      "replication_np",
      "@dunnett",
      "@steel"
    ),
    no = c(
      "normality",
      "min_replicates",
      "min_replicates",
      "@no_valid_test",
      "@bonferroni_t",
      "@wilcoxon"
    ),
    reference = c(
      "EPA-821-R-02-013 Appendix B, section 4.2",
      "EPA-821-R-02-013 Appendix B, section 2.1",
      "EPA-821-R-02-013 Appendix B, section 3",
      "EPA-821-R-02-013 section 9.4.5.2",
      "EPA-821-R-02-013 Figure 2",
      "EPA-821-R-02-013 Figure 2"
    ),
    stringsAsFactors = FALSE
  )
}

#' Walk a flowchart, recording each branch point
#'
#' @param chart A chart data frame.
#' @param state A list carrying the data and anything the rules add to it.
#' @return A list with `decisions` (a data frame), `terminal` (the `@` node
#'   reached) and the final `state`.
#' @noRd
walk_flowchart <- function(chart, state) {
  node <- chart$node[1]
  rows <- list()
  step <- 0L

  while (!startsWith(node, "@")) {
    definition <- chart[chart$node == node, , drop = FALSE]
    if (!nrow(definition)) {
      chk::abort_chk("Unknown flowchart node: ", node, ".")
    }

    result <- apply_rule(definition$rule, state)
    state <- result$state
    step <- step + 1L

    rows[[step]] <- data.frame(
      step = step,
      node = node,
      question = definition$question,
      criterion = result$criterion,
      statistic = result$statistic,
      p_value = result$p_value,
      answer = result$answer,
      outcome = result$outcome,
      reference = definition$reference,
      stringsAsFactors = FALSE
    )

    node <- if (result$answer == "yes") definition$yes else definition$no
  }

  list(
    decisions = do.call(rbind, rows),
    terminal = node,
    state = state
  )
}

#' Evaluate one branch point
#'
#' @param rule The rule name from the chart.
#' @param state The traversal state.
#' @return A list with `answer`, `criterion`, `statistic`, `p_value`, `outcome`
#'   and the updated `state`.
#' @noRd
apply_rule <- function(rule, state) {
  switch(
    rule,
    transform = rule_transform(state),
    normality = rule_normality(state),
    variance = rule_variance(state),
    min_replicates = rule_min_replicates(state),
    balanced = rule_balanced(state),
    partial_two = rule_partial_two(state),
    partial_one = rule_partial_one(state),
    probit_fits = rule_probit_fits(state),
    sk_requirements = rule_sk_requirements(state),
    chk::abort_chk("Unknown flowchart rule: ", rule, ".")
  )
}

#' @noRd
rule_transform <- function(state) {
  quantal <- state$working$type == "quantal"

  if (quantal) {
    transformed <- state$working
    transformed$replicates$response <- arcsine_sqrt(
      state$working$replicates$proportion,
      n = state$working$replicates$n_exposed
    )
    # Downstream the transformed values are an ordinary measurement, so the
    # type is changed to match. The untransformed object is kept in
    # `state$data` for back-transforming the minimum significant difference.
    transformed$type <- "continuous"
    transformed$pooled <- pool_replicates(
      transformed$replicates,
      type = "continuous"
    )
    state$working <- transformed
    state$transform <- "arcsine_sqrt"
  } else {
    state$transform <- "none"
  }

  new_rule_result(
    answer = if (quantal) "yes" else "no",
    criterion = "A quantal endpoint is a proportion",
    outcome = if (quantal) {
      "arc sine square root transform applied"
    } else {
      "no transformation needed"
    },
    state = state
  )
}

#' @noRd
rule_normality <- function(state) {
  fit <- epa_normality(
    state$working,
    alpha_assumption = state$alpha_assumption
  )
  state$assumptions$normality <- fit

  new_rule_result(
    answer = if (fit$normal) "yes" else "no",
    criterion = paste0(
      fit$statistic_name,
      ", reject if ",
      if (is.na(fit$p_value)) {
        paste0("above the critical value ", signif(fit$critical, 4))
      } else {
        paste0("p <= ", state$alpha_assumption)
      }
    ),
    statistic = fit$statistic,
    p_value = fit$p_value,
    outcome = if (fit$normal) {
      "residuals consistent with normality"
    } else {
      "residuals not normally distributed"
    },
    state = state
  )
}

#' @noRd
rule_variance <- function(state) {
  fit <- epa_variance(
    state$working,
    alpha_assumption = state$alpha_assumption
  )
  state$assumptions$variance <- fit

  new_rule_result(
    answer = if (fit$homogeneous) "yes" else "no",
    criterion = paste0("Bartlett B, reject if p <= ", state$alpha_assumption),
    statistic = fit$statistic,
    p_value = fit$p_value,
    outcome = if (fit$homogeneous) {
      "variances not significantly different"
    } else {
      "variances heterogeneous"
    },
    state = state
  )
}

#' @noRd
rule_min_replicates <- function(state) {
  n <- state$working$pooled$n_rep
  enough <- all(n >= 4)

  new_rule_result(
    answer = if (enough) "yes" else "no",
    criterion = "at least 4 replicates at every concentration",
    statistic = min(n),
    outcome = if (enough) {
      paste0("smallest concentration has ", min(n), " replicates")
    } else {
      paste0(
        "only ",
        min(n),
        " replicates; no EPA-sanctioned test is available"
      )
    },
    state = state
  )
}

#' @noRd
rule_balanced <- function(state) {
  n <- state$working$pooled$n_rep
  balanced <- length(unique(n)) == 1L

  new_rule_result(
    answer = if (balanced) "yes" else "no",
    criterion = "equal replication at every concentration",
    outcome = paste0(
      paste(n, collapse = ", "),
      if (balanced) " replicates, balanced" else " replicates, unbalanced"
    ),
    state = state
  )
}

#' @noRd
new_rule_result <- function(
  answer,
  criterion,
  outcome,
  state,
  statistic = NA_real_,
  p_value = NA_real_
) {
  list(
    answer = answer,
    criterion = criterion,
    statistic = as.numeric(statistic),
    p_value = as.numeric(p_value),
    outcome = outcome,
    state = state
  )
}

#' The test each terminal names
#'
#' @return A named character vector mapping terminal to function name.
#' @noRd
flowchart_terminals <- function() {
  c(
    "@dunnett" = "dunnett",
    "@bonferroni_t" = "bonferroni_t",
    "@steel" = "steel",
    "@wilcoxon" = "wilcoxon_rank_sum"
  )
}

#' Extract the decision trail
#'
#' Returns the record of how a test was selected: every branch point of the EPA
#' flowchart, the statistic that answered it, the consequence, and the section
#' of the manual that justifies it.
#'
#' @param x An object carrying a decision trail, such as the result of
#'   [toxcalc()].
#' @param ... Unused.
#'
#' @return A data frame with columns `step`, `node`, `question`, `criterion`,
#'   `statistic`, `p_value`, `answer`, `outcome` and `reference`.
#'
#' @examples
#' decisions(toxcalc(fathead_c1, response = "weight"))
#'
#' @export
decisions <- function(x, ...) {
  UseMethod("decisions")
}

#' @describeIn decisions Decision trail from a `toxcalc` analysis.
#' @export
decisions.toxcalc <- function(x, ...) {
  chk::chk_unused(...)
  x$decisions
}
