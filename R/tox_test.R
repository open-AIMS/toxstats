#' Run the EPA hypothesis-testing analysis
#'
#' Walks the EPA flowchart, selects the hypothesis test the manual directs for
#' the data, runs it, and returns the endpoints together with a record of how
#' the test was chosen. This is the function that recreates what ToxCalc did.
#'
#' @details
#' The flowchart is Figure 2 of the chronic manual, which is also Figure 13 of
#' the acute manual. It branches on four things in order: whether the response
#' is a proportion, whether the pooled within-group residuals are normal,
#' whether the variances are homogeneous, and whether replication is equal. The
#' four terminals are [dunnett()], [bonferroni_t()], [steel()] and
#' [wilcoxon_rank_sum()].
#'
#' Two significance levels are involved and they differ. The assumption tests
#' use `alpha_assumption`, which section 9.4.6.1 sets at 0.01; the hypothesis
#' test uses `alpha`, set at 0.05 and one-sided.
#'
#' ## What is recorded
#'
#' The value of this function over calling the tests by hand is the decision
#' trail, returned by [decisions()]. Each row is a branch point: the question,
#' the criterion, the statistic, the answer, the consequence, and the manual
#' section that justifies it. `summary()` prints it above the endpoints.
#'
#' ## Exclusions
#'
#' Section 9.5.2 excludes some concentrations from the hypothesis test while
#' retaining them for point estimation. Two rules apply.
#'
#' A concentration at which every replicate showed a complete response, for
#' instance total mortality, is excluded automatically for a quantal endpoint.
#'
#' A concentration at which survival was significantly reduced is excluded from
#' the analysis of a **sublethal** endpoint such as growth or reproduction.
#' That requires a separate survival analysis, so it cannot be inferred here;
#' pass the concentrations to `exclude`. The Appendix E example does exactly
#' this, dropping the 50 per cent concentration from the reproduction analysis
#' on the strength of the Appendix G survival result.
#'
#' ## The percent minimum significant difference override
#'
#' Section 10.2.8.2.5 states that a concentration shall not be declared toxic
#' if its relative difference from the control is less than the lower bound
#' published for the test method. This overrides the hypothesis test. It is
#' applied only when `pmsd_bounds` is supplied, and each concentration it
#' reverses is added to the decision trail.
#'
#' ## Forcing a test
#'
#' Passing `test` runs that test whatever the flowchart selected. When the two
#' disagree a warning is signalled, `overridden` is set, and a row recording
#' the override is added to the trail, so a forced analysis can never be
#' printed as though the manual had chosen it.
#'
#' @inheritParams toxstats_params
#' @param test Optional. Name of a test to run instead of the one the
#'   flowchart selects: one of `"dunnett"`, `"bonferroni_t"`, `"steel"` or
#'   `"wilcoxon_rank_sum"`.
#' @param exclude Optional numeric vector of concentrations to exclude from the
#'   hypothesis test, over and above those excluded automatically.
#' @param pmsd_bounds Optional. Passed to [pmsd()]; supplying it also enables
#'   the section 10.2.8.2.5 override.
#' @param branch `"hypothesis"` (the default) gives the NOEC, LOEC, MSD and
#'   PMSD. `"both"` adds a point estimate alongside them, which is what a
#'   laboratory report usually contains: a median lethal concentration through
#'   [lc50()] for a quantal endpoint, or an inhibition concentration through
#'   [icp()] for a continuous one.
#'
#'   There is no `"point"` option, because [lc50()] and [icp()] already provide
#'   that directly and are easier to find. The default is not `"both"` because
#'   the point branch runs a bootstrap for a continuous endpoint, which takes
#'   appreciably longer, so it is asked for rather than assumed.
#' @param p Percentages for the point branch, or `NULL` for the default of
#'   whichever method applies.
#' @param nboot,seed Passed to [icp()] when the point branch applies to a
#'   continuous endpoint.
#'
#' @return An object of class `tox_test`, a list with elements `data`,
#'   `working`, `decisions`, `assumptions`, `comparison`, `msd`, `pmsd`,
#'   `noec`, `loec`, `excluded`, `overridden`, `transform`, `flowchart`,
#'   `alpha` and `alpha_assumption`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Figure 2 and sections 9.4 to 9.6.
#'
#' @seealso [decisions()] for the trail, [msd()] and [pmsd()] for sensitivity.
#'
#' @examples
#' # Appendix C growth data: normal, homogeneous and balanced, so the
#' # flowchart selects Dunnett's procedure.
#' fit <- tox_test(fathead_c1, response = "weight")
#' summary(fit)
#'
#' @export
tox_test <- function(
  x,
  ...,
  alpha = 0.05,
  alpha_assumption = 0.01,
  test = NULL,
  exclude = NULL,
  pmsd_bounds = NULL,
  branch = c("hypothesis", "both"),
  p = NULL,
  nboot = 200,
  seed = NULL
) {
  branch <- match.arg(branch)
  chk::chk_number(alpha)
  chk::chk_range(alpha, c(0, 1))
  chk::chk_number(alpha_assumption)
  chk::chk_range(alpha_assumption, c(0, 1))
  chk::chk_null_or(test, vld = chk::vld_string)
  chk::chk_null_or(exclude, vld = chk::vld_numeric)
  if (!is.null(test)) {
    # Williams' test is not a flowchart terminal, but it can be forced, which
    # is how ToxCalc's one non-EPA test is reached through the driver.
    chk::chk_subset(test, c(unname(flowchart_terminals()), "williams"))
  }

  x <- as_tox_data(x, ...)
  excluded <- excluded_concentrations(x, exclude)
  analysis <- drop_concentrations(x, excluded$conc)

  walked <- walk_flowchart(
    flowchart_hypothesis_multi(),
    list(
      data = analysis,
      working = analysis,
      alpha = alpha,
      alpha_assumption = alpha_assumption,
      assumptions = list()
    )
  )

  if (walked$terminal == "@no_valid_test") {
    chk::abort_chk(
      "The residuals are not normally distributed, or the variances are not ",
      "homogeneous, and at least one concentration has fewer than four ",
      "replicates. Section 9.4.5.2 rules out the non-parametric tests, so no ",
      "EPA-sanctioned test is available for this design."
    )
  }

  selected <- unname(flowchart_terminals()[walked$terminal])
  chosen <- selected
  overridden <- FALSE
  trail <- walked$decisions

  if (!is.null(test) && test != selected) {
    overridden <- TRUE
    chosen <- test
    warning(
      "The flowchart selected `",
      selected,
      "` but `test = \"",
      test,
      "\"` was given. The result is not the analysis the EPA manual ",
      "prescribes for these data.",
      call. = FALSE
    )
    trail <- rbind(
      trail,
      data.frame(
        step = nrow(trail) + 1L,
        node = "override",
        question = "Was the selected test overridden?",
        criterion = "user supplied `test`",
        statistic = NA_real_,
        p_value = NA_real_,
        answer = "yes",
        outcome = paste0(
          "user forced ",
          test,
          "; the flowchart selected ",
          selected
        ),
        reference = "Not an EPA decision",
        stringsAsFactors = FALSE
      )
    )
  }

  comparison <- do.call(chosen, list(walked$state$working, alpha = alpha))

  sensitivity <- pmsd_or_null(walked$state$working, alpha, pmsd_bounds)
  override <- apply_pmsd_override(comparison, walked$state$working, sensitivity)
  comparison <- override$comparison
  if (nrow(override$rows)) {
    override$rows$step <- nrow(trail) + seq_len(nrow(override$rows))
    trail <- rbind(trail, override$rows)
  }

  trail <- rbind(
    trail,
    data.frame(
      step = nrow(trail) + 1L,
      node = "test",
      question = "Which test was run?",
      criterion = NA_character_,
      statistic = NA_real_,
      p_value = NA_real_,
      answer = chosen,
      outcome = comparison$method,
      reference = comparison$reference,
      stringsAsFactors = FALSE
    )
  )
  rownames(trail) <- NULL

  point <- if (branch == "hypothesis") {
    NULL
  } else {
    point_branch(x, p = p, alpha = alpha, nboot = nboot, seed = seed)
  }

  structure(
    list(
      data = x,
      working = walked$state$working,
      point = point,
      branch = branch,
      decisions = trail,
      assumptions = walked$state$assumptions,
      comparison = comparison,
      msd = if (is.null(sensitivity)) NULL else sensitivity$msd,
      pmsd = sensitivity,
      noec = comparison$noec,
      loec = comparison$loec,
      excluded = excluded,
      overridden = overridden,
      selected = selected,
      transform = walked$state$transform,
      flowchart = "EPA-821-R-02-013 Figure 2",
      alpha = alpha,
      alpha_assumption = alpha_assumption,
      call = match.call()
    ),
    class = "tox_test"
  )
}

#' Concentrations excluded from the hypothesis test
#'
#' @param x A `tox_data` object.
#' @param exclude User-supplied concentrations to exclude.
#' @return A data frame with columns `conc` and `reason`.
#' @noRd
excluded_concentrations <- function(x, exclude) {
  rows <- list()

  if (x$type == "quantal") {
    complete <- x$pooled$conc[x$pooled$proportion == 1]
    complete <- setdiff(complete, x$control)
    if (length(complete)) {
      rows[[length(rows) + 1L]] <- data.frame(
        conc = complete,
        reason = "complete response in every replicate (section 9.5.2)",
        stringsAsFactors = FALSE
      )
    }
  }

  if (!is.null(exclude) && length(exclude)) {
    chk::chk_subset(exclude, x$pooled$conc)
    rows[[length(rows) + 1L]] <- data.frame(
      conc = exclude,
      reason = "excluded by the caller (section 9.5.2)",
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) {
    return(data.frame(
      conc = numeric(0),
      reason = character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  out[!duplicated(out$conc), , drop = FALSE]
}

#' Remove concentrations from a tox_data object
#'
#' @param x A `tox_data` object.
#' @param conc Concentrations to drop.
#' @return A `tox_data` object.
#' @noRd
drop_concentrations <- function(x, conc) {
  if (!length(conc)) {
    return(x)
  }
  if (x$control %in% conc) {
    chk::abort_chk("The control concentration cannot be excluded.")
  }

  keep <- !x$replicates$conc %in% conc
  x$replicates <- x$replicates[keep, , drop = FALSE]
  rownames(x$replicates) <- NULL
  x$pooled <- pool_replicates(x$replicates, type = x$type)
  x
}

#' Sensitivity measures, or NULL when no bounds were requested
#'
#' @noRd
pmsd_or_null <- function(x, alpha, bounds) {
  # The minimum significant difference is defined parametrically, and section
  # 10.2.8 asks for it even when the flowchart chose a non-parametric test.
  test <- if (length(unique(x$pooled$n_rep)) == 1L) "dunnett" else "bonferroni"
  pmsd(x, alpha = alpha, test = test, bounds = bounds)
}

#' Apply the section 10.2.8.2.5 lower-bound override
#'
#' @param comparison A `tox_comparison`.
#' @param x The working `tox_data`.
#' @param sensitivity The `tox_pmsd` result, or NULL.
#' @return A list with the possibly amended `comparison` and the decision
#'   `rows` to append.
#' @noRd
apply_pmsd_override <- function(comparison, x, sensitivity) {
  empty <- decision_rows(character(0))
  if (is.null(sensitivity) || is.null(sensitivity$bounds)) {
    return(list(comparison = comparison, rows = empty))
  }

  control_mean <- x$pooled$mean[x$pooled$conc == x$control]
  relative <- 100 *
    abs(control_mean - comparison$comparisons$mean) /
    abs(control_mean)
  reverse <- comparison$comparisons$significant &
    relative < sensitivity$bounds[1]

  if (!any(reverse)) {
    return(list(comparison = comparison, rows = empty))
  }

  comparison$comparisons$significant[reverse] <- FALSE
  comparison$comparisons$pmsd_override <- reverse
  endpoints <- derive_noec_loec(
    comparison$comparisons$conc,
    comparison$comparisons$significant
  )
  comparison$noec <- endpoints$noec
  comparison$loec <- endpoints$loec
  comparison$monotone <- endpoints$monotone

  rows <- data.frame(
    step = NA_integer_,
    node = "pmsd_override",
    question = paste0(
      "Is the relative difference at ",
      comparison$comparisons$conc[reverse],
      " below the lower PMSD bound?"
    ),
    criterion = paste0(
      "relative difference < ",
      sensitivity$bounds[1],
      " per cent"
    ),
    statistic = relative[reverse],
    p_value = NA_real_,
    answer = "yes",
    outcome = "significance reversed; the concentration is not declared toxic",
    reference = "EPA-821-R-02-013 section 10.2.8.2.5",
    stringsAsFactors = FALSE
  )

  list(comparison = comparison, rows = rows)
}

#' An empty decision trail with the right columns
#'
#' @noRd
decision_rows <- function(node) {
  data.frame(
    step = integer(0),
    node = node,
    question = character(0),
    criterion = character(0),
    statistic = numeric(0),
    p_value = numeric(0),
    answer = character(0),
    outcome = character(0),
    reference = character(0),
    stringsAsFactors = FALSE
  )
}

#' @export
print.tox_test <- function(x, ...) {
  cat("EPA WET hypothesis test\n")
  cat("  Flowchart: ", x$flowchart, "\n", sep = "")
  cat("  Selected:  ", x$comparison$method, "\n", sep = "")
  if (x$overridden) {
    cat("  OVERRIDDEN: the flowchart selected ", x$selected, "\n", sep = "")
  }
  cat("\n")
  print_endpoints(x)
  cat("\nUse summary() for the decision trail.\n")
  invisible(x)
}

#' @describeIn tox_test Print the decision trail above the endpoints.
#'
#' @param object A `tox_test` object.
#'
#' @export
summary.tox_test <- function(object, ...) {
  chk::chk_unused(...)

  cat("EPA WET hypothesis test\n")
  cat("Flowchart: ", object$flowchart, "\n\n", sep = "")

  if (nrow(object$excluded)) {
    cat("Excluded from the hypothesis test:\n")
    for (i in seq_len(nrow(object$excluded))) {
      cat(
        "  ",
        object$excluded$conc[i],
        " -- ",
        object$excluded$reason[i],
        "\n",
        sep = ""
      )
    }
    cat("\n")
  }

  trail <- object$decisions
  for (i in seq_len(nrow(trail))) {
    cat(
      "  ",
      formatC(trail$step[i], width = 2),
      "  ",
      trail$question[i],
      "\n",
      sep = ""
    )
    detail <- trail$outcome[i]
    if (!is.na(trail$statistic[i])) {
      detail <- paste0(
        signif(trail$statistic[i], 4),
        if (!is.na(trail$p_value[i])) {
          paste0(", p = ", signif(trail$p_value[i], 3))
        } else {
          ""
        },
        " -> ",
        detail
      )
    }
    cat("      ", detail, "\n", sep = "")
    cat("      (", trail$reference[i], ")\n", sep = "")
  }

  cat("\n")
  print_endpoints(object)
  invisible(object)
}

#' @noRd
print_endpoints <- function(x) {
  conc <- x$comparison$comparisons$conc
  cat(
    "  NOEC ",
    format_endpoint(x$noec, conc, "noec"),
    "    LOEC ",
    format_endpoint(x$loec, conc, "loec"),
    "\n",
    sep = ""
  )
  if (!is.null(x$pmsd)) {
    cat(
      "  MSD  ",
      signif(x$pmsd$msd$msd[1], 4),
      "    PMSD ",
      round(x$pmsd$pmsd[1], 1),
      " per cent",
      sep = ""
    )
    if (!is.null(x$pmsd$bounds)) {
      cat(
        " (EPA bounds ",
        x$pmsd$bounds[1],
        " to ",
        x$pmsd$bounds[2],
        ": ",
        x$pmsd$status[1],
        ")",
        sep = ""
      )
    }
    cat("\n")
  }
  point <- point_estimates(x)
  if (!is.null(point)) {
    cat(
      "
"
    )
    out <- point[, c("endpoint", "estimate", "lower", "upper")]
    out[] <- lapply(out, function(column) {
      if (is.numeric(column)) signif(column, 5) else column
    })
    print(out, row.names = FALSE)
  }
  if (!x$comparison$monotone) {
    cat(
      "
  Note: the pattern of significance is not monotone. Section 9.6.5.1
",
      "  advises that such results be used with caution.
",
      sep = ""
    )
  }
  invisible(x)
}

#' @describeIn tox_test Return every endpoint as one tidy row.
#'
#' @param row.names Unused.
#' @param optional Unused.
#'
#' @export
as.data.frame.tox_test <- function(x, row.names = NULL, optional = FALSE, ...) {
  chk::chk_unused(...)

  hypothesis <- data.frame(
    endpoint = c("NOEC", "LOEC", "MSD", "PMSD"),
    value = c(
      x$noec,
      x$loec,
      if (is.null(x$pmsd)) NA_real_ else unname(x$pmsd$msd$msd[1]),
      if (is.null(x$pmsd)) NA_real_ else unname(x$pmsd$pmsd[1])
    ),
    lower = NA_real_,
    upper = NA_real_,
    method = x$comparison$test,
    stringsAsFactors = FALSE
  )

  point <- point_estimates(x)
  if (!is.null(point)) {
    hypothesis <- rbind(
      hypothesis,
      data.frame(
        endpoint = point$endpoint,
        value = point$estimate,
        lower = point$lower,
        upper = point$upper,
        method = if (inherits(x$point, "tox_lc50")) {
          x$point$selected
        } else {
          "icp"
        },
        stringsAsFactors = FALSE
      )
    )
  }

  hypothesis$selected_by_flowchart <- x$selected
  hypothesis$overridden <- x$overridden
  hypothesis$transform <- x$transform
  hypothesis$conc_units <- x$data$conc_units
  rownames(hypothesis) <- NULL
  hypothesis
}

#' Run the point-estimation branch appropriate to the endpoint
#'
#' A quantal endpoint gets a median lethal concentration from the Figure 6
#' chart; a continuous one gets an inhibition concentration by linear
#' interpolation. Either way the **original** data is used, with nothing
#' excluded: section 9.5.2 retains for point estimation the concentrations it
#' drops from the no-observed-effect concentration.
#'
#' @param x The unmodified `tox_data` object.
#' @param p Percentages to estimate, or NULL for the default of each method.
#' @param alpha Significance level, used only by the probit fit check.
#' @param nboot,seed Passed to [icp()] for a continuous endpoint.
#' @return A `tox_lc50` or a `tox_estimate`.
#' @noRd
point_branch <- function(x, p, alpha, nboot, seed) {
  if (x$type == "quantal") {
    lc50(x, p = if (is.null(p)) c(1, 50) else p, alpha = alpha)
  } else {
    icp(x, p = if (is.null(p)) 25 else p, nboot = nboot, seed = seed)
  }
}

#' The point-estimation results as a data frame, or an empty frame
#'
#' @param x A `tox_test` object.
#' @return A data frame with the same columns `endpoint_frame()` produces.
#' @noRd
point_estimates <- function(x) {
  if (is.null(x$point)) {
    return(NULL)
  }
  if (inherits(x$point, "tox_lc50")) {
    x$point$estimate$estimates
  } else {
    x$point$estimates
  }
}
