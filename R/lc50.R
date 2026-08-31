#' The EPA point-estimation flowchart
#'
#' Figure 6 of the acute manual, which chooses between the four methods of
#' estimating a median lethal concentration. Held as data for the same reasons
#' as the hypothesis-testing chart; see [decisions()].
#'
#' @noRd
flowchart_lc50 <- function() {
  data.frame(
    node = c("partial_two", "probit_fits", "partial_one", "sk_requirements"),
    question = c(
      "Are there two or more partial responses?",
      "Is the probit model appropriate?",
      "Is there at least one partial response?",
      "Does the response run from zero at the lowest concentration to complete at the highest?"
    ),
    rule = c(
      "partial_two",
      "probit_fits",
      "partial_one",
      "sk_requirements"
    ),
    yes = c("probit_fits", "@probit", "sk_requirements", "@spearman_karber"),
    no = c(
      "partial_one",
      "sk_requirements",
      "@graphical",
      "@trimmed_spearman_karber"
    ),
    reference = rep("EPA-821-R-02-012 Figure 6", 4),
    stringsAsFactors = FALSE
  )
}

#' @noRd
lc50_terminals <- function() {
  c(
    "@probit" = "probit_lc",
    "@spearman_karber" = "spearman_karber",
    "@trimmed_spearman_karber" = "trimmed_spearman_karber",
    "@graphical" = "graphical_lc50"
  )
}

#' @noRd
rule_partial_two <- function(state) {
  # Counted on the smoothed, adjusted proportions, not the raw ones. On the
  # manual's own Spearman-Karber example the raw data show two partial
  # responses but the smoothed adjusted data show one, and only the latter
  # routes to the method the manual actually uses.
  partial <- sum(state$parts$adjusted > 0 & state$parts$adjusted < 1)
  new_rule_result(
    answer = if (partial >= 2) "yes" else "no",
    criterion = "at least two concentrations with a partial response",
    statistic = partial,
    outcome = paste0(partial, " partial response(s)"),
    state = state
  )
}

#' @noRd
rule_partial_one <- function(state) {
  # Counted on the smoothed, adjusted proportions, not the raw ones. On the
  # manual's own Spearman-Karber example the raw data show two partial
  # responses but the smoothed adjusted data show one, and only the latter
  # routes to the method the manual actually uses.
  partial <- sum(state$parts$adjusted > 0 & state$parts$adjusted < 1)
  new_rule_result(
    answer = if (partial >= 1) "yes" else "no",
    criterion = "at least one concentration with a partial response",
    statistic = partial,
    outcome = if (partial >= 1) {
      "a model can be fitted"
    } else {
      "all-or-nothing data; only the graphical method applies"
    },
    state = state
  )
}

#' @noRd
rule_probit_fits <- function(state) {
  fit <- try(
    probit_lc(state$data, p = 50, alpha = state$alpha),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(new_rule_result(
      answer = "no",
      criterion = "the probit model must be fittable",
      outcome = "the probit model could not be fitted",
      state = state
    ))
  }

  adequate <- fit$df > 0 && fit$chisq <= fit$chisq_critical
  state$probit <- fit

  new_rule_result(
    answer = if (adequate) "yes" else "no",
    criterion = paste0(
      "Pearson chi-square for heterogeneity on ",
      fit$df,
      " df, reject if above ",
      signif(fit$chisq_critical, 4)
    ),
    statistic = fit$chisq,
    outcome = if (adequate) {
      "the probit model fits"
    } else {
      "heterogeneity is significant; the probit model does not fit"
    },
    state = state
  )
}

#' @noRd
rule_sk_requirements <- function(state) {
  adjusted <- state$parts$adjusted
  k <- length(adjusted)
  met <- adjusted[1] == 0 && adjusted[k] == 1

  new_rule_result(
    answer = if (met) "yes" else "no",
    criterion = "adjusted response 0 at the lowest and 1 at the highest",
    outcome = paste0(
      "adjusted response runs ",
      signif(adjusted[1], 3),
      " to ",
      signif(adjusted[k], 3)
    ),
    state = state
  )
}

#' Estimate a median lethal concentration
#'
#' Walks the EPA point-estimation flowchart, selects the method the manual
#' directs for the data, runs it, and returns the estimate together with a
#' record of how the method was chosen.
#'
#' @details
#' The flowchart is Figure 6 of the acute manual. It branches on how many
#' concentrations produced a partial response, on whether the probit model
#' fits, and on whether the response runs the full way from none to complete.
#' The four terminals are [probit_lc()], [spearman_karber()],
#' [trimmed_spearman_karber()] and [graphical_lc50()].
#'
#' Point estimation uses all the data. Unlike the hypothesis-testing branch, no
#' concentration is excluded for complete response: section 9.5.2 of the
#' chronic manual is explicit that concentrations dropped from the
#' no-observed-effect concentration are retained here.
#'
#' @inheritParams toxstats_params
#' @param method Optional. Name of a method to use instead of the one the
#'   flowchart selects: one of `"probit_lc"`, `"spearman_karber"`,
#'   `"trimmed_spearman_karber"` or `"graphical_lc50"`. Forcing a method warns
#'   and is recorded, as in [tox_test()].
#'
#' @return An object of class `tox_lc50`, a list with elements `estimate`
#'   (a `tox_estimate`), `decisions`, `selected`, `overridden` and `data`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-012, Figure 6 and section 11.2.
#'
#' @examples
#' # Acute Table 20, probit column: two partial responses and an adequate fit,
#' # so the flowchart selects the probit method.
#' fit <- lc50(
#'   acute_table20,
#'   response = "probit", n_exposed = "exposed", type = "quantal"
#' )
#' summary(fit)
#'
#' @export
lc50 <- function(
  x,
  ...,
  p = c(1, 50),
  ci_level = 0.95,
  alpha = 0.05,
  method = NULL
) {
  chk::chk_null_or(method, vld = chk::vld_string)
  if (!is.null(method)) {
    chk::chk_subset(method, unname(lc50_terminals()))
  }

  x <- as_tox_data(x, ...)
  parts <- quantal_prep(x)

  walked <- walk_flowchart(
    flowchart_lc50(),
    list(data = x, parts = parts, alpha = alpha, ci_level = ci_level)
  )

  selected <- unname(lc50_terminals()[walked$terminal])
  chosen <- selected
  overridden <- FALSE
  trail <- walked$decisions

  if (!is.null(method) && method != selected) {
    overridden <- TRUE
    chosen <- method
    warning(
      "The flowchart selected `",
      selected,
      "` but `method = \"",
      method,
      "\"` was given. The result is not the analysis the EPA manual ",
      "prescribes for these data.",
      call. = FALSE
    )
    trail <- rbind(
      trail,
      data.frame(
        step = nrow(trail) + 1L,
        node = "override",
        question = "Was the selected method overridden?",
        criterion = "user supplied `method`",
        statistic = NA_real_,
        p_value = NA_real_,
        answer = "yes",
        outcome = paste0(
          "user forced ",
          method,
          "; the flowchart selected ",
          selected
        ),
        reference = "Not an EPA decision",
        stringsAsFactors = FALSE
      )
    )
  }

  arguments <- list(x)
  if (chosen == "probit_lc") {
    arguments <- c(
      arguments,
      list(p = p, ci_level = ci_level, alpha = alpha)
    )
  } else if (chosen != "graphical_lc50") {
    arguments <- c(arguments, list(ci_level = ci_level))
  }
  estimate <- do.call(chosen, arguments)

  trail <- rbind(
    trail,
    data.frame(
      step = nrow(trail) + 1L,
      node = "method",
      question = "Which method was used?",
      criterion = NA_character_,
      statistic = NA_real_,
      p_value = NA_real_,
      answer = chosen,
      outcome = estimate$method,
      reference = estimate$reference,
      stringsAsFactors = FALSE
    )
  )
  rownames(trail) <- NULL

  structure(
    list(
      data = x,
      estimate = estimate,
      decisions = trail,
      selected = selected,
      overridden = overridden,
      flowchart = "EPA-821-R-02-012 Figure 6",
      call = match.call()
    ),
    class = "tox_lc50"
  )
}

#' @describeIn decisions Decision trail from an [lc50()] analysis.
#' @export
decisions.tox_lc50 <- function(x, ...) {
  chk::chk_unused(...)
  x$decisions
}

#' @export
print.tox_lc50 <- function(x, ...) {
  cat("EPA WET point estimate\n")
  cat("  Flowchart: ", x$flowchart, "\n", sep = "")
  cat("  Selected:  ", x$estimate$method, "\n", sep = "")
  if (x$overridden) {
    cat("  OVERRIDDEN: the flowchart selected ", x$selected, "\n", sep = "")
  }
  cat("\n")
  out <- x$estimate$estimates
  numeric_columns <- vapply(out, is.numeric, logical(1))
  out[numeric_columns] <- lapply(out[numeric_columns], signif, digits = 5)
  print(out, row.names = FALSE)
  cat("\nUse summary() for the decision trail.\n")
  invisible(x)
}

#' @describeIn lc50 Print the decision trail above the estimate.
#'
#' @param object A `tox_lc50` object.
#'
#' @export
summary.tox_lc50 <- function(object, ...) {
  chk::chk_unused(...)

  cat("EPA WET point estimate\n")
  cat("Flowchart: ", object$flowchart, "\n\n", sep = "")

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
      detail <- paste0(signif(trail$statistic[i], 4), " -> ", detail)
    }
    cat("      ", detail, "\n", sep = "")
    cat("      (", trail$reference[i], ")\n", sep = "")
  }

  cat("\n")
  out <- object$estimate$estimates
  numeric_columns <- vapply(out, is.numeric, logical(1))
  out[numeric_columns] <- lapply(out[numeric_columns], signif, digits = 5)
  print(out, row.names = FALSE)
  invisible(object)
}

#' @describeIn lc50 Return the estimates as a data frame.
#'
#' @param row.names Unused.
#' @param optional Unused.
#'
#' @export
as.data.frame.tox_lc50 <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...
) {
  chk::chk_unused(...)
  out <- x$estimate$estimates
  out$method <- x$estimate$method
  out$selected_by_flowchart <- x$selected
  out$overridden <- x$overridden
  out$conc_units <- x$data$conc_units
  out
}
