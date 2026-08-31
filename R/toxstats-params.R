#' Shared parameter documentation
#'
#' Documentation-only object holding the `@param` tags shared across the
#' package. Every exported function inherits from it with
#' `@inheritParams toxstats_params`, so each argument is described in exactly
#' one place and cannot drift between topics.
#'
#' @param x A `tox_data()` object, or a data frame that can be coerced to
#'   one by passing it through `tox_data()`.
#' @param data A data frame with one row per replicate test chamber.
#' @param conc Name of the column holding the exposure concentration. A string.
#' @param response Name of the column holding the response. For
#'   `type = "continuous"` this is the measured value; for `type = "quantal"`
#'   it is the number of organisms affected. A string.
#' @param n_exposed Name of the column holding the number of organisms exposed
#'   in each replicate. Required when `type = "quantal"`, otherwise `NULL`.
#' @param replicate Name of an optional column labelling replicates within a
#'   concentration, or `NULL`.
#' @param type The response type. Either `"continuous"` (growth, reproduction,
#'   biomass) or `"quantal"` (survival, fertilisation, hatching).
#' @param direction The expected direction of the concentration-response
#'   relationship for the measured response. Either `"decreasing"` (the usual
#'   case, where toxicity reduces the response) or `"increasing"`.
#' @param control The concentration treated as the control. Defaults to `0`.
#' @param conc_units Units of `conc`, used only for labelling output. A string.
#' @param response_units Units of `response`, used only for labelling output.
#'   A string, or `NULL`.
#' @param alpha Significance level for the hypothesis tests themselves. The EPA
#'   manuals specify 0.05, one-sided.
#' @param alpha_assumption Significance level for the assumption tests
#'   (normality and homogeneity of variance). The EPA manuals specify 0.01.
#' @param p The effect percentage or percentages to estimate, on a 1-99 scale.
#'   For example `25` requests an IC25, and `50` an LC50.
#' @param ci_level Confidence level for interval estimates. Defaults to 0.95.
#' @param nboot Number of bootstrap resamples used for inhibition
#'   concentration confidence limits.
#' @param seed Optional integer seed, set before any resampling so that results
#'   are reproducible. `NULL` leaves the random number stream untouched.
#' @param ... Additional arguments passed to `tox_data()` when `x` is a
#'   bare data frame, and otherwise unused.
#'
#' @name toxstats_params
#' @keywords internal
NULL
