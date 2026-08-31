#' Prepare toxicity test data for analysis
#'
#' Validates replicate-level concentration-response data and puts it into the
#' single form every other function in the package accepts. One row per
#' replicate test chamber covers both kinds of endpoint the EPA WET manuals
#' recognise.
#'
#' @details
#' The manuals distinguish two response types.
#'
#' A **continuous** endpoint is a measurement, such as growth in milligrams or
#' the number of young per female. `response` holds the measured value and
#' `n_exposed` is not used.
#'
#' A **quantal** endpoint is one where each organism either responds or does
#' not, such as survival, hatching or fertilisation. `response` holds the
#' number of organisms **affected** in the replicate and `n_exposed` the number
#' exposed. Both are required, because the manuals analyse quantal data twice:
#' the hypothesis-testing flowchart works on replicate-level proportions after
#' an arc sine square root transformation, while the point-estimation flowchart
#' works on counts pooled across replicates within a concentration. Recording
#' the number exposed is what makes both possible from one input.
#'
#' Rows with a missing response are dropped, and the number dropped is recorded
#' in the returned object. This matters for the linear interpolation method,
#' whose bootstrap resamples within a concentration and so depends on how many
#' replicates each concentration actually has.
#'
#' @inheritParams toxstats_params
#'
#' @return An object of class `tox_data`, a list with elements:
#'   \describe{
#'     \item{`replicates`}{a data frame with one row per replicate, holding
#'       `conc`, `replicate`, `response`, and for quantal data `n_exposed` and
#'       `proportion`.}
#'     \item{`pooled`}{a data frame with one row per concentration, holding
#'       `n_rep`, `mean`, `sd` and `var` of the analysis response, and for
#'       quantal data the totals `n_exposed`, `n_affected` and the pooled
#'       `proportion`.}
#'     \item{`type`, `direction`, `control`, `conc_units`, `response_units`}{as
#'       supplied.}
#'     \item{`n_dropped`}{the number of rows removed for a missing response.}
#'   }
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, sections 9.4 to 9.6.
#'
#' @examples
#' growth <- data.frame(
#'   conc = rep(c(0, 32, 64, 128, 256), each = 4),
#'   growth = c(
#'     0.711, 0.662, 0.054, 0.785,
#'     0.517, 0.501, 0.723, 0.760,
#'     0.602, 0.669, 0.694, 0.706,
#'     0.348, 0.400, 0.041, 0.512,
#'     0.216, 0.277, 0.328, 0.347
#'   )
#' )
#' tox_data(growth, response = "growth")
#'
#' @export
tox_data <- function(
  data,
  conc = "conc",
  response = "response",
  n_exposed = NULL,
  replicate = NULL,
  type = c("continuous", "quantal"),
  direction = c("decreasing", "increasing"),
  control = 0,
  conc_units = "%",
  response_units = NULL
) {
  type <- match.arg(type)
  direction <- match.arg(direction)

  chk::chk_data(data)
  chk::chk_string(conc)
  chk::chk_string(response)
  chk::chk_subset(conc, names(data))
  chk::chk_subset(response, names(data))
  chk::chk_null_or(n_exposed, vld = chk::vld_string)
  chk::chk_null_or(replicate, vld = chk::vld_string)
  chk::chk_number(control)
  chk::chk_string(conc_units)
  chk::chk_null_or(response_units, vld = chk::vld_string)

  if (!nrow(data)) {
    chk::abort_chk("`data` must have at least one row.")
  }
  if (!is.null(n_exposed)) {
    chk::chk_subset(n_exposed, names(data))
  }
  if (!is.null(replicate)) {
    chk::chk_subset(replicate, names(data))
  }
  if (type == "quantal" && is.null(n_exposed)) {
    chk::abort_chk(
      "`n_exposed` must be supplied when `type` is \"quantal\", because the ",
      "point-estimation methods work on counts pooled within a concentration."
    )
  }

  out <- data.frame(
    conc = data[[conc]],
    response = data[[response]],
    stringsAsFactors = FALSE
  )
  chk::chk_numeric(out$conc)
  chk::chk_not_any_na(out$conc)
  chk::chk_gte(min(out$conc), 0)
  chk::chk_numeric(out$response)

  out$replicate <- if (is.null(replicate)) {
    NA_character_
  } else {
    as.character(data[[replicate]])
  }
  out$n_exposed <- if (is.null(n_exposed)) {
    NA_real_
  } else {
    as.numeric(data[[n_exposed]])
  }

  keep <- !is.na(out$response)
  n_dropped <- sum(!keep)
  out <- out[keep, , drop = FALSE]
  if (!nrow(out)) {
    chk::abort_chk("`data` has no rows with a non-missing response.")
  }

  if (type == "quantal") {
    chk::chk_not_any_na(out$n_exposed)
    if (any(out$n_exposed <= 0)) {
      chk::abort_chk("`n_exposed` must be all positive.")
    }
    if (any(out$response < 0)) {
      chk::abort_chk(
        "`response` must be non-negative for quantal data; it is the number ",
        "of organisms affected."
      )
    }
    if (any(out$response > out$n_exposed)) {
      chk::abort_chk(
        "`response` must not exceed `n_exposed`; more organisms cannot be ",
        "affected than were exposed."
      )
    }
    out$proportion <- out$response / out$n_exposed
  }

  if (!control %in% out$conc) {
    chk::abort_chk(
      "`control` (",
      control,
      ") must be one of the concentrations present ",
      "in the data: ",
      paste(sort(unique(out$conc)), collapse = ", "),
      "."
    )
  }

  # Replicate labels are only ever used for reporting, so a missing label is
  # filled in per concentration rather than treated as an error.
  if (is.null(replicate)) {
    out$replicate <- as.character(stats::ave(
      out$conc,
      out$conc,
      FUN = seq_along
    ))
  }

  out <- out[order(out$conc), , drop = FALSE]
  rownames(out) <- NULL

  structure(
    list(
      replicates = out,
      pooled = pool_replicates(out, type = type),
      type = type,
      direction = direction,
      control = control,
      conc_units = conc_units,
      response_units = response_units,
      n_dropped = n_dropped,
      call = match.call()
    ),
    class = "tox_data"
  )
}

#' Summarise replicates within each concentration
#'
#' @param x The validated replicate data frame.
#' @param type Either "continuous" or "quantal".
#' @return A data frame with one row per concentration.
#' @noRd
pool_replicates <- function(x, type) {
  # The analysis response is the measured value for a continuous endpoint and
  # the replicate proportion for a quantal one, because that is what the
  # hypothesis-testing flowchart operates on in each case.
  value <- if (type == "quantal") x$proportion else x$response
  split_value <- split(value, x$conc)

  out <- data.frame(
    conc = as.numeric(names(split_value)),
    n_rep = lengths(split_value),
    mean = vapply(split_value, mean, numeric(1)),
    sd = vapply(split_value, stats::sd, numeric(1)),
    var = vapply(split_value, stats::var, numeric(1)),
    stringsAsFactors = FALSE
  )

  if (type == "quantal") {
    exposed <- vapply(split(x$n_exposed, x$conc), sum, numeric(1))
    affected <- vapply(split(x$response, x$conc), sum, numeric(1))
    out$n_exposed <- exposed
    out$n_affected <- affected
    out$proportion <- affected / exposed
  }

  out <- out[order(out$conc), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @describeIn tox_data Print a compact description of the design.
#'
#' @param x A `tox_data` object.
#'
#' @export
print.tox_data <- function(x, ...) {
  cat("<tox_data>\n")
  cat(
    "  ",
    x$type,
    " response, expected to be ",
    x$direction,
    " with concentration\n",
    sep = ""
  )
  cat(
    "  ",
    nrow(x$pooled),
    " concentrations, ",
    nrow(x$replicates),
    " replicates",
    if (all(x$pooled$n_rep == x$pooled$n_rep[1])) {
      paste0(" (", x$pooled$n_rep[1], " per concentration, balanced)")
    } else {
      " (unbalanced)"
    },
    "\n",
    sep = ""
  )
  cat("  control at ", x$control, " ", x$conc_units, "\n", sep = "")
  if (x$n_dropped) {
    cat("  ", x$n_dropped, " row(s) dropped for a missing response\n", sep = "")
  }
  cat("\n")
  print(x$pooled)
  invisible(x)
}

#' @describeIn tox_data Return the per-concentration summary.
#'
#' @param object A `tox_data` object.
#' @param ... Unused, present for consistency with the generic.
#'
#' @export
summary.tox_data <- function(object, ...) {
  chk::chk_unused(...)
  object$pooled
}

#' @describeIn tox_data Return the validated replicate-level data.
#'
#' @param row.names Unused, present for consistency with the generic.
#' @param optional Unused, present for consistency with the generic.
#'
#' @export
as.data.frame.tox_data <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...
) {
  chk::chk_unused(...)
  x$replicates
}
