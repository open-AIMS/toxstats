#' Smooth a sequence to be monotone
#'
#' Replaces a sequence by the closest monotone sequence, by repeatedly pooling
#' adjacent values that violate the required direction and replacing them with
#' their weighted mean. The EPA manuals call this "smoothing" and require it
#' before the trimmed Spearman-Karber and linear interpolation methods, both of
#' which assume a monotone concentration-response relationship.
#'
#' @details
#' This is the pool-adjacent-violators algorithm (PAVA). The manuals describe
#' the operation as averaging a violating pair, but their own worked examples
#' pool runs of three and more values at once -- Appendix M pools four groups,
#' and Appendix K collapses five -- and the text acknowledges that "unusual
#' patterns in the deviations from monotonicity may require an additional step
#' of smoothing". PAVA is the algorithm that performs the pooling to
#' convergence, and it reproduces the published examples exactly; repeated
#' pairwise averaging does not, because it loses track of how many original
#' values each pooled block represents.
#'
#' Weights control how a pooled block is averaged. The default of equal weights
#' is what the EPA linear interpolation method requires, because it smooths the
#' concentration *means* rather than the underlying observations. Supply `w` to
#' weight by the number of observations behind each value, which is what
#' Williams' test requires for its isotonic dose means.
#'
#' @param y A numeric vector to smooth.
#' @param w Optional numeric vector of weights the same length as `y`, all
#'   positive. Defaults to equal weights.
#' @param direction The direction the smoothed sequence must follow.
#'   `"decreasing"` (the default) returns a non-increasing sequence, as
#'   expected for a response reduced by a toxicant. `"increasing"` returns a
#'   non-decreasing sequence, as expected for a proportion responding.
#'
#' @return A numeric vector the same length as `y`, monotone in the requested
#'   direction, with the same weighted mean as `y`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix L section 4 and Appendix M
#' section 7.
#'
#' Barlow RE, Bartholomew DJ, Bremner JM, Brunk HD (1972) *Statistical
#' Inference Under Order Restrictions*. Wiley, London.
#'
#' @examples
#' # Appendix K mortality proportions collapse to a single value
#' smooth_monotone(c(0.05, 0.00, 0.05, 0.00, 0.00), direction = "increasing")
#'
#' @export
smooth_monotone <- function(
  y,
  w = NULL,
  direction = c("decreasing", "increasing")
) {
  direction <- match.arg(direction)
  chk::chk_numeric(y)
  chk::chk_not_any_na(y)

  if (is.null(w)) {
    w <- rep(1, length(y))
  }
  chk::chk_numeric(w)
  chk::chk_not_any_na(w)
  if (length(w) != length(y)) {
    chk::abort_chk("`w` must be the same length as `y`.")
  }
  if (any(w <= 0)) {
    chk::abort_chk("`w` must be all positive.")
  }

  if (length(y) <= 1L) {
    return(y)
  }

  # PAVA fits a non-decreasing sequence. A non-increasing fit is the negated
  # fit of the negated input, which avoids maintaining two near-identical
  # implementations.
  if (direction == "decreasing") {
    -pava(-y, w)
  } else {
    pava(y, w)
  }
}

#' Pool adjacent violators, fitting a non-decreasing sequence
#'
#' @param y Numeric vector.
#' @param w Numeric vector of positive weights, same length as `y`.
#' @return Numeric vector the same length as `y`, non-decreasing.
#' @noRd
pava <- function(y, w) {
  n <- length(y)
  # Blocks are held on a stack: `value` is the weighted block mean, `weight`
  # the block's total weight and `size` the number of original elements it
  # spans. Merging on the stack is what lets a run of three or more violating
  # values pool to their common mean in one pass.
  value <- numeric(n)
  weight <- numeric(n)
  size <- integer(n)
  k <- 0L

  for (i in seq_len(n)) {
    k <- k + 1L
    value[k] <- y[i]
    weight[k] <- w[i]
    size[k] <- 1L

    while (k > 1L && value[k - 1L] > value[k]) {
      pooled <- weight[k - 1L] + weight[k]
      value[k - 1L] <-
        (weight[k - 1L] * value[k - 1L] + weight[k] * value[k]) / pooled
      weight[k - 1L] <- pooled
      size[k - 1L] <- size[k - 1L] + size[k]
      k <- k - 1L
    }
  }

  rep(value[seq_len(k)], size[seq_len(k)])
}

#' Is a sequence monotone?
#'
#' Reports whether a sequence already follows the required direction, and where
#' it first departs from it. Williams' test assumes monotonicity, and the EPA
#' linear interpolation method smooths towards it, so both need to know.
#'
#' @inheritParams smooth_monotone
#'
#' @return A list with elements `monotone` (a flag), `direction` (the direction
#'   tested) and `violations` (an integer vector of the positions at which the
#'   sequence departs from the required direction; empty when `monotone` is
#'   `TRUE`).
#'
#' @examples
#' is_monotone(c(10, 8, 9, 4), direction = "decreasing")
#'
#' @export
is_monotone <- function(y, direction = c("decreasing", "increasing")) {
  direction <- match.arg(direction)
  chk::chk_numeric(y)
  chk::chk_not_any_na(y)

  if (length(y) <= 1L) {
    return(list(
      monotone = TRUE,
      direction = direction,
      violations = integer(0)
    ))
  }

  differences <- diff(y)
  violations <- if (direction == "decreasing") {
    which(differences > 0) + 1L
  } else {
    which(differences < 0) + 1L
  }

  list(
    monotone = length(violations) == 0L,
    direction = direction,
    violations = violations
  )
}
