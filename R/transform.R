#' Arc sine square root transformation
#'
#' Transforms proportion data as specified in Appendix B, section 4.2 of the
#' EPA WET method manuals. The manuals require this transformation before the
#' parametric hypothesis tests are applied to a quantal endpoint such as
#' survival, because proportions have a variance that depends on the mean.
#'
#' @details
#' The transformation is `asin(sqrt(p))`, returned in **radians**, with two
#' endpoint adjustments that the manuals treat as part of the method:
#'
#' * an observed proportion of exactly 0 is replaced by `1 / (4 * n)`;
#' * an observed proportion of exactly 1 is replaced by `1 - 1 / (4 * n)`;
#'
#' where `n` is the number of organisms in that replicate. Without them a
#' replicate at 0 or 1 would contribute no variance, and the tests downstream
#' would understate the true variability. The adjustments follow Bartlett
#' (1937). Omitting them is the commonest way to fail to reproduce a published
#' EPA analysis, so `n` is required whenever any proportion is 0 or 1.
#'
#' @param p A numeric vector of proportions, each in `[0, 1]`.
#' @param n The number of organisms in each replicate, used only for the
#'   endpoint adjustments. Either a single number or a vector the same length
#'   as `p`. Required if any element of `p` is exactly 0 or exactly 1.
#'
#' @return A numeric vector of transformed values, in radians, the same length
#'   as `p`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix B, section 4.2.
#'
#' Bartlett MS (1937) Some examples of statistical methods of research in
#' agriculture and applied biology. *Supplement to the Journal of the Royal
#' Statistical Society* 4:137-183.
#'
#' @seealso [inv_arcsine_sqrt()] for the back-transformation.
#'
#' @examples
#' # Appendix B, section 4.2 worked values
#' arcsine_sqrt(0.60)
#' arcsine_sqrt(c(0, 1), n = 20)
#'
#' @export
arcsine_sqrt <- function(p, n = NULL) {
  chk::chk_numeric(p)
  chk::chk_not_any_na(p)
  chk::chk_range(p, c(0, 1))

  needs_n <- any(p == 0 | p == 1)
  if (needs_n) {
    if (is.null(n)) {
      chk::abort_chk(
        "`n` must be supplied when any element of `p` is exactly 0 or 1, ",
        "because the EPA endpoint adjustments are defined in terms of the ",
        "number of organisms in the replicate."
      )
    }
    chk::chk_numeric(n)
    chk::chk_gt(min(n), 0)
    if (length(n) != 1L && length(n) != length(p)) {
      chk::abort_chk("`n` must be length 1 or the same length as `p`.")
    }
  }

  n <- if (is.null(n)) rep(NA_real_, length(p)) else rep_len(n, length(p))

  adjusted <- p
  adjusted[p == 0] <- 1 / (4 * n[p == 0])
  adjusted[p == 1] <- 1 - 1 / (4 * n[p == 1])

  asin(sqrt(adjusted))
}

#' Back-transform an arc sine square root value
#'
#' Inverts [arcsine_sqrt()], returning a proportion from a value in radians.
#' The EPA manuals use this to express a minimum significant difference on the
#' original proportion scale (Appendix C, section 1.11.2), where it is more
#' readily interpreted than an angle.
#'
#' @details
#' The back-transformation is `sin(theta)^2`. It does not undo the endpoint
#' adjustments made by [arcsine_sqrt()]: a proportion of 0 transformed with
#' `n = 20` and then back-transformed returns `1 / 80`, not 0. That is intended
#' -- the adjusted value is the one the analysis actually used.
#'
#' @param theta A numeric vector of transformed values in radians, each in
#'   `[0, pi / 2]`.
#'
#' @return A numeric vector of proportions the same length as `theta`.
#'
#' @references
#' US EPA (2002) EPA-821-R-02-013, Appendix C, section 1.11.2.
#'
#' @seealso [arcsine_sqrt()].
#'
#' @examples
#' inv_arcsine_sqrt(arcsine_sqrt(0.60))
#'
#' @export
inv_arcsine_sqrt <- function(theta) {
  chk::chk_numeric(theta)
  chk::chk_not_any_na(theta)
  chk::chk_range(theta, c(0, pi / 2))

  sin(theta)^2
}

#' Abbott's correction for control response
#'
#' Adjusts observed proportions for the response already present in the
#' control, so that the adjusted values estimate the proportion responding to
#' the toxicant alone. The EPA point-estimation methods apply this before
#' estimating a median lethal concentration when control mortality is
#' non-zero.
#'
#' @details
#' The correction is
#'
#' \deqn{p' = \frac{p - p_0}{1 - p_0}}
#'
#' where \eqn{p_0} is the control proportion. It is undefined when
#' \eqn{p_0 = 1}, which is an error rather than a warning: a control in which
#' every organism responded carries no information about the toxicant.
#'
#' Where an observed proportion falls below the control, the raw correction is
#' negative. By default such values are clamped to 0, which is what the
#' point-estimation methods require. Set `clamp = FALSE` to see the raw values,
#' for instance when diagnosing a control problem.
#'
#' @param p A numeric vector of observed proportions, each in `[0, 1]`.
#' @param p_control The control proportion, a single number in `[0, 1)`.
#' @param clamp Should adjusted values outside `[0, 1]` be clamped to that
#'   range? A flag, defaulting to `TRUE`.
#'
#' @return A numeric vector of adjusted proportions the same length as `p`.
#'
#' @references
#' Abbott WS (1925) A method of computing the effectiveness of an insecticide.
#' *Journal of Economic Entomology* 18:265-267.
#'
#' @examples
#' abbott(c(0.02, 0.05, 0.80), p_control = 0.02)
#'
#' @export
abbott <- function(p, p_control, clamp = TRUE) {
  chk::chk_numeric(p)
  chk::chk_not_any_na(p)
  chk::chk_range(p, c(0, 1))
  chk::chk_number(p_control)
  chk::chk_range(p_control, c(0, 1))
  chk::chk_flag(clamp)

  if (p_control == 1) {
    chk::abort_chk(
      "`p_control` must be less than 1; a control in which every organism ",
      "responded carries no information about the toxicant."
    )
  }

  out <- (p - p_control) / (1 - p_control)
  if (clamp) {
    out <- pmin(pmax(out, 0), 1)
  }
  out
}
