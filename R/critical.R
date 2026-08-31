#' One-sided critical value for a many-to-one comparison
#'
#' @param test One of "dunnett", "bonferroni", "sidak".
#' @param alpha Familywise significance level.
#' @param df Residual degrees of freedom.
#' @param n_control Number of control replicates.
#' @param n_treat Number of replicates at each treatment concentration.
#' @return A single number.
#' @noRd
critical_value <- function(test, alpha, df, n_control, n_treat) {
  k <- length(n_treat)

  if (test == "dunnett") {
    return(dunnett_critical(alpha, n_control, n_treat, df))
  }

  per_comparison <- if (test == "bonferroni") {
    alpha / k
  } else {
    1 - (1 - alpha)^(1 / k)
  }
  stats::qt(1 - per_comparison, df)
}

#' Dunnett's one-sided critical value
#'
#' Computed by numerical integration rather than looked up in a table or
#' estimated by simulation, so that the same design always returns the same
#' value.
#'
#' The k comparisons all share the control, which gives the correlation matrix
#' a one-factor structure: `rho_ij = lambda_i lambda_j` with
#' `lambda_i = sqrt(n_i / (n_0 + n_i))`. Conditioning on the single common
#' factor and on the pooled standard deviation reduces the k-dimensional
#' multivariate t probability to a two-dimensional integral, which
#' `stats::integrate` evaluates deterministically.
#'
#' The alternative, `mvtnorm::qmvt()`, uses randomised quasi-Monte Carlo and
#' returns a different answer on every call -- values between 2.3552 and 2.3572
#' were observed for the Appendix C design. That is a poor property for an
#' analysis intended to be reproducible, and it is why this route was written.
#'
#' @param alpha Familywise one-sided significance level.
#' @param n_control Number of control replicates.
#' @param n_treat Integer vector of replicate counts per treatment.
#' @param df Residual degrees of freedom.
#' @return A single number.
#' @noRd
dunnett_critical <- function(alpha, n_control, n_treat, df) {
  # With one treatment there is no multiplicity, and the statistic is an
  # ordinary t. Taking this path also avoids integrating a degenerate case.
  if (length(n_treat) == 1L) {
    return(stats::qt(1 - alpha, df))
  }

  lambda <- sqrt(n_treat / (n_control + n_treat))
  stats::uniroot(
    function(d) dunnett_prob(d, lambda, df) - (1 - alpha),
    interval = c(0.5, 12),
    tol = 1e-9
  )$root
}

#' Probability that every one of k correlated t statistics is below d
#'
#' @param d The bound.
#' @param lambda One-factor loadings, `sqrt(n_i / (n_0 + n_i))`.
#' @param df Residual degrees of freedom.
#' @return A probability.
#' @noRd
dunnett_prob <- function(d, lambda, df) {
  # Inner integral: given the pooled scale s and the shared factor z, the k
  # comparisons are conditionally independent, so their joint probability is a
  # product of normal probabilities.
  conditional <- function(s) {
    vapply(
      s,
      function(si) {
        stats::integrate(
          function(z) {
            terms <- vapply(
              lambda,
              function(l) {
                stats::pnorm((l * z + d * si) / sqrt(1 - l^2))
              },
              numeric(length(z))
            )
            stats::dnorm(z) * apply(terms, 1, prod)
          },
          lower = -8.5,
          upper = 8.5,
          rel.tol = 1e-10
        )$value
      },
      numeric(1)
    )
  }

  stats::integrate(
    function(s) chi_scale_density(s, df) * conditional(s),
    lower = 0,
    upper = 8,
    rel.tol = 1e-10
  )$value
}

#' Density of s = sqrt(chi-squared_df / df)
#'
#' @param s Numeric vector of positive values.
#' @param df Degrees of freedom.
#' @return Numeric vector of densities.
#' @noRd
chi_scale_density <- function(s, df) {
  # Evaluated on the log scale, because (df/2)^(df/2) and gamma(df/2) both
  # overflow well before the degrees of freedom a large design could reach.
  out <- numeric(length(s))
  positive <- s > 0
  out[positive] <- exp(
    log(2) +
      (df / 2) * log(df / 2) -
      lgamma(df / 2) +
      (df - 1) * log(s[positive]) -
      df * s[positive]^2 / 2
  )
  out
}
