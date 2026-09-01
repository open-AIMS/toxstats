# Shared parameter documentation

Documentation-only object holding the `@param` tags shared across the
package. Every exported function inherits from it with
`@inheritParams toxstats_params`, so each argument is described in
exactly one place and cannot drift between topics.

## Arguments

- x:

  A
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  object, or a data frame that can be coerced to one by passing it
  through
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md).

- data:

  A data frame with one row per replicate test chamber.

- conc:

  Name of the column holding the exposure concentration. A string.

- response:

  Name of the column holding the response. For `type = "continuous"`
  this is the measured value; for `type = "quantal"` it is the number of
  organisms affected. A string.

- n_exposed:

  Name of the column holding the number of organisms exposed in each
  replicate. Required when `type = "quantal"`, otherwise `NULL`.

- replicate:

  Name of an optional column labelling replicates within a
  concentration, or `NULL`.

- type:

  The response type. Either `"continuous"` (growth, reproduction,
  biomass) or `"quantal"` (survival, fertilisation, hatching).

- direction:

  The expected direction of the concentration-response relationship for
  the measured response. Either `"decreasing"` (the usual case, where
  toxicity reduces the response) or `"increasing"`.

- control:

  The concentration treated as the control. Defaults to `0`.

- conc_units:

  Units of `conc`, used only for labelling output. A string.

- response_units:

  Units of `response`, used only for labelling output. A string, or
  `NULL`.

- alpha:

  Significance level for the hypothesis tests themselves. The EPA
  manuals specify 0.05, one-sided.

- alpha_assumption:

  Significance level for the assumption tests (normality and homogeneity
  of variance). The EPA manuals specify 0.01.

- p:

  The effect percentage or percentages to estimate, on a 1-99 scale. For
  example `25` requests an IC25, and `50` an LC50.

- ci_level:

  Confidence level for interval estimates. Defaults to 0.95.

- nboot:

  Number of bootstrap resamples used for inhibition concentration
  confidence limits.

- seed:

  Optional integer seed, set before any resampling so that results are
  reproducible. `NULL` leaves the random number stream untouched.

- ...:

  Additional arguments passed to
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  when `x` is a bare data frame, and otherwise unused.
