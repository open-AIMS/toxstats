# Estimate lethal concentrations by the probit method

Fits a probit regression of response on log10 concentration and reads
the requested lethal concentrations off it. The EPA flowchart selects
this method when at least two partial responses occur and the fit is
adequate.

## Usage

``` r
probit_lc(
  x,
  ...,
  p = c(1, 50),
  ci_level = 0.95,
  alpha = 0.05,
  heterogeneity = TRUE
)
```

## Arguments

- x:

  A
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  object, or a data frame that can be coerced to one by passing it
  through
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md).

- ...:

  Additional arguments passed to
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  when `x` is a bare data frame, and otherwise unused.

- p:

  The effect percentage or percentages to estimate, on a 1-99 scale. For
  example `25` requests an IC25, and `50` an LC50.

- ci_level:

  Confidence level for interval estimates. Defaults to 0.95.

- alpha:

  Significance level for the hypothesis tests themselves. The EPA
  manuals specify 0.05, one-sided.

- heterogeneity:

  Should the variance be inflated by `chi-square / df` when the
  heterogeneity chi-square is significant at `alpha`? A flag.

## Value

An object of class `tox_estimate`.

## Details

The model is `glm(family = binomial(link = "probit"))` on `log10(conc)`.
Where the control response is not zero the proportions are
Abbott-adjusted first, as the manual's own output does.

Adequacy of fit is judged by the **Pearson** chi-square statistic for
heterogeneity, compared with the chi-square distribution on `k - 2`
degrees of freedom. On the manual's worked example this gives 3.076
against the 3.076 printed, confirming that the statistic is Pearson's
and not the deviance, which is 3.859 for the same fit.

Confidence limits are Fieller's. That reproduces the manual's printed
limits exactly on the worked example: 18.787 to 27.846 for the LC50 and
4.147 to 10.959 for the LC1. The delta method does not, giving 19.04 to
27.47 and 5.17 to 12.14, so the choice is settled by evidence rather
than convention.

When the heterogeneity chi-square is significant the manual's practice
is to inflate the variance by `chi-square / df`. The worked example does
not exercise this, its chi-square being well below the tabular value, so
the behaviour is offered under `heterogeneity` and defaults to applying
it.

## References

US EPA (2002) EPA-821-R-02-012, section 11.2.5.

Finney DJ (1978) *Statistical Method in Biological Assay*, 3rd edition.

## Examples

``` r
# Acute Table 20; the manual prints LC50 22.872 (18.787, 27.846) and
# LC1 7.924 (4.147, 10.959), with a chi-square of 3.076.
probit_lc(
  acute_table20,
  response = "probit", n_exposed = "exposed", type = "quantal",
  p = c(1, 50)
)
#> Probit method
#>   EPA-821-R-02-012 section 11.2.5
#> 
#>  endpoint  p estimate   lower  upper                  ci_method
#>       LC1  1   7.9239  4.1468 10.959 Fieller on the log10 scale
#>      LC50 50  22.8720 18.7870 27.846 Fieller on the log10 scale
#> 
#>   Chi-square for heterogeneity: 3.076 on 3 df (critical 7.815)
```
