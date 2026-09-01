# Welch's t test against the control

Compares each concentration mean with the control mean without pooling
the variance, so it remains valid when the variances differ. Each
comparison carries its own degrees of freedom, from the
Welch-Satterthwaite approximation.

## Usage

``` r
welch_t(x, ..., alpha = 0.05, adjust = c("bonferroni", "sidak", "none"))
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

- alpha:

  Significance level for the hypothesis tests themselves. The EPA
  manuals specify 0.05, one-sided.

- adjust:

  How to adjust for multiplicity across the `k` comparisons.

## Value

An object of class `tox_comparison`, a list with elements `test`,
`method`, `comparisons` (a data frame with one row per treatment
concentration), `noec`, `loec`, `monotone`, `critical`, `sw`, `df`,
`alpha` and `reference`.

## Details

This is not on the EPA multi-concentration flowchart, which sends
heterogeneous variances to the non-parametric branch. It appears in the
single-concentration chart of the acute manual (Figure 12) and was
offered by ToxCalc, so it is provided here and flagged as a departure.

## References

US EPA (2002) EPA-821-R-02-012, Figure 12.

## Examples

``` r
welch_t(fathead_c1, response = "weight")
#> Welch's t test against the control, bonferroni adjusted (not an EPA multi-concentration method)
#>   EPA-821-R-02-012 Figure 12; not on the Figure 2 flowchart
#> 
#>  conc n   mean statistic    df p_value critical significant
#>    32 4 0.5752    1.9310 3.483 0.06807    3.777       FALSE
#>    64 4 0.6602    0.6861 5.448 0.26040    3.065       FALSE
#>   128 4 0.5650    1.9400 3.398 0.06847    3.835       FALSE
#>   256 4 0.4542    2.9640 3.230 0.02705    3.964       FALSE
#> 
#>   NOEC 256
#>   LOEC > 256
```
