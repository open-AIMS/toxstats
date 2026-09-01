# The t test with Bonferroni's adjustment

Compares each concentration mean with the control mean, dividing the
significance level by the number of comparisons. This is the parametric
test the EPA flowchart selects when the assumptions are met but
replication is unequal, which rules out Dunnett's procedure.

## Usage

``` r
bonferroni_t(x, ..., alpha = 0.05)
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

## Value

An object of class `tox_comparison`, a list with elements `test`,
`method`, `comparisons` (a data frame with one row per treatment
concentration), `noec`, `loec`, `monotone`, `critical`, `sw`, `df`,
`alpha` and `reference`.

## Details

The statistic is the same as in
[`dunnett()`](https://open-aims.github.io/toxstats/reference/dunnett.md);
only the critical value differs, being `qt(1 - alpha / k, df)` for `k`
treatment concentrations. This reproduces the manual's Table D.5
exactly: 2.510 for `alpha = 0.05` with 14 degrees of freedom and four
concentrations, and 2.241 in the row for infinite degrees of freedom
against a printed 2.242.

Bonferroni's adjustment sets an upper bound of `alpha` on the overall
error rate rather than fixing it there, so it is more conservative than
Dunnett's procedure and correspondingly less powerful.

## References

US EPA (2002) EPA-821-R-02-013, Appendix D.

## Examples

``` r
# Appendix D uses the Appendix C data with one replicate lost
lost <- fathead_c1[-19, ]
bonferroni_t(lost, response = "weight")
#> t test with Bonferroni's adjustment
#>   EPA-821-R-02-013 Appendix D
#> 
#>  conc n   mean statistic critical  p_value significant
#>    32 4 0.5752    1.6220     2.51 0.254200       FALSE
#>    64 4 0.6602    0.2703     2.51 1.000000       FALSE
#>   128 4 0.5650    1.7850     2.51 0.191900       FALSE
#>   256 3 0.4037    4.0280     2.51 0.002492        TRUE
#> 
#>   NOEC 128
#>   LOEC 256
```
