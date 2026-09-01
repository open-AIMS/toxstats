# The t test with Dunn-Sidak's adjustment

As
[`bonferroni_t()`](https://open-aims.github.io/toxstats/reference/bonferroni_t.md),
but using Dunn-Sidak's adjustment, which is exact when the comparisons
are independent and so slightly less conservative than Bonferroni's.

## Usage

``` r
dunn_sidak_t(x, ..., alpha = 0.05)
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

The critical value is `qt(1 - (1 - (1 - alpha)^(1 / k)), df)`. ToxCalc
offered this test; the EPA flowchart does not select it, so using it is
a documented departure from the manual.

## References

Sidak Z (1967) Rectangular confidence regions for the means of
multivariate normal distributions. *Journal of the American Statistical
Association* 62:626-633.

## Examples

``` r
dunn_sidak_t(fathead_c1, response = "weight")
#> t test with Dunn-Sidak's adjustment
#>   Not on the EPA flowchart; offered by ToxCalc
#> 
#>  conc n   mean statistic critical p_value significant
#>    32 4 0.5752    1.4860     2.48 0.28070       FALSE
#>    64 4 0.6602    0.2476     2.48 0.87370       FALSE
#>   128 4 0.5650    1.6350     2.48 0.22400       FALSE
#>   256 4 0.4542    3.2480     2.48 0.01077        TRUE
#> 
#>   NOEC 128
#>   LOEC 256
```
