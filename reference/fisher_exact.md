# Fisher's Exact Test

Compares the proportion responding at each concentration with the
control proportion, using the hypergeometric distribution. The EPA
manuals apply it to survival, whose outcome is counted rather than
measured.

## Usage

``` r
fisher_exact(x, ..., alpha = 0.05, adjust = c("none", "bonferroni"))
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

  How to adjust for multiplicity. `"none"` is the EPA default for this
  test.

## Value

An object of class `tox_comparison`, a list with elements `test`,
`method`, `comparisons` (a data frame with one row per treatment
concentration), `noec`, `loec`, `monotone`, `critical`, `sw`, `df`,
`alpha` and `reference`.

## Details

Requires quantal data, so `tox_data(type = "quantal")` with `n_exposed`
supplied.

**No multiplicity adjustment is applied by default.** Section 1 of
Appendix G states that because Fisher's Exact Test is itself
conservative, "a pair-wise comparison error rate of 0.05 is suggested
rather than an experiment-wise error rate". Set `adjust` to depart from
that.

The manual works the test by looking up a table of significance levels
(Table G.5). This package uses
[`stats::fisher.test()`](https://rdrr.io/r/stats/fisher.test.html),
which computes the same probability directly.

## References

US EPA (2002) EPA-821-R-02-013, Appendix G.

## Examples

``` r
# Appendix G worked example; the manual reports NOEC 12 and LOEC 25
fisher_exact(
  ceriodaphnia_g2,
  response = "dead",
  n_exposed = "exposed",
  type = "quantal"
)
#> Fisher's Exact Test
#>   EPA-821-R-02-013 Appendix G
#> 
#>  conc n_exposed n_affected proportion   p_value significant
#>     1        10          0          0 1.0000000       FALSE
#>     3        10          0          0 1.0000000       FALSE
#>     6        10          0          0 1.0000000       FALSE
#>    12        10          0          0 1.0000000       FALSE
#>    25        10         10          1 0.0001191        TRUE
#> 
#>   NOEC 12
#>   LOEC 25
```
