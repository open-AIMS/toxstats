# Wilcoxon Rank Sum Test with Bonferroni's adjustment

Compares each concentration with the control by ranking the two
together, dividing the significance level by the number of comparisons.
This is the test the EPA flowchart selects when the assumptions fail and
replication is unequal, which rules out Steel's test.

## Usage

``` r
wilcoxon_rank_sum(x, ..., alpha = 0.05)
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

The manual compares each rank sum with a tabulated critical value (Table
F.5). This package computes a p-value with
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) and
compares it with `alpha / k`, which avoids reproducing the table and
reproduces the manual's conclusion on its own worked example.

Where the data contain ties an exact p-value is unavailable and the
normal approximation with a continuity correction is used, which
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
signals with a warning. Ties are common in reproduction counts, so the
warning is suppressed here and the number of tied comparisons is
recorded in the result instead.

Bonferroni's adjustment bounds the overall error rate rather than fixing
it, so this test is more conservative than Steel's and correspondingly
less powerful. It requires at least four replicates per concentration
(section 9.4.5.2).

## References

US EPA (2002) EPA-821-R-02-013, Appendix F.

## Examples

``` r
# Appendix F uses the Appendix E data with two males removed, one from the
# control and one from the 12 per cent concentration.
repro <- ceriodaphnia_e1[ceriodaphnia_e1$conc < 50, ]
unequal <- repro[!(repro$conc == 0 & repro$replicate == "1") &
  !(repro$conc == 12 & repro$replicate == "5"), ]
wilcoxon_rank_sum(unequal, response = "young")
#> Wilcoxon Rank Sum Test with Bonferroni's adjustment
#>   EPA-821-R-02-013 Appendix F
#> 
#>  conc  n  mean rank_sum statistic   p_value significant
#>     3 10 18.20       79        24 4.277e-02       FALSE
#>     6 10 19.80       57         2 6.495e-05        TRUE
#>    12  9 21.89       58        13 6.993e-03        TRUE
#>    25 10  8.90       55         0 1.083e-05        TRUE
#> 
#>   NOEC 3
#>   LOEC 6
```
