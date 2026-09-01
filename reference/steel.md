# Steel's Many-One Rank Test

Compares each concentration with the control by ranking the two
together, holding the overall error rate at `alpha` across all the
comparisons. This is the test the EPA flowchart selects when the
residuals are not normal, or the variances are not homogeneous, and
replication is equal.

## Usage

``` r
steel(
  x,
  ...,
  alpha = 0.05,
  method = c("asymptotic", "simulated", "exact"),
  nsim = 10000
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

- alpha:

  Significance level for the hypothesis tests themselves. The EPA
  manuals specify 0.05, one-sided.

- method:

  How the p-value is obtained: `"asymptotic"` (the default),
  `"simulated"`, or `"exact"`. Passed to
  [`kSamples::Steel.test()`](https://rdrr.io/pkg/kSamples/man/Steel.test.html).

- nsim:

  Number of simulations when `method = "simulated"`.

## Value

An object of class `tox_comparison`, a list with elements `test`,
`method`, `comparisons` (a data frame with one row per treatment
concentration), `noec`, `loec`, `monotone`, `critical`, `sw`, `df`,
`alpha` and `reference`.

## Details

For each control and concentration pair the observations are combined
and ranked, and the ranks belonging to the concentration are summed. A
small rank sum indicates a response below the control.

The manual compares that rank sum with a tabulated critical value (Table
E.5). This package instead obtains a multiplicity-adjusted p-value from
[`kSamples::Steel.test()`](https://rdrr.io/pkg/kSamples/man/Steel.test.html),
which evaluates the joint null distribution of all `k` comparisons
directly. That avoids reproducing a third-party table, extends beyond
the tabulated designs, and reproduces the manual's conclusion on its own
worked example, including the borderline comparison whose rank sum
equals the tabulated critical value exactly.

Steel's test requires equal replication. The manual directs unequal
replication to
[`wilcoxon_rank_sum()`](https://open-aims.github.io/toxstats/reference/wilcoxon_rank_sum.md),
and this function signals an error rather than returning a value the
manual would not accept. It also requires at least four replicates per
concentration (section 9.4.5.2).

## References

US EPA (2002) EPA-821-R-02-013, Appendix E.

Steel RGD (1959) A multiple comparison rank sum test: treatments versus
control. *Biometrics* 15:560-572.

## Examples

``` r
# Appendix E worked example. The 50 per cent concentration is excluded from
# the reproduction analysis because survival there was significantly
# reduced; see section 9.5.2 and the Appendix G example.
repro <- ceriodaphnia_e1[ceriodaphnia_e1$conc < 50, ]
steel(repro, response = "young")
#> Steel's Many-One Rank Test (asymptotic p-values)
#>   EPA-821-R-02-013 Appendix E
#> 
#>  conc  n mean rank_sum statistic   p_value significant
#>     3 10 18.2     84.0    -1.597 0.1645000       FALSE
#>     6 10 19.8     63.5    -3.156 0.0032070        TRUE
#>    12 10 21.7     76.0    -2.205 0.0476900        TRUE
#>    25 10  8.9     55.0    -3.802 0.0003019        TRUE
#> 
#>   NOEC 3
#>   LOEC 6
```
