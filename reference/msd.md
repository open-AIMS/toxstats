# Minimum significant difference

The smallest difference from the control that the test could have
detected as significant. Reporting it alongside a no-observed-effect
concentration is what distinguishes a test that found no effect from a
test that could not have found one.

## Usage

``` r
msd(x, ..., alpha = 0.05, test = c("dunnett", "bonferroni", "sidak"))
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

- test:

  The multiple-comparison procedure whose critical value is used.
  `"dunnett"` is the EPA method for a balanced design; `"bonferroni"`
  and `"sidak"` correspond to the t test with Bonferroni's or
  Dunn-Sidak's adjustment, used when replication is unequal.

## Value

An object of class `tox_msd`, a list with elements `msd` (named by
concentration), `critical`, `sw`, `df`, `n_control`, `n`, `test`,
`alpha` and `balanced`.

## Details

The minimum significant difference is

\$\$MSD = d\\ S_w \sqrt{1/n_1 + 1/n_i}\$\$

where \\d\\ is the one-sided critical value of the multiple-comparison
procedure, \\S_w\\ is the square root of the within mean square from a
one-way analysis of variance, \\n_1\\ is the number of control
replicates and \\n_i\\ the number at concentration \\i\\ (Appendix C,
section 1.10).

The critical value for Dunnett's procedure is computed from the
multivariate t distribution by numerical integration. Because all `k`
comparisons share the control, the correlation matrix has the one-factor
form \\\rho\_{ij} = \lambda_i \lambda_j\\ with \\\lambda_i = \sqrt{n_i /
(n_1 + n_i)}\\, and conditioning on that single factor reduces the
problem to a two-dimensional integral.

The manuals instead reproduce a table credited to Miller (1981).
Computation is used here for three reasons: it agrees with the table,
returning 2.3561 against a printed 2.36 on the Appendix C example; it
extends to unbalanced designs the table does not cover; and it avoids
reproducing a third-party table whose permission does not transfer.

Integration is used in preference to
[`mvtnorm::qmvt()`](https://rdrr.io/pkg/mvtnorm/man/qmvt.html), which is
a randomised quasi-Monte Carlo method and returns a slightly different
answer on every call. The same data must give the same critical value
every time.

When the design is unbalanced the difference detectable is not the same
at every concentration, so a value is returned for each.

## References

US EPA (2002) EPA-821-R-02-013, Appendix C, section 1.10.

## See also

[`pmsd()`](https://open-aims.github.io/toxstats/reference/pmsd.md) to
express this as a percentage of the control mean.

## Examples

``` r
# Appendix C worked example; the manual prints MSD = 0.162
msd(fathead_c1, response = "weight")
#> Minimum significant difference
#>   EPA-821-R-02-013 Appendix C, section 1.10
#> 
#>   dunnett critical value = 2.3561 (one-sided, alpha = 0.05, df = 15)
#>   within mean square root Sw = 0.0971
#> 
#>        32        64       128       256 
#> 0.1617744 0.1617744 0.1617744 0.1617744 
```
