# Dunnett's procedure

Compares each concentration mean with the control mean, holding the
overall error rate at `alpha` across all the comparisons. This is the
parametric test the EPA flowchart selects when the residuals are normal,
the variances are homogeneous, and every concentration has the same
number of replicates.

## Usage

``` r
dunnett(x, ..., alpha = 0.05)
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

The statistic for concentration `i` is

\$\$t_i = \frac{\bar{Y}\_1 - \bar{Y}\_i}{S_w \sqrt{1/n_1 + 1/n_i}}\$\$

with `Sw` the square root of the within mean square from a one-way
analysis of variance. It is compared with a one-sided critical value
computed as described in
[`msd()`](https://open-aims.github.io/toxstats/reference/msd.md); the
manual reads the same value from a table credited to Miller (1981).

Dunnett's procedure requires equal replication. When replication is
unequal the manual directs the analysis to
[`bonferroni_t()`](https://open-aims.github.io/toxstats/reference/bonferroni_t.md)
instead, and this function signals an error rather than silently
returning a value the manual would not accept.

## References

US EPA (2002) EPA-821-R-02-013, Appendix C.

Dunnett CW (1955) A multiple comparison procedure for comparing several
treatments with a control. *Journal of the American Statistical
Association* 50:1096-1121.

## Examples

``` r
# Appendix C worked example; the manual reports NOEC 128 and LOEC 256
dunnett(fathead_c1, response = "weight")
#> Dunnett's procedure
#>   EPA-821-R-02-013 Appendix C
#> 
#>  conc n   mean statistic critical  p_value significant
#>    32 4 0.5752    1.4860    2.356 0.207900       FALSE
#>    64 4 0.6602    0.2476    2.356 0.710700       FALSE
#>   128 4 0.5650    1.6350    2.356 0.167200       FALSE
#>   256 4 0.4542    3.2480    2.356 0.009069        TRUE
#> 
#>   NOEC 128
#>   LOEC 256
```
