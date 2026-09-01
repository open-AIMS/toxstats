# Test the normality assumption

Tests whether the pooled within-group centred residuals are normally
distributed, the first branch point of the EPA hypothesis-testing
flowchart. The parametric tests downstream – Dunnett's procedure and the
t test with Bonferroni's adjustment – assume normality; failing this
test sends the analysis to the non-parametric branch.

## Usage

``` r
epa_normality(
  x,
  ...,
  method = c("auto", "shapiro_wilk", "kolmogorov"),
  alpha_assumption = 0.01
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

- method:

  Which statistic to use. `"auto"` (the default) follows the manuals:
  Shapiro-Wilk for 50 or fewer observations, Kolmogorov above that.

- alpha_assumption:

  Significance level for the assumption tests (normality and homogeneity
  of variance). The EPA manuals specify 0.01.

## Value

An object of class `tox_htest`, a list with elements `method`,
`statistic`, `statistic_name`, `p_value`, `critical`, `alpha`, `n`,
`normal`, `conclusion` and `reference`.

## Details

The manuals are specific about what is tested. Each observation is
centred by subtracting the mean of its own concentration, and the
resulting residuals are pooled across all concentrations and tested
together (Appendix B, section 2.3). Testing the raw values, or testing
each concentration separately, gives a different answer.

Two statistics are provided, and `"auto"` picks between them as the
manuals direct: Shapiro-Wilk for 50 or fewer observations, and the
Kolmogorov "D" statistic above that.

The Shapiro-Wilk statistic is computed by
[`stats::shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html),
which uses Royston's AS R94 algorithm. The manuals instead tabulate
Conover's coefficients and quantiles, a hand-calculation aid from before
this was available in software. The two agree closely: on the manual's
own Appendix B worked example, Royston returns `W = 0.9601` where the
manual prints `0.959`, and the two agree exactly once the manual's
rounding of the concentration means is applied. Royston's algorithm is
used here because it reproduces the published result, is not restricted
to the tabulated sample sizes, and returns a p-value rather than a
single fixed cutoff.

The Kolmogorov statistic is Stephens' (1974) modified form,
`D* = D (sqrt(n) - 0.01 + 0.85 / sqrt(n))`, compared against the
critical values in the manuals' Table B.11. It has no p-value, only a
decision at the tabulated alpha levels of 0.010, 0.025, 0.050, 0.100 and
0.150.

## References

US EPA (2002) EPA-821-R-02-013, Appendix B, sections 2.1 and 2.10.

Royston P (1995) A remark on Algorithm AS 181: the W test for normality.
*Applied Statistics* 44:547-551.

Stephens MA (1974) EDF statistics for goodness of fit and some
comparisons. *Journal of the American Statistical Association*
69:730-737.

## Examples

``` r
# Appendix B worked example; the manual prints W = 0.959 and concludes
# the data are normally distributed.
epa_normality(fathead_b1, response = "weight")
#> Shapiro-Wilk test for normality (Royston AS R94)
#>   EPA-821-R-02-013 Appendix B, section 2.1
#> 
#>   W = 0.96007, p = 0.5452
#>   alpha = 0.01
#>   The residuals are consistent with a normal distribution.
```
