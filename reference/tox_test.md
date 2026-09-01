# Run the EPA hypothesis-testing analysis

Walks the EPA flowchart, selects the hypothesis test the manual directs
for the data, runs it, and returns the endpoints together with a record
of how the test was chosen. This is the function that recreates what
ToxCalc did.

## Usage

``` r
tox_test(
  x,
  ...,
  alpha = 0.05,
  alpha_assumption = 0.01,
  test = NULL,
  exclude = NULL,
  pmsd_bounds = NULL,
  branch = c("hypothesis", "both"),
  p = NULL,
  nboot = 200,
  seed = NULL
)

# S3 method for class 'tox_test'
summary(object, ...)

# S3 method for class 'tox_test'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
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

- alpha_assumption:

  Significance level for the assumption tests (normality and homogeneity
  of variance). The EPA manuals specify 0.01.

- test:

  Optional. Name of a test to run instead of the one the flowchart
  selects: one of `"dunnett"`, `"bonferroni_t"`, `"steel"` or
  `"wilcoxon_rank_sum"`.

- exclude:

  Optional numeric vector of concentrations to exclude from the
  hypothesis test, over and above those excluded automatically.

- pmsd_bounds:

  Optional. Passed to
  [`pmsd()`](https://open-aims.github.io/toxstats/reference/pmsd.md);
  supplying it also enables the section 10.2.8.2.5 override.

- branch:

  `"hypothesis"` (the default) gives the NOEC, LOEC, MSD and PMSD.
  `"both"` adds a point estimate alongside them, which is what a
  laboratory report usually contains: a median lethal concentration
  through
  [`lc50()`](https://open-aims.github.io/toxstats/reference/lc50.md) for
  a quantal endpoint, or an inhibition concentration through
  [`icp()`](https://open-aims.github.io/toxstats/reference/icp.md) for a
  continuous one.

  There is no `"point"` option, because
  [`lc50()`](https://open-aims.github.io/toxstats/reference/lc50.md) and
  [`icp()`](https://open-aims.github.io/toxstats/reference/icp.md)
  already provide that directly and are easier to find. The default is
  not `"both"` because the point branch runs a bootstrap for a
  continuous endpoint, which takes appreciably longer, so it is asked
  for rather than assumed.

- p:

  Percentages for the point branch, or `NULL` for the default of
  whichever method applies.

- nboot, seed:

  Passed to
  [`icp()`](https://open-aims.github.io/toxstats/reference/icp.md) when
  the point branch applies to a continuous endpoint.

- object:

  A `tox_test` object.

- row.names:

  Unused.

- optional:

  Unused.

## Value

An object of class `tox_test`, a list with elements `data`, `working`,
`decisions`, `assumptions`, `comparison`, `msd`, `pmsd`, `noec`, `loec`,
`excluded`, `overridden`, `transform`, `flowchart`, `alpha` and
`alpha_assumption`.

## Details

The flowchart is Figure 2 of the chronic manual, which is also Figure 13
of the acute manual. It branches on four things in order: whether the
response is a proportion, whether the pooled within-group residuals are
normal, whether the variances are homogeneous, and whether replication
is equal. The four terminals are
[`dunnett()`](https://open-aims.github.io/toxstats/reference/dunnett.md),
[`bonferroni_t()`](https://open-aims.github.io/toxstats/reference/bonferroni_t.md),
[`steel()`](https://open-aims.github.io/toxstats/reference/steel.md) and
[`wilcoxon_rank_sum()`](https://open-aims.github.io/toxstats/reference/wilcoxon_rank_sum.md).

Two significance levels are involved and they differ. The assumption
tests use `alpha_assumption`, which section 9.4.6.1 sets at 0.01; the
hypothesis test uses `alpha`, set at 0.05 and one-sided.

### What is recorded

The value of this function over calling the tests by hand is the
decision trail, returned by
[`decisions()`](https://open-aims.github.io/toxstats/reference/decisions.md).
Each row is a branch point: the question, the criterion, the statistic,
the answer, the consequence, and the manual section that justifies it.
[`summary()`](https://rdrr.io/r/base/summary.html) prints it above the
endpoints.

### Exclusions

Section 9.5.2 excludes some concentrations from the hypothesis test
while retaining them for point estimation. Two rules apply.

A concentration at which every replicate showed a complete response, for
instance total mortality, is excluded automatically for a quantal
endpoint.

A concentration at which survival was significantly reduced is excluded
from the analysis of a **sublethal** endpoint such as growth or
reproduction. That requires a separate survival analysis, so it cannot
be inferred here; pass the concentrations to `exclude`. The Appendix E
example does exactly this, dropping the 50 per cent concentration from
the reproduction analysis on the strength of the Appendix G survival
result.

### The percent minimum significant difference override

Section 10.2.8.2.5 states that a concentration shall not be declared
toxic if its relative difference from the control is less than the lower
bound published for the test method. This overrides the hypothesis test.
It is applied only when `pmsd_bounds` is supplied, and each
concentration it reverses is added to the decision trail.

### Forcing a test

Passing `test` runs that test whatever the flowchart selected. When the
two disagree a warning is signalled, `overridden` is set, and a row
recording the override is added to the trail, so a forced analysis can
never be printed as though the manual had chosen it.

## Methods (by generic)

- `summary(tox_test)`: Print the decision trail above the endpoints.

- `as.data.frame(tox_test)`: Return every endpoint as one tidy row.

## References

US EPA (2002) EPA-821-R-02-013, Figure 2 and sections 9.4 to 9.6.

## See also

[`decisions()`](https://open-aims.github.io/toxstats/reference/decisions.md)
for the trail,
[`msd()`](https://open-aims.github.io/toxstats/reference/msd.md) and
[`pmsd()`](https://open-aims.github.io/toxstats/reference/pmsd.md) for
sensitivity.

## Examples

``` r
# Appendix C growth data: normal, homogeneous and balanced, so the
# flowchart selects Dunnett's procedure.
fit <- tox_test(fathead_c1, response = "weight")
summary(fit)
#> EPA WET hypothesis test
#> Flowchart: EPA-821-R-02-013 Figure 2
#> 
#>    1  Is the response a proportion requiring transformation?
#>       no transformation needed
#>       (EPA-821-R-02-013 Appendix B, section 4.2)
#>    2  Are the pooled within-group residuals normally distributed?
#>       0.9507, p = 0.378 -> residuals consistent with normality
#>       (EPA-821-R-02-013 Appendix B, section 2.1)
#>    3  Are the variances homogeneous across concentrations?
#>       7.856, p = 0.097 -> variances not significantly different
#>       (EPA-821-R-02-013 Appendix B, section 3)
#>    4  Is replication equal across all concentrations?
#>       4, 4, 4, 4, 4 replicates, balanced
#>       (EPA-821-R-02-013 Figure 2)
#>    5  Which test was run?
#>       Dunnett's procedure
#>       (EPA-821-R-02-013 Appendix C)
#> 
#>   NOEC 128    LOEC 256
#>   MSD  0.1618    PMSD 23.9 per cent
```
