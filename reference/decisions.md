# Extract the decision trail

Returns the record of how a test was selected: every branch point of the
EPA flowchart, the statistic that answered it, the consequence, and the
section of the manual that justifies it.

## Usage

``` r
decisions(x, ...)

# S3 method for class 'tox_test'
decisions(x, ...)

# S3 method for class 'tox_lc50'
decisions(x, ...)
```

## Arguments

- x:

  An object carrying a decision trail, such as the result of
  [`tox_test()`](https://open-aims.github.io/toxstats/reference/tox_test.md).

- ...:

  Unused.

## Value

A data frame with columns `step`, `node`, `question`, `criterion`,
`statistic`, `p_value`, `answer`, `outcome` and `reference`.

## Methods (by class)

- `decisions(tox_test)`: Decision trail from a `tox_test` analysis.

- `decisions(tox_lc50)`: Decision trail from an
  [`lc50()`](https://open-aims.github.io/toxstats/reference/lc50.md)
  analysis.

## Examples

``` r
decisions(tox_test(fathead_c1, response = "weight"))
#>   step          node
#> 1    1     transform
#> 2    2     normality
#> 3    3      variance
#> 4    4 replication_p
#> 5    5          test
#>                                                      question
#> 1      Is the response a proportion requiring transformation?
#> 2 Are the pooled within-group residuals normally distributed?
#> 3        Are the variances homogeneous across concentrations?
#> 4             Is replication equal across all concentrations?
#> 5                                         Which test was run?
#>                                  criterion statistic    p_value  answer
#> 1       A quantal endpoint is a proportion        NA         NA      no
#> 2                   W, reject if p <= 0.01 0.9506998 0.37783763     yes
#> 3          Bartlett B, reject if p <= 0.01 7.8564131 0.09698185     yes
#> 4 equal replication at every concentration        NA         NA     yes
#> 5                                     <NA>        NA         NA dunnett
#>                                 outcome
#> 1              no transformation needed
#> 2   residuals consistent with normality
#> 3 variances not significantly different
#> 4    4, 4, 4, 4, 4 replicates, balanced
#> 5                   Dunnett's procedure
#>                                  reference
#> 1 EPA-821-R-02-013 Appendix B, section 4.2
#> 2 EPA-821-R-02-013 Appendix B, section 2.1
#> 3   EPA-821-R-02-013 Appendix B, section 3
#> 4                EPA-821-R-02-013 Figure 2
#> 5              EPA-821-R-02-013 Appendix C
```
