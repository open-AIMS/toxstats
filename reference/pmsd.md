# Percent minimum significant difference

The minimum significant difference expressed as a percentage of the
control mean, which is how the EPA method-variability criteria are
stated. A value outside the bounds published for a test method indicates
a test that was either unusually insensitive or unusually sensitive.

## Usage

``` r
pmsd(
  x,
  ...,
  alpha = 0.05,
  test = c("dunnett", "bonferroni", "sidak"),
  bounds = NULL
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

- test:

  The multiple-comparison procedure whose critical value is used.
  `"dunnett"` is the EPA method for a balanced design; `"bonferroni"`
  and `"sidak"` correspond to the t test with Bonferroni's or
  Dunn-Sidak's adjustment, used when replication is unequal.

- bounds:

  Optional. Either a row name of
  [epa_pmsd_bounds](https://open-aims.github.io/toxstats/reference/epa_pmsd_bounds.md)
  such as `"fathead_growth"`, or a length-two numeric vector giving the
  lower and upper bounds directly. `NULL` reports the value without
  comparison.

## Value

An object of class `tox_pmsd`, a list with elements `pmsd` (named by
concentration), `control_mean`, `bounds`, `status` and the underlying
`msd` object.

## Details

`PMSD = 100 * MSD / control mean`.

EPA publishes lower and upper bounds for three sublethal endpoints,
shipped here as
[epa_pmsd_bounds](https://open-aims.github.io/toxstats/reference/epa_pmsd_bounds.md).
Pass `bounds` to compare against one of them.

The manuals require the percent minimum significant difference to be
computed parametrically **even when the flowchart selected a
non-parametric test**, so this function does not consult the flowchart.

## References

US EPA (2002) EPA-821-R-02-013, section 10.2.8 and Table 6.

## Examples

``` r
# Appendix C worked example; the manual reports about 24 per cent
pmsd(fathead_c1, response = "weight", bounds = "fathead_growth")
#> Percent minimum significant difference
#>   EPA-821-R-02-013 section 10.2.8
#> 
#>   control mean = 0.6773
#> 
#>   32   64  128  256 
#> 23.9 23.9 23.9 23.9 
#> 
#>   EPA bounds: 12 to 30 per cent
#>       32       64      128      256 
#> "within" "within" "within" "within" 
```
