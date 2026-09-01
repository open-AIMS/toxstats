# Estimate the LC50 by the graphical method

Interpolates between the two concentrations that bracket 50 per cent
response, on a logarithmic concentration scale. The EPA flowchart
selects this method when there are no partial responses, so no model can
be fitted.

## Usage

``` r
graphical_lc50(x, ...)
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

## Value

An object of class `tox_estimate`.

## Details

The proportions are smoothed and adjusted for control response first,
then the estimate is read off the straight line joining the two
bracketing points on semi-logarithmic axes, which is linear
interpolation in `log10(conc)`.

No confidence interval is available. The manual gives none, and with no
partial response there is no information from which to construct one.

## References

US EPA (2002) EPA-821-R-02-012, section 11.2.2.

## Examples

``` r
# Acute Table 20, all-or-nothing column; the manual reads 35 per cent off
# the plot, and interpolation gives 35.4.
graphical_lc50(
  acute_table20,
  response = "graphical", n_exposed = "exposed", type = "quantal"
)
#> Graphical method
#>   EPA-821-R-02-012 section 11.2.2
#> 
#>  endpoint  p estimate lower upper ci_method
#>      LC50 50   35.355    NA    NA      <NA>
```
