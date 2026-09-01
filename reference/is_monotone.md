# Is a sequence monotone?

Reports whether a sequence already follows the required direction, and
where it first departs from it. Williams' test assumes monotonicity, and
the EPA linear interpolation method smooths towards it, so both need to
know.

## Usage

``` r
is_monotone(y, direction = c("decreasing", "increasing"))
```

## Arguments

- y:

  A numeric vector to smooth.

- direction:

  The direction the smoothed sequence must follow. `"decreasing"` (the
  default) returns a non-increasing sequence, as expected for a response
  reduced by a toxicant. `"increasing"` returns a non-decreasing
  sequence, as expected for a proportion responding.

## Value

A list with elements `monotone` (a flag), `direction` (the direction
tested) and `violations` (an integer vector of the positions at which
the sequence departs from the required direction; empty when `monotone`
is `TRUE`).

## Examples

``` r
is_monotone(c(10, 8, 9, 4), direction = "decreasing")
#> $monotone
#> [1] FALSE
#> 
#> $direction
#> [1] "decreasing"
#> 
#> $violations
#> [1] 3
#> 
```
