# Back-transform an arc sine square root value

Inverts
[`arcsine_sqrt()`](https://open-aims.github.io/toxstats/reference/arcsine_sqrt.md),
returning a proportion from a value in radians. The EPA manuals use this
to express a minimum significant difference on the original proportion
scale (Appendix C, section 1.11.2), where it is more readily interpreted
than an angle.

## Usage

``` r
inv_arcsine_sqrt(theta)
```

## Arguments

- theta:

  A numeric vector of transformed values in radians, each in
  `[0, pi / 2]`.

## Value

A numeric vector of proportions the same length as `theta`.

## Details

The back-transformation is `sin(theta)^2`. It does not undo the endpoint
adjustments made by
[`arcsine_sqrt()`](https://open-aims.github.io/toxstats/reference/arcsine_sqrt.md):
a proportion of 0 transformed with `n = 20` and then back-transformed
returns `1 / 80`, not 0. That is intended – the adjusted value is the
one the analysis actually used.

## References

US EPA (2002) EPA-821-R-02-013, Appendix C, section 1.11.2.

## See also

[`arcsine_sqrt()`](https://open-aims.github.io/toxstats/reference/arcsine_sqrt.md).

## Examples

``` r
inv_arcsine_sqrt(arcsine_sqrt(0.60))
#> [1] 0.6
```
