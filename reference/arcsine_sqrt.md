# Arc sine square root transformation

Transforms proportion data as specified in Appendix B, section 4.2 of
the EPA WET method manuals. The manuals require this transformation
before the parametric hypothesis tests are applied to a quantal endpoint
such as survival, because proportions have a variance that depends on
the mean.

## Usage

``` r
arcsine_sqrt(p, n = NULL)
```

## Arguments

- p:

  A numeric vector of proportions, each in `[0, 1]`.

- n:

  The number of organisms in each replicate, used only for the endpoint
  adjustments. Either a single number or a vector the same length as
  `p`. Required if any element of `p` is exactly 0 or exactly 1.

## Value

A numeric vector of transformed values, in radians, the same length as
`p`.

## Details

The transformation is `asin(sqrt(p))`, returned in **radians**, with two
endpoint adjustments that the manuals treat as part of the method:

- an observed proportion of exactly 0 is replaced by `1 / (4 * n)`;

- an observed proportion of exactly 1 is replaced by `1 - 1 / (4 * n)`;

where `n` is the number of organisms in that replicate. Without them a
replicate at 0 or 1 would contribute no variance, and the tests
downstream would understate the true variability. The adjustments follow
Bartlett (1937). Omitting them is the commonest way to fail to reproduce
a published EPA analysis, so `n` is required whenever any proportion is
0 or 1.

## References

US EPA (2002) EPA-821-R-02-013, Appendix B, section 4.2.

Bartlett MS (1937) Some examples of statistical methods of research in
agriculture and applied biology. *Supplement to the Journal of the Royal
Statistical Society* 4:137-183.

## See also

[`inv_arcsine_sqrt()`](https://open-aims.github.io/toxstats/reference/inv_arcsine_sqrt.md)
for the back-transformation.

## Examples

``` r
# Appendix B, section 4.2 worked values
arcsine_sqrt(0.60)
#> [1] 0.8860771
arcsine_sqrt(c(0, 1), n = 20)
#> [1] 0.1120376 1.4587587
```
