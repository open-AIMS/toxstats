# Abbott's correction for control response

Adjusts observed proportions for the response already present in the
control, so that the adjusted values estimate the proportion responding
to the toxicant alone. The EPA point-estimation methods apply this
before estimating a median lethal concentration when control mortality
is non-zero.

## Usage

``` r
abbott(p, p_control, clamp = TRUE)
```

## Arguments

- p:

  A numeric vector of observed proportions, each in `[0, 1]`.

- p_control:

  The control proportion, a single number in `[0, 1)`.

- clamp:

  Should adjusted values outside `[0, 1]` be clamped to that range? A
  flag, defaulting to `TRUE`.

## Value

A numeric vector of adjusted proportions the same length as `p`.

## Details

The correction is

\$\$p' = \frac{p - p_0}{1 - p_0}\$\$

where \\p_0\\ is the control proportion. It is undefined when \\p_0 =
1\\, which is an error rather than a warning: a control in which every
organism responded carries no information about the toxicant.

Where an observed proportion falls below the control, the raw correction
is negative. By default such values are clamped to 0, which is what the
point-estimation methods require. Set `clamp = FALSE` to see the raw
values, for instance when diagnosing a control problem.

## References

Abbott WS (1925) A method of computing the effectiveness of an
insecticide. *Journal of Economic Entomology* 18:265-267.

## Examples

``` r
abbott(c(0.02, 0.05, 0.80), p_control = 0.02)
#> [1] 0.00000000 0.03061224 0.79591837
```
