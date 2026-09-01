# Estimate the LC50 by the Spearman-Karber method

Estimates the mean of the distribution of the log10 tolerance, which is
the median when that distribution is symmetric. The EPA flowchart
selects this method when partial responses occur but the probit model
does not fit, and the response runs from zero at the lowest
concentration to complete at the highest.

## Usage

``` r
spearman_karber(x, ..., ci_level = 0.95)
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

- ci_level:

  Confidence level for interval estimates. Defaults to 0.95.

## Value

An object of class `tox_estimate`.

## Details

After smoothing and Abbott adjustment,

\$\$m = \sum\_{i=1}^{k-1} (p^a\_{i+1} - p^a_i)\\(X_i + X\_{i+1})/2\$\$

with \\X_i\\ the log10 concentration, and

\$\$V(m) = \sum\_{i=2}^{k-1} \frac{p^a_i (1 - p^a_i)(X\_{i+1} -
X\_{i-1})^2}{4(n_i - 1)}\$\$

The interval is `m ± 2 sqrt(V(m))` on the log scale, back-transformed.
The manual writes the multiplier as 2.0 rather than 1.96; that rounding
is kept because it is what reproduces the printed interval.

The method requires the adjusted proportion to be zero at the lowest
concentration and one at the highest. When it is not, the flowchart
directs the analysis to
[`trimmed_spearman_karber()`](https://open-aims.github.io/toxstats/reference/trimmed_spearman_karber.md).

## References

US EPA (2002) EPA-821-R-02-012, section 11.2.3.

Finney DJ (1978) *Statistical Method in Biological Assay*, 3rd edition.
Charles Griffin, London.

## Examples

``` r
# Acute Table 20; the manual gives m = 1.656527, V(m) = 0.0010977 and an
# LC50 of 45.3 per cent with limits 38.9 and 52.8.
spearman_karber(
  acute_table20,
  response = "spearman_karber", n_exposed = "exposed", type = "quantal"
)
#> Spearman-Karber method
#>   EPA-821-R-02-012 section 11.2.3
#> 
#>  endpoint  p estimate  lower  upper                 ci_method
#>      LC50 50   45.344 38.928 52.817 normal on the log10 scale
```
