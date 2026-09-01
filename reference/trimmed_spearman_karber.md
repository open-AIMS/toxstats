# Estimate the LC50 by the trimmed Spearman-Karber method

Estimates a trimmed mean of the log10 tolerance distribution, so that
the estimate can be formed when the response does not run all the way
from zero to complete. The EPA flowchart selects this method when the
probit model does not fit and the untrimmed Spearman-Karber requirements
are not met.

## Usage

``` r
trimmed_spearman_karber(x, ..., trim = NULL, ci_level = 0.95)
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

- trim:

  The proportion to trim from each tail, or `NULL` for the automatic
  trim defined above.

- ci_level:

  Confidence level for interval estimates. Defaults to 0.95.

## Value

An object of class `tox_estimate`.

## Details

### The trim

Section 11.2.4.3 step 4 defines it as

\$\$\mathrm{trim} = \max(p^a_1,\\ 1 - p^a_k)\$\$

where the proportions are the smoothed, Abbott-adjusted ones. The
ordering matters: smoothing first, then the Abbott adjustment, then the
trim. Computing it from the raw proportions gives a different answer.
This definition reproduces the 20.51 per cent printed in the acute
manual and the 20.41 per cent printed in Appendix K of the chronic
manual.

`trim = NULL`, the default, is the automatic trim ToxCalc advertised.

### The estimate

The concentrations at which the adjusted response crosses the trim and
its complement are found by interpolation on the log10 scale, the
interior points are rescaled to run from zero to one, and the
Spearman-Karber formula is applied to that set. On the acute Table 20
data this gives 77.1105 against the 77.11 the manual's program prints.

### The interval

The manual delegates the interval to a program whose source it does not
supply, and it states neither the variance formula nor the multiplier.
Two choices together recover the printed limits, and both were
established by testing candidates against the printed output rather than
assumed:

- the variance is the delta method applied analytically to the
  estimator, **holding the trim fixed** even when it was chosen from the
  data. This is the variance of Hamilton et al. (1977). Differentiating
  through the choice of trim as well was tried, and does not reproduce
  the published limits.

- the multiplier is 2.0, the same rounding of 1.96 that section 11.2.3.3
  uses for the untrimmed method, rather than 1.96 itself.

On the acute Table 20 data this gives 69.74 to 85.26, which is what the
manual prints, to every digit printed.

The derivation is set out with the internal `tsk_variance()`. It is
derived here rather than transcribed: Hamilton's corrected expression is
not retrievable from any public source found, and the one implementation
of it that could be located, in the `ecotoxicology` package, is under a
licence this package cannot absorb. That implementation is used instead
as an independent check, and the two agree to machine precision on every
configuration but one; see the note on `V6` in the package vignette.

## References

US EPA (2002) EPA-821-R-02-012, section 11.2.4.

Hamilton MA, Russo RC, Thurston RV (1977) Trimmed Spearman-Karber method
for estimating median lethal concentrations in toxicity bioassays.
*Environmental Science & Technology* 11:714-719, with the correction at
12:417 (1978).

## Examples

``` r
# Acute Table 20; the manual's program prints a trim of 20.51 per cent and
# an LC50 of 77.11.
trimmed_spearman_karber(
  acute_table20,
  response = "trimmed", n_exposed = "exposed", type = "quantal"
)
#> Trimmed Spearman-Karber method (automatic trim 20.51 per cent)
#>   EPA-821-R-02-012 section 11.2.4
#> 
#>  endpoint  p estimate  lower  upper                            ci_method
#>      LC50 50   77.111 69.743 85.256 Hamilton variance on the log10 scale
```
