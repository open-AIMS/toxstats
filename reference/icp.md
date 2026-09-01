# Estimate an inhibition concentration by linear interpolation

Estimates the concentration causing a given percentage reduction in a
sublethal response such as growth or reproduction, by interpolating
between the two concentrations whose smoothed mean responses bracket it.
This is the EPA linear interpolation method, implemented in the manual
by the ICPIN program.

## Usage

``` r
icp(x, ..., p = 25, nboot = 200, seed = NULL, ci_level = 0.95)
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

- p:

  The percentage reductions to estimate, on a 1 to 99 scale.

- nboot:

  Number of bootstrap resamples used for inhibition concentration
  confidence limits.

- seed:

  Optional integer seed, set before any resampling so that results are
  reproducible. `NULL` leaves the random number stream untouched.

- ci_level:

  Confidence level for interval estimates. Defaults to 0.95.

## Value

An object of class `tox_estimate`, with an additional `boot` element
holding the resampled estimates, the count that could not be computed,
and the seed used.

## Details

### The estimate

The concentration means are smoothed to be monotone non-increasing, then

\$\$IC_p = C_J + \[M_1 (1 - p/100) - M_J\]\frac{C\_{J+1} -
C_J}{M\_{J+1} - M_J}\$\$

where \\C_J\\ and \\C\_{J+1}\\ bracket the target response and \\M_1\\
is the smoothed control mean. On the Appendix M worked example this
gives 8.5716 for the IC25 against the 8.5715 the manual's program
prints, and 10.893 for the IC50 against 10.89.

Smoothing is by pool-adjacent-violators with **equal weight on each
concentration mean**, not weighted by the number of replicates. That is
what reproduces the manual's smoothed means of 28.75 across the control
and the three lowest concentrations. With unequal replication the manual
does not say what it does, and this choice is a documented assumption.

Where the target response falls outside the range of the smoothed means
the estimate is reported as an inequality: greater than the highest
concentration, or less than the lowest.

### The interval

Ordinary interval methods do not apply to an interpolated estimate, so
the manual uses a bootstrap. Replicate values are resampled **with
replacement within each concentration**, the control included, and the
means are recomputed, **re-smoothed** and re-interpolated on every
iteration.

The re-smoothing is the detail most easily missed. The manual's step
list places smoothing before resampling, which reads as though it
happens once. Omitting it from the loop does not merely widen the
interval; it biases the estimate badly. On the Appendix M data the
resampled estimates then average about 10.4 against a point estimate of
8.57, with a standard deviation of 0.49 rather than 0.14, because
unsmoothed resampled means are not monotone and the interpolation
selects the wrong bracketing pair.

Limits are the empirical order statistics the manual describes, "the
second smallest and second largest" of 80 resamples, generalised to
`floor((nboot + 1) * (1 - ci_level) / 2)` in from each end.

**The interval is not exactly reproducible from the manual.** Its
printed limits of 8.3112 and 9.0418 come from 80 resamples drawn with
the seed -641671986 by a Turbo Pascal generator, so no other program can
return them. The estimate, the smoothed means and the bootstrap standard
deviation are reproducible, and are what the test suite checks.

The default of 200 resamples departs from the manual's 80. The manual
itself warns that "confidence limits based on the empirical quantiles of
a bootstrap distribution of 80 samples may be unstable", and permits up
to 1000.

### Expanded limits

The manual reports a second, wider interval when there are fewer than
seven replicates per concentration. Its definition appears only in the
ICPIN version 2.0 program documentation, which is not publicly
retrievable, so it is **not** implemented. A warning is raised when the
design would have triggered it.

## References

US EPA (2002) EPA-821-R-02-013, Appendix M.

Norberg-King TJ (1993) *A linear interpolation method for sublethal
toxicity: the inhibition concentration (ICp) approach*, version 2.0.
EPA/600/M-91/037.

Efron B (1982) *The Jackknife, the Bootstrap, and Other Resampling
Plans*. Society for Industrial and Applied Mathematics, Philadelphia.

## Examples

``` r
# Appendix M worked example; the manual gives IC25 = 8.5715 and
# IC50 = 10.89.
icp(ceriodaphnia_m1, response = "young", p = c(25, 50), seed = 42)
#> Linear interpolation method
#>   EPA-821-R-02-013 Appendix M
#> 
#>  endpoint  p estimate   lower   upper                  ci_method bound
#>      IC25 25   8.5715  8.3009  8.8863 bootstrap order statistics  <NA>
#>      IC50 50  10.8930 10.3520 11.5230 bootstrap order statistics  <NA>
#>  boot_mean boot_sd n_undefined
#>     8.5779 0.14997           0
#>    10.9060 0.29994           0
```
