# Ceriodaphnia dubia mortality, EPA Table G.2

Survival data from a *Ceriodaphnia dubia* survival and reproduction
test, used for the Fisher's Exact Test worked example in Appendix G of
the EPA chronic freshwater manual. One row per concentration, already
pooled across replicates, as the manual presents it.

## Usage

``` r
ceriodaphnia_g2
```

## Format

A data frame with 6 rows and 3 columns:

- conc:

  effluent concentration, per cent

- dead:

  number of adults dead at the end of the test

- exposed:

  number of live adults at the start of the test

## Source

US EPA (2002) EPA-821-R-02-013, Table G.2.

## Details

`exposed` is the number of live adults at the beginning of the test,
which is nine in the control rather than ten.

The manual concludes that only the 25 per cent concentration differs
significantly from the control, giving a NOEC of 12 per cent and a LOEC
of 25 per cent for survival. That result is what excludes the 50 per
cent concentration from the reproduction analysis of
[ceriodaphnia_e1](https://open-aims.github.io/toxstats/reference/ceriodaphnia_e1.md).

Appendix G applies a pairwise error rate of 0.05 rather than an
experiment-wise one, because Fisher's Exact Test is itself conservative.
