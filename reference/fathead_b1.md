# Fathead minnow larval growth, EPA Table B.1

Growth data from a Fathead Minnow Larval Survival and Growth Test, used
for the Shapiro-Wilk and Bartlett worked examples in Appendix B of the
EPA chronic freshwater manual.

## Usage

``` r
fathead_b1
```

## Format

A data frame with 20 rows and 3 columns:

- conc:

  sodium pentachlorophenate concentration, micrograms per litre

- replicate:

  replicate label, A to D

- weight:

  larval dry weight, milligrams

## Source

US EPA (2002) EPA-821-R-02-013, Table B.1.

## Details

The manual reports the following results for these data, which the
package test suite checks against:

- concentration means 0.714, 0.674, 0.677, 0.624 and 0.580 mg;

- Shapiro-Wilk `W = 0.959`, against a critical value of 0.868 at
  `alpha = 0.01` with 20 observations, so the data are normally
  distributed;

- Bartlett `B = 7.691`, against a critical value of 13.277, so the
  variances are not different.

The printed Bartlett statistic is not reproducible from these data,
which give 6.836. See
[`epa_variance()`](https://open-aims.github.io/toxstats/reference/epa_variance.md)
for why.

Note that Appendix C states its data are "the same data used in
Appendices B and D". They are not; see
[fathead_c1](https://open-aims.github.io/toxstats/reference/fathead_c1.md).
