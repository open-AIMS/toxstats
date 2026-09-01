# Ceriodaphnia dubia reproduction, EPA Table B.7

Reproduction data from a *Ceriodaphnia dubia* Survival and Reproduction
Test, used for the Kolmogorov "D" worked example in Appendix B of the
EPA chronic freshwater manual. With 60 observations the manual applies
the Kolmogorov statistic rather than Shapiro-Wilk.

## Usage

``` r
ceriodaphnia_b7
```

## Format

A data frame with 60 rows and 3 columns:

- conc:

  effluent concentration, per cent

- replicate:

  replicate label, 1 to 10

- young:

  number of young produced per female

## Source

US EPA (2002) EPA-821-R-02-013, Table B.7.

## Details

The manual reports `D+ = 0.0525`, `D- = 0.0597`, `D = 0.0597` and a
modified statistic `D* = 0.4684`, against a critical value of 1.035 at
`alpha = 0.01`, so the data are normally distributed.
