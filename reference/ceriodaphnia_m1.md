# Ceriodaphnia dubia reproduction, EPA Table M.1

Reproduction data used for the linear interpolation worked example in
Appendix M of the EPA chronic freshwater manual.

## Usage

``` r
ceriodaphnia_m1
```

## Format

A data frame with 60 rows and 3 columns:

- conc:

  effluent concentration, per cent

- replicate:

  replicate label, 1 to 10

- young:

  number of young produced

## Source

US EPA (2002) EPA-821-R-02-013, Table M.1.

## Details

This shares four columns with
[ceriodaphnia_b7](https://open-aims.github.io/toxstats/reference/ceriodaphnia_b7.md)
but is not the same dataset. Table M.1 has no counterpart to Table B.7's
12.5 per cent column, its own 12.5 per cent column holds what B.7
records at 25 per cent, and it adds a concentration at which
reproduction stopped entirely.

The manual reports, for these data:

- concentration means 22.4, 26.3, 34.6, 31.7, 9.4 and 0;

- smoothed means of 28.75 across the control and the three lowest
  concentrations, then 9.4 and 0;

- `IC25 = 8.5715` and `IC50 = 10.89`.

The ICPIN program output also prints standard deviations of 6.931,
8.001, 4.835, 2.946, 3.893 and 0, a bootstrap mean of 8.5891 with a
standard deviation of 0.1831 from 80 resamples, and confidence limits of
8.3112 and 9.0418. Those limits were drawn with the seed -641671986 by a
Turbo Pascal generator and cannot be reproduced by any other program;
see [`icp()`](https://open-aims.github.io/toxstats/reference/icp.md).
