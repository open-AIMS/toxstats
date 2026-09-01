# EPA variability criteria for sublethal endpoints

Lower and upper bounds on the percent minimum significant difference for
the three sublethal endpoints EPA publishes criteria for. A test whose
percent minimum significant difference falls outside these bounds was
either unusually insensitive or unusually sensitive.

## Usage

``` r
epa_pmsd_bounds
```

## Format

A data frame with 3 rows and 6 columns:

- id:

  short identifier, accepted by the `bounds` argument of
  [`pmsd()`](https://open-aims.github.io/toxstats/reference/pmsd.md)

- method:

  EPA method number

- test:

  full test method name

- endpoint:

  the sublethal endpoint the bounds apply to

- lower:

  lower percent minimum significant difference bound

- upper:

  upper percent minimum significant difference bound

## Source

US EPA (2002) EPA-821-R-02-013, Table 6.

## Details

EPA derived the bounds from the 10th and 90th percentiles of the percent
minimum significant difference recorded in its own WET Interlaboratory
Variability Study, so this is an EPA-generated table rather than a
reproduction of a third-party one.

Section 10.2.8.2.5 of the manual attaches a decision rule to the lower
bound: a concentration is not to be declared toxic if its relative
difference from the control is less than the lower bound, whatever the
hypothesis test concluded.
