# Fathead minnow larval growth, EPA Table C.1

Growth data from a Fathead Minnow Larval Survival and Growth Test, used
for the Dunnett's procedure worked example in Appendix C of the EPA
chronic freshwater manual.

## Usage

``` r
fathead_c1
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

US EPA (2002) EPA-821-R-02-013, Table C.1.

## Details

Despite the manual stating that Appendix C uses "the same data used in
Appendices B and D", Table C.1 differs from Table B.1. Both are shipped
separately for that reason.

The manual reports the following results for these data:

- concentration means 0.677, 0.575, 0.660, 0.565 and 0.454 mg;

- within mean square 0.0094, so `Sw = 0.097` on 15 degrees of freedom;

- `t` values 1.487, 0.248, 1.633 and 3.251;

- a one-sided Dunnett critical value of 2.36;

- `MSD = 0.162`, about 24 per cent of the control mean;

- NOEC 128 and LOEC 256 micrograms per litre.

Section 1.11 of Appendix C describes the minimum significant difference
as 0.087 mg. That is a typographical error: 0.162 is the value the same
section computes, and 24 per cent of the 0.677 control mean is 0.162.
