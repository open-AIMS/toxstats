# Ceriodaphnia dubia reproduction, EPA Table E.1

Reproduction data from a *Ceriodaphnia dubia* seven-day chronic test,
used for the Steel's Many-One Rank Test worked example in Appendix E of
the EPA chronic freshwater manual, and, with two values removed, for the
Wilcoxon Rank Sum worked example in Appendix F.

## Usage

``` r
ceriodaphnia_e1
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

US EPA (2002) EPA-821-R-02-013, Table E.1.

## Details

The 50 per cent concentration is included here but is **excluded from
the reproduction analysis** in the manual, because survival at that
concentration was significantly reduced (section 9.5.2; see
[ceriodaphnia_g2](https://open-aims.github.io/toxstats/reference/ceriodaphnia_g2.md)
and
[`fisher_exact()`](https://open-aims.github.io/toxstats/reference/fisher_exact.md)).
Subset to `conc < 50` to reproduce the worked example.

Appendix F uses the same data with two males presumed to have occurred,
one in the control and one at 12 per cent, giving unequal replication.
Remove the control replicate 1 and the 12 per cent replicate 5 to
reproduce it.

The manual reports the following for the Appendix E example:

- rank sums of 84, 64, 76 and 55 against a critical rank sum of 76;

- NOEC 3 per cent and LOEC 6 per cent.

The printed rank sum of 64 at 6 per cent is a slip. The ranks the manual
itself lists in Table E.3 for that concentration sum to 63.5. The
conclusion is unaffected, both being at or below the critical value.

For the Appendix F example the manual reports rank sums of 79, 57, 58
and 55, critical values of 72 for the ten-replicate concentrations and
60 for the nine-replicate one, and the same NOEC and LOEC.
