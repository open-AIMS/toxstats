# LC50 worked-example mortality data, EPA acute Table 20

Mortality counts used for the four median lethal concentration worked
examples in section 11.2 of the EPA acute manual. Twenty organisms were
exposed in the control and at every concentration.

## Usage

``` r
acute_table20
```

## Format

A data frame with 6 rows and 6 columns:

- conc:

  effluent concentration, per cent

- exposed:

  number of organisms exposed

- graphical:

  number dead, for the graphical method example

- spearman_karber:

  number dead, for the Spearman-Karber example

- trimmed:

  number dead, for the trimmed Spearman-Karber example

- probit:

  number dead, for the probit example

## Source

US EPA (2002) EPA-821-R-02-012, Table 20.

## Details

Each method column is a different set of counts, chosen so that the
Figure 6 flowchart routes it to that method. Together they exercise all
four terminals of the chart, which is how
[`lc50()`](https://open-aims.github.io/toxstats/reference/lc50.md) is
validated.

The published results are:

- `graphical`: LC50 read off the plot as 35 per cent; interpolating on
  the logarithmic scale gives 35.36. No confidence interval.

- `spearman_karber`: `m = 1.656527`, `V(m) = 0.0010977`, LC50 45.3 per
  cent with limits 38.9 and 52.8.

- `trimmed`: an automatic trim of 20.51 per cent and an LC50 of 77.11
  with limits 69.74 and 85.26. The limits here differ slightly; see
  [`trimmed_spearman_karber()`](https://open-aims.github.io/toxstats/reference/trimmed_spearman_karber.md).

- `probit`: a Pearson chi-square for heterogeneity of 3.076 on 3 degrees
  of freedom against a tabular 7.815, LC50 22.872 with limits 18.787 and
  27.846, and LC1 7.924 with limits 4.147 and 10.959.

The probit column has no control mortality, which is why the manual's
printed output shows its Abbott-adjusted proportions equal to the
observed ones.
