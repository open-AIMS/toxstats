
<!-- README.md is generated from README.Rmd. Please edit that file -->

# toxcalc

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![license](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-lightgrey.svg)](https://choosealicense.com/)
<!-- badges: end -->

## Overview

`toxcalc` is an open-source implementation of the statistical methods
specified in the US EPA Whole Effluent Toxicity (WET) method manuals. It
recreates the capability of ToxCalc v5.0 (Tidepool Scientific),
commercial software that has been withdrawn — the vendor product page
now returns a 404.

The package provides both branches of the EPA analysis:

- the **hypothesis-testing** flowchart, giving the no-observed-effect
  concentration (NOEC), the lowest-observed-effect concentration (LOEC),
  the minimum significant difference (MSD) and the percent MSD (PMSD);
- the **point-estimation** flowchart, giving the median lethal
  concentration (LC50) by the graphical, Spearman-Karber, trimmed
  Spearman-Karber or probit method, and the inhibition concentration
  (ICp) by linear interpolation with bootstrap confidence limits.

## Why this exists

The EPA manuals specify not only the statistical tests but the decision
rules that select between them. Those rules are the part that is hardest
to reproduce by hand and easiest to leave undocumented.

Every `toxcalc` analysis therefore carries an **audit trail**: a table
recording each branch point of the flowchart — the question asked, the
statistic, the answer, the consequence, and the section of the EPA
manual that justifies it.

``` r
summary(tox_test(fathead_c1, response = "weight",
                pmsd_bounds = "fathead_growth"))
#> EPA WET hypothesis test
#> Flowchart: EPA-821-R-02-013 Figure 2
#> 
#>    1  Is the response a proportion requiring transformation?
#>       no transformation needed
#>       (EPA-821-R-02-013 Appendix B, section 4.2)
#>    2  Are the pooled within-group residuals normally distributed?
#>       0.9507, p = 0.378 -> residuals consistent with normality
#>       (EPA-821-R-02-013 Appendix B, section 2.1)
#>    3  Are the variances homogeneous across concentrations?
#>       7.856, p = 0.097 -> variances not significantly different
#>       (EPA-821-R-02-013 Appendix B, section 3)
#>    4  Is replication equal across all concentrations?
#>       4, 4, 4, 4, 4 replicates, balanced
#>       (EPA-821-R-02-013 Figure 2)
#>    5  Which test was run?
#>       Dunnett's procedure
#>       (EPA-821-R-02-013 Appendix C)
#> 
#>   NOEC 128    LOEC 256
#>   MSD  0.1618    PMSD 23.9 per cent (EPA bounds 12 to 30: within)
```

Those are the growth data from Appendix C of the chronic manual, and the
flowchart arrives at Dunnett’s procedure, the NOEC of 128 and the LOEC
of 256 that the manual itself reports.

## Running both branches

A laboratory report usually carries a no-effect concentration and a
point estimate together. `branch = "both"` returns them from one call,
and `as.data.frame()` gives one tidy row per endpoint.

``` r
fit <- tox_test(ceriodaphnia_m1, response = "young",
               branch = "both", seed = 42)
as.data.frame(fit)[, c("endpoint", "value", "lower", "upper", "method")]
#>   endpoint     value    lower    upper method
#> 1     NOEC  6.250000       NA       NA  steel
#> 2     LOEC 12.500000       NA       NA  steel
#> 3      MSD  5.274541       NA       NA  steel
#> 4     PMSD 23.547058       NA       NA  steel
#> 5     IC25  8.571544 8.300892 8.886289    icp
```

The point estimate uses all the data, including any concentration
excluded from the no-effect concentration, which is what section 9.5.2
requires.

## Point estimation

`lc50()` walks the Figure 6 chart of the acute manual. Each of the four
columns of its Table 20 was constructed to exercise one of the four
methods, and each routes to the right one:

``` r
for (column in c("graphical", "spearman_karber", "trimmed", "probit")) {
  chosen <- lc50(acute_table20, response = column,
                 n_exposed = "exposed", type = "quantal")
  cat(format(column, width = 16), "->", chosen$selected, "\n")
}
#> graphical        -> graphical_lc50 
#> spearman_karber  -> spearman_karber 
#> trimmed          -> trimmed_spearman_karber 
#> probit           -> probit_lc
```

## Fidelity

The target is faithfulness to the EPA method manuals, validated against
the worked examples printed in them. It is **not** a bit-for-bit clone
of ToxCalc; no ToxCalc output was available for comparison.

Almost every published value reproduces exactly. The exceptions are
documented rather than smoothed over, and each is reconciled by a test:

- **seven printed values do not follow from the manuals’ own data**, all
  through hand-calculation rounding or arithmetic slips, and one dataset
  is not the dataset its own text says it is;
- **the ICp interval** cannot be matched, because the program that
  produced it names the random seed it used and nothing else can
  reproduce that draw.

That is now the only one. The trimmed Spearman-Karber interval was the
other, the manual having delegated it to a program whose variance
formula it never states; it is reproduced exactly since the formula was
rederived. The vignette records what it took.

The vignette `vignettes/recreating-toxcalc.qmd` records every point at
which the source material was ambiguous, internally inconsistent, or at
odds with modern practice, and what was decided in each case, with the
evidence.

`williams()` implements the one ToxCalc feature that is not an EPA
method. Neither manual gives a worked example or a table of critical
values for it, so the values are simulated and the implementation is
validated against mathematical identities instead. It is labelled as an
extension throughout, and neither flowchart will ever select it.

## Status

**Experimental, under active development.** Both branches are complete
and `R CMD check` is clean. Still to come: the marine manual’s
endpoints, designs with more than one control, and the remaining
alternative tests.

See `notes/TOXSTATS-human.md` for the plan and
`notes/TOXSTATS-claude.md` for the specification.

## Installation

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("beckyfisher/toxstats")
```

## Standing and licence

This package is **not** produced, endorsed, validated or reviewed by the
US Environmental Protection Agency, and it is not affiliated with
Tidepool Scientific or with the ToxCalc software it recreates. It
implements the methods described in the EPA method manuals, which are US
Government works, and it is offered without warranty of any kind. A
laboratory remains responsible for satisfying itself, and its regulator,
that an analysis is fit for the use it is put to.

No table of critical values that the EPA manuals reproduce from a
third-party source is included here. Those values are computed instead;
the reasoning is in the vignette.

## References

- US EPA (2002) *Methods for Measuring the Acute Toxicity of Effluents
  and Receiving Waters to Freshwater and Marine Organisms*, 5th ed.
  EPA-821-R-02-012.
- US EPA (2002) *Short-term Methods for Estimating the Chronic Toxicity
  of Effluents and Receiving Waters to Freshwater Organisms*, 4th ed.
  EPA-821-R-02-013.
- US EPA (2002) *Short-term Methods for Estimating the Chronic Toxicity
  of Effluents and Receiving Waters to Marine and Estuarine Organisms*,
  3rd ed. EPA-821-R-02-014.
- Hamilton MA, Russo RC, Thurston RV (1977) Trimmed Spearman-Karber
  method for estimating median lethal concentrations in toxicity
  bioassays. *Environmental Science & Technology* 11:714-719.
- Norberg-King TJ (1993) *A linear interpolation method for sublethal
  toxicity: the inhibition concentration (ICp) approach*, version 2.0.
  EPA/600/M-91/037.
- Williams DA (1971) A test for differences between treatment means when
  several dose levels are compared with a zero dose control.
  *Biometrics* 27:103-117.
