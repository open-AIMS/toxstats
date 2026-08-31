# toxcalc

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![license](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-lightgrey.svg)](https://choosealicense.com/)
<!-- badges: end -->

## Overview

`toxcalc` is an open-source implementation of the statistical methods specified
in the US EPA Whole Effluent Toxicity (WET) method manuals. It recreates the
capability of ToxCalc v5.0 (Tidepool Scientific), commercial software that has
been withdrawn — the vendor product page now returns a 404.

The package provides both branches of the EPA analysis:

- the **hypothesis-testing** flowchart, giving the no-observed-effect
  concentration (NOEC), the lowest-observed-effect concentration (LOEC), the
  minimum significant difference (MSD) and the percent MSD (PMSD);
- the **point-estimation** flowchart, giving the median lethal concentration
  (LC50) by the graphical, Spearman-Karber, trimmed Spearman-Karber or probit
  method, and the inhibition concentration (ICp) by linear interpolation with
  bootstrap confidence limits.

## Why this exists

The EPA manuals specify not only the statistical tests but the decision rules
that select between them. Those rules are the part that is hardest to reproduce
by hand and easiest to leave undocumented.

Every `toxcalc` analysis therefore carries an **audit trail**: a table recording
each branch point of the flowchart — the question asked, the statistic, the
answer, the consequence, and the section of the EPA manual that justifies it.
The printed output shows the reasoning above the result, so a reviewer can check
the analysis against the manual line by line.

Individually, some of the required tests already exist in R — Steel's Many-One
Rank Test in `kSamples`, the assumption tests in `stats`. Others exist nowhere:
Williams' test, the trimmed Spearman-Karber method with automatic trim, the EPA
arc sine square root transformation with its endpoint adjustments, and the EPA
linear-interpolation ICp method are all implemented here from the primary
sources.

## Status

**Experimental, under active development.** The hypothesis-testing branch is
complete and usable. `toxcalc()` walks the EPA flowchart, selects the test the
manual directs, runs it, and returns the NOEC, LOEC, MSD and PMSD together with
the decision trail that justifies the selection.

```r
library(toxcalc)
summary(toxcalc(fathead_c1, response = "weight",
                pmsd_bounds = "fathead_growth"))
```

```
EPA WET hypothesis test
Flowchart: EPA-821-R-02-013 Figure 2

   1  Is the response a proportion requiring transformation?
      no transformation needed
      (EPA-821-R-02-013 Appendix B, section 4.2)
   2  Are the pooled within-group residuals normally distributed?
      0.9507, p = 0.378 -> residuals consistent with normality
      (EPA-821-R-02-013 Appendix B, section 2.1)
   3  Are the variances homogeneous across concentrations?
      7.856, p = 0.097 -> variances not significantly different
      (EPA-821-R-02-013 Appendix B, section 3)
   4  Is replication equal across all concentrations?
      4, 4, 4, 4, 4 replicates, balanced
      (EPA-821-R-02-013 Figure 2)
   5  Which test was run?
      Dunnett's procedure
      (EPA-821-R-02-013 Appendix C)

  NOEC 128    LOEC 256
  MSD  0.1618    PMSD 23.9 per cent (EPA bounds 12 to 30: within)
```

Both of the manual's multi-concentration worked examples run end to end and the
flowchart independently selects the same test the manual selects: Dunnett's
procedure for the Appendix C growth data, and Steel's Many-One Rank Test for
the Appendix E reproduction data.

The point-estimation branch is complete too. `lc50()` walks the Figure 6 chart
and selects between the graphical, Spearman-Karber, trimmed Spearman-Karber and
probit methods. Each of the four columns of the acute manual's Table 20 was
constructed to exercise one of those methods, and each routes to the right one
and reproduces the published estimate:

| Column | Method selected | Estimate | Manual |
|---|---|---|---|
| `graphical` | Graphical | 35.36 | 35 (read off a plot) |
| `spearman_karber` | Spearman-Karber | 45.3 (38.9, 52.8) | 45.3 (38.9, 52.8) |
| `trimmed` | Trimmed Spearman-Karber | 77.11 (69.6, 85.4) | 77.11 (69.74, 85.26) |
| `probit` | Probit | 22.872 (18.787, 27.846) | 22.872 (18.787, 27.846) |

The trimmed Spearman-Karber interval is the one quantity in the package that
does not reproduce a published EPA value. The manual delegates it to a program
whose formula it does not state, and Hamilton's corrected variance expression is
not publicly retrievable, so the delta method is used instead. The trim itself
and the point estimate both reproduce exactly. See the vignette.

`icp()` completes the sublethal side, estimating an inhibition concentration by
linear interpolation with bootstrap limits. It reproduces the Appendix M
example exactly: smoothed means of 28.75 across the control and the three
lowest concentrations, `IC25 = 8.5715` and `IC50 = 10.893`.

`williams()` implements the one ToxCalc feature that is not an EPA method.
Neither manual gives a worked example or a table of critical values for it, so
the values are simulated and the implementation is validated against
mathematical identities instead. It is labelled as an extension throughout, and
neither flowchart will ever select it.

`R CMD check` is clean at every phase.

Still to come: the marine manual's endpoints, designs with more than one
control, and the remaining alternative tests.

The vignette `vignettes/recreating-toxcalc.qmd` records every point at which the
source material was ambiguous, internally inconsistent, or at odds with modern
practice, and what was decided in each case, with the evidence. Seven of the
manual's printed values have so far been found not to follow from its own data;
each is reconciled by a test rather than quietly worked around.

See `notes/TOXCALC-human.md` for the plan and `notes/TOXCALC-claude.md` for the
specification.

## Installation

```r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("beckyfisher/toxcalc")
```

## Fidelity

The target is faithfulness to the EPA method manuals, validated against the
worked examples printed in them. It is **not** a bit-for-bit clone of ToxCalc;
no ToxCalc output was available for comparison. Where a manual is ambiguous, the
choice made is recorded in `notes/TOXCALC-claude.md` rather than left implicit.

## References

- US EPA (2002) *Methods for Measuring the Acute Toxicity of Effluents and
  Receiving Waters to Freshwater and Marine Organisms*, 5th ed.
  EPA-821-R-02-012.
- US EPA (2002) *Short-term Methods for Estimating the Chronic Toxicity of
  Effluents and Receiving Waters to Freshwater Organisms*, 4th ed.
  EPA-821-R-02-013.
- US EPA (2002) *Short-term Methods for Estimating the Chronic Toxicity of
  Effluents and Receiving Waters to Marine and Estuarine Organisms*, 3rd ed.
  EPA-821-R-02-014.
- Hamilton MA, Russo RC, Thurston RV (1977) Trimmed Spearman-Karber method for
  estimating median lethal concentrations in toxicity bioassays.
  *Environmental Science & Technology* 11:714-719.
- Norberg-King TJ (1993) *A linear interpolation method for sublethal toxicity:
  the inhibition concentration (ICp) approach*, version 2.0. EPA/600/M-91/037.
- Williams DA (1971) A test for differences between treatment means when several
  dose levels are compared with a zero dose control. *Biometrics* 27:103-117.
