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

**Experimental, under active development.** Phases 0 to 2 are complete: the data
layer, the EPA transformations and monotone smoothing, the normality and
homogeneity-of-variance tests, and the minimum significant difference. The
package passes `R CMD check` cleanly at every phase. The individual hypothesis
tests, the flowchart engine and the point-estimation methods are still to come.

The vignette `vignettes/recreating-toxcalc.qmd` records every point at which the
source material was ambiguous, internally inconsistent, or at odds with modern
practice, and what was decided in each case, with the evidence. Four of the
manuals' printed statistics have so far been found not to follow from their own
data; each is reconciled by a test rather than quietly worked around.

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
