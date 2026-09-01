# Package index

## Start here

The drivers. Each walks an EPA flowchart, selects the method the manual
directs, and records why.

- [`tox_test()`](https://open-aims.github.io/toxstats/reference/tox_test.md)
  [`summary(`*`<tox_test>`*`)`](https://open-aims.github.io/toxstats/reference/tox_test.md)
  [`as.data.frame(`*`<tox_test>`*`)`](https://open-aims.github.io/toxstats/reference/tox_test.md)
  : Run the EPA hypothesis-testing analysis
- [`lc50()`](https://open-aims.github.io/toxstats/reference/lc50.md)
  [`summary(`*`<tox_lc50>`*`)`](https://open-aims.github.io/toxstats/reference/lc50.md)
  [`as.data.frame(`*`<tox_lc50>`*`)`](https://open-aims.github.io/toxstats/reference/lc50.md)
  : Estimate a median lethal concentration
- [`decisions()`](https://open-aims.github.io/toxstats/reference/decisions.md)
  : Extract the decision trail

## Preparing data

The input contract and the transformations every method begins with.

- [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  [`print(`*`<tox_data>`*`)`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  [`summary(`*`<tox_data>`*`)`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  [`as.data.frame(`*`<tox_data>`*`)`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  : Prepare toxicity test data for analysis
- [`arcsine_sqrt()`](https://open-aims.github.io/toxstats/reference/arcsine_sqrt.md)
  : Arc sine square root transformation
- [`inv_arcsine_sqrt()`](https://open-aims.github.io/toxstats/reference/inv_arcsine_sqrt.md)
  : Back-transform an arc sine square root value
- [`abbott()`](https://open-aims.github.io/toxstats/reference/abbott.md)
  : Abbott's correction for control response
- [`smooth_monotone()`](https://open-aims.github.io/toxstats/reference/smooth_monotone.md)
  : Smooth a sequence to be monotone
- [`is_monotone()`](https://open-aims.github.io/toxstats/reference/is_monotone.md)
  : Is a sequence monotone?

## Assumptions and sensitivity

The tests that decide which branch of the hypothesis-testing flowchart
is taken, and the measures of how sensitive the test was.

- [`epa_normality()`](https://open-aims.github.io/toxstats/reference/epa_normality.md)
  : Test the normality assumption
- [`epa_variance()`](https://open-aims.github.io/toxstats/reference/epa_variance.md)
  : Test the homogeneity of variance assumption
- [`msd()`](https://open-aims.github.io/toxstats/reference/msd.md) :
  Minimum significant difference
- [`pmsd()`](https://open-aims.github.io/toxstats/reference/pmsd.md) :
  Percent minimum significant difference

## Hypothesis tests

The tests the flowchart can select.
[`williams()`](https://open-aims.github.io/toxstats/reference/williams.md)
is an extension rather than an EPA method and is never selected
automatically.

- [`dunnett()`](https://open-aims.github.io/toxstats/reference/dunnett.md)
  : Dunnett's procedure
- [`bonferroni_t()`](https://open-aims.github.io/toxstats/reference/bonferroni_t.md)
  : The t test with Bonferroni's adjustment
- [`dunn_sidak_t()`](https://open-aims.github.io/toxstats/reference/dunn_sidak_t.md)
  : The t test with Dunn-Sidak's adjustment
- [`welch_t()`](https://open-aims.github.io/toxstats/reference/welch_t.md)
  : Welch's t test against the control
- [`steel()`](https://open-aims.github.io/toxstats/reference/steel.md) :
  Steel's Many-One Rank Test
- [`wilcoxon_rank_sum()`](https://open-aims.github.io/toxstats/reference/wilcoxon_rank_sum.md)
  : Wilcoxon Rank Sum Test with Bonferroni's adjustment
- [`fisher_exact()`](https://open-aims.github.io/toxstats/reference/fisher_exact.md)
  : Fisher's Exact Test
- [`williams()`](https://open-aims.github.io/toxstats/reference/williams.md)
  : Williams' test

## Point estimation

Median lethal concentrations from the Figure 6 flowchart, and inhibition
concentrations by linear interpolation.

- [`graphical_lc50()`](https://open-aims.github.io/toxstats/reference/graphical_lc50.md)
  : Estimate the LC50 by the graphical method
- [`spearman_karber()`](https://open-aims.github.io/toxstats/reference/spearman_karber.md)
  : Estimate the LC50 by the Spearman-Karber method
- [`trimmed_spearman_karber()`](https://open-aims.github.io/toxstats/reference/trimmed_spearman_karber.md)
  : Estimate the LC50 by the trimmed Spearman-Karber method
- [`probit_lc()`](https://open-aims.github.io/toxstats/reference/probit_lc.md)
  : Estimate lethal concentrations by the probit method
- [`icp()`](https://open-aims.github.io/toxstats/reference/icp.md) :
  Estimate an inhibition concentration by linear interpolation

## EPA worked-example data

The datasets the manuals use in their worked examples, each documented
with the results the manual publishes for it.

- [`fathead_b1`](https://open-aims.github.io/toxstats/reference/fathead_b1.md)
  : Fathead minnow larval growth, EPA Table B.1
- [`fathead_c1`](https://open-aims.github.io/toxstats/reference/fathead_c1.md)
  : Fathead minnow larval growth, EPA Table C.1
- [`ceriodaphnia_b7`](https://open-aims.github.io/toxstats/reference/ceriodaphnia_b7.md)
  : Ceriodaphnia dubia reproduction, EPA Table B.7
- [`ceriodaphnia_e1`](https://open-aims.github.io/toxstats/reference/ceriodaphnia_e1.md)
  : Ceriodaphnia dubia reproduction, EPA Table E.1
- [`ceriodaphnia_g2`](https://open-aims.github.io/toxstats/reference/ceriodaphnia_g2.md)
  : Ceriodaphnia dubia mortality, EPA Table G.2
- [`ceriodaphnia_m1`](https://open-aims.github.io/toxstats/reference/ceriodaphnia_m1.md)
  : Ceriodaphnia dubia reproduction, EPA Table M.1
- [`acute_table20`](https://open-aims.github.io/toxstats/reference/acute_table20.md)
  : LC50 worked-example mortality data, EPA acute Table 20
- [`epa_pmsd_bounds`](https://open-aims.github.io/toxstats/reference/epa_pmsd_bounds.md)
  : EPA variability criteria for sublethal endpoints
