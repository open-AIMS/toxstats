#' Fathead minnow larval growth, EPA Table B.1
#'
#' Growth data from a Fathead Minnow Larval Survival and Growth Test, used for
#' the Shapiro-Wilk and Bartlett worked examples in Appendix B of the EPA
#' chronic freshwater manual.
#'
#' @details
#' The manual reports the following results for these data, which the package
#' test suite checks against:
#'
#' * concentration means 0.714, 0.674, 0.677, 0.624 and 0.580 mg;
#' * Shapiro-Wilk `W = 0.959`, against a critical value of 0.868 at
#'   `alpha = 0.01` with 20 observations, so the data are normally distributed;
#' * Bartlett `B = 7.691`, against a critical value of 13.277, so the variances
#'   are not different.
#'
#' The printed Bartlett statistic is not reproducible from these data, which
#' give 6.836. See [epa_variance()] for why.
#'
#' Note that Appendix C states its data are "the same data used in Appendices B
#' and D". They are not; see [fathead_c1].
#'
#' @format A data frame with 20 rows and 3 columns:
#' \describe{
#'   \item{conc}{sodium pentachlorophenate concentration, micrograms per litre}
#'   \item{replicate}{replicate label, A to D}
#'   \item{weight}{larval dry weight, milligrams}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-013, Table B.1.
"fathead_b1"

#' Fathead minnow larval growth, EPA Table C.1
#'
#' Growth data from a Fathead Minnow Larval Survival and Growth Test, used for
#' the Dunnett's procedure worked example in Appendix C of the EPA chronic
#' freshwater manual.
#'
#' @details
#' Despite the manual stating that Appendix C uses "the same data used in
#' Appendices B and D", Table C.1 differs from Table B.1. Both are shipped
#' separately for that reason.
#'
#' The manual reports the following results for these data:
#'
#' * concentration means 0.677, 0.575, 0.660, 0.565 and 0.454 mg;
#' * within mean square 0.0094, so `Sw = 0.097` on 15 degrees of freedom;
#' * `t` values 1.487, 0.248, 1.633 and 3.251;
#' * a one-sided Dunnett critical value of 2.36;
#' * `MSD = 0.162`, about 24 per cent of the control mean;
#' * NOEC 128 and LOEC 256 micrograms per litre.
#'
#' Section 1.11 of Appendix C describes the minimum significant difference as
#' 0.087 mg. That is a typographical error: 0.162 is the value the same section
#' computes, and 24 per cent of the 0.677 control mean is 0.162.
#'
#' @format A data frame with 20 rows and 3 columns:
#' \describe{
#'   \item{conc}{sodium pentachlorophenate concentration, micrograms per litre}
#'   \item{replicate}{replicate label, A to D}
#'   \item{weight}{larval dry weight, milligrams}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-013, Table C.1.
"fathead_c1"

#' Ceriodaphnia dubia reproduction, EPA Table B.7
#'
#' Reproduction data from a *Ceriodaphnia dubia* Survival and Reproduction
#' Test, used for the Kolmogorov "D" worked example in Appendix B of the EPA
#' chronic freshwater manual. With 60 observations the manual applies the
#' Kolmogorov statistic rather than Shapiro-Wilk.
#'
#' @details
#' The manual reports `D+ = 0.0525`, `D- = 0.0597`, `D = 0.0597` and a modified
#' statistic `D* = 0.4684`, against a critical value of 1.035 at
#' `alpha = 0.01`, so the data are normally distributed.
#'
#' @format A data frame with 60 rows and 3 columns:
#' \describe{
#'   \item{conc}{effluent concentration, per cent}
#'   \item{replicate}{replicate label, 1 to 10}
#'   \item{young}{number of young produced per female}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-013, Table B.7.
"ceriodaphnia_b7"

#' EPA variability criteria for sublethal endpoints
#'
#' Lower and upper bounds on the percent minimum significant difference for the
#' three sublethal endpoints EPA publishes criteria for. A test whose percent
#' minimum significant difference falls outside these bounds was either
#' unusually insensitive or unusually sensitive.
#'
#' @details
#' EPA derived the bounds from the 10th and 90th percentiles of the percent
#' minimum significant difference recorded in its own WET Interlaboratory
#' Variability Study, so this is an EPA-generated table rather than a
#' reproduction of a third-party one.
#'
#' Section 10.2.8.2.5 of the manual attaches a decision rule to the lower
#' bound: a concentration is not to be declared toxic if its relative
#' difference from the control is less than the lower bound, whatever the
#' hypothesis test concluded.
#'
#' @format A data frame with 3 rows and 6 columns:
#' \describe{
#'   \item{id}{short identifier, accepted by the `bounds` argument of [pmsd()]}
#'   \item{method}{EPA method number}
#'   \item{test}{full test method name}
#'   \item{endpoint}{the sublethal endpoint the bounds apply to}
#'   \item{lower}{lower percent minimum significant difference bound}
#'   \item{upper}{upper percent minimum significant difference bound}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-013, Table 6.
"epa_pmsd_bounds"

#' Ceriodaphnia dubia reproduction, EPA Table E.1
#'
#' Reproduction data from a *Ceriodaphnia dubia* seven-day chronic test, used
#' for the Steel's Many-One Rank Test worked example in Appendix E of the EPA
#' chronic freshwater manual, and, with two values removed, for the Wilcoxon
#' Rank Sum worked example in Appendix F.
#'
#' @details
#' The 50 per cent concentration is included here but is **excluded from the
#' reproduction analysis** in the manual, because survival at that
#' concentration was significantly reduced (section 9.5.2; see
#' [ceriodaphnia_g2] and [fisher_exact()]). Subset to `conc < 50` to reproduce
#' the worked example.
#'
#' Appendix F uses the same data with two males presumed to have occurred, one
#' in the control and one at 12 per cent, giving unequal replication. Remove
#' the control replicate 1 and the 12 per cent replicate 5 to reproduce it.
#'
#' The manual reports the following for the Appendix E example:
#'
#' * rank sums of 84, 64, 76 and 55 against a critical rank sum of 76;
#' * NOEC 3 per cent and LOEC 6 per cent.
#'
#' The printed rank sum of 64 at 6 per cent is a slip. The ranks the manual
#' itself lists in Table E.3 for that concentration sum to 63.5. The conclusion
#' is unaffected, both being at or below the critical value.
#'
#' For the Appendix F example the manual reports rank sums of 79, 57, 58 and
#' 55, critical values of 72 for the ten-replicate concentrations and 60 for
#' the nine-replicate one, and the same NOEC and LOEC.
#'
#' @format A data frame with 60 rows and 3 columns:
#' \describe{
#'   \item{conc}{effluent concentration, per cent}
#'   \item{replicate}{replicate label, 1 to 10}
#'   \item{young}{number of young produced}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-013, Table E.1.
"ceriodaphnia_e1"

#' Ceriodaphnia dubia mortality, EPA Table G.2
#'
#' Survival data from a *Ceriodaphnia dubia* survival and reproduction test,
#' used for the Fisher's Exact Test worked example in Appendix G of the EPA
#' chronic freshwater manual. One row per concentration, already pooled across
#' replicates, as the manual presents it.
#'
#' @details
#' `exposed` is the number of live adults at the beginning of the test, which
#' is nine in the control rather than ten.
#'
#' The manual concludes that only the 25 per cent concentration differs
#' significantly from the control, giving a NOEC of 12 per cent and a LOEC of
#' 25 per cent for survival. That result is what excludes the 50 per cent
#' concentration from the reproduction analysis of [ceriodaphnia_e1].
#'
#' Appendix G applies a pairwise error rate of 0.05 rather than an
#' experiment-wise one, because Fisher's Exact Test is itself conservative.
#'
#' @format A data frame with 6 rows and 3 columns:
#' \describe{
#'   \item{conc}{effluent concentration, per cent}
#'   \item{dead}{number of adults dead at the end of the test}
#'   \item{exposed}{number of live adults at the start of the test}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-013, Table G.2.
"ceriodaphnia_g2"

#' LC50 worked-example mortality data, EPA acute Table 20
#'
#' Mortality counts used for the four median lethal concentration worked
#' examples in section 11.2 of the EPA acute manual. Twenty organisms were
#' exposed in the control and at every concentration.
#'
#' @details
#' Each method column is a different set of counts, chosen so that the Figure 6
#' flowchart routes it to that method. Together they exercise all four
#' terminals of the chart, which is how [lc50()] is validated.
#'
#' The published results are:
#'
#' * `graphical`: LC50 read off the plot as 35 per cent; interpolating on the
#'   logarithmic scale gives 35.36. No confidence interval.
#' * `spearman_karber`: `m = 1.656527`, `V(m) = 0.0010977`, LC50 45.3 per cent
#'   with limits 38.9 and 52.8.
#' * `trimmed`: an automatic trim of 20.51 per cent and an LC50 of 77.11 with
#'   limits 69.74 and 85.26. The limits here differ slightly; see
#'   [trimmed_spearman_karber()].
#' * `probit`: a Pearson chi-square for heterogeneity of 3.076 on 3 degrees of
#'   freedom against a tabular 7.815, LC50 22.872 with limits 18.787 and
#'   27.846, and LC1 7.924 with limits 4.147 and 10.959.
#'
#' The probit column has no control mortality, which is why the manual's
#' printed output shows its Abbott-adjusted proportions equal to the observed
#' ones.
#'
#' @format A data frame with 6 rows and 6 columns:
#' \describe{
#'   \item{conc}{effluent concentration, per cent}
#'   \item{exposed}{number of organisms exposed}
#'   \item{graphical}{number dead, for the graphical method example}
#'   \item{spearman_karber}{number dead, for the Spearman-Karber example}
#'   \item{trimmed}{number dead, for the trimmed Spearman-Karber example}
#'   \item{probit}{number dead, for the probit example}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-012, Table 20.
"acute_table20"

#' Ceriodaphnia dubia reproduction, EPA Table M.1
#'
#' Reproduction data used for the linear interpolation worked example in
#' Appendix M of the EPA chronic freshwater manual.
#'
#' @details
#' This shares four columns with [ceriodaphnia_b7] but is not the same dataset.
#' Table M.1 has no counterpart to Table B.7's 12.5 per cent column, its own
#' 12.5 per cent column holds what B.7 records at 25 per cent, and it adds a
#' concentration at which reproduction stopped entirely.
#'
#' The manual reports, for these data:
#'
#' * concentration means 22.4, 26.3, 34.6, 31.7, 9.4 and 0;
#' * smoothed means of 28.75 across the control and the three lowest
#'   concentrations, then 9.4 and 0;
#' * `IC25 = 8.5715` and `IC50 = 10.89`.
#'
#' The ICPIN program output also prints standard deviations of 6.931, 8.001,
#' 4.835, 2.946, 3.893 and 0, a bootstrap mean of 8.5891 with a standard
#' deviation of 0.1831 from 80 resamples, and confidence limits of 8.3112 and
#' 9.0418. Those limits were drawn with the seed -641671986 by a Turbo Pascal
#' generator and cannot be reproduced by any other program; see [icp()].
#'
#' @format A data frame with 60 rows and 3 columns:
#' \describe{
#'   \item{conc}{effluent concentration, per cent}
#'   \item{replicate}{replicate label, 1 to 10}
#'   \item{young}{number of young produced}
#' }
#'
#' @source US EPA (2002) EPA-821-R-02-013, Table M.1.
"ceriodaphnia_m1"
