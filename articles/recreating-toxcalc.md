# Recreating ToxCalc: methods and decisions

## Purpose

This vignette records how `toxcalc` was built and, more importantly,
**every point at which the source material was ambiguous, internally
inconsistent, or at odds with modern practice, and what was decided in
each case**.

The package recreates ToxCalc v5.0 (Tidepool Scientific), commercial
software implementing the US EPA Whole Effluent Toxicity methods, which
has been withdrawn. No ToxCalc output was available for comparison, so
the specification is taken from the EPA method manuals themselves and
validated against the worked examples printed in them.

Every decision below is stated with the evidence that settled it. Where
the package departs from a figure printed in a manual, the departure is
demonstrated, not asserted: each is reconciled by a test in the package
test suite, so the discrepancy cannot later be mistaken for an
implementation error.

## Sources

- EPA-821-R-02-012, *Methods for Measuring the Acute Toxicity of
  Effluents and Receiving Waters to Freshwater and Marine Organisms*,
  5th edition, October
  2002. 
- EPA-821-R-02-013, *Short-term Methods for Estimating the Chronic
  Toxicity of Effluents and Receiving Waters to Freshwater Organisms*,
  4th edition, October
  2002. Section references below are to this manual unless stated
        otherwise.

Both are US Government works and are publicly available.

## The decisions

### 1. Statistical tables are computed, not transcribed

**Decision.** No table of critical values reproduced in the EPA manuals
from a third-party source is copied into this package. Critical values
are computed.

**Reasoning.** The manuals are US Government works, but several tables
*within* them are reproductions credited to other authors: Dunnett’s
Table C.5 to Miller (1981), the Shapiro-Wilk coefficients of Table B.4
to Conover (1980), and Williams’ critical values to *Biometrics*.
Permission to reproduce a table in a government manual does not transfer
to a third party redistributing it.

Two further considerations point the same way. Computation covers
designs the tables never did, in particular unbalanced replication. And
the reproduction in the source document is not always reliable: the
Dunnett table as printed contains values corrupted in transcription,
including 1.86 at 120 degrees of freedom with one treatment, where the
sequence requires a value near 1.66.

Only tables EPA generated itself are shipped: the Kolmogorov critical
values of Table B.11, and the variability criteria of Table 6, which EPA
derived from its own interlaboratory study.

### 2. Dunnett’s critical value is integrated, not simulated

**Decision.** The critical value is obtained by numerical integration,
using [`stats::integrate`](https://rdrr.io/r/stats/integrate.html) and
[`stats::uniroot`](https://rdrr.io/r/stats/uniroot.html). The obvious
alternative,
[`mvtnorm::qmvt()`](https://rdrr.io/pkg/mvtnorm/man/qmvt.html), is not
used, and `mvtnorm` is not a dependency.

**Reasoning.** `qmvt()` uses randomised quasi-Monte Carlo integration
and returns a different answer on every call. For the Appendix C design
it returned values between 2.3552 and 2.3572 across repeated calls. An
analysis submitted to a regulator must give the same answer every time
it is run.

Because all the comparisons share the same control, the correlation
matrix has a one-factor form: the correlation between comparisons `i`
and `j` is `lambda_i * lambda_j`, where
`lambda_i = sqrt(n_i / (n_control + n_i))`. Conditioning on that single
common factor, and on the pooled standard deviation, reduces the
multivariate t probability to a two-dimensional integral, which ordinary
quadrature evaluates to machine precision.

The result agrees with the manual’s table and is stable:

``` r

fit <- msd(fathead_c1, response = "weight")
fit$critical           # EPA Table C.5 prints 2.36
#> [1] 2.35614
replicate(3, msd(fathead_c1, response = "weight")$critical)
#> [1] 2.35614 2.35614 2.35614
```

### 3. Shapiro-Wilk follows Royston, not the manual’s tables

**Decision.** Normality is tested with
[`stats::shapiro.test()`](https://rdrr.io/r/stats/shapiro.test.html),
which uses Royston’s AS R94 algorithm, and the p-value is compared with
the assumption significance level. The manual’s procedure — looking up
Conover’s coefficients in Table B.4 and comparing the statistic with a
quantile in Table B.6 — is not implemented.

**Reasoning.** The two agree. On the manual’s own Appendix B worked
example Royston returns 0.9601 against a printed 0.959, and the
difference is entirely accounted for by the manual centring its
observations with concentration means rounded to three decimal places;
applying the same rounding gives 0.9594.

``` r

epa_normality(fathead_b1, response = "weight")
#> Shapiro-Wilk test for normality (Royston AS R94)
#>   EPA-821-R-02-013 Appendix B, section 2.1
#> 
#>   W = 0.96007, p = 0.5452
#>   alpha = 0.01
#>   The residuals are consistent with a normal distribution.
```

Royston’s algorithm is preferred because it is not limited to the
tabulated sample sizes, it returns a p-value rather than a decision at
one fixed cutoff, and it avoids reproducing a third-party table. The
manual’s tabulated route is a hand-calculation aid from before this was
available in software.

**Residual risk.** In a borderline case the two could send an analysis
down different branches of the flowchart. The package reports both the
statistic and the p-value so that such a case is visible.

### 4. Six printed values do not follow from the manual’s own data

The manuals were prepared for hand calculation, and several worked
examples carry intermediate values forward after rounding them.
Recomputing at full precision therefore does not always return the
printed figure. In every case found so far the conclusion is unchanged,
but the numbers differ, and anyone checking this package against the
manual needs to know which is which.

#### Bartlett’s test: 7.691 against 6.836

The Appendix B example prints `B = 7.691`. From the raw data the
statistic is 6.836. The manual computed it from its own printed table of
variances rounded to four decimal places, and then divided by a
correction factor also rounded, to 1.133.

``` r

epa_variance(fathead_b1, response = "weight")$statistic
#> [1] 6.836106

# reproducing the printed value from the manual's rounded intermediates
printed <- c(0.0018, 0.0020, 0.0001, 0.0059, 0.0037)
nu <- rep(3, 5)
(sum(nu) * log(sum(nu * printed) / sum(nu)) - sum(nu * log(printed))) / 1.133
#> [1] 7.690975
```

One of those rounded variances is itself wrong: the variance at 32 µg/L
is 0.002055, which the manual prints as 0.0020 rather than 0.0021.
Bartlett’s statistic depends on the logarithm of each variance, so a
small relative error in the smallest variance moves it appreciably.

**Decision.** Compute from the unrounded data. Both values are far below
the 13.277 critical value.

#### Kolmogorov D: 0.4684 against 0.4572

The Appendix B example prints `D* = 0.4684`. At full precision it is
0.4572. The manual standardises each residual and then looks the
probability up in a printed normal table, which requires rounding each
`z` to two decimal places first. Doing the same reproduces the printed
value:

``` r

epa_normality(ceriodaphnia_b7, response = "young")$statistic
#> [1] 0.4572353

residuals <- ceriodaphnia_b7$young -
  ave(ceriodaphnia_b7$young, ceriodaphnia_b7$conc,
      FUN = function(z) round(mean(z), 1))
n <- length(residuals)
p <- pnorm(round(sort(residuals) / sd(residuals), 2))
d <- max(max(seq_len(n) / n - p), max(p - (seq_len(n) - 1) / n))
d * (sqrt(n) - 0.01 + 0.85 / sqrt(n))
#> [1] 0.4682986
```

**Decision.** Compute at full precision. Both values are far below the
1.035 critical value.

#### Appendix C section 1.11: 0.087 against 0.162

The same section that computes `MSD = 0.162` then describes the minimum
detectable difference as 0.087 mg, “a decrease in growth of 24 % from
the control”. The control mean is 0.677, and 24 % of it is 0.162, so
0.162 is correct and 0.087 is a typographical error.

``` r

pmsd(fathead_c1, response = "weight", bounds = "fathead_growth")
#> Percent minimum significant difference
#>   EPA-821-R-02-013 section 10.2.8
#> 
#>   control mean = 0.6773
#> 
#>   32   64  128  256 
#> 23.9 23.9 23.9 23.9 
#> 
#>   EPA bounds: 12 to 30 per cent
#>       32       64      128      256 
#> "within" "within" "within" "within"
```

#### Appendix D Table D.3: a within mean square of 0.0029 against 0.0079

Table D.3 gives the within sum of squares as 0.111 on 14 degrees of
freedom and then the mean square as 0.0029. The quotient is 0.0079. The
manual’s own printed `t` values were computed from 0.0079, so the mean
square as printed is the error and everything downstream of it is right.

``` r

lost <- fathead_c1[-19, ]  # Appendix D loses one 256 ug/L replicate
fit <- bonferroni_t(lost, response = "weight")
c(sum_of_squares = fit$sw^2 * fit$df, mean_square = fit$sw^2)
#> sum_of_squares    mean_square 
#>    0.110724917    0.007908923
```

#### Appendix D Table D.4: a t value of 0.220 against 0.270

For the 64 µg/L comparison the manual prints 0.220. Its own
concentration means of 0.677 and 0.660 and its own pooled standard
deviation give 0.270, and the other three comparisons agree with the
manual to within rounding. Both values are far below the 2.510 critical
value.

``` r

fit$comparisons[, c("conc", "statistic")]
#>   conc statistic
#> 1   32 1.6220211
#> 2   64 0.2703369
#> 3  128 1.7850183
#> 4  256 4.0278434
```

#### Appendix E Table E.4: a rank sum of 64 against 63.5

The ranks the manual itself lists in Table E.3 for the 6 per cent
concentration are 3, 7.5, 1.5, 1.5, 11.5, 7.5, 4.5, 7.5, 11.5 and 7.5.
Those sum to 63.5, not the 64 printed in Table E.4. Both are at or below
the critical rank sum of 76, so the concentration is significant either
way.

``` r

repro <- ceriodaphnia_e1[ceriodaphnia_e1$conc < 50, ]
steel(repro, response = "young")$comparisons[, c("conc", "rank_sum")]
#>   conc rank_sum
#> 1    3     84.0
#> 2    6     63.5
#> 3   12     76.0
#> 4   25     55.0
```

### 5. Appendix C does not use the data it says it does

Appendix C states that its data are “the same data used in Appendices B
and D”. Table C.1 differs from Table B.1 at fifteen of twenty values.
The two are shipped as separate datasets, `fathead_b1` and `fathead_c1`,
and each is checked against the summary statistics printed beside it.

``` r

rbind(
  `Table B.1` = tapply(fathead_b1$weight, fathead_b1$conc, mean),
  `Table C.1` = tapply(fathead_c1$weight, fathead_c1$conc, mean)
)
#>                 0      32      64    128     256
#> Table B.1 0.71450 0.67375 0.67700 0.6235 0.58050
#> Table C.1 0.67725 0.57525 0.66025 0.5650 0.45425
```

Using the wrong one is not a small matter: Table B.1 gives a control
mean of 0.714 and a within mean square root of 0.052, against 0.677 and
0.097 for Table C.1, and the two lead to different no-observed-effect
concentrations.

### 6. Smoothing is full pool-adjacent-violators, not pairwise averaging

**Decision.**
[`smooth_monotone()`](https://open-aims.github.io/toxstats/reference/smooth_monotone.md)
implements the pool-adjacent-violators algorithm, with block sizes
tracked so that a run of several values pools to their common mean in
one step.

**Reasoning.** The manuals describe smoothing as averaging a violating
*pair*. Taken literally that is wrong, because their own worked examples
pool longer runs: the Appendix K mortality proportions collapse from
five values to a single 0.02, which repeated pairwise averaging does not
produce. The manuals acknowledge this obliquely, noting that “unusual
patterns in the deviations from monotonicity may require an additional
step of smoothing”.

``` r

smooth_monotone(c(0.05, 0.00, 0.05, 0.00, 0.00), direction = "increasing")
#> [1] 0.02 0.02 0.02 0.02 0.02
```

Weights are exposed because the two places the manuals use smoothing
need different ones. The linear interpolation method smooths
concentration means with equal weight, whereas the isotonic means of
Williams’ test are weighted by replication.

### 7. The arc sine transformation includes its endpoint adjustments

**Decision.**
[`arcsine_sqrt()`](https://open-aims.github.io/toxstats/reference/arcsine_sqrt.md)
substitutes `1 / (4 * n)` for an observed proportion of 0 and
`1 - 1 / (4 * n)` for a proportion of 1, where `n` is the number of
organisms in the replicate, and requires `n` whenever either occurs.

**Reasoning.** These adjustments are part of the method, not a
refinement. A replicate at 0 or 1 would otherwise contribute no variance
and the tests downstream would understate the true variability. They
reproduce the manual’s printed values exactly:

``` r

c(
  `RP = 0.60` = arcsine_sqrt(0.60),   # manual prints 0.8861
  `RP = 0`    = arcsine_sqrt(0, n = 20),  # manual prints 0.1120
  `RP = 1`    = arcsine_sqrt(1, n = 20)   # manual prints 1.4588
)
#> RP = 0.60    RP = 0    RP = 1 
#> 0.8860771 0.1120376 1.4587587
```

### 8. Williams’ test is an extension, not an EPA method

**Decision.** Williams’ test will be implemented, because ToxCalc
offered it, but it is documented as an extension and is not on the
flowchart the package walks by default.

**Reasoning.** The manuals do not include Williams’ test in the
hypothesis-testing flowchart of Figure 2. Section 9.4.1.2 mentions it
only as an alternative requiring additional assumptions, principally
that the concentration-response relationship is monotone. Presenting it
as an EPA method would misrepresent the manuals.

### 9. Two significance levels, not one

**Decision.** The package carries `alpha_assumption`, defaulting to
0.01, and `alpha`, defaulting to 0.05.

**Reasoning.** Section 9.4.6.1 sets the assumption tests at 0.01 and the
hypothesis tests at 0.05, one-sided. A single significance level
argument would be wrong for one or the other.

### 10. Dependencies are kept to what is genuinely needed

**Decision.** At runtime the package depends on `chk` and base `stats`
only. `MASS`, `multcomp`, `mvtnorm`, `nortest` and `boot` appear in
`Suggests`, where they are used in the test suite as independent checks
on methods implemented here from first principles.

**Reasoning.** Each was considered and rejected as a runtime dependency.
`mvtnorm` was displaced by the integration in decision 2. `multcomp`
would have been used for one quantile while pulling in four further
packages. `nortest` provides statistics the manuals do not specify.
[`MASS::dose.p`](https://rdrr.io/pkg/MASS/man/dose.p.html) supplies a
median lethal concentration but none of the Abbott adjustment,
heterogeneity test or Fieller limits the manuals also require. `boot`
cannot express the EPA bootstrap, which resamples within concentration,
re-smooths inside every iteration, and takes an order-statistic
interval.

`PMCMRplus`, the usual source for Williams’ and Steel’s tests, is not
used: it requires the external MPFR and GMP libraries, which makes it
impractical to install on Windows.

### 11. The flowchart is data, not control flow

**Decision.** Each EPA decision chart is stored as a data frame of
branch points, and a small engine walks it. The alternative, nesting
`if` statements inside the driver, was rejected.

**Reasoning.** Three things follow from holding the chart as data. Each
branch point can be tested on its own, without running an analysis end
to end. The traversal can be recorded without the recording being
tangled into the logic, which is what makes the audit trail trustworthy
rather than a narration written alongside the code. And a second chart —
the acute manual has three more — is added without touching the engine.

The chart carries the manual section for every branch point, so the
citation in the printed output is not a comment that can drift away from
the code:

``` r

chart <- toxstats:::flowchart_hypothesis_multi()
chart[, c("node", "yes", "no", "reference")]
#>             node            yes             no
#> 1      transform      normality      normality
#> 2      normality       variance min_replicates
#> 3       variance  replication_p min_replicates
#> 4 min_replicates replication_np @no_valid_test
#> 5  replication_p       @dunnett  @bonferroni_t
#> 6 replication_np         @steel      @wilcoxon
#>                                  reference
#> 1 EPA-821-R-02-013 Appendix B, section 4.2
#> 2 EPA-821-R-02-013 Appendix B, section 2.1
#> 3   EPA-821-R-02-013 Appendix B, section 3
#> 4         EPA-821-R-02-013 section 9.4.5.2
#> 5                EPA-821-R-02-013 Figure 2
#> 6                EPA-821-R-02-013 Figure 2
```

Two details of the chart are worth pointing out because they are easy to
get wrong when reading Figure 2 quickly. The normality branch goes
**directly** to the replication question when it fails: the variance
test is not reached at all. And there is an explicit `@no_valid_test`
terminal, because a design whose residuals are not normal and which has
fewer than four replicates has no test the manual sanctions — section
9.4.5.2 rules out the non-parametric branch. That case raises an error
rather than falling through to something plausible.

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

### 12. Two rules can overturn a test result

Both are in the manual and both are easy to miss, because neither is on
the flowchart.

**Exclusions, section 9.5.2.** A concentration at which every replicate
showed a complete response is excluded from the no-observed-effect
concentration, though it is retained for point estimation. That one is
applied automatically for a quantal endpoint. A concentration at which
*survival* was significantly reduced is excluded from the analysis of a
**sublethal** endpoint, which needs a separate survival analysis and so
cannot be inferred; it is passed to `exclude`. The Appendix E example
does exactly this, dropping the 50 per cent concentration from the
reproduction analysis on the strength of the Appendix G survival result.

**The lower PMSD bound, section 10.2.8.2.5.** A concentration shall
**not** be declared toxic if its relative difference from the control is
smaller than the lower bound published for the test method, whatever the
hypothesis test concluded. This reverses a significant result after the
fact. It is applied only when `pmsd_bounds` is given, and every
concentration it reverses gets its own row in the decision trail, so the
reversal is visible rather than silent.

### 13. Overriding the flowchart is possible but never quiet

**Decision.** `test` forces a particular test. When it disagrees with
the flowchart a warning is raised, `overridden` is set on the result, an
extra row is added to the decision trail, and
[`print()`](https://rdrr.io/r/base/print.html) reports it above the
endpoints.

**Reasoning.** A statistician sometimes has good reason to overrule a
borderline assumption test. What must not happen is that a forced
analysis is later read as the analysis the manual prescribes. Making the
override loud costs the user nothing and removes that risk.

``` r

fit <- withCallingHandlers(
  tox_test(fathead_c1, response = "weight", test = "steel"),
  warning = function(w) invokeRestart("muffleWarning")
)
c(selected_by_flowchart = fit$selected, actually_run = fit$comparison$test)
#> selected_by_flowchart          actually_run 
#>             "dunnett"               "steel"
decisions(fit)[decisions(fit)$node == "override", c("outcome", "reference")]
#>                                             outcome           reference
#> 5 user forced steel; the flowchart selected dunnett Not an EPA decision
```

### 14. The point-estimation methods, settled by evidence

Four of the manual’s ambiguities were resolved by testing candidate
constructions against the printed output rather than by choosing a
convention.

**The probit heterogeneity statistic is Pearson’s, not the deviance.**
Both are chi-square statistics on the same fit and the manual names
neither. For the worked example Pearson gives 3.0764 against a printed
3.076, and the deviance gives 3.8588, so the two are distinguishable and
the manual is using Pearson.

**The probit limits are Fieller’s, not the delta method.** Fieller
reproduces the manual’s printed limits to every digit — 18.787 to 27.846
for the LC50 and 4.147 to 10.959 for the LC1 — while the delta method
gives 19.04 to 27.47 and 5.17 to 12.14.

``` r

probit_lc(
  acute_table20,
  response = "probit", n_exposed = "exposed", type = "quantal",
  p = c(1, 50)
)
#> Probit method
#>   EPA-821-R-02-012 section 11.2.5
#> 
#>  endpoint  p estimate   lower  upper                  ci_method
#>       LC1  1   7.9239  4.1468 10.959 Fieller on the log10 scale
#>      LC50 50  22.8720 18.7870 27.846 Fieller on the log10 scale
#> 
#>   Chi-square for heterogeneity: 3.076 on 3 df (critical 7.815)
```

**The Spearman-Karber multiplier is 2.0.** Section 11.2.3.3 step 6
writes the interval as `m ± 2.0 sqrt(V(m))`. That is the manual’s
rounding of 1.96, and keeping it is what reproduces the printed 38.9 to
52.8. The exact normal quantile is used at other confidence levels.

**Partial responses are counted after smoothing and adjustment.** Figure
6 branches on how many concentrations show a partial response, without
saying whether that means the raw proportions or the smoothed, adjusted
ones. The manual’s own Spearman-Karber example settles it: the raw data
show two partial responses, which would route to the probit branch,
while the smoothed adjusted data show one, which routes to
Spearman-Karber, and Spearman-Karber is what the manual uses.

``` r

for (column in c("graphical", "spearman_karber", "trimmed", "probit")) {
  fit <- lc50(
    acute_table20,
    response = column, n_exposed = "exposed", type = "quantal"
  )
  cat(format(column, width = 16), "->", fit$selected, "\n")
}
#> graphical        -> graphical_lc50 
#> spearman_karber  -> spearman_karber 
#> trimmed          -> trimmed_spearman_karber 
#> probit           -> probit_lc
```

### 15. The one place the package does not reproduce the manual

**The trimmed Spearman-Karber confidence interval.**

The trim itself is settled. Section 11.2.4.3 step 4 defines it as
`max(p1, 1 - pk)` on the smoothed, Abbott-adjusted proportions, and that
reproduces both published trims exactly: 20.51 per cent in the acute
manual and 20.41 per cent in Appendix K of the chronic manual. The
ordering matters — smoothing, then the Abbott adjustment, then the trim
— and computing it from the raw proportions gives a different answer.

The point estimate is settled too: 77.1105 against the 77.11 the
manual’s program prints.

The interval took longer. The manual delegates it to a program whose
source is not available and states neither the variance formula nor the
multiplier. Hamilton et al. (1977) published a variance expression that
Hamilton et al. (1978) then corrected, and no retrievable source gives
the corrected form.

Three candidates were tried first, and all three were wrong. The plain
Spearman-Karber variance scaled by `1/(1-2a)^2` gives 70.8 to 84.0, too
narrow. A parametric bootstrap gives 68.9 to 86.3, too wide. A
**numerical** delta method, perturbing each observed proportion and
recomputing the estimate, gives 69.6 to 85.4 — close, and for a while
accepted as the best available.

It was wrong for two reasons, and finding them is what settled the
section.

**The trim must be held fixed.** The automatic trim is a function of the
data, so recomputing it inside the perturbation propagates the choice of
trim into the variance. That is arguably the more honest quantity, but
it is not what the manual’s program reports. Holding the trim fixed is
what Hamilton does.

**The multiplier is 2.0, not 1.96.** Section 11.2.3.3 writes the
untrimmed interval as `m ± 2.0 sqrt(V(m))`, and the trimmed method
inherits it. The two differ by about two per cent of the half-width,
which is the size of the residual discrepancy that had been attributed
to the variance.

With the derivative taken analytically instead of numerically, the trim
held fixed and the multiplier at 2.0, the package returns 69.74 to 85.26
— the manual’s limits, to every digit printed.

``` r

trimmed_spearman_karber(
  acute_table20,
  response = "trimmed", n_exposed = "exposed", type = "quantal"
)
#> Trimmed Spearman-Karber method (automatic trim 20.51 per cent)
#>   EPA-821-R-02-012 section 11.2.4
#> 
#>  endpoint  p estimate  lower  upper                            ci_method
#>      LC50 50   77.111 69.743 85.256 Hamilton variance on the log10 scale
```

### Deriving the variance rather than transcribing it

The derivation is set out in full beside the internal `tsk_variance()`.
In outline: summation by parts rewrites the Spearman-Karber sum so that
every dependence on the proportions is explicit, after which each
partial derivative is elementary. The interior proportions enter through
their own rescaled value; the four proportions bracketing the trim and
its complement enter through the interpolated endpoints.

This is written as one general expression rather than as a separate
formula for each possible number of interior concentrations, and that
choice earned itself. Hamilton’s variance is conventionally presented as
five or six named terms selected by how many concentrations fall inside
the trim, and the boundary cases are where such presentations go wrong.

### An independent check, and what it found

The `ecotoxicology` package is a translation of the EPA’s own Ecological
Exposure Research Division programs, discontinued in 1999, and it
carries Hamilton’s variance directly. It is the only external
implementation located, and it is used in the test suite as an oracle
rather than as a source: it is GPL (\>= 3), which this package cannot
absorb at GPL (\>= 2).

Across 1,253 randomly generated designs the two agree on the point
estimate to machine precision, and on the variance to machine precision
in every configuration except one. Where exactly one concentration falls
inside the trim, they disagree by up to 60 per cent of the standard
deviation.

That case is covered in `ecotoxicology` by a term named `V6`, and the
sign of its head contribution is opposite to the sign the same
contribution carries in `V2`, the term used when more than one
concentration falls inside. Flipping it brings all 216 such designs into
machine-precision agreement. A numerical derivative of this package’s
own estimator, which needs no external reference at all, confirms the
sign used here. The disagreement does not affect the EPA worked example,
which has no concentration inside the trim.

The comparison also found two defects in
[`ecotoxicology::TSK`](https://rdrr.io/pkg/ecotoxicology/man/TSK.html)
that prevent it being called on the manual’s own data: it formats a
non-integer trim with `%d`, and where the control response is not zero
it divides by `n[-1] - n[1]`, which is zero in a balanced design. Both
are confined to its control-handling branch, so the oracle test passes
it the treatments only.

Two smaller findings came out of the same comparison. The trimmed method
here had been taking the **last** concentration at or below `1 - a` as
the upper endpoint, where the correct convention is the **first**
concentration to reach it; the two differ only when the response
plateaus exactly on the boundary, which no EPA worked example does. And
`ecotoxicology` reports a heterogeneity chi-square of 3.076 for the
probit example, confirming from an independent implementation the
Pearson statistic chosen in section 14 over the deviance.

### 16. The inhibition concentration, and the step that is easy to miss

**Decision.** The bootstrap re-smooths the concentration means
**inside** every iteration.

**Reasoning.** The manual’s step list places smoothing at step 4 and
resampling at step 6, which reads as though the smoothing happens once.
Leaving it out of the loop is not a small error. Unsmoothed resampled
means are not monotone, so the interpolation selects the wrong
bracketing pair, and the resampled estimates drift away from the point
estimate altogether. On the Appendix M data they then average about 10.4
against a point estimate of 8.57, with a standard deviation of 0.49
rather than 0.14.

That is worth recording precisely because the original plan for this
package predicted the opposite — that omitting the re-smoothing would
make the interval too narrow. Measuring it showed a bias, not a
narrowing.

``` r

icp(ceriodaphnia_m1, response = "young", p = c(25, 50), seed = 42)
#> Linear interpolation method
#>   EPA-821-R-02-013 Appendix M
#> 
#>  endpoint  p estimate   lower   upper                  ci_method bound
#>      IC25 25   8.5715  8.3009  8.8863 bootstrap order statistics  <NA>
#>      IC50 50  10.8930 10.3520 11.5230 bootstrap order statistics  <NA>
#>  boot_mean boot_sd n_undefined
#>     8.5779 0.14997           0
#>    10.9060 0.29994           0
```

The estimates reproduce the manual exactly: 8.5715 for the IC25 against
a printed 8.5715, and 10.893 for the IC50 against 10.89. So do the
smoothed means, 28.75 across the control and the three lowest
concentrations, and the standard deviations the ICPIN output prints.

**Smoothing weights.** The means are smoothed with equal weight on each
concentration, not weighted by the number of replicates. That is what
reproduces the manual’s 28.75. With unequal replication the manual does
not say what it does, and this is a documented assumption rather than a
verified one.

**The interval cannot be reproduced, and not because of any
disagreement.** The manual’s limits of 8.3112 and 9.0418 come from 80
resamples drawn with the seed -641671986 by a Turbo Pascal generator,
which is stated in the output itself. No other program can return those
numbers. What can be checked is the shape of the interval: at 1000
resamples this package gives about 8.34 to 8.91, and the point estimate,
the smoothed means and the bootstrap standard deviation are all
reproducible under a seed.

The default of 200 resamples departs from the manual’s 80. The manual
itself warns that limits “based on the empirical quantiles of a
bootstrap distribution of 80 samples may be unstable”, and permits up to
1000.

### 17. Williams’ test: computed critical values, and validation without a reference

Williams’ test is the one method here with no EPA worked example and no
retrievable table, which makes it the hardest to validate and the
easiest to get quietly wrong.

**It is labelled as an extension, everywhere.** Section 9.4.1.2 of the
chronic manual mentions Williams’ test only as an alternative that
“requires additional assumptions”. Neither manual gives a worked
example, a critical value table, or a place for it on any flowchart. It
is in the package because ToxCalc offered it.
[`lc50()`](https://open-aims.github.io/toxstats/reference/lc50.md) and
[`tox_test()`](https://open-aims.github.io/toxstats/reference/tox_test.md)
never select it; reaching it always counts as an override and always
warns.

**Critical values are simulated.** Williams tabulated them in
*Biometrics*, so by decision 1 they are not transcribed. The null
distribution of the statistic is simulated instead for the design in
hand. That also removes the two errors the tables invite: the critical
value depends on the **position** of the concentration being tested, not
on the total number of concentrations, and the tables need interpolation
in `1/nu` rather than in `nu`.

The simulation is seeded by default, so the same design always returns
the same critical value, and it restores the caller’s random stream
afterwards so that a user’s own analysis is not shifted by having called
it.

**Validation without a reference.** With nothing published to compare
against, the implementation is pinned by three mathematical facts.

With one concentration there is no order restriction, so the statistic
is an ordinary one-sided `t` and the simulated critical value must equal
`qt(1 - alpha, nu)`. It does, to within Monte Carlo error, at every
degrees of freedom tested. That single identity is what establishes that
the simulation machinery is right; the rest rests on it.

Beyond one concentration the critical value must rise with position, and
must sit below Dunnett’s for the same design, because the order
restriction is information Dunnett’s procedure does not use. On the
Appendix C design Williams’ values run 1.76, 1.83, 1.88, 1.87 against
Dunnett’s 2.356.

``` r

fit <- suppressWarnings(
  williams(fathead_c1, response = "weight", nsim = 20000)
)
fit$comparisons[, c("conc", "mean", "isotonic", "statistic", "critical",
                    "tested", "significant")]
#>   conc    mean isotonic statistic critical tested significant
#> 1   32 0.57525  0.57525 1.4855647 1.736204  FALSE       FALSE
#> 2   64 0.66025  0.61775 0.8665794 1.826941  FALSE       FALSE
#> 3  128 0.56500  0.56500 1.6348494 1.887311   TRUE       FALSE
#> 4  256 0.45425  0.45425 3.2478522 1.847940   TRUE        TRUE
```

**Two behaviours worth pointing out in that output.** The `isotonic`
column shows the order restriction at work: the observed means rise from
0.575 at 32 to 0.660 at 64, and the isotonic step replaces both by their
average. And `tested` is `FALSE` for the two lowest concentrations,
because the procedure steps **down** from the highest and stops at the
first that is not significant. Everything below is declared not
significant without being tested, which is what makes this a trend test
rather than a set of pairwise comparisons.

The package warns when the observed means are not monotone, as they are
not here. The isotonic step will absorb the departure silently
otherwise, and the assumption that justifies using this test in the
first place deserves to be seen.

## Still open

Two of the four questions raised when this package was planned have
since been settled by evidence and are described in sections 14 and 15:
the probit interval, and the trimmed Spearman-Karber interval, the
latter after an initial answer that turned out to be wrong. These
remain.

- **The expanded interval for fewer than seven replicates**, used by the
  linear-interpolation method. It is defined only in the ICPIN version
  2.0 program documentation, which is not publicly retrievable. The
  current intention is to omit it and report the ordinary interval with
  a note.
- **Non-monotone patterns of significance.** Section 9.6.5.1 advises
  caution where a lower concentration is not significant between two
  that are, but prescribes no rule. The package reports the pattern as
  observed, sets the no-observed-effect concentration below the lowest
  significant one, and flags the result rather than forcing it monotone.
- **The heterogeneity correction to the probit variance.** The manual’s
  practice of inflating the variance by `chi-square / df` when
  heterogeneity is significant is implemented, but the worked example
  does not exercise it — its chi-square is well below the tabular value
  — so it is unverified against any published output.

## Checking this package against the manuals

Every claim above is a test. To run them:

``` r

devtools::test()
```

The datasets are rebuilt by `data-raw/DATASET.R`, which re-checks each
transcription against the summary statistics printed beside it in the
manual and fails if any value has drifted.
