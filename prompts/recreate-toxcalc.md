## Session: recreate-toxcalc
Date: 2026-08-31
Model: claude-opus-5 (1M context)

### Prompts and Responses

**User:** I want to recreate the toxcalc software

**Claude:** Established that the target is ToxCalc v5.0 (Tidepool Scientific), a
commercial statistical package for environmental toxicity testing implementing
the US EPA Whole Effluent Toxicity (WET) methods, cited in
`bayesnec/article/article.Rmd:414` as proprietary software providing limited
flexibility and no reproducible open-source workflow. The vendor product page
now returns 404 and redirects to CETIS, so the software is withdrawn and the
authoritative specification is the EPA method manuals.

Four scoping decisions were put to the user, who chose: R package only (no Shiny
or Quarto); the full EPA flowchart, both the hypothesis-testing and
point-estimation branches; fidelity to the EPA manuals rather than a
bit-for-bit ToxCalc clone; and validation fixtures transcribed from the public
EPA manuals, no ToxCalc outputs being available. Two further decisions were
taken: Williams' critical values from transcribed EPA tables with a computed
fallback, and bare function names matching the `toxval` convention.

An audit of the installed R library (642 packages, R 4.5.1) established what
must not be reimplemented — `stats` for the assumption tests, `kSamples` for
Steel's Many-One Rank Test and Jonckheere-Terpstra, `mvtnorm` for exact Dunnett
critical values, `MASS::dose.p` for probit — and what has no implementation
anywhere and must be written: Williams' test, trimmed Spearman-Karber with
auto-trim, ICp linear interpolation with bootstrap limits, Abbott's correction,
MSD and PMSD, the EPA arc sine square root transformation, and
pool-adjacent-violators smoothing. `PMCMRplus` and `DescTools` are not installed
and are impractical to add on Windows (they require MPFR and GMP system
libraries), so the design avoids them.

A design pass reading the operative sections of both manuals returned material
corrections to the initial plan, recorded in `notes/TOXCALC-claude.md`. The
plan was revised accordingly.

**Phase 0 (scaffold) implemented.** `DESCRIPTION` with a four-package runtime
dependency list (`chk`, `kSamples`, `mvtnorm`, `stats`); `R/toxcalc-package.R`
and `R/toxcalc-params.R` following the `toxval` house conventions; `air.toml`,
`.editorconfig`, `.gitattributes`, `LICENSE.md`, `_pkgdown.yml`,
`scripts/build.R` and the three GitHub Actions workflows copied from `toxval`;
`.Rbuildignore` extended with `^notes$`, `^superceded$` and `^air\.toml$`, the
entries missing from `toxval`; repo-level `CLAUDE.md`; `NEWS.md`;
`tests/testthat.R`; `toxcalc.Rproj`.

---
**User:** yes please continue. With the open decisions, I am not a user so if you can consult the internet for best practice or some other evidence on the right call please go ahead and implement a solution, just document the decision points clearly in a vignette describing the steps used in the re-creation.

**Claude:** Implemented Phases 1 and 2, resolving each open decision against
evidence rather than referring it back, and recording every decision in a new
Quarto vignette, `vignettes/recreating-toxcalc.qmd`.

**Phase 1.** `toxcalc_data()` validates replicate-level data into one long form
covering both continuous and quantal endpoints, retaining the number exposed so
quantal data can feed both branches of the EPA analysis. `arcsine_sqrt()` and
`inv_arcsine_sqrt()` implement the Appendix B 4.2 transformation with the
Bartlett endpoint adjustments, reproducing the printed 0.8861, 0.1120 and
1.4588. `smooth_monotone()` implements pool-adjacent-violators with block sizes
tracked, rather than the pairwise averaging the manuals describe, because the
Appendix K example pools five values at once and pairwise averaging cannot
reproduce it. `abbott()` and `is_monotone()` complete the primitives.

**Phase 2.** The EPA chronic freshwater manual was downloaded and its text
extracted, and the worked examples were read directly rather than recalled. Three
datasets were transcribed with `data-raw/DATASET.R`, which re-checks each against
the summary statistics printed beside it and fails on any drift.
`epa_normality()` applies Shapiro-Wilk or the Kolmogorov D statistic to the
pooled within-group centred residuals; `epa_variance()` applies Bartlett's test
with Levene and Fligner offered and flagged as non-EPA; `msd()` and `pmsd()`
give the minimum significant difference and its percentage of the control mean,
with the Table 6 bounds shipped as `epa_pmsd_bounds`.

**Decisions taken on evidence.**

*Shapiro-Wilk.* `stats::shapiro.test` (Royston AS R94) returns 0.9601 against
the manual's printed 0.959, and 0.9594 once the manual's own rounding of the
concentration means is applied. The agreement removed any reason to also
implement Conover's tabled route, which would have required transcribing two
third-party tables.

*Dunnett's critical value.* `mvtnorm::qmvt()` was found to be randomised,
returning values from 2.3552 to 2.3585 across calls on the same design. Since a
regulatory analysis must be reproducible, the value is instead obtained by
numerical integration: the comparisons share a control, so the correlation
matrix is one-factor and the multivariate t probability reduces to a
two-dimensional integral. This returns 2.35614 identically on every call,
against the manual's tabled 2.36, and removed `mvtnorm` from Imports. Runtime
dependencies are now `chk` and `stats` only.

**Four defects found in the manual, each confirmed numerically and reconciled by
a test.** Appendix C section 1.11 gives the MSD as 0.087 where the same section
computes 0.162; Appendix C states it uses the same data as Appendix B, which it
does not, so both tables are shipped separately; the printed Bartlett statistic
of 7.691 follows only from the manual's own variances rounded to four decimal
places, one of them misrounded, and a correction factor rounded to 1.133, where
the raw data give 6.836; and the printed Kolmogorov `D*` of 0.4684 follows only
from rounding each z to two decimal places for a printed normal table, where the
correct value is 0.4572. In every case the conclusion the manual draws is
unchanged.

The vignette is a Quarto document, built with the Quarto CLI bundled with
Positron; `quarto-dev/quarto-actions/setup@v2` was added to the R-CMD-check and
pkgdown workflows so CI can build it.

R CMD check clean at 0/0/0 throughout; 175 assertions passing.

---
**User:** Yes, continue

**Claude:** Implemented Phase 3, the seven hypothesis tests the EPA flowchart
can select, each validated against its worked example in the manual.

`R/parametric.R` adds `dunnett()`, `bonferroni_t()`, `dunn_sidak_t()` and
`welch_t()`; `R/nonparametric.R` adds `steel()`, `wilcoxon_rank_sum()` and
`fisher_exact()`; `R/comparisons.R` holds the shared `toxcalc_comparison`
class, the NOEC and LOEC derivation, and the print method. Two datasets were
added, `ceriodaphnia_e1` and `ceriodaphnia_g2`.

**Reproduction of the worked examples.** Appendix C: `t = 1.486, 0.248, 1.635,
3.248` against a printed 1.487, 0.248, 1.633, 3.251, critical value 2.3561
against a tabled 2.36, NOEC 128 and LOEC 256. Appendix D: critical value 2.510
exactly, NOEC 128 and LOEC 256. Appendix E: rank sums 84, 63.5, 76 and 55, NOEC
3 and LOEC 6. Appendix F: rank sums 79, 57, 58 and 55 exactly, NOEC 3 and LOEC
6. Appendix G: NOEC 12 and LOEC 25.

**Decisions taken.**

*Steel's test.* The manual compares each rank sum with a tabulated critical
value (Table E.5), which section 3 of the specification forbids transcribing.
`kSamples::Steel.test()` instead evaluates the joint null distribution of all
`k` comparisons and returns multiplicity-adjusted p-values. These reproduce the
manual's conclusion exactly, including the 12 per cent comparison whose rank sum
equals the tabulated critical value, where the adjusted p-value is 0.0477.

*Wilcoxon.* Same reasoning for Table F.5: the p-value comes from
`stats::wilcox.test()` and is compared with `alpha / k`. Reproduction counts tie
constantly, so an exact p-value is unavailable and the normal approximation is
used; the resulting warning is muffled deliberately and documented.

*Fisher's Exact Test.* No multiplicity adjustment is applied by default, because
Appendix G section 1 asks for a pairwise error rate rather than an
experiment-wise one on the grounds that the test is itself conservative.

*Bonferroni critical values.* EPA Table D.5 turns out to be exactly
`qt(1 - alpha / k, df)`, giving 2.510 for the worked example and 2.241 against a
printed 2.242 in the row for infinite degrees of freedom, so the table needs no
transcription either.

*Adjusted p-values.* Each parametric test reports a multiplicity-adjusted
p-value rather than the raw one, so that comparing it with alpha gives the same
answer as comparing the statistic with the familywise critical value. For
Dunnett's procedure the adjustment comes from the same integral that produces
the critical value.

*Refusing the wrong test.* `dunnett()` and `steel()` signal an error on an
unbalanced design rather than returning a value the manual would not accept, and
both non-parametric tests refuse fewer than four replicates per concentration
(section 9.4.5.2).

**Three further defects found in the manual**, bringing the total to seven, each
confirmed numerically and reconciled by a test. Appendix D Table D.3 gives the
within mean square as 0.0029 where its own sum of squares and degrees of freedom
give 0.0079, which is what its printed t values were computed from. Appendix D
Table D.4 gives `t = 0.220` at 64 µg/L where the manual's own means and pooled
standard deviation give 0.270. Appendix E Table E.4 gives the 6 per cent rank sum
as 64 where the ranks listed in its own Table E.3 sum to 63.5.

`scripts/build.R` now installs the package before checking, because Quarto
renders the vignette in a fresh session against the installed package rather
than the source tree.

R CMD check clean at 0/0/0; 254 assertions passing.

---
**User:** yes continue

**Claude:** Implemented Phase 4, the flowchart engine and the `toxcalc()`
driver. This is the part that actually recreates what ToxCalc did: selecting a
test automatically and recording why.

**The chart is data, not control flow.** `R/flowchart.R` holds the EPA Figure 2
chart as a data frame of branch points, with columns `node`, `question`,
`rule`, `yes`, `no` and `reference`, and `walk_flowchart()` traverses it. The
alternative, nesting `if` statements in the driver, was rejected for three
reasons: each branch point can be tested on its own; the traversal can be
recorded without the recording being tangled into the logic, which is what makes
the audit trail trustworthy rather than a narration maintained alongside the
code; and the three further charts in the acute manual can be added without
touching the engine. The manual section for each branch point lives in the chart
itself, so the citation in the printed output cannot drift from the code.

**Two details of Figure 2 that are easy to misread** are encoded explicitly. The
normality branch goes directly to the replication question when it fails, so the
variance test is not reached at all. And there is an `@no_valid_test` terminal,
because a design whose residuals are not normal and which has fewer than four
replicates has no test the manual sanctions, section 9.4.5.2 having ruled out
the non-parametric branch; that case raises an error rather than falling through
to something plausible.

**Both worked examples run end to end, and the chart independently selects the
same test the manual selects.** For the Appendix C growth data the residuals are
normal (W = 0.951, p = 0.378) and the variances homogeneous (B = 7.856,
p = 0.097) with balanced replication, so Dunnett's procedure is selected,
returning NOEC 128, LOEC 256, MSD 0.162 and PMSD 23.9 per cent, within the
Table 6 bounds. For the Appendix E reproduction data with the 50 per cent
concentration excluded, normality fails (W = 0.928, p = 0.0047) with ten
balanced replicates, so Steel's Many-One Rank Test is selected, returning
NOEC 3 and LOEC 6.

**Two rules that overturn a test result** were implemented, neither of which
appears on the flowchart. Section 9.5.2 excludes a concentration at which every
replicate showed a complete response, applied automatically for a quantal
endpoint, and excludes from a sublethal analysis any concentration at which
survival was significantly reduced, which needs a separate survival analysis and
so is passed through `exclude`. Section 10.2.8.2.5 states that a concentration
shall not be declared toxic if its relative difference from the control is below
the lower PMSD bound, whatever the test concluded; this reverses a significant
result after the fact, and every reversal gets its own row in the decision
trail.

**Overrides are possible but never quiet.** Passing `test` runs that test
whatever the chart selected; when the two disagree a warning is raised,
`overridden` is set, a row is added to the trail, and `print()` reports it above
the endpoints. A forced analysis can therefore never be read as the analysis the
manual prescribes.

The two summary outputs are covered by testthat snapshot tests, a departure from
`toxval`, which uses none, but appropriate for report-style output whose exact
wording is part of what is being delivered.

R CMD check clean at 0/0/0; 338 assertions passing.

---
**User:** yes continue

**Claude:** Implemented Phase 5, the LC50 point-estimation branch: the
graphical method, Spearman-Karber, trimmed Spearman-Karber with automatic trim,
probit, and the Figure 6 chart that selects between them. The acute manual
(EPA-821-R-02-012) was downloaded and its section 11.2 worked examples read
directly.

**All four columns of Table 20 route to the method they were built for and
reproduce the published estimate.** Graphical 35.36, where the manual reads 35
off a plot and the geometric mean of the two bracketing concentrations is
35.355. Spearman-Karber `m = 1.65652` against a printed 1.656527,
`V(m) = 0.001098` against 0.0010977, and 45.3 (38.9, 52.8) exactly. Trimmed
Spearman-Karber trim 20.51 per cent and LC50 77.11, both exact. Probit LC50
22.872 (18.787, 27.846) and LC1 7.924 (4.147, 10.959), all exact.

**Four ambiguities resolved by testing candidates against the printed output
rather than by choosing a convention.**

*The probit heterogeneity statistic is Pearson's, not the deviance.* Both are
chi-square statistics on the same fit and the manual names neither. Pearson
gives 3.0764 against a printed 3.076; the deviance gives 3.8588.

*The probit limits are Fieller's, not the delta method.* Fieller reproduces the
printed limits to every digit for both the LC1 and the LC50; the delta method
gives 19.04 to 27.47 and 5.17 to 12.14, matching neither.

*The Spearman-Karber multiplier is 2.0.* Section 11.2.3.3 step 6 writes it that
way, which is the manual's rounding of 1.96, and keeping it is what reproduces
the printed interval. The exact normal quantile is used at other confidence
levels.

*Partial responses are counted after smoothing and Abbott adjustment.* Figure 6
branches on how many concentrations show a partial response without saying which
proportions are meant. The manual's own Spearman-Karber example settles it: the
raw data show two partial responses, which would route to the probit branch,
while the smoothed adjusted data show one, which routes to Spearman-Karber, and
that is the method the manual uses. Counting the raw proportions misroutes that
example, which is how the ambiguity was found.

**The trimmed Spearman-Karber interval is the one quantity that does not
reproduce a published EPA value, and this is documented prominently rather than
buried.** The trim rule was verified exactly against both published trims,
20.51 per cent in the acute manual and 20.41 per cent in Appendix K of the
chronic manual, confirming both the definition `max(p1, 1 - pk)` and the
ordering: smooth, then Abbott-adjust, then trim. The point estimate is exact.
The interval is not, because the manual delegates it to a program whose formula
it does not state, and the corrected variance expression of Hamilton et al.
(1978) is not retrievable from any public source. Three candidates were tested
against the printed 69.74 to 85.26: the plain Spearman-Karber variance scaled by
`1/(1-2a)^2` gives 70.8 to 84.0, clearly too narrow; a parametric bootstrap
gives 68.9 to 86.3, too wide; the delta method applied to the estimator as
implemented gives 69.6 to 85.4. The delta method was chosen because it is exact
for the estimator actually implemented rather than adopted on faith, and is the
closest of the three.

Two implementation details worth recording. The probit model is fitted as
proportions with the group sizes as weights rather than as a two-column count
matrix, because the Abbott adjustment makes the implied counts non-integer.
And two glm warnings are muffled deliberately, the separation warning and the
non-integer successes warning, both of which fire on almost every EPA design
and say nothing actionable; separation is recorded on the result instead.

R CMD check clean at 0/0/0; 437 assertions passing.

---
**User:** continue with phase 6

**Claude:** Implemented Phase 6, `icp()`, the EPA linear interpolation method
for estimating an inhibition concentration, with bootstrap confidence limits.
Appendix M of the chronic manual and the ICPIN program output reproduced there
were read directly.

**The estimate reproduces the manual exactly.** The smoothed means come out at
28.75 across the control and the three lowest concentrations, then 9.4 and 0,
matching Table M.2. The concentration standard deviations match those the ICPIN
output prints. `IC25 = 8.5715` against a printed 8.5715 and `IC50 = 10.893`
against 10.89.

Smoothing uses pool-adjacent-violators with equal weight on each concentration
mean rather than weighting by replication, which is what produces 28.75. With
unequal replication the manual does not say what it does, so this is recorded as
a documented assumption rather than a verified one.

**A correction to the original plan.** The specification predicted that
omitting the re-smoothing from inside the bootstrap loop would make the interval
"far too narrow". Measuring it showed the opposite and worse: it biases the
estimate. On the Appendix M data the resampled estimates then average about 10.4
against a point estimate of 8.57, with a standard deviation of 0.49 rather than
0.14, because unsmoothed resampled means are not monotone and the interpolation
selects the wrong bracketing pair. The prediction was wrong in direction, not
just in degree, and the test now asserts the measured bias rather than the
predicted narrowing.

**The interval cannot be reproduced, and not because of any disagreement.** The
ICPIN output states the seed it used, -641671986, for 80 resamples drawn by a
Turbo Pascal generator. No other program can return its printed limits of 8.3112
and 9.0418. What is reproducible, and what the tests check, is the point
estimate, the smoothed means, the bootstrap standard deviation under a seed, and
the order-statistic rule the manual describes as "the second smallest and second
largest" of 80 resamples, generalised to `floor((nboot + 1) * (1 - ci_level)/2)`
in from each end. At 1000 resamples the interval comes out at about 8.34 to
8.91 against the manual's 8.31 to 9.04.

The default of 200 resamples departs from the manual's 80, on the manual's own
grounds: it warns that limits based on 80 samples "may be unstable" and permits
up to 1000.

The expanded interval the manual reports when there are fewer than seven
replicates per concentration is not implemented, its definition appearing only
in the ICPIN version 2.0 documentation, which is not publicly retrievable. A
warning is raised when a design would have triggered it, rather than the
omission being silent.

R CMD check clean at 0/0/0; 474 assertions passing.

---
**User:** Yes continue

**Claude:** Implemented Phase 7, `williams()`. This is the one ToxCalc feature
that is not an EPA method, and the hardest to validate: neither manual gives a
worked example, a critical value table, or a place for it on any flowchart, and
Williams' own tables in *Biometrics* are not retrievable.

**It is labelled as an extension everywhere.** Section 9.4.1.2 of the chronic
manual mentions the test only as an alternative that "requires additional
assumptions". The method string, the reference string, and the flowchart
behaviour all say so: neither `toxcalc()` nor `lc50()` will ever select it, and
forcing it always counts as an override and always warns.

**Critical values are simulated, not tabulated,** per the copyright decision
taken in Phase 2. That also removes the two errors the tables invite: the
critical value depends on the *position* of the concentration being tested,
not on the total number of concentrations, which is the commonest
implementation error; and the tables require interpolation in `1/nu` rather
than `nu`. One simulation serves every position, using common random numbers,
which is faster and keeps the critical values coherent in `i` rather than
slightly non-monotone through Monte Carlo noise. The simulation is seeded by
default so a given design always returns the same value, and it restores the
caller's random number stream on exit so that a user's own analysis cannot be
shifted by having called it.

**Validation without a reference.** With nothing published to compare against,
the implementation is pinned by three mathematical facts rather than by
assertion. With one concentration there is no order restriction, so the
statistic is an ordinary one-sided t and the simulated critical value must equal
`qt(1 - alpha, nu)`; it does at every degrees of freedom tested. Beyond one
concentration the critical value must rise with position, and must sit below
Dunnett's for the same design, because the order restriction is information
Dunnett's procedure does not use. On the Appendix C design the values run 1.76,
1.83, 1.88 and 1.87 against Dunnett's 2.356. The isotonic step and the
step-down procedure are tested directly.

Two behaviours worth recording. The isotonic step absorbs departures from
monotonicity, which on the Appendix C data means the means at 32 and 64 are
replaced by their average; the package warns when that happens, because the
assumption that justifies using the test at all should not be absorbed
silently. And the step-down means concentrations below the first
non-significant one are never tested, which the result records in a `tested`
column so the distinction between "tested and not significant" and "not tested"
is visible rather than implied.

A p-value is reported from the simulated null distribution, which is more
informative than reporting only which side of the critical value the statistic
fell.

One test expectation of mine was wrong on the first run: I asserted that the
weighted isotonic estimate of `c(4, 10)` with weights `c(1, 3)` was 5.5 when it
is 8.5. The implementation was right and the test was wrong.

R CMD check clean at 0/0/0; 514 assertions passing.

---
**User:** Yes go ahead and do phase 8

**Claude:** Implemented Phase 8: a driver that runs both branches, a tidy
output shape, an executable README and a grouped pkgdown reference index.

**Both branches from one call.** `toxcalc(branch = "both")` returns the NOEC,
LOEC, MSD and PMSD alongside a point estimate, choosing `lc50()` for a quantal
endpoint and `icp()` for a continuous one. That is what a laboratory report
usually contains. The point branch deliberately uses the **original** data,
including any concentration the hypothesis branch excluded, which is what
section 9.5.2 requires and is asserted in a test.

`as.data.frame()` now returns one tidy row per endpoint rather than one wide
row, which is the only shape that works once two branches can report together.

**Two design corrections came out of testing.**

A `branch = "point"` option was implemented first and then removed. It
duplicated `lc50()` and `icp()`, which are easier to find, and it left a
half-populated object carrying a NULL hypothesis test. Offering only
`"hypothesis"` and `"both"` is cleaner.

The same test exposed a genuine robustness defect: a design with a single
replicate at every concentration produced no within-concentration variation, and
`stats::shapiro.test` failed with "all 'x' values are identical". That now
raises a clear error explaining that such a design supports no hypothesis test
at all while the point-estimation methods, which work on pooled counts, handle
it. `acute_table20` is exactly such a design, so the case is real rather than
hypothetical.

**Documentation.** `README.Rmd` replaces the hand-written README, so the audit
trail, the tidy output and the Figure 6 routing shown there are executed rather
than pasted and cannot drift from the code. The pkgdown reference index is
grouped by branch of the EPA analysis, so it reads in the same order as the
manuals; `pkgdown::check_pkgdown()` reports no problems. The pkgdown site URL
was added to `DESCRIPTION` because pkgdown requires it there.

R CMD check clean at 0/0/0; 536 assertions passing.

---
**User:** let's use toxstats

**Claude:** Renamed the package from `toxcalc` to `toxstats`, and the interface
from `shinytoxcalc` to `shinytoxstats`.

**Why the original name had to go.** `toxcalc` collided with two things. The
Tidepool commercial product it recreates, whose trademark position could not be
cleared from public sources, and `USEPA/toxcalc`, an archived EPA repository for
an unrelated human-health dose-conversion calculator. Adding an R suffix was
considered and rejected: infringement turns on likelihood of confusion rather
than identity, and "ToxCalc for R" reads as an official port, which implies more
affiliation rather than less. It would also be inaccurate, since this is an
independent implementation from the EPA manuals rather than a port of anything.

**Why an intermediate name was also rejected.** `wetstats` was implemented and
then abandoned on the user's objection that it implies the package serves only
Whole Effluent Toxicity testing. That objection is correct. Dunnett's procedure,
Steel's test, probit, the trimmed Spearman-Karber method and linear
interpolation apply to any concentration-response toxicity data with a control.
The EPA WET manuals are where the methods and the selection rules are
*specified*, not the only setting they apply to. The Title and Description were
rewritten accordingly, and now lead with what the package computes rather than
with the regulatory context it came from.

`toxcalc()` became `tox_test()`, `toxcalc_data()` became `tox_data()`, and the
classes followed. Every reference to **ToxCalc the software** was deliberately
left alone: the package still describes itself as recreating it, which is
nominative fair use and belongs in the description rather than in the name.

Two renames in succession across nineteen R files, the test suite, the vignette
and both notes documents, with 536 and 101 assertions passing after each and
R CMD check clean on both packages. The test suites are what made a mechanical
rename of this size safe to attempt at all.

---
