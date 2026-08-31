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
