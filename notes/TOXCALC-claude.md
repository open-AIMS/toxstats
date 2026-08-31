# toxcalc — specification

Implementation specification for `toxcalc`. Written to be implemented from
directly. The collaborator-facing companion is `notes/TOXCALC-human.md`; every
decision is stated in full **here** and only summarised there.

Status: Phases 0 to 3 complete. Phases 4-9 specified below. Decisions taken
since this document was first written are recorded in the package vignette
`vignettes/recreating-toxcalc.qmd`, which is the collaborator-facing record.

---

## 1. Scope

Recreate the statistical capability of ToxCalc v5.0 (Tidepool Scientific), which
implemented the US EPA Whole Effluent Toxicity (WET) methods. Fidelity target is
the EPA method manuals, validated against the worked examples printed in them.
No ToxCalc outputs are available for comparison.

### Sources, and what was read

The operative sections of both primary manuals were extracted and read directly
rather than recalled: EPA-821-R-02-013 (chronic freshwater, 350 pp) and
EPA-821-R-02-012 (acute, 275 pp) — Figures 2, 6, 12 and 13, Appendices B, C, D,
E, G, K, L and M, and Sections 9.4-9.6 and 10.2.

- EPA-821-R-02-012, acute, 5th ed., October 2002.
- EPA-821-R-02-013, chronic freshwater, 4th ed., October 2002.
- EPA-821-R-02-014, chronic marine and estuarine (Phase 9).
- Williams (1971) *Biometrics* 27:103-117; Williams (1972) *Biometrics* 28:519-531.
- Hamilton, Russo & Thurston (1977) *Environ. Sci. Technol.* 11:714-719, with the
  1978 correction, 12:417.
- Norberg-King (1993), EPA/600/M-91/037, ICPIN v2.0.

---

## 2. Ground truth extracted from the manuals

These are the facts that drive the design. Several correct assumptions made
before the manuals were read; those corrections are marked.

### 2.1 The hypothesis-testing flowchart

Chronic Figure 2 and acute Figure 13 are the **same chart**:

```
response
  |
  +-- proportion? -> arc sine square root transform (App B 4.2)
  |
  +-- Shapiro-Wilk normality, alpha = 0.01,
  |     on POOLED WITHIN-GROUP CENTERED RESIDUALS      (App B 2.1)
  |     |
  |     +-- fail ------------------------------+
  |     |                                      |
  |     +-- pass -- Bartlett, alpha = 0.01     |
  |                   |                        |
  |                   +-- fail ----------------+
  |                   |                        |
  |                   |                        v
  |                   |             equal replicates?
  |                   |                yes -> Steel's Many-One Rank
  |                   |                no  -> Wilcoxon Rank Sum + Bonferroni
  |                   |
  |                   +-- pass -- equal replicates?
  |                                  yes -> Dunnett's
  |                                  no  -> t-test + Bonferroni
```

**Corrections to earlier assumptions.**

1. Shapiro-Wilk is applied to the **pooled within-group centered residuals**,
   not to each group separately and not to raw values.
2. Alpha is **0.01 for the assumption tests and 0.05 for the hypothesis tests**
   (9.4.6.1), one-sided. A single `alpha` argument is wrong; there are two.
3. The normality-fail arm goes **directly** to the replication question. It does
   not first test variance.
4. Nonparametric tests require **at least 4 replicates per concentration**
   (9.4.5.2). A design that fails normality with 3 replicates has **no
   EPA-sanctioned test** — this must be an explicit terminal, not a fall-through.
5. **Williams' test is not on this chart.** The manual names it only as an
   alternative "requiring additional assumptions" (9.4.1.2). It is in scope
   because ToxCalc offered it, but it is a documented extension, not the EPA
   method. This moves it out of the core phase.

### 2.2 The single-concentration / paired flowchart

Acute Figure 12 is a **different chart** and returns Pass/Fail, not a NOEC:

```
arc sine sqrt -> Shapiro-Wilk
                   fail -> Wilcoxon Rank Sum
                   pass -> F-test for equal variances
                             pass -> pooled t-test
                             fail -> Welch modified t-test
```

### 2.3 The LC50 flowchart

Acute Figure 6:

```
>= 2 partial mortalities?
  yes -> probit appropriate (chi-square heterogeneity not significant)?
           yes -> PROBIT
           no  -> fall through
  no  -> >= 1 partial mortality?
           no  -> GRAPHICAL METHOD
           yes -> fall through

fall through -> zero mortality at lowest conc AND 100% at highest?
                  yes -> SPEARMAN-KARBER
                  no  -> TRIMMED SPEARMAN-KARBER
```

**Correction.** The scope is larger than probit + TSK. The **graphical method**
and **plain untrimmed Spearman-Karber** are both required terminals of the EPA
chart and must be implemented.

### 2.4 Exclusion rules (9.5.2)

- Concentrations with significantly reduced survival are excluded from the
  **sublethal** hypothesis test.
- Concentrations at 100% mortality in all replicates are excluded from NOEC and
  LOEC, but **retained for point estimation**.
- Point estimation always uses an all-data approach.

### 2.5 PMSD and the lower-bound override (10.2.8)

`PMSD = 100 * MSD / control_mean`, computed **parametrically even when the
flowchart selected a nonparametric test**.

EPA Table 6 bounds (lower / upper), which are EPA-generated and therefore safe
to ship:

| Method | endpoint | lower | upper |
|---|---|---|---|
| Fathead minnow (1000.0) | growth | 12 | 30 |
| *Ceriodaphnia dubia* (1002.0) | reproduction | 13 | 47 |
| *Selenastrum* (1003.0) | growth | 9.1 | 29 |

**Rule 10.2.8.2.5, easy to miss:** a concentration **shall not be declared
toxic** if its relative difference from the control is less than the lower PMSD
bound. This is a post-hoc override of the significance result and must be
applied after the test, with its own row in the decision trail.

### 2.6 Manual errata, all confirmed numerically

**Every fixture must be checked for internal consistency before being frozen
into a test, not transcribed on faith.** Seven defects have been confirmed so
far, all in EPA-821-R-02-013. Each is reconciled by a test in the suite, so none
can later be mistaken for an implementation error.

1. **Appendix C section 1.11 gives the MSD as 0.087 mg.** The same section
   computes 0.162, and 24% of the 0.677 control mean is 0.162. The 0.087 is the
   typo. *Confirmed:* recomputing from Table C.1 gives `Sw = 0.0971`,
   `MSD = 0.1619`, `PMSD = 23.9%`.

2. **Appendix C states its data are "the same data used in Appendices B and
   D". They are not.** Table C.1 differs from Table B.1 at most values. Table
   B.1 has a control mean of 0.714 and `Sw = 0.052`; Table C.1 has 0.677 and
   0.097, and only Table C.1 reproduces the printed `t` values, critical value
   and MSD. Both are shipped separately as `fathead_b1` and `fathead_c1`.

3. **Bartlett `B = 7.691` (Appendix B 3.5) is not reproducible from the raw
   data, which give 6.836.** *Confirmed:* the manual computes from its own
   printed variances rounded to 4 dp, giving `sum(log) = -32.4771` exactly, and
   then divides by `C` rounded to 1.133. One of those variances is itself
   misrounded: 0.002055 is printed as 0.0020. Bartlett depends on the log of
   each variance, so the smallest variance dominates the error.

4. **Kolmogorov `D* = 0.4684` (Appendix B 2.10) is not reproducible; the
   correct value is 0.4572.** *Confirmed:* rounding each `z` to 2 dp, as
   required to use a printed normal table, reproduces `D = 0.0597` and
   `D* = 0.4683` against the printed 0.4684.

5. **Appendix D Table D.3 gives the within mean square as 0.0029.** The same
   table gives `SSW = 0.111` on 14 degrees of freedom, and 0.111/14 = 0.0079.
   *Confirmed:* the manual's own printed `t` values follow from 0.0079, so the
   printed mean square is the error and everything downstream of it is right.

6. **Appendix D Table D.4 gives `t = 0.220` at 64 ug/L; it is 0.270.**
   *Confirmed:* from the manual's own means (0.677, 0.660) and pooled standard
   deviation the statistic is 0.2704. The other three comparisons agree with
   the manual to within rounding.

7. **Appendix E Table E.4 gives the 6% rank sum as 64; it is 63.5.**
   *Confirmed:* the ranks the manual itself lists in Table E.3 for that
   concentration sum to 63.5. Both are at or below the critical 76.

In all seven cases the conclusion the manual draws is unchanged.

Two further printed values agree once the manual's own rounding is applied and
are therefore **not** errata: Shapiro-Wilk `W = 0.959` (Royston gives 0.9601
unrounded, 0.9594 with the manual's 3 dp concentration means) and the Dunnett
critical value 2.36 (integration gives 2.3561).

---

## 3. Statistical tables and third-party copyright

The EPA manuals are US Government works, but several tables **within** them are
reproduced from third-party sources under a permission that does not transfer:
Williams' critical values from *Biometrics*; Dunnett's Table C.5 credited to
Miller (1981); Steel's rank-sum tables; Conover's Shapiro-Wilk coefficients
(Table B.4) and quantiles (Table B.6).

**Decision: no third-party table is transcribed into this package.** Critical
values are computed instead:

| Need | Computed route | Why it is at least as good |
|---|---|---|
| Dunnett critical value | Numerical integration of the one-factor multivariate t (`stats::integrate` + `stats::uniroot`) | Exact and **deterministic**; covers the unbalanced case, which Table C.5 does not. `mvtnorm::qmvt()` was rejected: it is randomised and returned 2.3552 to 2.3585 across calls on the same design. The reproduced table in the source PDF also has visible OCR corruption. |
| Shapiro-Wilk | `stats::shapiro.test` (Royston AS R94) | Agrees with Conover's tabled `W` to about three decimals. |
| Williams' `t-bar` | Monte Carlo null quantile for the observed design | No transcription risk, no interpolation risk, and handles designs the published tables never covered. |
| Steel's Many-One Rank | `kSamples::Steel.test` | Implements the Appendix E procedure with an exact or simulated joint null, better than the tabled approximation. |

Only **EPA-generated** tables are shipped: Table B.11 (Kolmogorov D critical
values, 1.035 / 0.955 / 0.895 / 0.819 / 0.775) and Table 6 (PMSD bounds).

This **reverses** the earlier decision to transcribe the Williams appendix
tables. The reversal is on copyright grounds first and accuracy grounds second.

---

## 4. Data structures

### 4.1 Input

One long data frame, one row per replicate test chamber. This single shape
covers both endpoint families, because EPA analyses quantal data **twice** —
replicate-level proportions (arcsine-transformed) go through the hypothesis
flowchart, and concentration-pooled counts go into probit / SK / TSK. A
`(conc, affected, n_exposed)` frame supports both; a `(conc, response)` frame
supports only the first.

| column | continuous | quantal |
|---|---|---|
| `conc` | numeric, >= 0 | numeric, >= 0 |
| `response` | measured value | number affected (integer) |
| `n_exposed` | `NULL` | integer, organisms in that chamber |
| `replicate` | optional label | optional label |

```r
toxcalc_data(
  data,
  conc           = "conc",
  response       = "response",
  n_exposed      = NULL,
  replicate      = NULL,
  type           = c("continuous", "quantal"),
  direction      = c("decreasing", "increasing"),
  control        = 0,
  conc_units     = "%",
  response_units = NULL
)
```

Fields of the returned object:

- `replicates` — the validated input frame, ordered by `conc`
- `proportion` — `affected / n_exposed` per replicate (quantal only)
- `pooled` — per concentration: `n_rep`, `n_exposed`, `n_affected`, `mean`,
  `sd`, `var`

Validation, with `chk`:

```r
chk::chk_data(data)
chk::chk_string(conc); chk::chk_subset(conc, names(data))
chk::chk_numeric(data[[conc]]); chk::chk_gte(min(data[[conc]]), 0)
chk::chk_number(control); chk::chk_subset(control, unique(data[[conc]]))
chk::chk_range(alpha, c(0, 1)); chk::chk_range(p, c(1, 99))
```

### 4.2 The returned analysis object

```r
structure(
  list(
    data       = <toxcalc_data>,
    design     = list(n_conc, n_rep, balanced, control_level, control_response,
                      n_partial_mortality, monotone, min_reps),
    hypothesis = <toxcalc_hypothesis> | NULL,
    point      = <toxcalc_point> | NULL,
    decisions  = <data.frame>,
    flowchart  = "epa_hypothesis_multi",
    overridden = FALSE,
    alpha, alpha_assumption, call, version
  ),
  class = "toxcalc"
)
```

`toxcalc_hypothesis`:

```r
list(
  transform   = list(name, applied, fun, inverse, n),
  assumptions = list(normality = <htest>, variance = <htest>),
  test        = "steel",
  comparisons = data.frame(conc, n, mean, se, statistic, critical,
                           p_value, significant, direction),
  noec, loec, msd, pmsd,
  msd_untransformed, pmsd_untransformed,     # App C 1.11.2 back-transform
  pmsd_bound_status = "within" | "above_upper" | "below_lower" | NA,
  excluded    = data.frame(conc, reason),
  decisions   = <data.frame>
)
```

`toxcalc_point`:

```r
list(
  method    = "trimmed_spearman_karber",
  estimates = data.frame(endpoint, p, estimate, lower, upper, ci_method),
  working   = data.frame(conc, log10_conc, n, affected,
                         observed, smoothed, adjusted),
  trim      = 0.2051,            # TSK only
  fit       = <glm> | NULL,      # probit only
  heterogeneity = <htest> | NULL,
  boot      = list(replicates, n_na, seed),   # ICp only
  decisions = <data.frame>
)
```

**Design rule: no field is a printed string that cannot be recomputed.**
`print()` and `summary()` are pure renderers over these fields, so the printed
EPA-style report and the machine-readable `as.data.frame()` can never disagree.

---

## 5. The flowchart engine

### 5.1 Charts are data, not control flow

Encoding the decision trees as nested `if`/`else` would make them untestable and
would force a code change per manual. Each chart is an internal data frame with
columns `node`, `question`, `rule`, `yes`, `no`, `reference`. Nodes prefixed `@`
are terminals naming a fit function.

```r
walk_flowchart(chart, state)
# -> list(decisions = <df>, terminal = "@steel", state = state)
```

It calls `rule_<name>(state)` per node; a rule returns
`list(answer, statistic, p_value, note, state)`. The engine appends one row and
follows `yes` / `no`.

Charts shipped in v1: `flowchart_hypothesis_multi` (Fig. 2 / 13),
`flowchart_hypothesis_paired` (Fig. 12), `flowchart_lc50` (Fig. 6),
`flowchart_icp` (App M).

### 5.2 The audit trail

| column | meaning |
|---|---|
| `step` | traversal order |
| `node` | stable machine key, e.g. `"normality"` |
| `question` | human sentence |
| `criterion` | `"Shapiro-Wilk W, reject if p <= 0.01"` |
| `statistic` | numeric or `NA` |
| `p_value` | numeric or `NA` |
| `answer` | `"yes"` / `"no"` / a value |
| `outcome` | what the branch implies |
| `reference` | manual section |

Rendered by `print()` as:

```
EPA WET analysis -- chronic multi-concentration hypothesis test
Flowchart: EPA-821-R-02-013 Figure 2

  1  Response is a proportion               -> arc sine square root transform applied
     (EPA-821-R-02-013 Appendix B, 4.2)
  2  Normality (Shapiro-Wilk)  W = 0.846, p = 0.0043 <= 0.01   -> NOT normal
     (EPA-821-R-02-013 Appendix B, 2.1)
  3  Equal replicates?  4, 4, 4, 4          -> yes
     (EPA-821-R-02-013 9.6.1.3)
  4  Selected test: Steel's Many-One Rank Test
     (EPA-821-R-02-013 Appendix E)

NOEC 12.5 %    LOEC 25 %    MSD 0.162 (24.0 % of control)
PMSD 24.0 % -- within EPA Table 6 bounds for Method 1000.0 (12 - 30)
```

### 5.3 Overrides must be loud

If the user forces `test = "dunnett"` but the chart terminal was `@steel`, then
`overridden` is `TRUE`, a `warning()` fires, and an extra decision row records
`node = "override"`. A forced non-compliant analysis can never print as if it
were compliant.

---

## 6. The from-scratch algorithms

### 6.1 Trimmed Spearman-Karber — auto-trim rule **verified**

*Source:* Hamilton et al. (1977, corr. 1978); EPA-821-R-02-012 §11.2.3-11.2.4;
EPA-821-R-02-013 Appendix K.

Let `p_0..p_k` be observed proportion mortalities at control plus `k`
concentrations, `x_i = log10(conc_i)`.

1. **Smooth** to non-decreasing by pooling adjacent violators
   (`p_{i-1} = p_i = (p_i + p_{i-1})/2`, iterating), **including the control**.
2. **Abbott-adjust:** `p_i^a = (p_i^s - p_0^s) / (1 - p_0^s)`.
3. **Auto-trim:** `trim = max(p_1^a, 1 - p_k^a)`.
4. Interpolate `x_alpha` and `x_{1-alpha}` linearly in `x` where `p^a` crosses
   `alpha` and `1 - alpha`; renormalise interior proportions
   `q_i = (p_i^a - alpha)/(1 - 2*alpha)`; apply the SK formula to the augmented
   set `{(x_alpha, 0), interior (x_i, q_i), (x_{1-alpha}, 1)}`:
   `mu = sum over consecutive pairs of (q_{i+1} - q_i) * (x_i + x_{i+1}) / 2`,
   `LC50 = 10^mu`.
5. Variance: `sum p_i^a (1 - p_i^a) (x_{i+1} - x_{i-1})^2 / (4 (n_i - 1))` over
   interior concentrations, divided by `(1 - 2*alpha)^2`.
6. CI on the log10 scale, back-transformed.

**The trim rule was verified numerically against the manual's own output, not
assumed.** Chronic Appendix K (40 organisms per concentration; mortalities
2, 0, 2, 0, 0, 32): pooling collapses `0.05, 0.00, 0.05, 0.00, 0.00` to `0.02`;
Abbott gives `p_1^a..p_4^a = 0` and `p_5^a = (0.80 - 0.02)/0.98 = 0.79592`; so
`trim = max(0, 1 - 0.79592) = 0.20408`, and the manual prints **20.41%**. The
acute Table 20 example gives **20.51%** by the same route. This confirms both
the rule and the **ordering — smooth, then Abbott, then trim**.

**Most likely to be got wrong:**

- Computing the trim from raw or un-Abbott-adjusted proportions. The two
  fixtures above discriminate this immediately; write them first.
- Treating the control as a point on the tolerance curve. It is used only to set
  `p_0^s` for the Abbott adjustment.
- Variance denominator `n_i` vs `n_i - 1`, and CI multiplier 1.96 vs 2.0. EPA
  prints `m +/- 2*sqrt(V(m))` for plain SK; whether TSK uses 1.96 is not stated.
  **Resolve empirically against the fixtures; do not assume.** See §8.3.
- Which concentrations enter the variance sum: only interior untrimmed ones,
  using the interpolated endpoints as neighbours for the `(x_{i+1} - x_{i-1})`
  differences.

### 6.2 ICp linear interpolation with bootstrap

*Source:* Norberg-King (1993); EPA-821-R-02-013 Appendix M; Marcus & Holtzman
(1988); Efron (1982).

1. Per-concentration means `Ybar_i`.
2. **Smooth** to non-increasing by pooling adjacent means.
3. `ICp = C_J + [M_1 (1 - p/100) - M_J] * (C_{J+1} - C_J) / (M_{J+1} - M_J)`,
   where `C_J, C_{J+1}` bracket `M_1 (1 - p/100)`.
4. If `C_J` is the highest concentration, report `> C_J`; if extrapolating below
   the lowest, report `< C_1`.
5. Bootstrap: for `R` resamples (80-1000, multiples of 40; >= 250 recommended),
   resample **replicate values with replacement within each concentration
   group**, `n_i` draws per group, control included; recompute means,
   **re-smooth**, re-interpolate. SE is the sd of `ICp*`; the CI is from
   empirical order statistics.
6. **The reported point estimate is from the original smoothed means**, not the
   bootstrap mean.

**Most likely to be got wrong:**

- **Smoothing outside the bootstrap loop.** It must be inside, or the interval
  is far too narrow. EPA's step list (7.2) makes smoothing step 4 and resampling
  step 6, which reads as outside; the semantics require inside.
- **The smoothing algorithm.** EPA's prose describes a single forward sweep, but
  both worked examples require *repeated* pooling (Appendix M pools four groups;
  Appendix L pools four), and the manual hedges: "Unusual patterns in the
  deviations from monotonicity may require an additional step of smoothing."
  **Implement full PAVA with equal weights on the group means, not weighted by
  `n_i`** — this reproduces both examples. Expose
  `smooth = c("pava", "forward")`. With unequal replicates ICPIN's actual
  weighting is unverifiable; document the choice.
- **Resampling the wrong unit.** Not concentrations, not the pooled pool —
  replicate values *within* group.
- **Undefined `ICp*` in a resample** (target never bracketed). ICPIN's handling
  is undocumented. Record `NA`, report `n_na`, warn above 5%.
- **Quantile definition.** "Second smallest of 80" is an order statistic, not
  `quantile(type = 7)`. Default `ci_method = "epa_order"`.

### 6.3 Williams' test

*Source:* Williams (1971, 1972). **Not an EPA flowchart method** — in scope
because ToxCalc offered it.

1. One-way ANOVA over all `k+1` groups: `s^2 = MSE`, `nu = N - (k+1)`.
2. Isotonic (amalgamated) dose means under `mu_1 >= ... >= mu_k` (decreasing
   response), via the max-min formula
   `M_i = min_{u<=i} max_{v>=i} ( sum_{j=u..v} n_j Ybar_j / sum_{j=u..v} n_j )`,
   i.e. weighted PAVA. **The control is excluded from the order restriction.**
3. `tbar_i = (Ybar_0 - M_i) / (s * sqrt(1/n_i + 1/n_0))`.
4. Compare to `tbar(alpha; i, nu)`, indexed by **`i`, the position of the dose
   being tested**, not by `k`.
5. Unequal `n`: the Williams (1972) correction
   `tbar' = tbar - (beta/100) * (1 - n_0/n_i)`.
6. **Step down:** test `i = k` first; on significance decrement `i`; stop at the
   first non-significant dose and declare all lower doses non-significant. The
   NOEC is the highest non-significant dose.

**Most likely to be got wrong:**

- **Indexing the critical value by `k` instead of `i`.** The commonest
  implementation bug in this test.
- **Amalgamating the control with the doses.** Wrong; the control is a free
  parameter.
- **Interpolating in `nu` instead of `1/nu`.** Williams' tables are tabulated at
  `nu = 5..20, 24, 30, 40, 60, 120, inf`.
- **The unequal-`n` `beta` correction**, whose exact form varies between
  secondary sources. With PMCMRplus unavailable as an oracle, validate against
  Williams' own published numerical example.
- **Forgetting the monotonicity precondition.** Williams assumes it; run the
  monotonicity check and warn.

Per §3, critical values are obtained by Monte Carlo, not from the tables.

---

## 7. Dependencies

```
Imports:  chk (>= 0.10.0), kSamples (>= 1.2.10), stats
Suggests: boot, covr, knitr, MASS, multcomp, mvtnorm, nortest, quarto, rmarkdown,
          testthat (>= 3.0.0), withr
```

Imports are added to `DESCRIPTION` in the phase that first uses them, so
`R CMD check` stays clean at every phase. Phase 0 declares `stats` only.

**Why neither `multcomp` nor `mvtnorm`.** `multcomp` would be used for exactly one
thing, the multivariate-t Dunnett critical value and adjusted p-value, which is
`mvtnorm::qmvt()` and `pmvt()` in about ten lines. `mvtnorm` has no further
dependencies; `multcomp` pulls `survival`, `TH.data`, `sandwich` and
`codetools`. `multcomp` stays in Suggests as a test oracle.

**Why `MASS`, `nortest` and `boot` are Suggests only.**

- `MASS::dose.p` gives the LC50 and a delta-method SE from a probit glm, but EPA
  additionally needs the Abbott adjustment, the chi-square heterogeneity test,
  the heterogeneity variance factor and Fieller limits — none of which
  `dose.p` does. The LC50 itself is `-b0/b1`.
- `nortest` provides Shapiro-Francia, Anderson-Darling and Lilliefors, **none of
  which EPA specifies**. EPA specifies Shapiro-Wilk (in `stats`) and the
  Kolmogorov D with its own Table B.11 critical values, which are the Stephens
  (1974) modified statistic and do not match `nortest::lillie.test`.
- `boot` cannot express the ICPIN interval: the resample is stratified within
  group, the smoothing is inside the loop, and `boot.ci` cannot produce the
  "second smallest of 80" order-statistic rule. About 15 lines of base R.

`car` and `drc` are not needed. EPA WET never fits a parametric sigmoid, and
`stats::fligner.test` covers the nonparametric variance alternative.

---

## 8. Open decisions

These are **not to be guessed**. Each needs a recorded answer before the phase
that depends on it. Resolved items are kept here with the evidence that settled
them; the same reasoning appears, less densely, in
`vignettes/recreating-toxcalc.qmd`.

1. ~~**Shapiro-Wilk implementation**~~ **RESOLVED (Phase 2): Royston only.**
   `stats::shapiro.test` returns `W = 0.9601` on the Appendix B example against
   the printed 0.959, and 0.9594 once the manual's 3 dp rounding of the
   concentration means is applied. The agreement removed the reason to
   implement Conover's tabled route as well, and implementing it would have
   meant transcribing Tables B.4 and B.6, which §3 forbids. Both `W` and the
   p-value are reported so a borderline case is visible.

2. ~~**Dunnett critical values**~~ **RESOLVED (Phase 2): numerical
   integration, not `mvtnorm`.** `mvtnorm::qmvt()` is randomised quasi-Monte
   Carlo and returned values from 2.3552 to 2.3585 across calls on the
   Appendix C design; a regulatory analysis must be reproducible. The `k`
   comparisons share a control, so the correlation is one-factor
   (`rho_ij = lambda_i lambda_j`, `lambda_i = sqrt(n_i/(n_0+n_i))`) and
   conditioning on the common factor and the pooled scale reduces the problem
   to a two-dimensional integral, evaluated by `stats::integrate` and inverted
   by `stats::uniroot`. Returns 2.35614 against the tabled 2.36, identically on
   every call. `mvtnorm` moved from Imports to Suggests, where it is a test
   oracle. See `R/critical.R`.

3. **TSK confidence interval convention** (Phase 5). Multiplier (1.96 vs 2.0)
   and variance denominator (`n_i` vs `n_i - 1`) are not both determinable from
   the manual text. *Plan:* fit them against the two published fixtures. **If no
   combination reproduces both (77.11 / 69.74 / 85.26, and the Appendix K case),
   escalate rather than pick the closer one.**

4. **Probit specification** (Phase 5). Three coupled ambiguities: (a) Abbott-
   adjust the proportions then fit, versus fit a natural-response parameter
   after Finney; (b) Fieller versus delta-method limits — EPA's output
   (22.872, 18.787, 27.846) is nearly symmetric on the log scale, consistent
   with either when `g` is small; (c) whether the heterogeneity factor
   `chi2/df` is applied to the variance when the chi-square is significant.
   *Plan:* lock whichever combination reproduces the printed output exactly. If
   more than one does, the user chooses.

5. **ICPIN "expanded confidence interval" for fewer than 7 replicates**
   (Phase 6). Defined only in the ICPIN v2.0 program documentation, which is not
   publicly retrievable. *Proposed:* omit it, report only the original interval
   with a note.

6. **Test of Significant Toxicity (TST)**. EPA-833-R-10-003 (2010) supersedes
   parts of this framework in NPDES practice and post-dates ToxCalc v5.0.
   *Proposed:* out of scope for v1.

7. **Non-monotone significance patterns** (Phase 4). EPA 9.6.5.1 warns that if a
   lower concentration is non-significant between two significant ones "the
   results should be used with extreme caution", but prescribes no rule.
   *Proposed:* default `step_down = FALSE` — report the raw pattern, set the
   NOEC to the highest non-significant concentration *below the lowest
   significant one*, and add a warning row to the decision trail.

8. **Multiple controls** (Phase 9). The marine manual uses brine controls and,
   for some methods, both a dilution-water and a solvent control.
   `toxcalc_data(control = )` currently takes one value.

9. **Marine manual scope** (Phase 9). Sea urchin fertilisation is quantal per
   replicate; *Champia parvula* cystocarp counts have a different replicate
   structure.

---

## 9. Phases

| Phase | Content | Definition of done |
|---|---|---|
| **0** | Scaffold: DESCRIPTION, package/params files, air, editorconfig, gitattributes, LICENSE, NEWS, pkgdown, three workflows, build script, Rbuildignore, CLAUDE.md | `R CMD check` 0/0/0 with zero exports. **DONE.** |
| **1** | `toxcalc_data()` + methods; `arcsine_sqrt()` / `inv_arcsine_sqrt()` with the Bartlett (1937) endpoint adjustments; `smooth_monotone()` (PAVA, both directions, optional weights); `is_monotone()`; `abbott()` | **DONE.** Reproduces App B 4.2 (`RP = 0.60 -> 0.8861`; `RP = 0, n = 20 -> 0.1120`; `RP = 1 -> 1.4588`) and App K smoothing (`0.02`). PAVA, not pairwise averaging: the App K case pools five values at once. |
| **2** | `epa_normality()` (pooled centred within-group residuals; Royston, Kolmogorov D above n = 50 with Table B.11); `epa_variance()` (Bartlett, plus Levene and Fligner flagged non-EPA); `msd()`, `pmsd()`; deterministic Dunnett critical value in `R/critical.R`; the three EPA datasets and `epa_pmsd_bounds` | **DONE.** App B `W = 0.9601` vs printed 0.959; Bartlett 6.836 with the printed 7.691 reconciled; Kolmogorov `D* = 0.4572` with the printed 0.4684 reconciled; App C `Sw = 0.0971`, `t = 1.486, 0.248, 1.635, 3.248`, `d = 2.3561` vs tabled 2.36, `MSD = 0.1619`, `PMSD = 23.9%`. 175 assertions. |
| **3** | `dunnett()`, `bonferroni_t()`, `dunn_sidak_t()`, `welch_t()`, `steel()`, `wilcoxon_rank_sum()`, `fisher_exact()`, sharing the `toxcalc_comparison` class and the NOEC/LOEC derivation | **DONE.** App C: `t = 1.486, 0.248, 1.635, 3.248`, critical 2.3561, NOEC 128, LOEC 256. App D: `t = 1.622, 0.270, 1.785, 4.028`, critical 2.510, NOEC 128, LOEC 256. App E: rank sums 84, 63.5, 76, 55, NOEC 3, LOEC 6. App F: rank sums 79, 57, 58, 55, NOEC 3, LOEC 6. App G: NOEC 12, LOEC 25. 254 assertions. |
| **4** | `walk_flowchart()`, chart data, `decisions()`, `toxcalc()` hypothesis branch, print/summary with snapshot tests, exclusion rules (9.5.2), lower-PMSD override (10.2.8.2.5), override warning | End-to-end reproduction of the Section 12 fathead embryo-larval example and the App B/C growth example. **First releasable version.** |
| **5** | `graphical_lc50()`, `spearman_karber()`, `trimmed_spearman_karber()`, `probit_lc()`, and the Fig. 6 chart | Acute Table 20: Graphical 35%; SK `m = 1.656527`, `V(m) = 0.0010977`, LC50 45.3% (38.9, 52.8); TSK trim 20.51%, LC50 77.11 (69.74, 85.26); App K trim 20.41%, LC50 77.28; Probit chi-sq 3.076, LC50 22.872 (18.787, 27.846), LC1 7.924 (4.147, 10.959) |
| **6** | `icp()` with in-loop re-smoothing and the EPA order-statistic interval | App M `IC25 = 8.57%`, `IC50 = 10.89%`; reproducible under a seed; NA-replicate accounting reported |
| **7** | `williams_test()` with simulated critical values and a monotonicity precondition | Reproduces Williams (1971) published examples; documented as a ToxCalc feature and non-EPA extension |
| **8** | Top-level driver wiring both branches, `as.data.frame()`, README, pkgdown reference grouped by flowchart branch | Every worked example in both manuals runs through the single driver |
| **9** | Shirley-Williams, Dunn-Sidak, Welch many-one, Jonckheere-Terpstra monotonicity report, marine manual endpoints, multiple controls | — |

---

## 10. Rejected alternatives

- **Transcribing the EPA appendix tables of critical values.** Rejected on
  third-party copyright grounds (§3), and secondarily because the reproduced
  Dunnett table in the source PDF is visibly OCR-corrupted. Computation is
  cleaner and covers unbalanced designs the tables never did.
- **Separate `noec()` / `loec()` / `msd()` entry points as the primary API.**
  Rejected: all four endpoints come from **one** traversal, and separate entry
  points would each have to re-run the assumption tests and could disagree if
  arguments drifted. One traversal, one object, four accessors.
- **Nested `if`/`else` for the flowcharts.** Rejected: untestable in isolation,
  and a code change would be needed per manual. Charts are data (§5.1).
- **`multcomp` as a runtime Import.** Rejected: used for one quantile that
  `mvtnorm` provides directly, while pulling four further packages (§7).
- **`boot` for the ICp interval.** Rejected: cannot express in-loop re-smoothing
  or the order-statistic rule (§7).
- **`PMCMRplus` for Williams and Steel.** Rejected: not installed, and
  impractical to install on Windows (needs MPFR and GMP system libraries).
  `kSamples` covers Steel better; Williams is implemented here.
