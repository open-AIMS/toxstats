# toxcalc — what we are building and why

This document is for collaborators. It explains what the package does, why, in
what order it will be built, and what changes for someone analysing toxicity
data. The full specification, including algorithms and the evidence behind each
decision, is in `notes/TOXCALC-claude.md`; section references below point there.

---

## The problem

ToxCalc v5.0 was commercial software that carried out the statistical analysis
required by the US EPA Whole Effluent Toxicity (WET) method manuals. It has been
withdrawn. The vendor product page now returns a 404 and redirects to a
different product.

Laboratories that used it have no open replacement. An audit of 642 R packages
installed on this machine found no implementation of Williams' test, no
implementation of the trimmed Spearman-Karber method, and no implementation of
the EPA linear-interpolation ICp method. All three are required by the EPA
protocols.

The pieces that *do* exist are scattered. Normality and variance tests are in
base R, Steel's Many-One Rank Test is in `kSamples`, and Dunnett's test can be
assembled from `mvtnorm`. An analyst can call each by hand, but then the choice
of which test to run — which is the part the EPA manuals actually specify — is
undocumented and unreproducible.

## What we are building

An R package that walks the EPA decision charts and reports both the answer and
the reasoning.

The package takes replicate-level concentration-response data and returns the
EPA endpoints. From the hypothesis-testing chart it returns the NOEC (the
highest concentration at which no significant effect was detected), the LOEC
(the lowest at which one was), the MSD (the minimum significant difference, the
smallest difference from the control the test could have detected) and the PMSD
(that difference as a percentage of the control mean). From the point-estimation
chart it returns an LC50 or an ICp.

**The distinguishing feature is the audit trail.** Every analysis carries a
table recording each branch point: the question asked, the statistic, the
answer, the consequence, and the section of the EPA manual that justifies it.
Printed output shows the reasoning above the result, so a reviewer can check the
analysis line by line against the manual rather than taking the number on trust.

### What it will not do

It is an R package, with no graphical interface and no report templates. It
targets the EPA manuals, not ToxCalc's exact numbers; no ToxCalc output was
available to compare against. It does not implement the 2010 Test of Significant
Toxicity, which post-dates ToxCalc (specification §8.6).

## Three findings that changed the plan

Reading the manuals directly, rather than working from recollection, corrected
three things worth flagging.

**Williams' test is not an EPA method.** The manuals name it only as an
alternative "requiring additional assumptions". It is in the package because
ToxCalc offered it, but it is a documented extension rather than part of the EPA
chart. It therefore moves from an early phase to a late one, and will be
labelled as an extension in its documentation (specification §2.1, §6.3).

**No third-party statistical tables will be transcribed.** This reverses an
earlier decision. The EPA manuals are US Government works, but several tables
*inside* them are reproduced from other sources under a permission that does not
transfer to us: Williams' critical values from *Biometrics*, Dunnett's from
Miller (1981), and Conover's Shapiro-Wilk coefficients. Critical values are
computed instead. This is also more accurate: computation handles unbalanced
designs that the published tables never covered, and the Dunnett table
reproduced in the EPA PDF is visibly corrupted by optical character recognition
(specification §3).

Dunnett's critical value turned out to need care. The obvious tool,
`mvtnorm::qmvt()`, is a randomised method and returns a slightly different
answer every time it is called — between 2.3552 and 2.3585 on the example we
validate against. An analysis that goes to a regulator has to give the same
number every time, so the value is obtained by ordinary numerical integration
instead, which is exact and repeatable. `mvtnorm` is consequently not a
dependency.

**The point-estimation scope is larger than expected.** The EPA chart for the
LC50 has four terminal methods, not two. As well as probit analysis and the
trimmed Spearman-Karber method, it requires a graphical method and a plain
untrimmed Spearman-Karber method (specification §2.3).

**A note on where decisions are recorded.** Decisions taken during
implementation, and the evidence for each, are written up in the package
vignette `vignettes/recreating-toxcalc.qmd`, which ships with the package and
is the document to read when checking this work against the manuals.

## How the data goes in

One data frame, one row per replicate test chamber, for both kinds of endpoint.

For a continuous endpoint such as growth or reproduction, the columns are the
concentration and the measured value. For a quantal endpoint such as survival —
that is, one where each organism either responds or does not — the columns are
the concentration, the number of organisms affected, and the number exposed.

The reason both kinds share one shape is that EPA analyses quantal data twice.
Replicate-level proportions go through the hypothesis chart, and counts pooled
across replicates go into the point-estimation methods. Recording the number
exposed makes both possible from one input.

A typical call looks like this:

```r
result <- toxcalc(growth ~ conc, data = fathead, type = "continuous")
summary(result)
```

## Order of work

Nine phases. Each is separately testable, separately committable, and useful on
its own. Full definitions of done are in specification §9.

**Phase 0 — scaffold. Complete.** *Purpose:* a package that builds and checks cleanly
before any statistics are added, so later phases never debug the build and the
analysis at the same time. *Done:* `R CMD check` reports zero errors, warnings
and notes. **Complete.**

**Phase 1 — data handling and primitives. Complete.** *Purpose:* the input contract, plus
the two transformations everything else depends on — the EPA arc sine square
root transformation and monotone smoothing. *Why it is needed:* the EPA arcsine
has a specific adjustment when an observed proportion is exactly 0 or exactly 1,
and the smoothing step is shared by three later methods. Neither exists in any
CRAN package, so this phase is useful even alone. *Done:* reproduces four
worked transformations printed in the manuals.

**Phase 2 — assumption tests and sensitivity. Complete.** *Purpose:* the normality and
equal-variance tests that decide which branch of the chart is taken, and the MSD
and PMSD that describe how sensitive the test was. *Done:* reproduces the
manuals' printed statistics for both, and the Appendix C minimum significant
difference.

**Phase 3 — the individual hypothesis tests. Complete.** *Purpose:* the seven tests the
chart can select, each callable and testable on its own. *Done:* each reproduces
its worked example in the manuals.

**Phase 4 — the chart engine and the audit trail. Complete.** *Purpose:* the actual
recreation of what ToxCalc did — selecting a test automatically and recording
why. *Done:* two complete worked examples run end to end and select the same
test the manual selects. **This is the first version worth releasing.**

**Phase 5 — point estimation for survival data. Complete.** *Purpose:* the LC50 branch,
with all four EPA methods and the chart that chooses between them. *Done:*
reproduces five worked examples, including the trimmed Spearman-Karber trims of
20.51% and 20.41% printed in the two manuals.

**Phase 6 — inhibition concentrations. Complete.** *Purpose:* the ICp method for sublethal
endpoints, with bootstrap confidence limits. *Done:* reproduces the Appendix M
IC25 of 8.57% and IC50 of 10.89%.

**Phase 7 — Williams' test. Complete.** *Purpose:* the ToxCalc feature that is not an EPA
method, kept separate so its status is unambiguous. *Done:* reproduces Williams'
own published examples.

**Phase 8 — the top-level driver and documentation. Complete.** *Purpose:* one function
that runs either branch, plus the README and website. *Done:* every worked
example in both manuals runs through the single driver.

**Phase 9 — extensions. Not started, and not needed.** Phases 0 to 8 leave a
package that does the job. Everything below is a self-contained addition, and
none depends on another. Specification §11 sets out each one with its scope,
its effort and how you would know it was finished. In the order I would do
them:

1. **The acute manual's single-concentration chart.** A second
   hypothesis-testing chart, for a test at one concentration against a control,
   returning Pass or Fail rather than a NOEC. Three of its four terminals
   already exist, so this is the smallest piece of real work remaining.
2. **The marine and estuarine manual.** Mostly reading. Its charts are believed
   to match the freshwater ones; its endpoints do not all have the same
   replicate structure. That manual has not been downloaded.
3. **Designs with more than one control.** The marine manual uses brine and
   solvent controls. This is not an implementation detail: what a second
   control *means* for the comparison has to be read from the manual, not
   chosen.
4. **A formal monotonicity test** to accompany the warning Williams' test
   already gives. Small; the dependency is already present.
5. **Shirley-Williams.** Neither ToxCalc nor the manuals include it, so there
   is no reason from the brief to add it at all. Recorded only because the
   original design listed it.

The 2010 Test of Significant Toxicity is **deliberately out of scope**: it
post-dates ToxCalc, so adding it would be a different project rather than a
continuation of this one.

## What had to be decided, and what still does

Nine points in the manuals were ambiguous enough that guessing would have been
wrong. Six have since been settled by testing candidate answers against the
manuals' own printed output. Each is written up with its evidence in the
package vignette, `vignettes/recreating-toxcalc.qmd`, which is the document to
read when checking this work.

**Settled.** R's Shapiro-Wilk is used, because it reproduces the manual's
printed statistic. Dunnett's critical value is computed by integration rather
than taken from `mvtnorm`, because `mvtnorm` is randomised and gives a slightly
different answer every call. The probit heterogeneity statistic is Pearson's,
not the deviance, and its limits are Fieller's, not the delta method — in both
cases because that is what reproduces the printed output and the alternative
does not. The Spearman-Karber interval uses the manual's own multiplier of 2.0.
Partial responses are counted after smoothing and adjustment, because counting
the raw ones sends the manual's own example to the wrong method.

**Still open.** Three, none of which blocks anything.

The **expanded confidence interval** the linear-interpolation method reports
when there are fewer than seven replicates is defined only in program
documentation that is not publicly retrievable. It is omitted, and a design
that would have triggered it warns.

A **non-monotone pattern of significance** — a middle concentration not
significant between two that are — draws advice to use caution from the manual
but no rule. The pattern is reported as observed, the NOEC set below the lowest
significant concentration, and the result flagged.

The **heterogeneity correction to the probit variance** is implemented as the
manual describes, but the manual's own example does not exercise it, so it is
unverified against any published output.

One further quantity, the **trimmed Spearman-Karber confidence interval**,
could not be settled: the manual delegates it to a program whose formula it
never states. The delta method is used and the departure is documented. The
trim and the point estimate both reproduce exactly.

## Dependencies

Three packages at runtime: `chk` for argument checking, `kSamples` for Steel's
test, and base `stats`.

Five more are used only in the test suite, as independent checks on methods
implemented here from first principles: `MASS`, `multcomp`, `mvtnorm`,
`nortest` and `boot`. Keeping them out of the runtime list means a user
installs three packages, not eight, and it forces each method to be implemented
rather than delegated. The reasoning for each is in specification §7.

`mvtnorm` was originally going to supply Dunnett's critical value and is in the
suggested list rather than the required one for the reason given above: it is
randomised, and the value is computed by integration instead.

`PMCMRplus`, which would otherwise be the obvious source for Williams' test, is
not used. It is not installed here and needs external MPFR and GMP libraries,
which makes it impractical to install on Windows.

## A caution about the manuals

The manuals contain arithmetic errors. Appendix C computes a minimum significant
difference of 0.162 and then, a paragraph later, describes it as 0.087 — but 24%
of the stated control mean of 0.677 is 0.162, so the second figure is the typo.

Every published value used as a test fixture will therefore be checked for
internal consistency before it is frozen into the test suite, rather than
transcribed on trust. This is budgeted into Phases 2 and 5.
