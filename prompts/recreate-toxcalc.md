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
