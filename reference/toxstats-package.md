# toxstats: Standard Endpoint Statistics for Toxicity Test Data

Computes the standard endpoints of a concentration-response toxicity
test: the no- and lowest-observed-effect concentrations, the minimum
significant difference and its percentage of the control, median lethal
concentrations, and inhibition concentrations. The method used is
selected by the decision rules set out in the US EPA Whole Effluent
Toxicity method manuals, and every analysis carries a record of which
test was selected and which section of the manual justifies the choice.
Those rules are written for a regulatory setting but apply to any
concentration-response toxicity data with a control. Recreates the
capability of the withdrawn ToxCalc software (Tidepool Scientific), with
which this package is not affiliated. Written with generative AI, and
intended for testing and validation only: it must not be used to derive
toxicity estimates for regulatory submission, compliance reporting, or
any other official purpose.

## Intended use

**This package was written with generative AI, and it is for testing and
validation purposes only. It must not be used to derive toxicity
estimates for regulatory submission, compliance reporting, or any other
official purpose.**

Every method is checked against the worked examples printed in the EPA
method manuals, and those checks are in the test suite rather than
asserted in prose. That establishes less than it may appear to.
Agreement with a printed example shows that one path through one method
reproduces one published number; it is not independent verification of
the package, it does not cover the paths no manual exercises, and it is
not the review a regulator would expect of software used to produce a
submitted result.

Two further limits are worth stating plainly. The package is
experimental and its interface may change. And it is not, and does not
claim to be, a validated replacement for the withdrawn ToxCalc software:
no ToxCalc output was available to compare against, so nothing here has
been shown to match what that program produced.

Use it to check an analysis done another way, to explore how the EPA
decision rules behave, or to develop and test other software. For a
result that will be relied on, obtain it from software your regulator
accepts.

## See also

Useful links:

- <https://open-aims.github.io/toxstats>

- <https://github.com/open-AIMS/toxstats>

- Report bugs at <https://github.com/open-AIMS/toxstats/issues>

## Author

**Maintainer**: Rebecca Fisher <r.fisher@aims.gov.au>

Authors:

- Rebecca Fisher <r.fisher@aims.gov.au>
