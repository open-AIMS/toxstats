#' @keywords internal
#'
#' @section Intended use:
#'
#' **This package was written with generative AI, and it is for testing and
#' validation purposes only. It must not be used to derive toxicity estimates
#' for regulatory submission, compliance reporting, or any other official
#' purpose.**
#'
#' Every method is checked against the worked examples printed in the EPA
#' method manuals, and those checks are in the test suite rather than asserted
#' in prose. That establishes less than it may appear to. Agreement with a
#' printed example shows that one path through one method reproduces one
#' published number; it is not independent verification of the package, it does
#' not cover the paths no manual exercises, and it is not the review a
#' regulator would expect of software used to produce a submitted result.
#'
#' Two further limits are worth stating plainly. The package is experimental
#' and its interface may change. And it is not, and does not claim to be, a
#' validated replacement for the withdrawn ToxCalc software: no ToxCalc output
#' was available to compare against, so nothing here has been shown to match
#' what that program produced.
#'
#' Use it to check an analysis done another way, to explore how the EPA
#' decision rules behave, or to develop and test other software. For a result
#' that will be relied on, obtain it from software your regulator accepts.
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats ave sd var
## usethis namespace: end

## mockable bindings: start
## mockable bindings: end
NULL
