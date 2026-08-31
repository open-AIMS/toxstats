# Builds the EPA worked-example datasets shipped in data/.
#
# Every value here is transcribed from the printed tables of
# EPA-821-R-02-013, "Short-term Methods for Estimating the Chronic Toxicity of
# Effluents and Receiving Waters to Freshwater Organisms", 4th edition,
# October 2002. Each dataset is checked against the summary statistics the
# manual prints alongside it; the checks at the foot of this file fail loudly
# if a transcription error is introduced.
#
# Note that the manual's Appendix C states its data are "the same data used in
# Appendices B and D". They are not: Table C.1 differs from Table B.1. Both are
# shipped, separately, for that reason.

# Table B.1 -- fathead minnow larval growth, used for the Shapiro-Wilk and
# Bartlett worked examples in Appendix B.
fathead_b1 <- data.frame(
  conc = rep(c(0, 32, 64, 128, 256), each = 4),
  replicate = rep(c("A", "B", "C", "D"), times = 5),
  weight = c(
    0.711,
    0.662,
    0.718,
    0.767,
    0.646,
    0.626,
    0.723,
    0.700,
    0.669,
    0.669,
    0.694,
    0.676,
    0.629,
    0.680,
    0.513,
    0.672,
    0.650,
    0.558,
    0.606,
    0.508
  )
)

# Table C.1 -- fathead minnow larval growth, used for the Dunnett worked
# example in Appendix C.
fathead_c1 <- data.frame(
  conc = rep(c(0, 32, 64, 128, 256), each = 4),
  replicate = rep(c("A", "B", "C", "D"), times = 5),
  weight = c(
    0.711,
    0.662,
    0.646,
    0.690,
    0.517,
    0.501,
    0.723,
    0.560,
    0.602,
    0.669,
    0.694,
    0.676,
    0.566,
    0.612,
    0.410,
    0.672,
    0.455,
    0.502,
    0.606,
    0.254
  )
)

# Table B.7 -- Ceriodaphnia dubia reproduction, used for the Kolmogorov "D"
# worked example in Appendix B. Sixty observations, so Shapiro-Wilk is not the
# method the manual applies here.
ceriodaphnia_b7 <- data.frame(
  conc = rep(c(0, 1.56, 3.12, 6.25, 12.5, 25.0), each = 10),
  replicate = rep(as.character(1:10), times = 6),
  young = c(
    27,
    30,
    29,
    31,
    16,
    15,
    18,
    17,
    14,
    27,
    32,
    35,
    32,
    26,
    18,
    29,
    27,
    16,
    35,
    13,
    39,
    30,
    33,
    33,
    36,
    33,
    33,
    27,
    38,
    44,
    27,
    34,
    36,
    34,
    31,
    27,
    33,
    31,
    33,
    31,
    19,
    25,
    26,
    17,
    16,
    21,
    23,
    15,
    18,
    10,
    10,
    13,
    7,
    7,
    7,
    10,
    10,
    16,
    12,
    2
  )
)

# Transcription checks against the summary statistics the manual prints ------

# Compare an unrounded statistic with the value the manual prints, allowing
# half a unit in the last printed place. Comparing round(x, digits) instead
# would fail on an exact half such as the Table B.1 control mean of 0.7145,
# which the manual prints as 0.714 while R's round() returns 0.715.
check_printed <- function(x, value, group, fun, expected, digits, what) {
  observed <- as.numeric(tapply(x[[value]], x[[group]], fun))
  tolerance <- 0.5 * 10^(-digits) + 1e-12
  bad <- abs(observed - expected) > tolerance
  if (any(bad)) {
    stop(
      what,
      " does not match the printed table.\n",
      "  transcribed: ",
      paste(signif(observed, 6), collapse = " "),
      "\n",
      "  manual:      ",
      paste(expected, collapse = " "),
      "\n",
      "  differing:   ",
      paste(which(bad), collapse = " "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Table B.1 prints these means.
check_printed(
  fathead_b1,
  "weight",
  "conc",
  mean,
  c(0.714, 0.674, 0.677, 0.624, 0.580),
  3,
  "Table B.1 means"
)

# Table B.1 also prints the within-concentration variances as
# 0.0018 0.0020 0.0001 0.0059 0.0037. The second of these is a rounding error
# in the manual: the variance at 32 ug/L is 0.002055, which rounds to 0.0021.
# The error matters, because it is what produces the manual's printed Bartlett
# statistic of 7.691 in Appendix B section 3.5 -- recomputing from the raw data
# gives 6.836. The check below therefore uses the manual's own printed values,
# with that one entry corrected, so the transcription is still verified.
check_printed(
  fathead_b1,
  "weight",
  "conc",
  var,
  c(0.0018, 0.0021, 0.0001, 0.0059, 0.0037),
  4,
  "Table B.1 variances"
)

# Table C.1 prints these means and totals.
check_printed(
  fathead_c1,
  "weight",
  "conc",
  mean,
  c(0.677, 0.575, 0.660, 0.565, 0.454),
  3,
  "Table C.1 means"
)
stopifnot(isTRUE(all.equal(
  round(as.numeric(tapply(fathead_c1$weight, fathead_c1$conc, sum)), 3),
  c(2.709, 2.301, 2.641, 2.260, 1.817)
)))

# Table B.7 prints these means.
check_printed(
  ceriodaphnia_b7,
  "young",
  "conc",
  mean,
  c(22.4, 26.3, 34.6, 31.7, 19.0, 9.4),
  1,
  "Table B.7 means"
)

usethis::use_data(fathead_b1, overwrite = TRUE)
usethis::use_data(fathead_c1, overwrite = TRUE)
usethis::use_data(ceriodaphnia_b7, overwrite = TRUE)

# Table 6 -- EPA variability criteria (PMSD bounds) for sublethal hypothesis
# testing endpoints submitted under NPDES permits. Derived by EPA from the
# 10th and 90th percentiles of its own WET Interlaboratory Variability Study,
# so this is an EPA-generated table rather than a third-party reproduction.
epa_pmsd_bounds <- data.frame(
  id = c("fathead_growth", "ceriodaphnia_reproduction", "selenastrum_growth"),
  method = c("1000.0", "1002.0", "1003.0"),
  test = c(
    "Fathead Minnow Larval Survival and Growth Test",
    "Ceriodaphnia dubia Survival and Reproduction Test",
    "Selenastrum capricornutum Growth Test"
  ),
  endpoint = c("growth", "reproduction", "growth"),
  lower = c(12, 13, 9.1),
  upper = c(30, 47, 29),
  stringsAsFactors = FALSE
)

usethis::use_data(epa_pmsd_bounds, overwrite = TRUE)
