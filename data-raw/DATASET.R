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

# Table E.1 -- Ceriodaphnia dubia reproduction, used for the Steel's Many-One
# Rank Test worked example in Appendix E, and (with two values removed) for
# the Wilcoxon Rank Sum worked example in Appendix F.
#
# The 50 per cent concentration is retained here but is excluded from the
# reproduction analysis in the manual, because survival at that concentration
# was significantly reduced (section 9.5.2, and the Appendix G example below).
ceriodaphnia_e1 <- data.frame(
  conc = rep(c(0, 3, 6, 12, 25, 50), each = 10),
  replicate = rep(as.character(1:10), times = 6),
  young = c(
    20,
    26,
    26,
    23,
    24,
    27,
    26,
    23,
    27,
    24,
    13,
    15,
    14,
    13,
    23,
    26,
    0,
    25,
    26,
    27,
    18,
    22,
    13,
    13,
    23,
    22,
    20,
    22,
    23,
    22,
    14,
    22,
    20,
    23,
    20,
    23,
    25,
    24,
    25,
    21,
    9,
    0,
    9,
    7,
    6,
    10,
    12,
    14,
    9,
    13,
    rep(0, 10)
  )
)

# Table G.2 -- Ceriodaphnia dubia mortality, used for the Fisher's Exact Test
# worked example in Appendix G. `exposed` is the number of live adults at the
# start of the test, which is nine in the control.
ceriodaphnia_g2 <- data.frame(
  conc = c(0, 1, 3, 6, 12, 25),
  dead = c(1, 0, 0, 0, 0, 10),
  exposed = c(9, 10, 10, 10, 10, 10)
)

# Table E.4 prints rank sums of 84, 64, 76 and 55. The 6 per cent value is a
# slip: the ranks the manual itself lists in Table E.3 for that concentration
# are 3, 7.5, 1.5, 1.5, 11.5, 7.5, 4.5, 7.5, 11.5 and 7.5, which sum to 63.5.
# The conclusion is unaffected, both being at or below the critical 76.
rank_sum_check <- function(control, treatment) {
  combined <- rank(c(control, treatment))
  sum(combined[-seq_along(control)])
}
repro <- split(
  ceriodaphnia_e1$young[ceriodaphnia_e1$conc < 50],
  ceriodaphnia_e1$conc[ceriodaphnia_e1$conc < 50]
)
stopifnot(isTRUE(all.equal(
  vapply(repro[-1], function(v) rank_sum_check(repro[[1]], v), numeric(1)),
  c(`3` = 84, `6` = 63.5, `12` = 76, `25` = 55)
)))

# Table G.2 totals.
stopifnot(
  sum(ceriodaphnia_g2$dead) == 11,
  sum(ceriodaphnia_g2$exposed) == 59
)

usethis::use_data(ceriodaphnia_e1, overwrite = TRUE)
usethis::use_data(ceriodaphnia_g2, overwrite = TRUE)
