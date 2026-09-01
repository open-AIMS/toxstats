# Changelog

## toxstats 0.0.0.9000

- [`trimmed_spearman_karber()`](https://open-aims.github.io/toxstats/reference/trimmed_spearman_karber.md)
  now reproduces the confidence limits printed in EPA-821-R-02-012 Table
  20 exactly, where it previously differed in the third significant
  figure. The variance is Hamilton’s, derived analytically and holding
  the trim fixed, and the multiplier is the manual’s 2.0 rather than
  1.96. This was the last quantity in the hypothesis-testing and
  point-estimation branches known to depart from a published EPA value
  for a reason other than an error in the manual.

- [`trimmed_spearman_karber()`](https://open-aims.github.io/toxstats/reference/trimmed_spearman_karber.md)
  takes the first concentration reaching `1 - a` as the upper endpoint
  of the trimmed range, rather than the last concentration at or below
  it. The two differ only where the response plateaus exactly on the
  boundary, in which case the estimate was previously too high.

- `ecotoxicology` added to `Suggests`, as a test oracle for the trimmed
  Spearman-Karber method.

- Initial package scaffold.
