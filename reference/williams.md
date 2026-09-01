# Williams' test

Compares each concentration with the control while assuming the
underlying concentration-response relationship is monotone. Using that
assumption makes the test more powerful than Dunnett's procedure when it
holds, and misleading when it does not.

## Usage

``` r
williams(x, ..., alpha = 0.05, nsim = 20000, seed = 1L)
```

## Arguments

- x:

  A
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  object, or a data frame that can be coerced to one by passing it
  through
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md).

- ...:

  Additional arguments passed to
  [`tox_data()`](https://open-aims.github.io/toxstats/reference/tox_data.md)
  when `x` is a bare data frame, and otherwise unused.

- alpha:

  Significance level for the hypothesis tests themselves. The EPA
  manuals specify 0.05, one-sided.

- nsim:

  Number of simulations used for each critical value.

- seed:

  Integer seed for the simulation. The random number stream is restored
  afterwards, so calling this function does not disturb any other
  analysis.

## Value

An object of class `tox_comparison`, a list with elements `test`,
`method`, `comparisons` (a data frame with one row per treatment
concentration), `noec`, `loec`, `monotone`, `critical`, `sw`, `df`,
`alpha` and `reference`.

## Details

### This is not an EPA method

Williams' test is **not** on the EPA flowchart. Section 9.4.1.2 of the
chronic manual mentions it only as an alternative that "requires
additional assumptions", and neither manual gives a worked example or a
table of critical values for it. It is provided here because ToxCalc
offered it, and every result is labelled as an extension so it cannot be
mistaken for the analysis the manual prescribes.

### The method

The concentration means are replaced by their isotonic estimates under
the restriction that the response is monotone in concentration. **The
control is not part of that restriction**; it is estimated freely. The
statistic for concentration `i` is then

\$\$\bar{t}\_i = \frac{\bar{Y}\_0 - M_i}{s\sqrt{1/n_i + 1/n_0}}\$\$

with `s` the square root of the within mean square from the full
analysis of variance, and `M_i` the isotonic estimate at concentration
`i` formed from concentrations 1 to `i` only.

The test proceeds **downwards** from the highest concentration. Once a
concentration is not significant, every lower one is declared not
significant without being tested. That is what makes the test a trend
procedure rather than a set of pairwise comparisons, and it is recorded
in the `tested` column of the result.

### Critical values are simulated

Williams tabulated his critical values, and those tables are reproduced
from *Biometrics*, so they are not transcribed here (see the package
vignette). They are obtained instead by simulating the null distribution
of \\\bar{t}\_i\\ for the design in hand: the concentration means and
the pooled variance are drawn under the hypothesis of no effect, the
isotonic estimate is formed, and the upper `alpha` quantile is taken.

Simulation also removes the two things that go wrong with the tables.
The critical value depends on `i`, the position of the concentration
being tested, not on the total number of concentrations, which is the
commonest implementation error. And the tables are tabulated only at
particular degrees of freedom and require interpolation in `1/nu`, which
simulation avoids entirely.

`seed` defaults to a fixed value so that the same design always returns
the same critical value. The Monte Carlo standard error is reported
alongside it; raise `nsim` to reduce it.

### Validation

With no EPA worked example and no retrievable table, this implementation
is checked against mathematical identities rather than published output.
With a single concentration there is no order restriction, so the
critical value must equal `qt(1 - alpha, nu)`, and the simulation
reproduces it. Beyond one concentration the critical value must fall
below Dunnett's for the same design, because the order restriction is
additional information, and it must increase with `i`. All three are
asserted in the test suite.

## References

Williams DA (1971) A test for differences between treatment means when
several dose levels are compared with a zero dose control. *Biometrics*
27:103-117.

Williams DA (1972) The comparison of several dose levels with a zero
dose control. *Biometrics* 28:519-531.

## Examples

``` r
# Not an EPA method; compare with dunnett() on the same data.
williams(fathead_c1, response = "weight", nsim = 2000)
#> Warning: Williams' test assumes a monotone concentration-response relationship. The observed means depart from it at position(s) 2. The isotonic step will absorb the departure, but the assumption should be considered before the result is used.
#> Williams' test (not an EPA method; 2,000 simulations)
#>   Williams (1971) Biometrics 27:103-117; not on the EPA-821-R-02-013 Figure 2 flowchart
#> 
#>  conc n   mean isotonic statistic critical   mc_se p_value tested significant
#>    32 4 0.5752   0.5752    1.4860    1.896 0.05930  0.0945  FALSE       FALSE
#>    64 4 0.6602   0.6178    0.8666    1.926 0.05340  0.2600  FALSE       FALSE
#>   128 4 0.5650   0.5650    1.6350    1.956 0.05227  0.0815   TRUE       FALSE
#>   256 4 0.4542   0.4542    3.2480    1.917 0.05084  0.0045   TRUE        TRUE
#> 
#>   NOEC 128
#>   LOEC 256
```
