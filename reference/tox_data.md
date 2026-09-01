# Prepare toxicity test data for analysis

Validates replicate-level concentration-response data and puts it into
the single form every other function in the package accepts. One row per
replicate test chamber covers both kinds of endpoint the EPA WET manuals
recognise.

## Usage

``` r
tox_data(
  data,
  conc = "conc",
  response = "response",
  n_exposed = NULL,
  replicate = NULL,
  type = c("continuous", "quantal"),
  direction = c("decreasing", "increasing"),
  control = 0,
  conc_units = "%",
  response_units = NULL
)

# S3 method for class 'tox_data'
print(x, ...)

# S3 method for class 'tox_data'
summary(object, ...)

# S3 method for class 'tox_data'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- data:

  A data frame with one row per replicate test chamber.

- conc:

  Name of the column holding the exposure concentration. A string.

- response:

  Name of the column holding the response. For `type = "continuous"`
  this is the measured value; for `type = "quantal"` it is the number of
  organisms affected. A string.

- n_exposed:

  Name of the column holding the number of organisms exposed in each
  replicate. Required when `type = "quantal"`, otherwise `NULL`.

- replicate:

  Name of an optional column labelling replicates within a
  concentration, or `NULL`.

- type:

  The response type. Either `"continuous"` (growth, reproduction,
  biomass) or `"quantal"` (survival, fertilisation, hatching).

- direction:

  The expected direction of the concentration-response relationship for
  the measured response. Either `"decreasing"` (the usual case, where
  toxicity reduces the response) or `"increasing"`.

- control:

  The concentration treated as the control. Defaults to `0`.

- conc_units:

  Units of `conc`, used only for labelling output. A string.

- response_units:

  Units of `response`, used only for labelling output. A string, or
  `NULL`.

- x:

  A `tox_data` object.

- ...:

  Unused, present for consistency with the generic.

- object:

  A `tox_data` object.

- row.names:

  Unused, present for consistency with the generic.

- optional:

  Unused, present for consistency with the generic.

## Value

An object of class `tox_data`, a list with elements:

- `replicates`:

  a data frame with one row per replicate, holding `conc`, `replicate`,
  `response`, and for quantal data `n_exposed` and `proportion`.

- `pooled`:

  a data frame with one row per concentration, holding `n_rep`, `mean`,
  `sd` and `var` of the analysis response, and for quantal data the
  totals `n_exposed`, `n_affected` and the pooled `proportion`.

- `type`, `direction`, `control`, `conc_units`, `response_units`:

  as supplied.

- `n_dropped`:

  the number of rows removed for a missing response.

## Details

The manuals distinguish two response types.

A **continuous** endpoint is a measurement, such as growth in milligrams
or the number of young per female. `response` holds the measured value
and `n_exposed` is not used.

A **quantal** endpoint is one where each organism either responds or
does not, such as survival, hatching or fertilisation. `response` holds
the number of organisms **affected** in the replicate and `n_exposed`
the number exposed. Both are required, because the manuals analyse
quantal data twice: the hypothesis-testing flowchart works on
replicate-level proportions after an arc sine square root
transformation, while the point-estimation flowchart works on counts
pooled across replicates within a concentration. Recording the number
exposed is what makes both possible from one input.

Rows with a missing response are dropped, and the number dropped is
recorded in the returned object. This matters for the linear
interpolation method, whose bootstrap resamples within a concentration
and so depends on how many replicates each concentration actually has.

## Methods (by generic)

- `print(tox_data)`: Print a compact description of the design.

- `summary(tox_data)`: Return the per-concentration summary.

- `as.data.frame(tox_data)`: Return the validated replicate-level data.

## References

US EPA (2002) EPA-821-R-02-013, sections 9.4 to 9.6.

## Examples

``` r
growth <- data.frame(
  conc = rep(c(0, 32, 64, 128, 256), each = 4),
  growth = c(
    0.711, 0.662, 0.054, 0.785,
    0.517, 0.501, 0.723, 0.760,
    0.602, 0.669, 0.694, 0.706,
    0.348, 0.400, 0.041, 0.512,
    0.216, 0.277, 0.328, 0.347
  )
)
tox_data(growth, response = "growth")
#> <tox_data>
#>   continuous response, expected to be decreasing with concentration
#>   5 concentrations, 20 replicates (4 per concentration, balanced)
#>   control at 0 %
#> 
#>   conc n_rep    mean         sd         var
#> 1    0     4 0.55300 0.33648675 0.113223333
#> 2   32     4 0.62525 0.13523899 0.018289583
#> 3   64     4 0.66775 0.04646414 0.002158917
#> 4  128     4 0.32525 0.20147684 0.040592917
#> 5  256     4 0.29200 0.05865720 0.003440667
```
