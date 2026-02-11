# Simple feature summaries

Simple function to calculate feature summaries. These include counts and
percentages of samples in which a chromatographic peak is present for
each feature and counts and percentages of samples in which more than
one chromatographic peak was annotated to the feature. Also relative
standard deviations (RSD) are calculated for the integrated peak areas
per feature across samples. For `perSampleCounts = TRUE` also the
individual chromatographic peak counts per sample are returned.

## Usage

``` r
featureSummary(
  x,
  group,
  perSampleCounts = FALSE,
  method = "maxint",
  skipFilled = TRUE
)
```

## Arguments

- x:

  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with correspondence results.

- group:

  `numeric`, `logical`, `character` or `factor` with the same length
  than `x` has samples to aggregate counts by the groups defined in
  `group`.

- perSampleCounts:

  `logical(1)` whether feature wise individual peak counts per sample
  should be returned too.

- method:

  `character` passed to the
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  function. See respective help page for more information.

- skipFilled:

  `logical(1)` whether filled-in peaks should be excluded (default) or
  included in the summary calculation.

## Value

`matrix` with one row per feature and columns:

- `"count"`: the total number of samples in which a peak was found.

- `"perc"`: the percentage of samples in which a peak was found.

- `"multi_count"`: the total number of samples in which more than one
  peak was assigned to the feature.

- `"multi_perc"`: the percentage of those samples in which a peak was
  found, that have also multiple peaks annotated to the feature.
  Example: for a feature, at least one peak was detected in 50 samples.
  In 5 of them 2 peaks were assigned to the feature. `"multi_perc"` is
  in this case 10%.

- `"rsd"`: relative standard deviation (coefficient of variation) of the
  integrated peak area of the feature's peaks.

- The same 4 columns are repeated for each unique element (level) in
  `group` if `group` was provided.

If `perSampleCounts = TRUE` also one column for each sample is returned
with the peak counts per sample.

## Author

Johannes Rainer
