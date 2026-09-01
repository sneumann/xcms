# Feature to chromatographic peak mapping

During the correspondence step in the preprocessing, chromatographic
peaks get assigned (grouped) to features. The abundances of these
resulting LC-MS features are supposed to represent signal from the same
ion across all analyzed samples. Depending on the correspondence
analysis method used, multiple chromatographic peaks (also eventually
from the **same** sample) are assigned to a feature. This mapping
between features and chromatographic peaks is (for
[XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
and
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object) stored in the `"peakidx"` column of the
[`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
data frame. Alternatively, the mapping can be extracted from an *xcms*
result object using the functions:

- `featureChromPeaks()`: returns a two-column `data.frame` with the IDs
  of the features and the IDs of the associated chromatographic peaks.
  Each row in this `data.frame` represents the mapping of one
  chromatographic peak with one feature. The order of the features in
  the `data.frame` matches the order of the features in
  [`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md).

- `featurePeakidx()`: returns a named `list` of `integer` indices of the
  rows in the
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  matrix that are assigned to a feature. The names of the `list` are the
  feature IDs. The length and order of the `list` matches the number of
  rows and order of features in
  [`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md).

## Usage

``` r
featureChromPeaks(object, ...)

featurePeakidx(object, ...)

# S4 method for class 'XcmsResult'
featurePeakidx(object, msLevel = integer())

# S4 method for class 'XcmsResult'
featureChromPeaks(object, msLevel = integer())

# S4 method for class 'XcmsExperimentHdf5'
featureChromPeaks(object, msLevel = integer())

# S4 method for class 'XcmsExperimentHdf5'
featurePeakidx(object, msLevel = integer())
```

## Arguments

- object:

  An *xcms* result object with correspondence analysis results being
  present.

- ...:

  Optional parameters. Currently ignored.

- msLevel:

  Optional `integer` to restrict to features from a certain MS level.

## Value

See description above.

## Examples

``` r

## Load preprocessing results
library(MsExperiment)
xmse <- loadXcmsData()

## Get the mapping between features and chromatographic peaks
map <- featureChromPeaks(xmse)

head(map)
#>   feature_id chrom_peak_id
#> 1      FT001        CP0511
#> 2      FT001        CP1261
#> 3      FT001        CP2957
#> 4      FT001        CP3129
#> 5      FT001        CP3447
#> 6      FT001        CP3536

## Column `"feature_id"` contains the IDs for the features defined in
## `featureDefinitions()`
featureDefinitions(xmse) |> head()
#>       mzmed mzmin mzmax    rtmed    rtmin    rtmax npeaks KO WT      peakidx
#> FT001 200.1 200.1 200.1 2902.634 2882.603 2922.664      2  2  0 458, 116....
#> FT002 205.0 205.0 205.0 2789.901 2782.955 2796.531      8  4  4 44, 443,....
#> FT003 206.0 206.0 206.0 2789.405 2781.389 2794.219      7  3  4 29, 430,....
#> FT004 207.1 207.1 207.1 2718.560 2714.047 2727.347      7  4  3 16, 420,....
#> FT005 233.0 233.0 233.1 3023.579 3015.145 3043.959      7  3  4 69, 959,....
#> FT006 241.1 241.1 241.2 3683.299 3661.586 3695.886      8  3  4 276, 284....
#>       ms_level
#> FT001        1
#> FT002        1
#> FT003        1
#> FT004        1
#> FT005        1
#> FT006        1

## Column `"chrom_peak_id"` contains the IDs of the chromatographic peaks
chromPeaks(xmse) |> head()
#>           mz mzmin mzmax       rt    rtmin    rtmax     into     intb  maxo sn
#> CP0001 594.0 594.0 594.0 2607.809 2587.465 2643.803 161042.2 146073.3  7850 11
#> CP0002 577.0 577.0 577.0 2610.939 2587.465 2632.848 136105.2 128067.9  6215 11
#> CP0003 307.0 307.0 307.0 2625.024 2598.419 2651.628 284782.4 264907.0 16872 20
#> CP0004 302.0 302.0 302.0 2623.459 2601.549 2646.933 687146.6 669778.1 30552 43
#> CP0005 370.1 370.1 370.1 2679.797 2650.063 2706.592 449284.6 417225.3 25672 17
#> CP0006 427.0 427.0 427.0 2681.362 2650.063 2690.804 283334.7 263943.2 11025 13
#>        sample
#> CP0001      1
#> CP0002      1
#> CP0003      1
#> CP0004      1
#> CP0005      1
#> CP0006      1

## Alternatively, get the mapping as a `list` of `integer` indices
featurePeakidx(xmse) |> head()
#> $FT001
#> [1]  458 1161 2677 2849 3167 3256 3369 3530
#> 
#> $FT002
#> [1]   44  443  947 1155 1404 1798 2157 2389
#> 
#> $FT003
#> [1]   29  430 1145 1388 1786 2145 2376 2850
#> 
#> $FT004
#> [1]   16  420  930 1135 1376 1775 2361 3370
#> 
#> $FT005
#> [1]   69  959 1168 1451 1814 2183 2417 2763
#> 
#> $FT006
#> [1]  276  284 1067 1313 1654 1993 2308 2611 2764
#> 
```
