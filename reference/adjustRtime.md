# Alignment: Retention time correction methods.

The `adjustRtime` method(s) perform retention time correction
(alignment) between chromatograms of different samples/dataset.
Alignment is performed by default on MS level 1 data. Retention times of
spectra from other MS levels, if present, are subsequently adjusted
based on the adjusted retention times of the MS1 spectra. Note that
calling `adjustRtime` on a *xcms* result object will remove any
eventually present previous alignment results as well as any
correspondence analysis results. To run a second round of alignment, raw
retention times need to be replaced with adjusted ones using the
[`applyAdjustedRtime()`](https://sneumann.github.io/xcms/reference/applyAdjustedRtime.md)
function.

The alignment method can be specified (and configured) using a dedicated
`param` argument.

Supported `param` objects are:

- `ObiwarpParam`: performs retention time adjustment based on the full
  m/z - rt data using the *obiwarp* method (Prince (2006)). It is based
  on the [original code](http://obi-warp.sourceforge.net) but supports
  in addition alignment of multiple samples by aligning each against a
  *center* sample. The alignment is performed directly on the
  *profile-matrix* and can hence be performed independently of the peak
  detection or peak grouping.

- `PeakGroupsParam`: performs retention time correction based on the
  alignment of features defined in all/most samples (corresponding to
  *house keeping compounds* or marker compounds) (Smith 2006). First the
  retention time deviation of these features is described by fitting
  either a polynomial (`smooth = "loess"`) or a linear
  (`smooth = "linear"`) function to the data points. These are then
  subsequently used to adjust the retention time of each spectrum in
  each sample (even from spectra of MS levels different than MS 1).
  Since the function is based on features (i.e. chromatographic peaks
  grouped across samples) a initial correspondence analysis has to be
  performed **before** using the
  [`groupChromPeaks()`](https://sneumann.github.io/xcms/reference/groupChromPeaks.md)
  function. Alternatively, it is also possible to manually define a
  `numeric` matrix with retention times of markers in each samples that
  should be used for alignment. Such a `matrix` can be passed to the
  alignment function using the `peakGroupsMatrix` parameter of the
  `PeakGroupsParam` parameter object. By default the
  `adjustRtimePeakGroups` function is used to define this `matrix`. This
  function identifies peak groups (features) for alignment in `object`
  based on the parameters defined in `param`. See also
  [`do_adjustRtime_peakGroups()`](https://sneumann.github.io/xcms/reference/do_adjustRtime_peakGroups.md)
  for the core API function.

- `LamaParama`: This function performs retention time correction by
  aligning chromatographic data to an external reference dataset
  (concept and initial implementation by Carl Brunius). The process
  involves identifying and aligning peaks within the experimental
  chromatographic data, represented as an `XcmsExperiment` object, to a
  predefined set of landmark features called "lamas". These landmark
  features are characterized by their mass-to-charge ratio (m/z) and
  retention time. see
  [`LamaParama()`](https://sneumann.github.io/xcms/reference/LamaParama.md)
  for more information on the method.

## Usage

``` r
adjustRtime(object, param, ...)

adjustRtimePeakGroups(object, param, ...)

# S4 method for class 'MsExperiment,ObiwarpParam'
adjustRtime(object, param, chunkSize = 2L, BPPARAM = bpparam())

# S4 method for class 'MsExperiment,PeakGroupsParam'
adjustRtime(object, param, msLevel = 1L, ...)

PeakGroupsParam(
  minFraction = 0.9,
  extraPeaks = 1,
  smooth = "loess",
  span = 0.2,
  family = "gaussian",
  peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
  subset = integer(),
  subsetAdjust = c("average", "previous")
)

ObiwarpParam(
  binSize = 1,
  centerSample = integer(),
  response = 1L,
  distFun = "cor_opt",
  gapInit = numeric(),
  gapExtend = numeric(),
  factorDiag = 2,
  factorGap = 1,
  localAlignment = FALSE,
  initPenalty = 0,
  subset = integer(),
  subsetAdjust = c("average", "previous"),
  rtimeDifferenceThreshold = 5
)

# S4 method for class 'OnDiskMSnExp,ObiwarpParam'
adjustRtime(object, param, msLevel = 1L)

# S4 method for class 'ObiwarpParam'
binSize(object) <- value

# S4 method for class 'XCMSnExp,PeakGroupsParam'
adjustRtime(object, param, msLevel = 1L)

# S4 method for class 'XCMSnExp,ObiwarpParam'
adjustRtime(object, param, msLevel = 1L)
```

## Arguments

- object:

  For `adjustRtime`: an
  [`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html),
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md),
  [`MsExperiment::MsExperiment()`](https://rdrr.io/pkg/MsExperiment/man/MsExperiment.html)
  or
  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  object.

- param:

  The parameter object defining the alignment method (and its setting).

- ...:

  ignored.

- chunkSize:

  For `adjustRtime` if `object` is either an `MsExperiment` or
  `XcmsExperiment`: `integer(1)` defining the number of files (samples)
  that should be loaded into memory and processed at the same time.
  Alignment is then performed in parallel (per sample) on this subset of
  loaded data. This setting thus allows to balance between memory demand
  and speed (due to parallel processing). Because parallel processing
  can only performed on the subset of data currently loaded into memory
  in each iteration, the value for `chunkSize` should match the defined
  parallel setting setup. Using a parallel processing setup using 4 CPUs
  (separate processes) but using
  `chunkSize = `1`will not perform any parallel processing, as only the data from one sample is loaded in memory at a time. On the other hand, setting`chunkSize\`
  to the total number of samples in an experiment will load the full MS
  data into memory and will thus in most settings cause an out-of-memory
  error.

- BPPARAM:

  parallel processing setup. Defaults to `BPPARAM = bpparam()`. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details.

- msLevel:

  For `adjustRtime`: `integer(1)` defining the MS level on which the
  alignment should be performed.

- minFraction:

  For `PeakGroupsParam`: `numeric(1)` between 0 and 1 defining the
  minimum required proportion of samples in which peaks for the peak
  group were identified. Peak groups passing this criteria will be
  aligned across samples and retention times of individual spectra will
  be adjusted based on this alignment. For `minFraction = 1` the peak
  group has to contain peaks in all samples of the experiment. Note that
  if `subset` is provided, the specified fraction is relative to the
  defined subset of samples and not to the total number of samples
  within the experiment (i.e., a peak has to be present in the specified
  proportion of subset samples).

- extraPeaks:

  For `PeakGroupsParam`: `numeric(1)` defining the maximal number of
  additional peaks for all samples to be assigned to a peak group
  (feature) for retention time correction. For a data set with 6
  samples, `extraPeaks = 1` uses all peak groups with a total peak count
  `<= 6 + 1`. The total peak count is the total number of peaks being
  assigned to a peak group and considers also multiple peaks within a
  sample that are assigned to the group. This parameter is ignored for
  `adjustRtime()` on an
  [`XcmsExperimentHdf5()`](https://sneumann.github.io/xcms/reference/XcmsExperimentHdf5.md).

- smooth:

  For `PeakGroupsParam`: `character(1)` defining the function to be used
  to interpolate corrected retention times for all peak groups. Can be
  either `"loess"` or `"linear"`.

- span:

  For `PeakGroupsParam`: `numeric(1)` defining the degree of smoothing
  (if `smooth = "loess"`). This parameter is passed to the internal call
  to [`stats::loess()`](https://rdrr.io/r/stats/loess.html).

- family:

  For `PeakGroupsParam`: `character(1)` defining the method for loess
  smoothing. Allowed values are `"gaussian"` and `"symmetric"`. See
  [`stats::loess()`](https://rdrr.io/r/stats/loess.html) for more
  information.

- peakGroupsMatrix:

  For `PeakGroupsParam`: optional `matrix` of (raw) retention times for
  the (marker) peak groups on which the alignment should be performed.
  Each column represents a sample, each row a feature/peak group. The
  `adjustRtimePeakGroups` method is used by default to determine this
  matrix on the provided `object`.

- subset:

  For `ObiwarpParam` and `PeakGroupsParam`: `integer` with the indices
  of samples within the experiment on which the alignment models should
  be estimated. Samples not part of the subset are adjusted based on the
  closest subset sample. See *Subset-based alignment* section for
  details.

- subsetAdjust:

  For `ObiwarpParam` and `PeakGroupsParam`: `character(1)` specifying
  the method with which non-subset samples should be adjusted. Supported
  options are `"previous"` and `"average"` (default). See *Subset-based
  alignment* section for details.

- binSize:

  `numeric(1)` defining the bin size (in mz dimension) to be used for
  the *profile matrix* generation. See `step` parameter in
  [profile-matrix](https://sneumann.github.io/xcms/reference/profMat-xcmsSet.md)
  documentation for more details.

- centerSample:

  `integer(1)` defining the index of the center sample in the
  experiment. It defaults to
  `floor(median(1:length(fileNames(object))))`. Note that if `subset` is
  used, the index passed with `centerSample` is within these subset
  samples.

- response:

  For `ObiwarpParam`: `numeric(1)` defining the *responsiveness* of
  warping with `response = 0` giving linear warping on start and end
  points and `response = 100` warping using all bijective anchors.

- distFun:

  For `ObiwarpParam`: `character(1)` defining the distance function to
  be used. Allowed values are `"cor"` (Pearson's correlation),
  `"cor_opt"` (calculate only 10% diagonal band of distance matrix;
  better runtime), `"cov"` (covariance), `"prd"` (product) and `"euc"`
  (Euclidian distance). The default value is `distFun = "cor_opt"`.

- gapInit:

  For `ObiwarpParam`: `numeric(1)` defining the penalty for gap opening.
  The default value for depends on the value of `distFun`:
  `distFun = "cor"` and `distFun = "cor_opt"` it is `0.3`, for
  `distFun = "cov"` and `distFun = "prd"` `0.0` and for
  `distFun = "euc"` `0.9`.

- gapExtend:

  For `ObiwarpParam`: `numeric(1)` defining the penalty for gap
  enlargement. The default value for `gapExtend` depends on the value of
  `distFun`: for `distFun = "cor"` and `distFun = "cor_opt"` it is
  `2.4`, `distFun = "cov"` `11.7`, for `distFun = "euc"` `1.8` and for
  `distFun = "prd"` `7.8`.

- factorDiag:

  For `ObiwarpParam`: `numeric(1)` defining the local weight applied to
  diagonal moves in the alignment.

- factorGap:

  For `ObiwarpParam`: `numeric(1)` defining the local weight for gap
  moves in the alignment.

- localAlignment:

  For `ObiwarpParam`: `logical(1)` whether a local alignment should be
  performed instead of the default global alignment.

- initPenalty:

  For `ObiwarpParam`: `numeric(1)` defining the penalty for initiating
  an alignment (for local alignment only).

- rtimeDifferenceThreshold:

  For `ObiwarpParam`: `numeric(1)` defining the threshold to identify a
  *gap* in the sequence of retention times of (MS1) spectra of a
  sample/file. A gap is defined if the difference in retention times
  between consecutive spectra is `> rtimeDifferenceThreshold` of the
  median observed difference or retenion times of that data sample/file.
  Spectra with an retention time after such a *gap* will not be
  adjusted. The default for this parameter is
  `rtimeDifferenceThreshold = 5`. For Waters data with lockmass scans or
  LC-MS/MS data this might however be a too low threshold and it should
  be increased. See also [issue
  \#739](https://github.com/sneumann/xcms/issues/739).

- value:

  For all assignment methods: the value to set/replace.

## Value

`adjustRtime` on an `OnDiskMSnExp` or `XCMSnExp` object will return an
`XCMSnExp` object with the alignment results.

`adjustRtime` on an `MsExperiment` or `XcmsExperiment` will return an
`XcmsExperiment` with the adjusted retention times stored in an new
*spectra variable* `rtime_adjusted` in the object's `spectra`.

`ObiwarpParam`, `PeakGroupsParam` and `LamaParama` return the respective
parameter object.

`adjustRtimeGroups` returns a `matrix` with the retention times of
*marker* features in each sample (each row one feature, each row one
sample).

## Subset-based alignment

All alignment methods allow to perform the retention time correction on
a user-selected subset of samples (e.g. QC samples) after which all
samples not part of that subset will be adjusted based on the adjusted
retention times of the *closest* subset sample (close in terms of index
within `object` and hence possibly injection index). It is thus
suggested to load MS data files in the order in which their samples were
injected in the measurement run(s).

How the non-subset samples are adjusted depends also on the parameter
`subsetAdjust`: with `subsetAdjust = "previous"`, each non-subset sample
is adjusted based on the closest *previous* subset sample which results
in most cases with adjusted retention times of the non-subset sample
being identical to the subset sample on which the adjustment bases. The
second, default, option is `subsetAdjust = "average"` in which case each
non subset sample is adjusted based on the average retention time
adjustment from the previous and following subset sample. For the
average, a weighted mean is used with weights being the inverse of the
distance of the non-subset sample to the subset samples used for
alignment.

See also section *Alignment of experiments including blanks* in the
*xcms* vignette for more details.

## References

Prince, J. T., and Marcotte, E. M. (2006) "Chromatographic Alignment of
ESI-LC-MS Proteomic Data Sets by Ordered Bijective Interpolated Warping"
*Anal. Chem.*, 78 (17), 6140-6152. doi:
[10.1021/ac0605344](https://doi.org/10.1021/ac0605344)

Smith, C.A., Want, E.J., O'Maille, G., Abagyan, R. and Siuzdak, G.
(2006). "XCMS: Processing Mass Spectrometry Data for Metabolite
Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
*Anal. Chem.* 78:779-787. doi:
[10.1021/ac051437y](https://doi.org/10.1021/ac051437y)

## See also

[`plotAdjustedRtime()`](https://sneumann.github.io/xcms/reference/plotAdjustedRtime.md)
for visualization of alignment results.

## Author

Colin Smith, Johannes Rainer, Philippine Louail, Carl Brunius
