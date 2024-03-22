# xcms 4.1

## Changes in version 4.1.12

- Implementation of the `LamaParama` class and method for the `adjustRtime()`
 function. Allowing alignment of a dataset based on landmarks (lamas) from an
 external reference dataset.
- Implementation of related user-level function `matchLamasChromPeaks()`,
  `summarizeMatchLama()` and `plot(LamaParama)` which allows for evaluation of
  matching between lamas and chromPeaks.

## Changes in version 4.1.11

- Clean up of required and suggested packages and namespace imports.
- Re-creation of bundled data objects.

## Changes in version 4.1.10

- Ensure backward compatibility for parameter objects that gained additional
  slots.

## Changes in version 4.1.9

- Fix bug in `filterFeatures,PercentMissingFilter`.

## Changes in version 4.1.8

- Fixing issue #716: edit of `.empty_chrom_peaks` function so an `sn` column is
  returned. Fixes extracting and plotting of peaks after using
  `manualChromPeaks`

## Changes in version 4.1.7

- Implementation of `filterFeatures` function with `filter` parameters:
  `RsdFilter`, `DratioFilter`, `PercentMissingFilter`, `BlankFlag`. They can be
  used ot filter features from `XcmsResult` and `SummarizedExperiment` objects.
- Addition of a section in the main xcms vignette to describe how to use it.

## Changes in version 4.1.6

- Import `filterSpectra` from `MsExperiment`.
- Import `breaks_ppm` from `MsCoreUtils`.
- Update `featureArea` function to consider all chromatographic peaks per
  feature, not only the one with the highest intensity. As a consequence,
  returned m/z and rt ranges might be higher which has an influence in
  `featureChromatograms`, EIC-based feature grouping and, to a lesser extent
  also in gap-filling. Related documentation was updated.
- Improve performance of the `featureArea` function (and related of the
  `PeakAreaParam`-based gap filling).
- Add parameter `ppm` to `PeakDensityParam` to enable peak-density-based
  correspondence throgh m/z-dependent bins along the m/z.

## Changes in version 4.1.5

- Improve performance of the `chromatogram` call for `XcmsExperiment` objects.
- Remove internal (not exported) normalization functions. These have been
  transferred to the MetaboCoreUtils package.
- Support subsetting of `XcmsExperiment` with negative indices.

## Changes in version 4.1.4

- Rename variable `data` in the vignette to `faahko`.
- Fix issue in adjustRtime resulting in corrupt processHistory.
- Add support to perform peakGroups alignment using pre-defined anchor peak
  matrix (i.e., the numeric matrix with retention times of anchor peaks in
  the samples that can be used to align these samples).
- Fix errors related to invalid `Chromatogram` objects extracted from xcms
  results: ensure MS level in `chromPeaksMatrix` is `integer`.
- Fix definition of anchor peaks for peakGroups alignment with subset
  (issue #702).
- Add `filterMsLevel` method for `MsExperiment` and `XcmsExperiment`.
- Ensure chunk-wise processing of Spectra (introduced with version 1.13.2) is
  disabled when xcms is using its own chunk-wise processing.

## Changes in version 4.1.3

- Add parameter `verboseBetaColumns` to `CentWaveParam` to enable calculation
  of additional peak quality metrics comparing the EIC to an idealized bell
  curve.

## Changes in version 4.1.2

- Add a `param =` to generic function `storeResults`: `PlainTextParam` to save
  an `XcmsExperiment` or `MsExperiment` object as colleciton of plain text
  files.

## Changes in version 4.1.1

- Add method `storeResults` and one of its `param =`: `RDataParam` to save an
  `XcmsExperiment` object as an .RData file.


# xcms 3.99

## Changes in version 3.99.6

- Add method to coerce a `XcmsExperiment` to a `xcmsSet` (issue #696).
- Support providing only `mz` or `rt` also for `chromatogram,MsExperiment`.

## Changes in version 3.99.5

- Only `mz` or `rt` need to be provided for `chromatogram`.

## Changes in version 3.99.4

- Add `chromPeakChromatograms` function to extract (EIC) chromatograms for
  chromatographic peaks.

## Changes in version 3.99.3

- Small fixes in the *direct injection* vignette.
- Add parameter `isolationWindowTargetMz` to the `chromatogram` function for
  `MsExperiment` and `XcmsExperiment` to ensure MS2 chromatographic data is
  extracted from the MS2 spectra containing fragments of the compound of
  interest.

## Changes in version 3.99.2

- Add the `xmse` data set representing an `XcmsExperiment` object.
- Update the *compounding* vignette to use the new objects.
- Add `loadXcmsData` to load test data objects (and fix/update paths).
- Add `groupFeatures` methods for `XcmsExperiment`.
- Fix issue in `featureArea` for `XcmsExperiment`.
- Update main vignette to use and describe the new data objects.
- Add `findChromPeaksIsolationWindow` method for `MsExperiment` and
  `XcmsExperiment`.
- Make `reconstructChromPeakSpectra` a method.
- Add `reconstructChromPeakSpectra` implementation for `XcmsExperiment`.
- Add `filterIsolationWindow` for `MsExperiment` and `XcmsExperiment` to filter
  spectra (and eventually chromatographic peaks) based on the isolation window.
- Update the LC-MS/MS vignette adding also an example how to deisotope SWATH
  MS2 spectra.

## Changes in version 3.99.1

- `featureSummary` and `overlappingFeatures` gain support for `XcmsExperiment`.
- Fix in `featureChromatograms` to ensure a valid object is returned.

## Changes in version 3.99.0

- Add `XcmsExperiment` and support for `MsExperiment`/`Spectra`: add all
  functionality for a full xcms processing on a `MsExperiment` object.
- Fix issue in `refineChromPeaks` with `MergeNeighboringPeaksParam` where a
  wrong apex position was considered in the evaluation whether candidate peaks
  should be merged (would only happen for merging of > 2 candidate peaks).
- Re-write the `reconstructChromPeakSpectra` for DIA data analysis to fix an
  issue with chromatographic peaks in overlapping SWATH isolation windows and
  generally to improve performance.


# xcms 3.21

## Changes in version 3.21.5

- Fix issue in `chromatogram` after filtering a result object (issue #511).

## Changes in version 3.21.4

- Move multtest from Suggests to Imports in dependencies

## Changes in version 3.21.3

- Only fixes in the long running tests

## Changes in version 3.21.1

- Fix error with `fillChromPeaks` on sparse data (many empty spectra) and peak
  detection performed with `MatchedFilterParam` (issue #653).
- Update to newer function names in the `rgl` package (issue #654).


# xcms 3.19

## Changes in version 3.19.2

- Update/expand documentation for the `firstBaselineCheck` parameter of
  centWave.

## Changes in version 3.19.1

- Update documentation to reference updates in MassSpecWavelet package.


# xcms 3.17

## Changes in version 3.17.6

- Rewrite code to subset features and chromatographic peaks. This results in a
  perfomance improvement for `filterFile` and similar functions.
- Add parameter `expandMz` to `featureChromatograms`
  https://github.com/sneumann/xcms/issues/612.

## Changes in version 3.17.5

- Change the way the m/z value for a chromatographic peak is determined by
  centWave: if a ROI contains more than one peak for one scan (spectrum) an
  intensity-weighted m/z is reported for that scan. The m/z of the
  chromatographic peak is then calculated based on these reported m/z values for
  each scan (spectrum). In the original version the mean m/z for a scan was
  reported instead. As a result, m/z values of chromatographic peaks are now
  slightly different but are expected to be more accurate. See
  https://github.com/sneumann/xcms/issues/590 for more details.

## Changes in version 3.17.4

- Add `transformIntensity` method.
- Fix issue when calling `chromPeakSpectra` or `featureSpectra` on an object
  that contains also files with only MS1 spectra
  (https://github.com/sneumann/xcms/issues/603).

## Changes in version 3.17.2

- Use mzML instead of mzData files in testing and vignettes,
  since mzR drop mzData reading and msdata package will drop mzData files as well

## Changes in version 3.17.1

- Fix bug in feature grouping by EIC correlation that would return a
  non-symmetric similarity matrix.
- Fix error message from issue [584](https://github.com/sneumann/xcms/issues/584).


# xcms 3.15

## Changes in version 3.15.5

- Disable testing on windows i386, providing some speedup
- Disable parallel processing on Windows, causing an issue in testthat on BioC build check

## Changes in version 3.15.4

- Fix in `plot` with `type = "XIC"` to plot an empty plot if no data is present.
- Skip re-indexing of peaks to features if not necessary. This results in
  performance improvements for MS1 only data.

## Changes in version 3.15.3

- Add `manualFeatures` allowing to manually define and add features to an
  `XCMSnExp` object.
- Add `plotChromatogramsOverlay` function to support plotting of multiple EICs
  from the same sample into the same plot (eventually stacked).
- Add feature grouping by EIC similarity: `EicSimilarityParam`.
- Import `compareChromatograms` from `MSnbase`.
- Add feature grouping by similar retention time: `SimilarRtimeParams.
- Add feature grouping by similarity of feature abundances across samples:
  `AbundanceSimilarityParam`.
- Add feature grouping methodology based on `MsFeatures`.

## Changes in version 3.15.2

- Fix LC-MS/MS vignette.

## Changes in version 3.15.1

- Compatibility fix for nls() in R >= 4.1, contributed by Rick Helmus.


# xcms 3.13

## Changes in version 3.13.8

- Fix plotQC() for XCMSnExp objects

## Changes in version 3.13.7

- Add `featureArea` function to extract the m/z-rt region for features.
- Fix `featureSpectra` function.
- Re-add the LC-MS/MS vignette.
- Feature: plotQC() supports XCMSnExp objects now

## Changes in version 3.13.6

- Fix issue #545: skip second centWave run with CentWavePredIsoParam in regions
  of interest with undefined peak boundaries/scan ranges.
- Temporarily remove the LC-MS/MS vignette (until MsBackendMgf is added to
  Bioconductor).

## Changes in version 3.13.5

- Add `filterChromPeaks` method to filter chromatographic peaks in a
  `XChromatogram` or `XChromatograms` object.
- Add `filterChromPeaks` method for `XCMSnExp` (issue #541).
- Support return of `Spectra` objects by `chromPeakSpectra`, `featureSpectra`
  and `reconstructChromPeakSpectra`.
- Support extraction of MS1 spectra with `chromPeakSpectra`.
- Support extraction of the spectrum with the largest total signal or largest
  base peak signal in `chromPeakSpectra`.
- Add support for extraction of spectra for selected/individual peaks/features
  using the `peaks` and `features` parameter in `chromPeakSpectra` and
  `featureSpectra`, respectively.

## Changes in version 3.13.4

- Import `Param` object from `ProtGenerics`.
- Import `filterIntensity`, `normalize` and `alignRt` for `Chromatogram` and
  `MChromatograms` from `MSnbase`.

## Changes in version 3.13.3

- `align,Chromatogram` gains new method `"none"` which will only keep values
  with identical retention times. For `method = "matchRtime"` the (much faster)
  matching function `closest` from the `MsCoreUtils` package is used.
- Method `correlate,Chromatogram` gains parameter `useIntensitiesAbove` to
  perform the correlation only with values larger than this threshold
  (avoiding thus high correlation because of many 0-values).
- Add method `filterIntensity,Chromatogram` that allows to filter a chromatogram
  object keeping only data points with an intensity above a user provided
  threshold.

## Changes in version 3.13.2

- Add new function `manualChromPeaks` allowing to manually add and integrate
  chromatographic peaks.

## Changes in version 3.13.1

- Support subsetting of `XChromatograms` with `drop = FALSE`.


# xcms 3.11

## Changes in version 3.11.8

- Disable parallel processing in vignettes.

## Changes in version 3.11.7

- More efficient splitting data per file especially for larger data sets.
- Disable parallel processing in examples.

## Changes in version 3.11.6

- Add `FilterIntensityParam` to filter chromatographic peaks on intensity
  (issue #502).
- Add `estimatePrecursorIntensity` function to determine the precursor intensity
  for MS2 spectra from the neighboring MS1 spectra.

## Changes in version 3.11.4

- Change from `Spectra` and `Chromatograms` to `MSpectra` and `MChromatograms`
  from MSnbase version >= 2.15.3.

## Changes in version 3.11.3

- `reconstructChromPeakSpectra`: report also polarity and `precusorIntensity`.
- `reconstructChromPeakSpectra`: ensure a retention time is reported for
  reconstructed MS2 spectra (issue #485).
- Change default for `expandRt` to `0` in `reconstructChromPeakSpectra`.
- Fix error in `refineChromPeaks,MergeNeighboringPeaksParam` if no peaks found
  to be merged.

## Changes in version 3.11.2

- Add `fillChromPeaks,ChromPeakAreaParam` to base the area from which missing
  peak data should be filled-in on the actually detected chromatographic peaks
  of a feature.
- Potential fix for issue #481: function should no longer throw an error because
  retention times are of length 0.
- More efficient splitting of processing which should increase the speed of
  the findChromPeaks, refineChromPeaks, reconstructChromPeakSpectra and
  chromPeakSpectra calls.

## Changes in version 3.11.1

- Fix issue #471: conversion from `XCMSnExp` to `xcmsSet` looses phenodata
  (thanks to Andris Jankevics for reporting and providing a solution).
- Add `normalize` method for `Chromatogram` and `Chromatograms` objects.
- `featureChromatograms` gets new parameter `n` and `value` to extract EICs
  only from the top n samples with highest intensities.
- `filterFile` gets new parameter `keepFeatures` to support retaining
  correspondence results even if a data set is filtered by file.
- Export the virtual `Param` class.
- Add filterColumnsIntensityAbove method for Chromatograms object that allows
  to select columns (samples) of an Chromatograms object for which intensities
  of its chromatographic data are higher than a threshold.
- Add removeIntensity method for Chromatogram, Chromatograms, XChromatogram
  and XChromatograms objects allowing to *remove* intensities based on different
  criteria.
- Add correlate method for Chromatograms allowing to correlate multiple
  chromatograms with each other.


# xcms 3.9

## Changes in version 3.9.4

- Fix issue in centWave which skips peak detection depending on minimum
  peakwidth (issue #445): add parameter `extendLengthMSW` in `CentWaveParam`.
  Thanks to William Kumler for contributing the fix.
- Tentatively reduce memory requirements in `fillChromPeaks`.
- Fix issue #467 for fillPeaks() of an xcmsSet converted from an XCMSnSet

## Changes in version 3.9.3

- Move multtest from Imports to Suggests to avoid duplicated method definition
  for plot (issue #459).
- Add support for peak filling from MS level > 1 to fillChromPeaks.
- featureValues gains parameter msLevel to extract feature values for features
  of all, or from a specific MS level.
- refineChromPeaks supports different MS levels.
- Added support to perform correspondence analysis on MS level > 1 and add the
  respective results to already present feature definitions.
- hasChromPeaks and hasFeatures gain parameter msLevel to check for presence of
  chromatographic peaks or features from a specific MS level.

## Changes in version 3.9.2

- Fix featureChromatograms and chromatograms on a XCMSnExp object with features:
  features can be duplicated across rows (EICs).
- findChromPeaks: add parameter `add` to allow several rounds of peak detections
  on the same object.
- Small performance enhancement in fillChromPeaks.
- Better support for MS > 1 data in fillChromPeaks: skip MS level 2 spectra for
  filling in.
- Add refineChromPeaks for XChromatogram and XChromatograms objects.
- Add groupOverlaps function to group arbitrary ranges.
- Add quantify,XCMSnExp object to quantify an XCMSnExp into a
  SummarizedExperiment.
- Fine-tune MergeNeighboringPeaks peak refinement method: the average of the
  3 data points between candidate peaks is used to evaluate whether the peaks
  should be merged making the approach more robust against outliers. In
  addition, an ion chromatogram for candidate peaks is extracted with an m/z
  range expanded depending on the expandMz and ppm setting ensuring that low
  intensity data points between candidate peaks are not missed out (because
  their m/z might be slightly shifted on ToF instruments). The mzmin and mzmax
  of the merged peak represents also the minimum and maximum m/z of all data
  points in that extracted ion chromatogram.

## Changes in version 3.9.1

- Fix problem of not shown/plotted peak positions in plotChromPeakSpectra
  for experiments in which peaks were not detected in the first sample(s).
- Add method *from_to* to missing value imputation method `imputeRowMinRand`.
- Show warning in findChromPeaks if empty spectra are detected.
- Add refineChromPeaks method and CleanPeaksParam class to allow removal of
  chromatographic peaks exceeding a user-definable maximal peak width.
- Add MergeNeighboringPeaksParam for refineChromPeaks to allow merging of
  chromatographic peaks close in m/z and retention time with a signal between
  them higher than a certain threshold (issue #414).
- Fix misspelled parameter `mzd` in LC-MS/MS vignette.


# xcms 3.7

## Changes in version 3.7.5

- Remove xcmsMSn vignette (based on old xcms).

## Changes in version 3.7.4

- mzClust correspondence analysis: check and fix missing values in column mz of
  the peaks matrix (issue #416).

## Changes in version 3.7.3

- plot type = "XIC" on an XCMSnExp object will draw rectangles indicating the
  identified chromatographic peaks.
- Add a vignette describing LC-MS/MS data analysis with xcms.

## Changes in version 3.7.2

- Fix documentation (issue #401).
- Add support for SWATH data analysis.

## Changes in version 3.7.1

- Add correlate method for Chromatogram objects.
- Add parameter lwd to plotAdjustedRtime.
- Add align method for Chromatogram objects.
- Add findChromPeaksIsolationWindow to enable chromatographic peak detection
  in isolation windows.
- Fix issue in chromPeakSpectra with method = "signal".
- chromPeakSpectra and featureSpectra return now MS2 spectra with an precursor
  m/z >= mzmin, <= mzmax and retention time >= rtmin, <= rtmax.
- Improve performance of chromPeakSpectra and featureSpectra.


# xcms 3.5

## Changes in version 3.5.5

- Add dirname and dirname<- methods for OnDiskMSnExp to change the path to the
  raw data files.
- Add section "Subset-based alignment" to the xcms vignette to describe the
  alignment possibility to perform alignments based on a subset of samples
  (e.g. QC samples).

## Changes in version 3.5.4

- Fix problem in featureChromatograms with include = "feature_only" that could
  return a non-valid object.
- Ensure that XCMSnExp objects are updated if necessary in all analysis methods.

## Changes in version 3.5.3

- Fix unit tests.

## Changes in version 3.5.2

- Small changes in fillChromPeaks,XCMSnExp to reduce memory demand.
- Fix issue #359.
- Fix issue #360: rawEIC skipped last scan/spectrum if rtrange was provided.
- filterMsLevel keeps now chromatographic peaks and feature definitions from the
  specified MS levels (issue #362).
- Fix bug in `xcmsRaw` that leads to a netCDF error message (issue #363).
- Add parameter msLevel to chromPeaks for XCMSnExp objects.
- Add chromPeakData to allow adding arbitrary annotation to chromatographic
  peaks.
- Change default of parameter value in featureValues from value = "index" to
  value = "into".
- Add parameter isFilledColumn to chromPeaks allowing the old behaviour to
  include the is_filled column in the chromatographic peak matrix.

## Changes in version 3.5.1

- Fix issue #349.
- Add updateObject function for XCMSnExp objects (issue #347).
- Add dropFilledChromPeaks methods for XChromatogram and XChromatograms objects.
- Add parameter filled = FALSE to chromatogram and featureChromatograms
  functions.
- Fix matchedFilter peak detection problems with empty spectra (issue #325).
- featureChromatograms extracts by default only chromatographic peaks associated
  with a feature.
- chromatogram,XCMSnExp extracts an XChromatogram containing also
  chromatographic peaks and feature definitions.
- Add featureValues method for XChromatograms objects (issue #336).
- Add correspondence analysis (peak grouping) for chromatographic data (for now
  only with PeakDensity method; issue #336).
- Add featureDefinitions slot to XChromatograms object and related accessor
  methods.
- Add subset alignment option subsetAdjust = "average" to adjust left-out
  samples (blanks or simply non-subset samples) based on an interpolation from
  the results of the previous and subsequent subset sample.
- Add parameter subsetAdjust to PeakGroupsParam allowing to switch between
  different methods to adjust samples left out in the alignment process.
- Alignment based on a sample subset for the peak groups method (issue #335):
  sample subset can be defined with the subset parameter, samples not included
  in the subset will be aligned based on the adjusted retention times of the
  closest sample in the subset.
- Add findChromPeaks,XChromatograms (issue #332).
- Add processHistory,XChromatograms.
- Add plot,XChromatograms method with automatic peak highlighting (issue #334).
- Add hasChromPeaks,XChromatograms method.
- Add XChromatograms class with constructor function and coercing method.
- Add hasChromPeaks,XChromatogram method.
- Add filterRt,XChromatogram, filterMz,XChromatogram.
- Add plot,XChromatogram function supporting of highlighting/drawing identified
  chromatographic peaks.
- findChromPeaks,Chromatogram returns an XChromatogram object (issue #329).
- Add chromPeaks,XChromatogram (issue #329).
- Add XChromatogram object (issue #329).
- Fix higlightChromPeaks with type = "polygon": peak filling represents now the
  full detected peak and is no longer cut by the provided rt.
- Add argument peakIds to highlightChromPeaks allowing to specify the IDs of
  peaks to be highlighted.
- Add example on clustering of base peak chromatograms to the vignette (issue
  #328).
- Small fix in the vignette (issue #327).
- Add parameter groupval to exportMetaboAnalyst (issue #296).
- Fix bug in show,XCMSnExp that would throw an error if no process history is
  present.


# xcms 3.3

## Changes in version 3.3.6

- Add type = "polygon" to highlightChromPeaks allowing to fill the actual
  signal area of identified chromatographic peaks.

## Changes in version 3.3.5

- Performance enhancement of the chromPeakSpectra and featureSpectra functions.

## Changes in version 3.3.4

- Add featureChromatograms to extract ion chromatograms for each feature.
- Add hasFilledChromPeaks function.
- Add argument skipFilled to the featureSummary function.

## Changes in version 3.3.3

- Add chromPeakSpectra and featureSpectra functions to extract MS2 spectra
  for chromatographic peaks and features, respectively (issue #321).
- Fix profMat to handle also data files with empty spectra (issue #312).
- Add argument ylim to plotAdjustedRtime (issue #314).
- Add imputeRowMin and imputeRowMinRand, two simple missing value imputation
  helper functions.
- Fix additional problem mentioned in issue #301 with obiwarp retention time
  correction if some spectra have m/z values of `NA`.
- Fix issue #300 avoiding chromatographic peaks with rtmin > rtmax.
- Fixes for issues #291, #296.
- Add parameter 'missing' to diffreport allowing to replace NA with arbitrary
  numbers.
- Add exportMetaboAnalyst function to export the feature matrix in MetaboAnalyst
  format.
- Add parameter missing to featureValues allowing to specify how to handle/
  report missing values.
- The chromPeaks matrix has now rownames to uniquely identify chromatographic
  peaks in an experiment. Chromatographic peak IDs start with "CP" followed by
  a number.

## Changes in version 3.3.2

- Add writeMSData method for XCMSnExp allowing to write mzML/mzXML files with
  adjusted retention times (issue #294).
- Fix profEIC call for single-scan-peak (pull request #287 from @trljcl).
- Fix centWave avoiding that the same peak is reported multiple times if
  fitgauss = TRUE is used (issue #284).
- featureSummary reports also RSD (relative standard deviations) of features
  across samples (issue #286).
- Add parameters fixedMz and fixedRt to FillChromPeaksParam that allow to
  increase the features' m/z and rt widths by a constant factor.
- Add option "sum" to featureValues' method parameter allowing to sum the
  intensities of peaks that are assigned to the same feature in a file/sample.

## Changes in version 3.3.1

- Add overlappingFeatures function to identify overlapping or close features.
- Add support for type = "apex_within" for featureDefinitions.
- Fix a bug in fillChromPeaks that would return the integrated signal being Inf.
- Fix for issue #267: error in fillChromPeaks when the retention time of the
  peaks are outside of the retention time range of certain files.
- New featureSummary function to calculate basic feature summaries (number of
  samples in which peaks were found etc).
- Parameter 'type' added to plotChromPeakDensity and 'whichPeaks' to
  highlightChromPeaks. Both parameters are passed to the 'type' argument
  of chromPeaks.
- Parameter 'type' in chromPeaks gets additional option "apex_within" to return
  chromatographic peaks that have their apex within the defined rt and/or m/z
  range.
- Add functions rla and rowRla to calculate RLA (relative log abundances).
- Add peaksWithMatchedFilter to perform peak detection in chromatographic
  (MRM/SRM) data (issues #277 and #278).
- Add peaksWithCentWave to perform centWave peak detection in chromatographic
  (MRM/SRM) data (issue #279).
- Add findChromPeaks,Chromatogram methods for CentWaveParam and
  MatchedFilterParam (issue #280).


# xcms 3.1

## Changes in version 3.1.3

- Fix misplaced parenthesis in the check for multiple spectra in
  findChromPeaks,OnDiskMSnExp,MSWParam. Thanks to @RonanDaly (PR #276).
- Update link to correct metlin page in diffreport result (issue #204).

## Changes in version 3.1.2

- Add filterFeatureDefinitions function.
- Fix #273: better error message in case not a single feature could be defined
  by groupChromPeaks.

## Changes in version 3.1.1

- Reading raw files using xcmsSet or xcmsRaw uses now the automatic file type
  detection feature from mzR.
- c function to concatenate XCMSnExp objects.
- groupnames method for XCMSnExp objects (issue #250).
- Fix #237: findPeaks.MSW was not throwing an error if applied to multi-spectrum
  MS file.
- Fix #249: quantile call in adjustRtime PeakGroups without na.rm = TRUE.
- Fix #259


# xcms 2.99

## Changes in version 2.99.10

- Fix #230: Failing vignettes on Windows.

## Changes in version 2.99.9

- Chromatographic peak detection uses adjusted retention times on an aligned
  XCMSnExp object (issue #213, #208).
- New parameter msLevel for processHistory,XCMSnExp.
- New parameter keepAdjustedRtime for filterMsLevel,XCMSnExp, dropChromPeaks,
  XCMSnExp and dropFeatureDefinitions,XCMSnExp.
- Add parameter msLevel to chromatogram,XCMSnExp method (issue #205).
- Obiwarp alignment is now performed on one MS level and adjustment is applied
  to all MS levels (issue #214).
- Add function plotMsData to plot intensity against retention time and m/z
  against retention time for a MS slice in one sample.
- Add argument msLevel = 1L to extractMsData method (issue #223).
- New applyAdjustedRtime function to consolidate the alignment results, i.e.
  replace the raw retention times in the XCMSnExp with the adjusted retention
  times.
- [,XCMSnExp method gains argument keepAdjustedRtime to allow keeping adjusted
  retention times in the sub-setting.
- Implement spectrapply,XCMSnExp to ensure returned results use adjusted
  retention times (if present).
- [[,XCMSnExp method returns a Spectrum object with adjusted retention time, if
  the XCMSnExp contains adjusted retention times.
- Argument 'sampleGroups' is mandatory for 'PeakDensityParam' (issue #228).
- Fix #191: Excessive memory use in fillPeaks.
- Fix #220: peaks matrix is missing column "sample" if no peaks were found in
  the first sample.
- Fix #222: findChromPeaks does not return an XCMSnExp object filtered to a
  single MS level despite peak detection is performed on a single level.
- Fix problem in plotMsData causing wrong colors to be used to label the data
  points.

## Changes in version 2.99.8

- Replace xcmsMSn Rnw with Rmd vignette to fix Windows build errors.

## Changes in version 2.99.7

- Fix #201: Warnings: 'readMSData2' is deprecated, thanks to L. Gatto.
- Merge with BioC git after transition

## Changes in version 2.99.6

- calibrate,XCMSnExp method that allows to calibrate chromatographic peaks.
- Export phenoDataFromPaths function (issue $195).
- Add arguments mz and rt to featureDefinitions method allowing to extract
  features within the specified ranges.
- Increase n for the density function call in group density-based correspondence
  by 2.
- Replace xcmsDirect.Rnw with rmarkdown-based vignette using the new user
  interface.
- issue #196: removed the unnecessary requirement for same-dimension profile
  matrices in adjustRtime,XCMSnExp,ObiwarpParam.
- issue #194: fixes in retcor.obiwarp: 1) subset raw data if scanrange != NULL.
  2) if the mz range of the two files to be aligned differ, expand them
  correctly. Depending on the profStep and the mz values/ranges the matrices
  were not expanded correctly.
- Potential problems in the plotChromPeakDensity function.

## Changes in version 2.99.5

- Re-enable sleep parameter in findPeaks.centWave and findPeaks.matchedFilter.

## Changes in version 2.99.4

- Add plotChromPeaks function to plot the definition (rt and mz range) of
  detected chromatographic peaks of one file into the mz-rt plane.
- Add plotChromPeakImage function to plot the number of detected peaks along
  the retention time axis per file as an image plot.
- Move Chromatogram class and functionality to the MSnbase package
- Add argument msLevel to the findChromPeaks method to allow (chromatographic)
  peak detection also on MS level > 1.
- Polarity information was not read from mzXML files (issue #192).

## Changes in version 2.99.3

- issue #188: determine file type from file content if file ending not known.

## Changes in version 2.99.2

- issue #181: problem when isCentroided,Spectrum method returns NA because of
  too few peaks in a spectrum. Fixed by checking in such cases all spectra in
  the file.
- issue #184: add parameter sleep to do_groupChromPeaks_density function to be
  backwards compatible with the old group.density code.

## Changes in version 2.99.1

- extractMsData to extract raw MS data as a data.frame (issue #120).
- issue #175: an error is now thrown if no peak group was identified for peak
  group retention time correction.
- issue #178: scanrange was collapsed when the adjusted range was reported
  (pull request by Jan Stanstrup).
- issue #180: error when both parameters method and smooth are provided in the
  retcor method.

## Changes in version 2.99.0

- plotChromatogram and highlightChromPeaks functions.
- plotChromPeakDensity function.
- clean method for Chromatogram classes.
- Change default for ppm parameter in chromPeaks method to 0.
- extractChromatograms supports extraction of multiple rt and mz ranges.
- New parameter missing for extractChromatograms allowing to specify the
  intensity value to be used for rts for which no signal is available within
  the mz range.
- extractChromatograms returns Chromatograms of length equal to the number of
  scans within the specified rt range, even if no signals are measured
  (intensity values are NA).


# xcms 1.53

## Changes in version 1.53.1

- Increase parameter n for the density call in the peak density correspondence
  method. This enables to separate neighboring peaks using small n (issue #161).
  Thanks to Jan Stanstrup.

# xcms 1.51

## Changes in version 1.51.11

- Parameter "filled" for featureValues (issue #157).
- Parameters "rt" and "mz" in chromPeaks method allowing to extract
  chromatographic peaks from the specified ranges (issue #156).
- Fixed possible memory problem in obiwarp (issue #159).
- Update getPeaks to use non-deprecated API (issue #163).

## Changes in version 1.51.10

- filterRt for Chromatogram class (issue #142).
- adjustRtimePeakGroups function (issue #147).
- adjustRtime,XCMSnExp,PeakGroupsParam and do_adjustRtime_peakGroups support
  use of pre-defined matrix to perform alignment (issue #153).
- plotAdjustedRtime to visualize alignment results (issue #141).
- featureDefinitions and featureValues return DataFrame and matrix with rownames
  corresponding to arbitrary feature IDs (issue #148).
- New peakGroupsMatrix slot for PeakGroupsParam class (issue #153).
- Issue #146: ensure adjusted retention times returned by the peakGroups method
  to be in the same order than the raw retention times.

## Changes in version 1.51.9

- fillChromPeaks, dropFilledChromPeaks methods and FillChromPeaksParam class.
- featureValues method.
- Extended new_functionality vignette.
- Change default backend for reading mzML files to pwiz.
- Issue #135: fix peak signal integration for centWave.
- Issue #139: problem with expand.mz and expand.rt in fillPeaks.chrom.
- Issue #137: Error in findChromPeaks if no peaks are found.

## Changes in version 1.51.8

- Add Chromatogram class and extractChromatograms method.
- Issue #118: failing unit test on Windows build machine.
- Issue #133: error with c() and xcmsSet without peaks.
- Issue #134: xcmsSet constructor endless loop.

## Changes in version 1.51.7

- Major renaming of methods and classes to follow the naming convention:
  - chromatographic peak (chromPeak): the peaks identified in rt dimension.
  - feature: mz-rt feature, being the grouped chromatographic peaks within and
    across samples.
- Issue #127: failing unit test on Windows build machine.

## Changes in version 1.51.6

- groupFeatures and adjustRtime methods for XCMSnExp objects.
- New Param classes for groupFeatures and adjustRtime analysis methods:
  FeatureDensityParam, MzClustParam, NearestFeaturesParam, FeatureGroupsParam
  and ObiwarpParam.
- Issue #124 (filterRt,XCMSnExp returned empty object).

## Changes in version 1.51.5

- MsFeatureData and XCMSnExp objects.
- features, features<-, adjustedRtime, adjustedRtime<-, featureGroups,
  featureGroups<-, hasAlignedFeatures, hasAdjustedRtime and hasDetectedFeatures
  methods.
- dropFeatures, dropFeatureGroups and dropAdjustedRtime methods.
- filterMz, filterRt, filterFile etc implemented.
- mz, intensity and rtime methods for XCMSnExp allowing to return values grouped
  by sample.
- Issue #99 (rtrange outside of retention time range in getEIC,xcmsSet).
- Issue #101 (xcmsRaw function returns NULL if mslevel = 1 is specified).
- Issue #102 (centWave returns empty matrix if scales not OK). Thanks to
  J. Stanstrup.
- Issue #91 (warning instead of error if no peaks in ROI). Thanks to J. Stanstrup.

## Changes in version 1.51.4

- added deepCopy to avoid corrupting the original object, thanks to
  J. Stanstrup, closes #93

## Changes in version 1.51.3

- binYonX binning function.
- imputeLinInterpol function providing linear interpolation of missing values.
- breaks_on_binSize and breaks_on_nBins functions to calculate breaks defining
  bins.
- New vignette "new_functionality.Rmd" describing new and modified functionality
  in xcms.
- Add do_detectFeatures_matchedFilter function.
- Add do_detectFeatures_centWave function.
- Add do_detectFeatures_centWaveWithPredIsoROIs function and unit test.
- Implement a new data import function.
- Add do_detectFeatures_MSW function and unit test.
- Argument stopOnError in xcmsSet function that allows to perform feature
  detection on all files without stopping on errors.
- Method showError for xcmsSet objects that list all errors during feature
  detection (if stopOnError = FALSE in the xcmsSet function).
- [ method to subset xcmsRaw objects by scans.
- profMat method to extract/create the profile matrix from/for an xcmsRaw.
- Add new detectFeatures methods for MSnExp and OnDiskMSnExp objects from the
  MSnbase package.
- Add new CentWaveParam, MatchedFilterParam, MassifquantParam, MSWParam and
  CentWavePredIsoParam parameter class to perform method dispatch in the
  detectFeatures method.
- retcor.obiwarp uses the new binning methods for profile matrix generation.
- scanrange,xcmsRaw reports always a scanrange of 1 and length(object@scantime).
- scanrange,xcmsSet reports the scanrange eventually specified by the user in
  the xcmsSet function.
- Fixed bug in rawMat (issue #58).
- Fix issue #60: findPeaks.massifquant always returns a xcmsPeaks object.

## Changes in version 1.51.2

- As suggested by Jan Stanstrup, do not error if a centWave ROI
  contains no data, closes #90

## Changes in version 1.51.1

- Fix incorrrect indexing getEIC function reported by Will Edmands, closes #92


# xcms 1.49

## Changes in version 1.49.7

- Fix documentation warnings.

## Changes in version 1.49.6

- Peak Picking function findPeaks.centWaveWithPredictedIsotopeROIs() and
  findPeaks.addPredictedIsotopeFeatures(), which allow more sensitive detection
  of isotope features.

## Changes in version 1.49.5

- Some documentation updates.
- Preparation for a new binning function

## Changes in version 1.49.4

- Fix getXcmsRaw that would prevent retention time correction to be applied
  (issue #44 reported by Aleksandr).

## Changes in version 1.49.3

- updateObject method for xcmsSet.
- xcms uses now BiocParallel for parallel processing. All other parallel
  processing functions have been deprecated.
- Added missing package imports.
- Fix bug in fillPeaksChromPar referencing a non-existing variables i and
  object.
- Fix bug in group.nearest: variable scoreList was mis-spelled (coreList).
- Remove all DUP = FALSE from the .C calls as they are ignored anyways.
- Re-organization of class, function and method definitions in R-files.
- Use roxygen2 to manage the DESCRIPTION's collate field.

## Changes in version 1.49.2

- Initial support for exporint mzTab format. Since Changes are
  still to be expected, xcms:::writeMzTab() is not yet exported.

## Changes in version 1.49.1

- The raw CDF/mzXML/mzData/mzML is assumed to have scans sorted by m/z.
  Instead of throwing an "m/z sort assumption violated !" error,
  the data is re-read and on-demand sorted by m/z.


# xcms 1.47

## Changes in version 1.47.3

- Disable parallel processing in unit tests causing a timeout
  on BioC build machines

## Changes in version 1.47.2

- Fix problem in getEIC on xcmsSet objects reported by Alan Smith in issue #7 and
  add a RUnit test case to test for this (test.issue7 in runit.getEIC.R).
- Changed some unnecessary warnings into messages.

## Changes in version 1.47.2

- Disabled parallel processing in unit tests
- migrate dependencies from ncdf -> ncdf4


# xcms 1.45

## Changes in version 1.45.7

- Disabled Rmpi support and usage on Windows

## Changes in version 1.45.6

- J. Rainer implemented a [ method that allows to subset an xcmsSet.
- Fixed a problem in split.xcmsSet that did not split the phenoData properly.
  Added some details to the documentation of xcmsSet-class.

## Changes in version 1.45.5

- The sampclass method for xcmsSet will now return the content of the
  column "class" from the data.frame in the phenoData slot, or if not
  present, the interaction of all factors (columns) of that data.frame.
- The sampclass<- method replaces the content of the "class" column in
  the phenoData data.frame. If a data.frame is submitted, the interaction
  of its columns is calculated and stored into the "class" column.
- Fixed a bug that resulted in a cryptic error message
  when no input files are available to the xcmsSet function.

## Changes in version 1.45.4

- Fixed a bug in the levelplot method for xcmsSet.

## Changes in version 1.45.3

- xcmsSet now allows phenoData to be an AnnotatedDataFrame.
- new slots for xcmsRaw:
  - mslevel: store the mslevel parameter submitted to xcmsRaw.
  - scanrange: store the scanrange parameter submitted to xcmsRaw.
- new slots for xcmsSet:
  - mslevel: stores the mslevel argument from the xcmsSet method.
  - scanrange: to keep track of the scanrange argument of the xcmsSet method.
- new methods for xcmsRaw:
  - levelplot: similar to the image method, plots m/z vs RT with color coded
    intensities.
  - mslevel: returns the value for the .mslevel slot. For downstream
    compatibility, this method returns NULL if the object does not have the same
    named slot.
  - profinfo: same functionality as the profinfo method for xcmsSet.
  - scanrange: returns the value for the scanrange slot. For downstream
    compatibility, this method returns NULL if the object does not have the same
    named slot.
- new methods for xcmsSet:
  - getXcmsRaw: returns a xcmsRaw object for one or more files in the xcmsSet,
    eventually applying retention time correction etc.
  - levelplot: similar to the image method, plots m/z vs RT with color coded
    intensities. Allows in addition to highlight identified peaks.
  - mslevel: returns the value for the mslevel slot. For downstream
    compatibility, this method returns NULL if the object does not have the same
    named slot.
  - profMethod: same functionality as the profMethod method of xcmsRaw.
  - profStep: same functionality as the profStep method of xcmsRaw.
  - scanrange: returns the value for the scanrange slot. For downstream
    compatibility, this method returns NULL if the object does not have the same
    named slot.
- show method for xcmsSet updated to display also informations about the mslevel
  and scanrange.
- Elaborated some documentation entries.
- rtrange and mzrange for xcmsRaw method plotEIC use by default the full RT and
  m/z range.
- Added arguments "lty" and "add" to plotEIC method for xcmsRaw.
- getEIC without specifying mzrange returns the ion chromatogram for the full
  m/z range (i.e. the base peak chromatogram).
- Checking if phenoData is a data.frame or AnnotatedDataFrame and throw an error
  otherwise.
- xcmsSet getEIC method for water Lock mass corrected files for a subset of
  files did not evaluate whether the specified files were corrected.

## Changes in version 1.45.2

- The xcms split() function now accepts factors that are shorter than the number
  of samples in the xcmsSet, following more closely the standard split()
  behaviour

## Changes in version 1.45.1

- plotrt now allows col to be a vector of color definition,
  same as the plots for retcor methods.
- Added $ method to access phenoData columns in a eSet/ExpressionSet like
  manner.
- Allow to use the "parallel" package for parallel processing of the functions
  xcmsSet and fillPeaks.chrom.
- Thanks to J. Rainer!


# xcms 1.43

## Changes in version 1.43.3

- Give a more verbose error message when file not found

## Changes in version 1.43.2

- Use ProtGenerics, adapted xcms peaks()

## Changes in version 1.43.1

- function plotQC() for plotting various QC plots on RT and m/z


# xcms 1.41

## Changes in version 1.41.1

- fix sampclass generation from phenoData if some combinations of factors don't
  exist
- disable parallel code in manpages to avoid issues on BioC windows build farm
  machines


# xcms 1.39

## Changes in version 1.39.6

- Massifquant reports the maximum intensity for each isotope trace (peak). This
  is useful for interactive parameter optimization.
- Major memory reduction in parallel fillPeaks() thanks to Jan Stanstrup. Now
  using an environment to mirror gvals to each list item in the very large
  argList.

## Changes in version 1.39.4

- Fixed write.cdf(), which had an intensity offset of +1, added a unit test

## Changes in version 1.39.3

- New R-devel check unload better. Lingering ramp code removed, import from
  mzR. Cleaned up other errors in package check.

## Changes in version 1.39.1

- Updated doubleMatrix c code to allow for larger profile matrixes
- Moved inst/doc to vignettes


# xcms 1.37

## Changes in version 1.37.6

- Introducing write.mzQuantML(xcmsSet) to export the peak list and grouped
  matrix to the PSI format mzQuantML (see http://www.psidev.info/mzquantml)
- Add Brigham Young University to LICENSE file for copyright purposes.
- Add copyright information display when running findPeaks.massifquant()
  within xcmsRaw.R
- Clean and update documentation for findPeaks.massifquant-methods.Rd
- Remove unused parameters in findKalmanROIs() within xcmsRaw.R

## Changes in version 1.37.5

- fixed bug in retcor.obiwarp where the scanrange of the first sample would be
  checked instead of the center sample

## Changes in version 1.37.4

- Skip t-test in diffreport() if one class has less than 2 samples.

## Changes in version 1.37.3

- fixed bug in patternVsRowScore (group.nearest) that was introduced by the
  modifications in rev 65169 and caused features to be aligned that were far
  outside the given m/z and retention time windows.

## Changes in version 1.37.1

- fixed fillPeaks, which 1) dropped non-standard columns and 2) failed if
  nothing to do, based on patches by Tony Larson.

## Changes in version 1.37.1

- Introducing msn2xcmsRaw, to allow findPeaks() on MS2 and MSn data


# xcms 1.35

## Changes in version 1.35.7

- fixed indexing bug in group.nearest, which under certain circumstances caused
  all peaks in the first sample to be ignored (reported by Tony Larson)

## Changes in version 1.35.6

- Obiwarp retention time alignment error-ed if scanrange was used as a parameter
  setting during xcmsSet/peak detection The method now tries to automatically
  find the set scanrange and uses this range for alignment.

## Changes in version 1.35.4

- Introducing parallel fillPeaks
- Replace snow requirement with minimum R version 2.14.0

## Changes in version 1.35.3

- if group.density was used with very low minfrac settings (< 0.5) it did not
  return all feature groups, but only those that include features from at least
  50% of samples in a group. This limitation was removed.

## Changes in version 1.35.2

- Behind the scenes xcms now uses the xcmsSource class to read raw data.  This
  allows e.g. to write a class that pulls raw data from e.g. a database
- massifquant: simplified logic structure of Tracker::claimDataIdx resolved
  failure on new test case.
- massifquant: reporting features data structure compatible with multiple sample
  comparison within XCMS.

## Changes in version 1.35.1

- The mzData export is now much faster and uses less memory


# xcms 1.33

## Changes in version 1.33.16

- diffreport and plotEIC have a new parameter mzdec, with is the number of
  decimal places of the m/z values in the EIC plot title

## Changes in version 1.33.16

- Lock mass gap filler now works with netCDF lock mass function file to find the
  exact times of the scans and works with the newer Waters MS instruments.

## Changes in version 1.33.15

- scanrage is now honoured in xcmsSet, also when in parallel mode

## Changes in version 1.33.14

- scanrage is now honoured in xcmsRaw, and consequently also in
  xcmsSet(matchedFilter), where previously it was ignored.

## Changes in version 1.33.13

- write.cdf() has been fixed to write files AMDIS can read

## Changes in version 1.33.12

- write.mzData adds Polarity to the file if available

## Changes in version 1.33.11

- centWave uses a new method to estimate local noise which improves detection of
  closely spaced peaks
- group.mzClust was failing when result had one peak

# xcms 1.32 and before

For more details and all changes before May 2012 please see the (now
discontinued) CHANGELOG in the source package (inst/ folder).

## CHANGED BEHAVIOUR since Version 1.32:

Other Changes since Version 1.32:
- improved mzData writing, now includes MSn spectra and less verbose.
- improved netCDF writing, but not yet good enough for AMDIS

## CHANGED BEHAVIOUR since Version 1.14:

- centWave may report a smaller set of peaks, due to a small bug
  in the ROI algorithm some features with mass deviation > ppm were retained.

Other Changes since Version 1.14:

- New method for grouping: an algorithm inspired by mzMine
  group(method="nearest") has been implemented. It is slower
  than group(method="density"). It can individually group
  close-eluting peaks of very similar mass

- New method for retention time correction:
  The retcor(method="obiwarp") algorithm operates on the raw data,
  and thus allows to correct runs without well-behaving
  peak groups, or without peak picking at all.

- fillPeaks(method="MSW") is now also available
  for direct infusion spectra. The findPeaks(method="MSW")
  now returns several intensities, and correctly reports
  mzmin and mzmax for peaks.

- centWave now uses dynamic memory allocation, needs much less memory,
  and these BUF related errors should be a thing of the past.

- centWave gains an optional argument "noise",
  which is useful for data that was centroided without any intensity threshold,
  centroids with intensity < "noise" are omitted from ROI detection

- the fillPeaks() methods now remember which was
  an observed, and which was a "filled" peak.

- For direct infusion spectra diffreport() now shows
  the raw peak shapes, and also indicated "real" and "filled" peaks.

- xcmsRaw can now filter for positive/negative spectra,
  if the file includes both polarities. xcmsSet() can pass
  the polarity to contain positive/negative peaks only.
