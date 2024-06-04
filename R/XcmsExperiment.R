#' @title Next Generation `xcms` Result Object
#'
#' @aliases XcmsExperiment-class show,XcmsExperiment-method filterChromPeaks
#'
#' @description
#'
#' The `XcmsExperiment` is a data container for `xcms` preprocessing results
#' (i.e. results from chromatographic peak detection, alignment and
#' correspondence analysis).
#'
#' It provides the same functionality than the [XCMSnExp] object, but uses the
#' more advanced and modern MS infrastructure provided by the `MsExperiment`
#' and `Spectra` Bioconductor packages. With this comes a higher flexibility on
#' how and where to store the data.
#'
#' Documentation of the various functions for `XcmsExperiment` objects are
#' grouped by topic and provided in the sections below.
#'
#' The default `xcms` workflow is to perform
#'
#' - chromatographic peak detection using [findChromPeaks()]
#'
#' - optionally refine identified chromatographic peaks using
#'   [refineChromPeaks()]
#'
#' - perform an alignment (retention time adjustment) using [adjustRtime()].
#'   Depending on the method used this requires to run a correspondence
#'   analysis first
#'
#' - perform a correspondence analysis using the [groupChromPeaks()] function
#'   to group chromatographic peaks across samples to define the LC-MS
#'   features.
#'
#' - optionally perform a gap-filling to *rescue* signal in samples in which
#'   no chromatographic peak was identified and hence a missing value would
#'   be reported. This can be performed using the [fillChromPeaks()] function.
#'
#'
#' @section Subsetting and filtering:
#'
#' - `[`: subset an `XcmsExperiment` by **sample** (parameter `i`). Subsetting
#'   will by default drop correspondence results (as subsetting by samples will
#'   obviously affect the feature definition) and alignment results (adjusted
#'   retention times) while identified chromatographic peaks (for the selected
#'   samples) will be retained. Which preprocessing results should be
#'   kept or dropped can also be configured with optional parameters
#'   `keepChromPeaks` (by default `TRUE`), `keepAdjustedRtime` (by default
#'   `FALSE`) and `keepFeatures` (by default `FALSE`).
#'
#' - `filterChromPeaks`: filter chromatographic peaks of an `XcmsExperiment`
#'   keeping only those specified with parameter `keep`. Returns the
#'   `XcmsExperiment` with the filtered data. Chromatographic peaks to
#'   retain can be specified either by providing their index in the
#'   `chromPeaks` matrix, their ID (rowname in `chromPeaks`) or with a
#'   `logical` vector with the same length than number of rows of
#'   `chromPeaks`. Assignment of chromatographic peaks are updated to
#'   eventually present feature definitions after filtering.
#'
#' - `filterFeatureDefinitions`: filter feature definitions of an
#'   `XcmsExperiment` keeping only those defined with parameter `features`,
#'   which can be a `logical` of length equal to the number of features,
#'   an `integer` with the index of the features in
#'   `featureDefinitions(object)` to keep or a `character` with the feature
#'   IDs (i.e. row names in `featureDefinitions(object)`).
#'
#' - `filterFile`: filter an `XcmsExperiment` (or `MsExperiment`) by *file*
#'   (sample). The index of the samples to which the data should be subsetted
#'   can be specified with parameter `file`. The sole purpose of this function
#'   is to provide backward compatibility with the `MSnbase` package. Wherever
#'   possible, the `[` function should be used instead for any sample-based
#'   subsetting. Parameters `keepChromPeaks`, `keepAdjustedRtime` and
#'   `keepChromPeaks` can be passed using `...`.
#'   Note also that in contrast to `[`, `filterFile` does not support subsetting
#'   in arbitrary order.
#'
#' - `filterIsolationWindow`: filter the **spectra** within an `MsExperiment`
#'   or `XcmsExperiment` object keeping only those with an isolation window
#'   containing the specified m/z (i.e., keeping spectra with an
#'   `"isolationWindowLowerMz"` smaller than the user-provided `mz` and an
#'   `"isolationWindowUpperMz"` larger than `mz`). For an `XcmsExperiment` also
#'   all chromatographic peaks (and subsequently also features) are removed for
#'   which the range of their `"isolationWindowLowerMz"` and
#'   `"isolationWindowUpperMz"` (columns in `chromPeakData`) do not contain
#'   the user provided `mz`.
#'
#' - `filterMsLevel`: filter the data of the `XcmsExperiment` or `MsExperiment`
#'   to keep only data of the MS level(s) specified with parameter `msLevel.`.
#'
#' - `filterMz`, `filterMzRange`: filter the spectra within an
#'   `XcmsExperiment` or `MsExperiment` to the specified m/z range (parameter
#'   `mz`). For `XcmsExperiment` also identified chromatographic peaks and
#'   features are filtered keeping only those that are within the specified
#'   m/z range (i.e. for which the m/z of the peak apex is within the m/z
#'   range). Parameter `msLevels.` allows to restrict the filtering to
#'   only specified MS levels. By default data from all MS levels are
#'   filtered.
#'
#' - `filterRt`: filter an `XcmsExperiment` keeping only data within the
#'   specified retention time range (parameter `rt`). This function will keep
#'   all preprocessing results present within the retention time range: all
#'   identified chromatographic peaks with the retention time of the apex
#'   position within the retention time range `rt` are retained along, if
#'   present, with the associated features.
#'   Parameter `msLevel.` is currently ignored, i.e. filtering will always
#'   performed on **all** MS levels of the object.
#'
#' @section Functionality related to chromatographic peaks:
#'
#' - `chromatogram`: extract chromatographic data from a data set. Parameters
#'   `mz` and `rt` allow to define specific m/z - retention time regions to
#'   extract the data from (to e.g. for extracted ion chromatograms EICs).
#'   Both parameters are expected to be numerical two-column matrices with
#'   the first column defining the lower and the second the upper margin.
#'   Each row can define a separate m/z - retention time region. Currently
#'   the function returns a [MChromatograms()] object for `object` being a
#'   `MsExperiment` or, for `object` being an `XcmsExperiment`, either a
#'   `MChromatograms` or [XChromatograms()] depending on parameter
#'   `return.type` (can be either `"MChromatograms"` or `"XChromatograms"`).
#'   For the latter also chromatographic peaks detected within the provided
#'   m/z and retention times are returned. Parameter `chromPeaks` allows
#'   to specify which chromatographic peaks should be reported. See
#'   documentation on the `chromPeaks` parameter for more information.
#'   If the `XcmsExperiment` contains correspondence results, also the
#'   associated feature definitions will be included in the returned
#'   `XChromatograms`. By default the function returns chromatograms from MS1
#'   data, but by setting parameter `msLevel = 2L` it is possible to e.g.
#'   extract also MS2 chromatograms. For `msLevel` other than 1 it is in
#'   addition important to also specify the `isolationWindowTargetMz` for which
#'   MS2 data should be extracted (e.g. for SWATH data MS2 spectra are created
#'   for different m/z isolation windows and the `isolationWindowTargetMz`
#'   parameter allows to define from which of these the MS2 chromatogram
#'   should be extracted.
#'   Note that in future more efficient data structures for chromatographic
#'   data will be available as well.
#'
#' - `chromPeaks`: returns a `numeric` matrix with the identified
#'   chromatographic peaks. Each row represents a chromatographic peak
#'   identified in one sample (file). The number of columns depends on the
#'   peak detection algorithm (see [findChromPeaks()]) but most methods return
#'   the following columns: `"mz"` (intensity-weighted mean of the m/z values
#'   of all mass peaks included in the chromatographic peak), `"mzmin"` (
#'   smallest m/z value of any mass peak in the chromatographic peak), `"mzmax"`
#'   (largest m/z value of any mass peak in the chromatographic peak), `"rt"`
#'   (retention time of the peak apex), `"rtmin"` (retention time of the first
#'   scan/mass peak of the chromatographic peak), `"rtmax"` (retention time of
#'   the last scan/mass peak of the chromatographic peak), `"into"` (integrated
#'   intensity of the chromatographic peak), `"maxo"` (maximal intensity of any
#'   mass peak of the chromatographic peak), `"sample"` (index of the sample
#'   in `object` in which the peak was identified). Parameters `rt`, `mz`,
#'   `ppm`, `msLevel` and `type` allow to extract subsets of identified
#'   chromatographic peaks from the `object`. See parameter description below
#'   for details.
#'
#' - `chromPeakData`: returns a `DataFrame` with potential additional
#'   *annotations* for the identified chromatographic peaks. Each row in this
#'   `DataFrame` corresponds to a row (same index and row name) in the
#'   `chromPeaks` matrix. The default *annotations* are `"ms_level"` (the MS
#'   level in which the peak was identified) and `"is_filled"` (whether the
#'   chromatographic peak was *detected* (by `findChromPeaks`) or *filled-in*
#'   (by `fillChromPeaks`).
#'
#' - `chromPeakSpectra`: extract MS spectra for identified chromatographic
#'   peaks. This can be either all (full scan) MS1 spectra with retention
#'   times between the retention time range of a chromatographic peak, all
#'   MS2 spectra (if present) with a retention time within the retention
#'   time range of a (MS1) chromatographic peak and a precursor m/z within
#'   the m/z range of the chromatographic peak or single, selected spectra
#'   depending on their total signal or highest signal. Parameter `msLevel`
#'   allows to define from which MS level spectra should be extracted,
#'   parameter `method` allows to define if all or selected spectra should
#'   be returned. See [chromPeakSpectra()] for details.
#'
#' - `dropChromPeaks`: removes (all) chromatographic peak detection results
#'   from `object`. This will also remove any correspondence results (i.e.
#'   features) and eventually present adjusted retention times from the object
#'   if the alignment was performed **after** the peak detection.
#'   Alignment results (adjusted retention times) can be retained if parameter
#'   `keepAdjustedRtime` is set to `TRUE`.
#'
#' - `dropFilledChromPeaks`: removes chromatographic peaks added by gap filling
#'   with `fillChromPeaks`.
#'
#' - `fillChromPeaks`: perform *gap filling* to integrate signal missing
#'   values in samples in which no chromatographic peak was found. This
#'   depends on correspondence results, hence `groupChromPeaks` needs to be
#'   called first. For details and options see [fillChromPeaks()].
#'
#' - `findChromPeaks`: perform chromatographic peak detection. See
#'   [findChromPeaks()] for details.
#'
#' - `hasChromPeaks`: whether the object contains peak detection results.
#'   Parameter `msLevel` allows to check whether peak detection results are
#'   available for the specified MS level(s).
#'
#' - `hasFilledChromPeaks`: whether gap-filling results (i.e., filled-in
#'   chromatographic peaks) are present.
#'
#' - `manualChromPeaks`: *manually* add chromatographic peaks by defining
#'   their m/z and retention time ranges. See [manualChromPeaks()] for
#'   details and examples.
#'
#' - `plotChromPeakImage`: show the *density* of identified chromatographic
#'   peaks per file along the retention time. See [plotChromPeakImage()] for
#'   details.
#'
#' - `plotChromPeaks`: indicate identified chromatographic peaks from one
#'   sample in the RT-m/z space. See [plotChromPeaks()] for details.
#'
#' - `plotPrecursorIons`: general visualization of precursor ions of
#'   LC-MS/MS data. See [plotPrecursorIons()] for details.
#'
#' - `refineChromPeaks`: *refines* identified chromatographic peaks in `object`.
#'   See [refineChromPeaks()] for details.
#'
#' @section Functionality related to alignment:
#'
#' - `adjustedRtime`: extract adjusted retention times. This is just an
#'   alias for `rtime(object, adjusted = TRUE)`.
#'
#' - `adjustRtime`: performs retention time adjustment (alignment) of the data.
#'   See [adjustRtime()] for details.
#'
#' - `applyAdjustedRtime`: replaces the original (raw) retention times with the
#'   adjusted ones. See [applyAdjustedRtime()] for more information.
#'
#' - `dropAdjustedRtime`: drops alignment results (adjusted retention time) from
#'   the result object. This also reverts the retention times of identified
#'   chromatographic peaks if present in the result object. Note that any
#'   results from a correspondence analysis (i.e. feature definitions) will be
#'   dropped too (if the correspondence analysis was performed **after** the
#'   alignment). This can be overruled with `keepAdjustedRtime = TRUE`.
#'
#' - `hasAdjustedRtime`: whether alignment was performed on the object (i.e.,
#'   the object contains alignment results).
#'
#' - `plotAdjustedRtime`: plot the alignment results; see [plotAdjustedRtime()]
#'   for more information.
#'
#' @section Functionality related to correspondence analysis:
#'
#' - `dropFeatureDefinitions`: removes any correspondence analysis results from
#'   `object` as well as any filled-in chromatographic peaks. By default
#'   (with parameter `keepAdjustedRtime = FALSE`) also all alignment results
#'   will be removed if alignment was performed **after** the correspondence
#'   analysis. This can be overruled with `keepAdjustedRtime = TRUE`.
#'
#' - `featureArea`: returns a `matrix` with columns `"mzmin"`, `"mzmax"`,
#'   `"rtmin"` and `"rtmax"` with the m/z and retention time range for each
#'   feature (row) in `object`. By default these represent the minimal m/z
#'   and retention times as well as maximal m/z and retention times for
#'   all chromatographic peaks assigned to that feature. Parameter
#'   `features` allows to extract these values for selected features only.
#'   Parameters `mzmin`, `mzmax`, `rtmin` and `rtmax` allow to define
#'   the function to calculate the reported `"mzmin"`, `"mzmax"`, `"rtmin"`
#'   and `"rtmax"` values.
#'
#' - `featureChromatograms`: extract ion chromatograms (EICs) for each
#'   feature in `object`. See [featureChromatograms()] for more details.
#'
#' - `featureDefinitions`: returns a `data.frame` with feature definitions or
#'   an empty `data.frame` if no correspondence analysis results are present.
#'   Parameters `msLevel`, `mz`, `ppm` and `rt` allow to define subsets of
#'   feature definitions that should be returned with the parameter `type`
#'   defining how these parameters should be used to subset the returned
#'   `data.frame`. See parameter descriptions for details.
#'
#' - `featureSpectra`: returns a [Spectra()] or `List` of `Spectra` with
#'   (MS1 or MS2) spectra associated to each feature. See [featureSpectra()]
#'   for more details and available parameters.
#'
#' - `featuresSummary`: calculate a simple summary on features. See
#'   [featureSummary()] for details.
#'
#' - `groupChromPeaks`: performs the correspondence analysis (i.e., grouping
#'   of chromatographic peaks into LC-MS *features*). See [groupChromPeaks()]
#'   for details.
#'
#' - `hasFeatures`: whether correspondence analysis results are presentin in
#'   `object`. The optional parameter `msLevel` allows to define the  MS
#'   level(s) for which it should be determined if feature definitions are
#'   available.
#'
#' - `overlappingFeatures`: identify features that overlapping or close in
#'   m/z - rt dimension. See [overlappingFeatures()] for more information.
#'
#' @section Extracting data and results from an `XcmsExperiment`:
#'
#' Preprocessing results can be extracted using the following functions:
#'
#' - `chromPeaks`: extract identified chromatographic peaks. See section on
#'   chromatographic peak detection for details.
#'
#' - `featureDefinitions`: extract the definition of *features* (chromatographic
#'   peaks grouped across samples). See section on correspondence analysis for
#'   details.
#'
#' - `featureValues`: extract a `matrix` of *values* for features from each
#'   sample (file). Rows are features, columns samples. Which *value* should be
#'   returned can be defined with parameter `value`, which can be any column of
#'   the `chromPeaks` matrix. By default (`value = "into"`) the integrated
#'   chromatographic peak intensities are returned. With parameter `msLevel` it
#'   is possible to extract values for features from certain MS levels.
#'   During correspondence analysis, more than one chromatographic peak per
#'   sample can be assigned to the same feature (e.g. if they are very close in
#'   retention time). Parameter `method` allows to define the strategy to deal
#'   with such cases: `method = "medret"`: report the value from the
#'   chromatographic peak with the apex position closest to the feautre's
#'   median retention time. `method = "maxint"`: report the value from the
#'   chromatographic peak with the largest signal (parameter `intensity` allows
#'   to define the column in `chromPeaks` that should be selected; defaults to
#'   `intensity = "into"). `method = "sum"`: sum the values for all
#'   chromatographic peaks assigned to the feature in the same sample.
#'
#' - `quantify`: extract the correspondence analysis results as a
#'   [SummarizedExperiment()]. The feature *values* are used as `assay` in
#'   the returned `SummarizedExperiment`, `rowData` contains the
#'   `featureDefinitions` (without column `"peakidx"`) and `colData` the
#'   `sampleData` of `object`. Additional parameters to the `featureValues`
#'   function (that is used to extract the feature value matrix) can be
#'   passed *via* `...`.
#'
#' @section Visualization:
#'
#' - `plot`: plot for each file the position of individual peaks in the m/z -
#'   retention time space (with color-coded intensity) and a base peak
#'   chromatogram. This function should ideally be called only on a data subset
#'   (i.e. after using `filterRt` and `filterMz` to restrict to a region of
#'   interest). Parameter `msLevel` allows to define from which MS level the
#'   plot should be created. If `x` is a `XcmsExperiment` with available
#'   identified chromatographic peaks, also the region defining the peaks
#'   are indicated with a rectangle. Parameter `peakCol` allows to define the
#'   color of the border for these rectangles.
#'
#' - `plotAdjustedRtime`: plot the alignment results; see [plotAdjustedRtime()]
#'   for more information.
#'
#' - `plotChromPeakImage`: show the *density* of identified chromatographic
#'   peaks per file along the retention time. See [plotChromPeakImage()] for
#'   details.
#'
#' - `plotChromPeaks`: indicate identified chromatographic peaks from one
#'   sample in the RT-m/z space. See [plotChromPeaks()] for details.
#'
#' @section General functionality and functions for backward compatibility:
#'
#' - `uniqueMsLevels`: returns the unique MS levels of the spectra in `object`.
#'
#' The functions listed below ensure compatibility with the *older*
#' [XCMSnExp()] xcms result object.
#'
#' - `fileNames`: returns the original data file names for the spectra data.
#'   Ideally, the `dataOrigin` or `dataStorage` spectra variables from the
#'   object's `spectra` should be used instead.
#'
#' - `fromFile`: returns the file (sample) index for each spectrum within
#'   `object`. Generally, subsetting by sample using the `[` is the preferred
#'    way to get spectra from a specific sample.
#'
#' - `polarity`: returns the polarity information for each spectrum in
#'   `object`.
#'
#' - `processHistory`: returns a `list` with [ProcessHistory] *process history*
#'   objects that contain also the parameter object used for the different
#'   processings. Optional parameter `type` allows to query for specific
#'   processing steps.
#'
#' - `rtime`: extract retention times of the **spectra** from the
#'   `MsExperiment` or `XcmsExperiment` object. It is thus a shortcut for
#'   `rtime(spectra(object))` which would be the preferred way to extract
#'   retention times from an `MsExperiment`. The `rtime` method for
#'   `XcmsExperiment` has an additional parameter `adjusted` which allows to
#'   define whether adjusted retention times (if present - `adjusted = TRUE`)
#'   or *raw* retention times (`adjusted = FALSE`) should be returned. By
#'   default adjusted retention times are returned if available.
#'
#' @section Differences compared to the [XCMSnExp()] object:
#'
#' - Subsetting by `[` supports arbitrary ordering.
#'
#' @param adjusted For `rtime,XcmsExperiment`: whether adjusted or *raw*
#'     retention times should be returned. The default is to return adjusted
#'     retention times, if available.
#'
#' @param aggregationFun For `chromatogram`: `character(1)` defining the
#'     function that should be used to *aggregate* intensities for retention
#'     time (i.e. each spectrum) along the specified m/z range (parameter
#'     `mz`). Defaults to `aggregationFun = "sum"` and hence all intensities
#'     will be summed up. Alternatively, use `aggregationFun = "max"` to use
#'     the maximal intensity per m/z range to create a base peak
#'     chromatogram (BPC).
#'
#' @param BPPARAM For `chromatogram`: parallel processing setup. Defaults
#'     to `BPPARAM = bpparam()`. See [bpparam()] for more information.
#'
#' @param chromPeaks For `chromatogram`: `character(1)` defining which
#'     chromatographic peaks should be returned. Can be either
#'     `chromPeaks = "apex_within"` (default) to return all chromatographic
#'     peaks with the m/z and RT of their apex within the m/z and retention
#'     time window, `chromPeaks = "any"` for all chromatographic peaks that
#'     are overlapping with the m/z - retention time window or
#'     `chromPeaks = "none"` to not include any chromatographic peaks. See
#'     also parameter `type` below for additional information.
#'
#' @param chunkSize For `chromatogram`: `integer(1)` defining the number of
#'     files from which the data should be loaded at a time into memory.
#'     Defaults to `chunkSize = 2L`.
#'
#' @param drop For `[`: ignored.
#'
#' @param features For `filterFeatureDefinitions` and `featureArea`: `logical`,
#'     `integer` or `character` defining the features to keep or from which
#'     to extract the feature area, respectively. See function description
#'     for more information.
#'
#' @param file For `filterFile`: `integer` with the indices of the samples
#'     (files) to which the data should be subsetted.
#'
#' @param filled For `featureValues`: `logical(1)` specifying whether values
#'     for filled-in peaks should be reported. For `filled = TRUE` (the
#'     default) filled peak values are returned, otherwise `NA` is reported
#'     for the respective features in the samples in which no peak was
#'     detected.
#'
#' @param i For `[`: `integer` or `logical` defining the samples/files to
#'     subset.
#'
#' @param include For `chromatogram`: deprecated; use parameter `chromPeaks`
#'      instead.
#'
#' @param intensity For `featureValues`: `character(1)` specifying the name
#'     of the column in the `chromPeaks(objects)` matrix containing the
#'     intensity value of the peak that should be used for the conflict
#'     resolution if `method = "maxint"`.
#'
#' @param isFilledColumn For `chromPeaks`: `logical(1)` whether a column
#'     `"is_filled"` should be included in the returned `matrix` with the
#'     information whether a peak was detected or *only* filled-in. Note that
#'     this information is also provided in the `chromPeakData` data frame.
#'
#' @param isolationWindowTargetMz For `chromatogram`: `numeric` (of length
#'     equal to the number of rows of `rt` and `mz`) with the isolation window
#'     target m/z of the MS2 spectra from which the chromatgrom should be
#'     generated. For MS1 data (`msLevel = 1L`, the default), this parameter
#'     is ignored. See examples on `chromatogram` below for further
#'     information.
#'
#' @param j For `[`: not supported.
#'
#' @param keep For `filterChromPeaks`: `logical`, `integer` or `character`
#'     specifying which chromatographic peaks to keep. If `logical` the
#'     length of `keep` needs to match the number of rows of `chromPeaks`.
#'     Alternatively, `keep` allows to specify the `index` (row) of peaks
#'     to keep or their ID (i.e. row name in `chromPeaks`).
#'
#' @param keepFeatures for most subsetting functions (`[`, `filterFile`):
#'     `logical(1)`: wheter eventually present feature definitions should
#'     be retained in the returned (filtered) object.
#'
#' @param keepAdjustedRtime `logical(1)`: whether adjusted retention times (if
#'     present) should be retained.
#'
#' @param method For `featureValues`: `character(1)` specifying the method to
#'     resolve multi-peak mappings within the same sample (correspondence
#'     analysis can assign more than one chromatographic peak within a sample
#'     to the same feature, e.g. if they are close in retention time). Options:
#'     `method = "medret"`: report the value for the chromatographic peak
#'     closest to the feature's median retention time.
#'     `method = "maxint"`: report the value for the chromatographic peak
#'     with the largest signal (parameter `intensity` allows to select the
#'     column in `chromPeaks` that should be used for *signal*).
#'     `method = "sum"`: sum the value for all chromatographic peaks in a
#'     sample assigned to the same feature. The default is `method = "medret"`.
#'     For `filterChromPeaks`: currently only `method = "keep"` is supported.
#'
#' @param missing For `featureValues`: default value for missing values.
#'     Allows to define the value that should be reported for a missing peak
#'     intensity. Defaults to `missing = NA_real_`.
#'
#' @param msLevel `integer` defining the MS level (or multiple MS level if the
#'     function supports it).
#'
#' @param msLevel. For `filterRt`: ignored. `filterRt` will always filter
#'     by retention times on all MS levels regardless of this parameter.
#'     For `chromatogram`: `integer` with the MS level from which the
#'     chromatogram(s) should be extracted. Has to be either of length 1 or
#'     length equal to the numer of rows of the parameters `mz` and `rt`
#'     defining the m/z and rt regions from which the chromatograms should
#'     be created. Defaults to `msLevel = 1L`.
#'     for `filterMsLevel`: `integer` defining the MS level(s) to which the
#'     data should be subset.
#'
#' @param mz For `chromPeaks` and `featureDefinitions`: `numeric(2)` optionally
#'     defining the m/z range for which chromatographic peaks or feature
#'     definitions should be returned. The full m/z range is used by default.
#'     For `chromatogram`: two-column numerical `matrix` with each row
#'     representing m/z range that should be aggregated into a chromatogram.
#'     If not provided the full m/z range of the data will be used (and hence
#'     a total ion chromatogram will be returned if `aggregationFun = "sum"`
#'     is used). For `filterIsolationWindow`: `numeric(1)` defining the m/z
#'     that should be contained within the spectra's isolation window.
#'
#' @param mzmax For `featureArea`: function to calculate the `"mzmax"` of
#'     a feature based on the `"mzmax"` values of the individual
#'     chromatographic peaks assigned to that feature. Defaults to
#'     `mzmax = max`.
#'
#' @param mzmin For `featureArea`: function to calculate the `"mzmin"` of
#'     a feature based on the `"mzmin"` values of the individual
#'     chromatographic peaks assigned to that feature. Defaults to
#'     `mzmin = min`.
#'
#' @param peakCol For `plot`: defines the border color of the rectangles
#'     indicating the identified chromatographic peaks. Only a single color
#'     is supported. Defaults to `peakCol = "#ff000060".
#'
#' @param ppm For `chromPeaks` and `featureDefinitions`: optional `numeric(1)`
#'     specifying the ppm by which the m/z range (defined by `mz` should be
#'     extended. For a value of `ppm = 10`, all peaks within `mz[1] - ppm / 1e6`
#'     and `mz[2] + ppm / 1e6` are returned.
#'
#' @param object An `XcmsExperiment` object.
#'
#' @param return.type For `chromPeakData`: `character(1)` defining the
#'     class of the returned object. Can be either `"DataFrame"` (the default)
#'     or `"data.frame"`. For `chromatogram`: `character(1)` defining the
#'     type of the returned object. Currently only
#'     `return.type = "MChromatograms"` is supported.
#'
#' @param rt For `chromPeaks` and `featureDefinitions`: `numeric(2)` defining
#'     the retention time range for which chromatographic peaks or features
#'     should be returned. The full range is used by default.
#'     For `chromatogram`: two column numerical `matrix` with each row
#'     representing the lower and upper retention time window(s) for the
#'     chromatograms. If not provided the full retention time range is used.
#'
#' @param rtmax For `featureArea`: function to calculate the `"rtmax"` of
#'     a feature based on the `"rtmax"` values of the individual
#'     chromatographic peaks assigned to that feature. Defaults to
#'     `rtmax = max`.
#'
#' @param rtmin For `featureArea`: function to calculate the `"rtmin"` of
#'     a feature based on the `"rtmin"` values of the individual
#'     chromatographic peaks assigned to that feature. Defaults to
#'     `rtmin = min`.
#'
#' @param type For `chromPeaks` and `featureDefinitions` and only if either
#'     `mz` and `rt` are defined too: `character(1)`: defining which peaks
#'     (or features) should be returned. For `type = "any"`: returns all
#'     chromatographic peaks or features also only partially overlapping any of
#'     the provided ranges. For `type = "within"`: returns only peaks or
#'     features completely within the region defined by `mz` and/or `rt`.
#'     For `type = "apex_within"`: returns peaks or features for which the m/z
#'     and retention time of the peak's apex is within the region defined by
#'     `mz` and/or `rt`.
#'     For `processHistory`: restrict returned processing steps to specific
#'     types. Use [processHistoryTypes()] to list all supported values.
#'
#' @param value For `featureValues`: `character(1)` defining which value should
#'     be reported for each feature in each sample. Can be any column of the
#'     `chromPeaks` matrix or `"index"` if simply the index of the assigned
#'     peak should be returned. Defaults to `value = "into"` thus the
#'     integrated peak area is reported.
#'
#' @param x An `XcmsExperiment` object.
#'
#' @param y For `plot`: should not be defined as it is not supported.
#'
#' @param ... Additional optional parameters. For `quantify`: any parameter
#'     for the `featureValues` call used to extract the feature value matrix.
#'
#' @name XcmsExperiment
#'
#' @author Johannes Rainer
#'
#' @exportClass XcmsExperiment
#'
#' @md
#'
#' @examples
#'
#' ## Creating a MsExperiment object representing the data from an LC-MS
#' ## experiment.
#' library(MsExperiment)
#'
#' ## Defining the raw data files
#' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'          system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'          system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#'
#' ## Defining a data frame with the sample characterization
#' df <- data.frame(mzML_file = basename(fls),
#'                 sample = c("ko15", "ko16", "ko18"))
#' ## Importing the data. This will initialize a `Spectra` object representing
#' ## the raw data and assign these to the individual samples.
#' mse <- readMsExperiment(spectraFiles = fls, sampleData = df)
#'
#' ## Extract a total ion chromatogram and base peak chromatogram
#' ## from the data
#' bpc <- chromatogram(mse, aggregationFun = "max")
#' tic <- chromatogram(mse)
#'
#' ## Plot them
#' par(mfrow = c(2, 1))
#' plot(bpc, main = "BPC")
#' plot(tic, main = "TIC")
#'
#' ## Extracting MS2 chromatographic data
#' ##
#' ## To show how MS2 chromatograms can be extracted we first load a DIA
#' ## (SWATH) data set.
#' mse_dia <- readMsExperiment(system.file("TripleTOF-SWATH",
#'     "PestMix1_SWATH.mzML", package = "msdata"))
#'
#' ## Extracting MS2 chromatogram requires also to specify the isolation
#' ## window from which to extract the data. Without that chromatograms
#' ## will be empty:
#' chr_ms2 <- chromatogram(mse_dia, msLevel = 2L)
#' intensity(chr_ms2[[1L]])
#'
#' ## First we list available isolation windows
#' table(isolationWindowTargetMz(spectra(mse_dia)))
#'
#' ## We can then extract the TIC of MS2 data for a specific isolation window
#' chr_ms2 <- chromatogram(mse_dia, msLevel = 2L,
#'     isolationWindowTargetMz = 244.05)
#' plot(chr_ms2)
#'
#' ####
#' ## Chromatographic peak detection
#'
#' ## Perform peak detection on the data using the centWave algorith. Note
#' ## that the parameters are chosen to reduce the run time of the example.
#' p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
#' xmse <- findChromPeaks(mse, param = p)
#' xmse
#'
#' ## Have a quick look at the identified chromatographic peaks
#' head(chromPeaks(xmse))
#'
#' ## Extract chromatographic peaks identified between 3000 and 3300 seconds
#' chromPeaks(xmse, rt = c(3000, 3300), type = "within")
#'
#' ## Extract ion chromatograms (EIC) for the first two chromatographic
#' ## peaks.
#' chrs <- chromatogram(xmse,
#'     mz = chromPeaks(xmse)[1:2, c("mzmin", "mzmax")],
#'     rt = chromPeaks(xmse)[1:2, c("rtmin", "rtmax")])
#'
#' ## An EIC for each sample and each of the two regions was extracted.
#' ## Identified chromatographic peaks in the defined regions are extracted
#' ## as well.
#' chrs
#'
#' ## Plot the EICs for the second defined region
#' plot(chrs[2, ])
#'
#' ## Subsetting the data to the results (and data) for the second sample
#' a <- xmse[2]
#' nrow(chromPeaks(xmse))
#' nrow(chromPeaks(a))
#'
#' ## Filtering the result by retention time: keeping all spectra and
#' ## chromatographic peaks within 3000 and 3500 seconds.
#' xmse_sub <- filterRt(xmse, rt = c(3000, 3500))
#' xmse_sub
#' nrow(chromPeaks(xmse_sub))
#'
#' ## Perform an initial feature grouping to allow alignment using the
#' ## peak groups method:
#' pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
#' xmse <- groupChromPeaks(xmse, param = pdp)
#'
#' ## Perform alignment using the peak groups method.
#' pgp <- PeakGroupsParam(span = 0.4)
#' xmse <- adjustRtime(xmse, param = pgp)
#'
#' ## Visualizing the alignment results
#' plotAdjustedRtime(xmse)
#'
#' ## Performing the final correspondence analysis
#' xmse <- groupChromPeaks(xmse, param = pdp)
#'
#' ## Show the definition of the first 6 features
#' featureDefinitions(xmse) |> head()
#'
#' ## Extract the feature values; show the results for the first 6 rows.
#' featureValues(xmse) |> head()
#'
#' ## The full results can also be extracted as a `SummarizedExperiment`
#' ## that would eventually simplify subsequent analyses with other packages.
#' ## Any additional parameters passed to the function are passed to the
#' ## `featureValues` function that is called to generate the feature value
#' ## matrix.
#' se <- quantify(xmse, method = "sum")
#'
#' ## EICs for all features can be extracted with the `featureChromatograms`
#' ## function. Note that, depending on the data set, extracting this for
#' ## all features might take some time. Below we extract EICs for the
#' ## first 10 features by providing the feature IDs.
#' chrs <- featureChromatograms(xmse,
#'     features = rownames(featureDefinitions(xmse))[1:10])
#' chrs
#'
#' plot(chrs[3, ])
NULL

.empty_chrom_peaks <- function(sample = TRUE) {
    cols <- c(.REQ_PEAKS_COLS, "maxo","sn")
    if (!sample)
        cols <- cols[cols != "sample"]
    matrix(numeric(), ncol = length(cols), nrow = 0,
           dimnames = list(character(), cols))
}

.empty_feature_definitions <- function() {
    res <- data.frame(mzmed = numeric(), mzmin = numeric(), mzmax = numeric(),
                      rtmed = numeric(), rtmin = numeric(), rtmax = numeric(),
                      peakidx = list(), ms_level = integer())
    res$peakidx <- list()
    res
}

setClass("XcmsExperiment",
         contains = "MsExperiment",
         slots = c(chromPeaks = "matrix",
                   chromPeakData = "data.frame",
                   featureDefinitions = "data.frame",
                   processHistory = "list"),
         prototype = prototype(
             chromPeaks = .empty_chrom_peaks(),
             chromPeakData = data.frame(ms_level = integer(),
                                        is_filled = logical()),
             featureDefinitions = .empty_feature_definitions(),
             processHistory = list()))

setClassUnion("XcmsResult", c("XcmsExperiment", "XCMSnExp"))

setValidity("XcmsExperiment", function(object) {
    msg <- c(.mse_valid_chrom_peaks(object@chromPeaks),
             .mse_valid_chrom_peak_data(object@chromPeakData),
             .mse_same_rownames(object@chromPeaks, object@chromPeakData),
             .mse_valid_feature_def(object@featureDefinitions))
    if (is.null(msg)) TRUE
    else msg
})

setMethod("show", "XcmsExperiment", function(object) {
    callNextMethod()
    cat(" xcms results:\n")
    if (hasChromPeaks(object))
        cat("  - chromatographic peaks:", nrow(object@chromPeaks),
            "in MS level(s):",
            paste(unique(object@chromPeakData$ms_level), collapse = ", "), "\n")
    if (hasAdjustedRtime(object))
        cat("  - adjusted retention times: mean absolute difference",
            format(mean(abs(rtime(spectra(object)) -
                           spectra(object)$rtime_adjusted)),
                  digits = 3), "seconds\n")
    if (hasFeatures(object))
        cat("  - correspondence results:", nrow(object@featureDefinitions),
            "features in MS level(s):",
            paste(unique(object@featureDefinitions$ms_level), collapse = ", "),
            "\n")
})

################################################################################
## Filtering and subsetting
################################################################################

#' @rdname XcmsExperiment
setMethod("[", "XcmsExperiment", function(x, i, j, ...) {
    if (!missing(j))
        stop("subsetting by j not supported")
    .subset_xcms_experiment(x, i = i, ...)
})

#' @rdname XcmsExperiment
setMethod(
    "filterIsolationWindow", "XcmsExperiment",
    function(object, mz = numeric()) {
        if (length(mz) > 1L)
            mz <- mz[1L]
        object <- filterSpectra(object, filterIsolationWindow, mz = mz)
        if (hasChromPeaks(object) && length(mz) &&
            all(c("isolationWindowLowerMz", "isolationWindowUpperMz") %in%
                colnames(object@chromPeakData))) {
            idx <- which(object@chromPeakData$isolationWindowLowerMz < mz &
                         object@chromPeakData$isolationWindowUpperMz > mz)
            object <- .filter_chrom_peaks(object, idx)
        }
        object
    })

#' @rdname XcmsExperiment
setMethod(
    "filterRt", "XcmsExperiment",
    function(object, rt, msLevel.) {
        if (missing(rt))
            return(object)
        rt <- range(rt)
        if (!missing(msLevel.))
            warning("Parameter 'msLevel.' currently ignored.", call. = FALSE)
        msLevel. <- uniqueMsLevels(object)
        if (hasChromPeaks(object))
            object <- .filter_chrom_peaks(
                object, base::which(between(.chromPeaks(object)[, "rt"], rt)))
        callNextMethod(object = object, rt = rt, msLevel. = msLevel.)
    })

#' @rdname XcmsExperiment
setMethod(
    "filterMzRange", "XcmsExperiment",
    function(object, mz = numeric(), msLevel. = uniqueMsLevels(object)) {
        if (missing(mz) || !length(mz))
            return(object)
        mz <- range(mz)
        if (hasChromPeaks(object)) {
            keep <- between(.chromPeaks(object)[, "mz"], mz)
            keep <- keep | (!.chromPeakData(object)$ms_level %in% msLevel.)
            object <- .filter_chrom_peaks(object, idx = base::which(keep))
        }
        callNextMethod(object = object, mz = mz, msLevel. = msLevel.)
    })

#' @rdname XcmsExperiment
setMethod(
    "filterMsLevel", "XcmsExperiment",
    function(object, msLevel. = uniqueMsLevels(object)) {
        if (!length(msLevel.))
            return(object)
        if (hasChromPeaks(object)) {
            keep <- .chromPeakData(object)$ms_level %in% msLevel.
            object <- .filter_chrom_peaks(object, idx = base::which(keep))
        }
        callNextMethod(object = object, msLevel. = msLevel.)
    })


################################################################################
## chromatographic peaks
################################################################################

#' @rdname findChromPeaks
setMethod(
    "findChromPeaks",
    signature(object = "MsExperiment", param = "Param"),
    function(object, param, msLevel = 1L, chunkSize = 2L, ...,
             BPPARAM = bpparam()) {
        if (length(msLevel) > 1)
            stop("Currently only peak detection in a single MS level is ",
                 "supported", call. = FALSE)
        if (chunkSize < 0) {
            res <- .mse_find_chrom_peaks(
                object, msLevel = msLevel, param = param,
                BPPARAM = BPPARAM)
        } else {
            res <- .mse_find_chrom_peaks_chunks(
                object, msLevel = msLevel, param = param,
                chunkSize = chunkSize, BPPARAM = BPPARAM)
        }
        ## Assign/define peak IDs.
        pkd <- data.frame(ms_level = rep(as.integer(msLevel), nrow(res)),
                          is_filled = rep(FALSE, nrow(res)))
        ph <- XProcessHistory(param = param,
                              type. = .PROCSTEP.PEAK.DETECTION,
                              fileIndex. = seq_along(object),
                              msLevel = msLevel)
        if (!is(object, "XcmsExperiment"))
            object <- as(object, "XcmsExperiment")
        object <- addProcessHistory(object, ph)
        object <- .mse_add_chrom_peaks(object, res, pkd)
        validObject(object)
        object
    })

#' @rdname findChromPeaks
setMethod(
    "findChromPeaks",
    signature(object = "XcmsExperiment", param = "Param"),
    function(object, param, msLevel = 1L, chunkSize = 2L, add = FALSE,
             ..., BPPARAM = bpparam()) {
        if (hasFeatures(object)) {
            message("Remove feature definitions")
            object <- dropFeatureDefinitions(
                object, keepAdjustedRtime = hasAdjustedRtime(object))
        }
        if (hasChromPeaks(object) && !add) {
            message("Remove previously identified chromatographic peaks")
            object <- dropChromPeaks(
                object, keepAdjustedRtime = hasAdjustedRtime(object))
        }
        callNextMethod()
    })

#' @rdname findChromPeaksIsolationWindow
setMethod(
    "findChromPeaksIsolationWindow", "MsExperiment",
    function(object, param, msLevel = 2L,
             isolationWindow = isolationWindowTargetMz(spectra(object)),
             chunkSize = 2L, ..., BPPARAM = bpparam()) {
        if (length(isolationWindow) != length(spectra(object)))
            stop("Length of 'isolationWindow' has to match the ",
                 "number of spectra in 'object'")
        if (all(is.na(isolationWindow)))
            stop("No non-missing values in 'isolationWindow'")
        object@spectra$isolationWindow <- isolationWindow
        if (hasAdjustedRtime(object)) {
            rts <- spectra(object)$rtime
            object@spectra$rtime <- object@spectra$rtime_adjusted
        }
        res <- lapply(.mse_split_spectra_variable(as(object, "MsExperiment"),
                                                  isolationWindow),
                      FUN = findChromPeaks, param = param, msLevel = msLevel,
                      chunkSize = chunkSize, BPPARAM = BPPARAM)
        pks <- lapply(res, .chromPeaks)
        lns <- lengths(pks)
        if (any(lns > 0)) {
            pks <- do.call(rbind, pks[lns > 0])
            pkd <- do.call(rbind, lapply(res[lns > 0], function(z) {
                p <- .chromPeakData(z)
                s <- z@spectra[1L]
                p$isolationWindow <- s$isolationWindow
                p$isolationWindowTargetMZ <- s$isolationWindowTargetMz
                p$isolationWindowLowerMz <- s$isolationWindowLowerMz
                p$isolationWindowUpperMz <- s$isolationWindowUpperMz
                p
            }))
        }
        xph <- XProcessHistory(param = param, date. = date(),
                               type. = .PROCSTEP.PEAK.DETECTION,
                               msLevel = msLevel)
        if (!is(object, "XcmsExperiment"))
            object <- as(object, "XcmsExperiment")
        object <- addProcessHistory(object, xph)
        if (hasAdjustedRtime(object))
            object@spectra$rtime <- rts
        object <- .mse_add_chrom_peaks(object, pks, pkd)
        validObject(object)
        object
    })

#' @rdname XcmsExperiment
setMethod("hasChromPeaks", "XcmsExperiment",
          function(object, msLevel = integer()) {
              if (length(msLevel))
                  any(object@chromPeakData$ms_level %in% msLevel)
              else as.logical(nrow(object@chromPeaks))
})

#' @rdname XcmsExperiment
setMethod(
    "dropChromPeaks", "XcmsExperiment",
    function(object, keepAdjustedRtime = FALSE) {
        if (hasChromPeaks(object)) {
            pt <- vapply(object@processHistory, processType, character(1))
            object@processHistory <- dropProcessHistoriesList(
                object@processHistory,
                type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                         .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                         .PROCSTEP.PEAK.REFINEMENT))
            object@chromPeaks <- .empty_chrom_peaks()
            object@chromPeakData <- data.frame(ms_level = integer(),
                                               is_filled = logical())
            if (hasAdjustedRtime(object) && !keepAdjustedRtime) {
                ## remove if alignment performed AFTER chrom peaks
                nom <- length(pt) + 1L
                idx_cp <- .match_last(.PROCSTEP.PEAK.DETECTION, pt,
                                      nomatch = nom)
                idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, pt,
                                      nomatch = nom)
                if (idx_al > idx_cp)
                    object <- dropAdjustedRtime(object)
            }
            if (hasFeatures(object))
                object <- dropFeatureDefinitions(
                    object, keepAdjustedRtime = keepAdjustedRtime)
        }
        object
    })

#' @rdname XcmsExperiment
setReplaceMethod("chromPeaks", "XcmsExperiment", function(object, value) {
    object@chromPeaks <- value
    object
})

#' @rdname XcmsExperiment
setMethod(
    "chromPeaks", "XcmsExperiment",
    function(object, rt = numeric(), mz = numeric(), ppm = 0,
             msLevel = integer(), type = c("any", "within", "apex_within"),
             isFilledColumn = FALSE) {
        type <- match.arg(type)
        pks <- object@chromPeaks
        if (isFilledColumn)
            pks <- cbind(
                pks, is_filled = as.numeric(object@chromPeakData$is_filled))
        pks[.index_chrom_peaks(object, rt = rt, mz = mz, ppm = ppm,
                               msLevel = msLevel, type = type), , drop = FALSE]
    })

#' @rdname XcmsExperiment
setReplaceMethod("chromPeakData", "XcmsExperiment", function(object, value) {
    object@chromPeakData <- value
    object
})

#' @rdname XcmsExperiment
setMethod(
    "chromPeakData", "XcmsExperiment",
    function(object, msLevel = integer(),
             return.type = c("DataFrame", "data.frame")) {
        return.type <- match.arg(return.type)
        if (return.type == "DataFrame")
            as(.chromPeakData(object, msLevel = msLevel), "DataFrame")
        else .chromPeakData(object, msLevel = msLevel)
    })

#' @rdname refineChromPeaks
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperiment", param = "CleanPeaksParam"),
    function(object, param = CleanPeaksParam(), msLevel = 1L) {
        if (!hasChromPeaks(object, msLevel = msLevel)) {
            warning("No chromatographic peaks for MS level ",
                    msLevel, " present", call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Removing feature definitions")
            object <- dropFeatureDefinitions(object)
        }
        validObject(param)
        rtw <- .chromPeaks(object)[, "rtmax"] - .chromPeaks(object)[, "rtmin"]
        keep_ms <- object@chromPeakData$ms_level %in% msLevel
        keep_rt <- rtw < param@maxPeakwidth & keep_ms
        keep <- which(keep_rt | !keep_ms)
        message("Removed ", nrow(.chromPeaks(object)) - length(keep), " of ",
                nrow(.chromPeaks(object)), " chromatographic peaks.")
        object@chromPeaks <- object@chromPeaks[keep, , drop = FALSE]
        object@chromPeakData <- object@chromPeakData[keep, , drop = FALSE]
        xph <- XProcessHistory(param = param, date. = date(),
                               type. = .PROCSTEP.PEAK.REFINEMENT,
                               fileIndex = seq_along(object),
                               msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

#' @rdname refineChromPeaks
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperiment", param = "MergeNeighboringPeaksParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        if (!hasChromPeaks(object, msLevel = msLevel)) {
            warning("No chromatographic peaks for MS level ",
                    msLevel, " present", call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Removing feature definitions")
            object <- dropFeatureDefinitions(object)
        }
        npks_orig <- nrow(.chromPeaks(object))
        validObject(param)
        res <- .xmse_apply_chunks(
            object, .xmse_merge_neighboring_peaks, msLevel = msLevel,
            expandRt = param@expandRt, expandMz = param@expandMz,
            ppm = param@ppm, minProp = param@minProp, BPPARAM = BPPARAM,
            keepAdjustedRtime = TRUE, ignoreHistory = TRUE,
            keepSampleIndex = FALSE, chunkSize = chunkSize)
        pks <- do.call(rbind, lapply(res, `[[`, 1L))
        pkd <- do.call(rbind.data.frame, c(lapply(res, `[[`, 2L),
                                           make.row.names = FALSE))
        npks <- unlist(lapply(res, `[[`, 3L), use.names = FALSE)
        pks[, "sample"] <- rep(seq_along(npks), npks)
        nas <- which(is.na(rownames(pks))) # merged peaks
        if (!any(colnames(pkd) == "merged"))
            pkd$merged <- FALSE
        pkd$merged[nas] <- TRUE
        ## Fix rownames
        maxi <- max(as.integer(sub("CP", "", rownames(object@chromPeaks))))
        rownames(pks)[nas] <- .featureIDs(length(nas), "CP", from = maxi + 1L)
        rownames(pkd) <- rownames(pks)
        ## Merge with existing peaks from **other** MS levels
        keep <- object@chromPeakData$ms_level != msLevel
        if (any(keep)) {
            object@chromPeaks <- rbind(object@chromPeaks[keep, ], pks)
            object@chromPeakData <- rbindFill(object@chromPeakData[keep, ], pkd)
        } else {
            object@chromPeaks <- pks
            object@chromPeakData <- pkd
        }
        message("Reduced from ", npks_orig, " to ", nrow(.chromPeaks(object)),
                " chromatographic peaks.")
        xph <- XProcessHistory(param = param, date. = date(),
                               type. = .PROCSTEP.PEAK.REFINEMENT,
                               fileIndex = seq_along(object),
                               msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

#' @rdname refineChromPeaks
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperiment", param = "FilterIntensityParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        if (!hasChromPeaks(object, msLevel = msLevel)) {
            warning("No chromatographic peaks for MS level ",
                    msLevel, " present", call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Removing feature definitions")
            object <- dropFeatureDefinitions(object)
        }
        npks_orig <- nrow(.chromPeaks(object))
        validObject(param)
        if (param@nValues == 1L) {
            if (!any(colnames(.chromPeaks(object)) == param@value))
                stop("Column '", param@value, "' not available.")
            keep <- .chromPeaks(object)[, param@value] >= param@threshold |
                .chromPeakData(object)$ms_level != msLevel
        } else
            keep <- unlist(.xmse_apply_chunks(
                object, .xmse_filter_peaks_intensities,nValues = param@nValues,
                threshold = param@threshold, msLevel = msLevel,
                keepAdjustedRtime = TRUE, ignoreHistory = TRUE,
                BPPARAM = BPPARAM, chunkSize = chunkSize), use.names = FALSE)
        object@chromPeaks <- object@chromPeaks[keep, , drop = FALSE]
        object@chromPeakData <- object@chromPeakData[keep, ]
        message("Reduced from ", npks_orig, " to ", nrow(.chromPeaks(object)),
                " chromatographic peaks.")
        xph <- XProcessHistory(param = param, date. = date(),
                               type. = .PROCSTEP.PEAK.REFINEMENT,
                               fileIndex = seq_along(object),
                               msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

#' @rdname manualChromPeaks
setMethod("manualChromPeaks", "MsExperiment",
          function(object, chromPeaks = matrix(numeric()),
                   samples = seq_along(object), msLevel = 1L,
                   chunkSize = 2L, BPPARAM = bpparam()) {
              manualChromPeaks(as(object, "XcmsExperiment"), chromPeaks,
                               samples, msLevel, chunkSize, BPPARAM)
          })

#' @rdname manualChromPeaks
setMethod(
    "manualChromPeaks", "XcmsExperiment",
    function(object, chromPeaks = matrix(numeric()),
             samples = seq_along(object), msLevel = 1L,
             chunkSize = 2L, BPPARAM = bpparam()) {
        if (length(msLevel) > 1L)
            stop("Can only add peaks from one MS level at a time.")
        if (is.data.frame(chromPeaks)) chromPeaks <- as.matrix(chromPeaks)
        if (!nrow(chromPeaks)) return(object)
        if (!all(c("mzmin", "mzmax", "rtmin", "rtmax") %in%
                 colnames(chromPeaks)))
            stop("'chromPeaks' lacks one or more of the required colums ",
                 "\"mzmin\", \"mzmax\", \"rtmin\" and \"rtmax\".")
        chromPeaks <- chromPeaks[, c("mzmin", "mzmax", "rtmin", "rtmax")]
        if (!all(samples %in% seq_along(object)))
            stop("'samples' out of bounds")
        if (hasFeatures(object))
            object <- dropFeatureDefinitions(object)
        pal <- lapply(samples, function(z) chromPeaks)
        names(pal) <- samples
        chunks <- split(samples, ceiling(seq_along(samples) / chunkSize))
        pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                               "total (:percent) in ",
                                               ":elapsed"),
                               total = length(chunks) + 1L, clear = FALSE)
        pb$tick(0)
        res <- lapply(chunks, function(z, ...) {
            pb$tick()
            .xmse_integrate_chrom_peaks(
                .subset_xcms_experiment(
                    object, i = z, keepAdjustedRtime = TRUE,
                    ignoreHistory = TRUE), pal = pal[as.character(z)],
                msLevel = msLevel, BPPARAM = BPPARAM)
        })
        res <- do.call(rbind, res)
        nr <- nrow(res)
        maxi <- max(
            0, as.integer(sub("CP", "", rownames(.chromPeaks(object)))))
        rownames(res) <- .featureIDs(nr, "CP", maxi + 1)
        pkd <- data.frame(ms_level = rep(as.integer(msLevel), nr),
                          is_filled = rep(FALSE, nr))
        rownames(pkd) <- rownames(res)
        pb$tick()
        object@chromPeakData <- rbindFill(object@chromPeakData, pkd)
        object@chromPeaks <- rbindFill(object@chromPeaks, res)
        validObject(object)
        object
    })


#' @rdname XcmsExperiment
setMethod(
    "filterChromPeaks", "XcmsExperiment",
    function(object, keep = rep(TRUE, nrow(.chromPeaks(object))),
             method = "keep", ...) {
        method <- match.arg(method)
        object <- switch(
            method,
            keep = {
                idx <- .i2index(keep, ids = rownames(.chromPeaks(object)),
                                name = "keep")
                .filter_chrom_peaks(object, idx)
            }
        )
        object
    })

#' @rdname chromPeakSpectra
setMethod(
    "chromPeakSpectra", "XcmsExperiment",
    function(object, method = c("all", "closest_rt", "closest_mz",
                                "largest_tic", "largest_bpi"),
             msLevel = 2L, expandRt = 0, expandMz = 0, ppm = 0,
             skipFilled = FALSE, peaks = character(),
             return.type = c("Spectra", "List"), BPPARAM = bpparam()) {
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        method <- match.arg(method)
        return.type <- match.arg(return.type)
        if (msLevel == 1L && method %in% c("closest_mz")) {
            warning("method = \"closest_mz\" is not supported for msLevel = 1.",
                    " Changing to method = \"all\".")
            method <- "all"
        }
        if (length(peaks))
            pkidx <- .i2index(peaks, rownames(.chromPeaks(object)), "peaks")
        else pkidx <- integer()
        res <- .mse_spectra_for_peaks(object, method, msLevel, expandRt,
                                      expandMz, ppm, skipFilled, pkidx,
                                      BPPARAM)
        if (!length(pkidx))
            peaks <- rownames(.chromPeaks(object))
        else peaks <- rownames(.chromPeaks(object))[pkidx]
        if (return.type == "Spectra")
            res <- res[as.matrix(findMatches(peaks, res$peak_id))[, 2L]]
        else
            as(split(res, factor(res$peak_id, levels = peaks)), "List")
    })

#' @rdname reconstructChromPeakSpectra
setMethod(
    "reconstructChromPeakSpectra", "XcmsExperiment",
    function(object, expandRt = 0, diffRt = 2, minCor = 0.8, intensity = "maxo",
             peakId = rownames(chromPeaks(object, msLevel = 1L)),
             BPPARAM = bpparam()) {
        if (!hasChromPeaks(object))
            stop("'object' does not contain chromatographic peaks")
        peakId <- intersect(peakId, rownames(chromPeaks(object, msLevel = 1L)))
        if (!length(peakId))
            stop("None of the IDs provided with 'peakId' match the ID of ",
                 "available chromatographic peaks")
        if (hasAdjustedRtime(object))
            object@spectra$rtime <- object@spectra$rtime_adjusted
        BPPARAM <- backendBpparam(object@spectra)
        ## If performance is bad we could alternatively split the object
        ## first and then apply on subsets. Unclear how to provide the fromFile
        res <- bplapply(seq_along(object), function(i) {
            .reconstruct_dia_ms2(
                .subset_xcms_experiment(
                    object, i, keepChromPeaks = TRUE, keepAdjustedRtime = TRUE,
                    keepFeatures = TRUE, keepSampleIndex = TRUE),
                expandRt = expandRt, diffRt = diffRt, minCor = minCor,
                column = intensity, peakId = peakId, fromFile = i)
        }, BPPARAM = BPPARAM)
        do.call(c, res)
    })

################################################################################
## alignment
################################################################################

#' @rdname adjustRtime
setMethod(
    "adjustRtime", signature(object = "MsExperiment", param = "ObiwarpParam"),
    function(object, param, chunkSize = 2L, BPPARAM = bpparam()) {
        msLevel <- 1L
        res <- .mse_obiwarp_chunks(
            object, param = param, chunkSize = chunkSize, BPPARAM = BPPARAM)
        ## Saving adjusted rtimes into $rtime_adjusted
        rt_adj <- rep(NA_real_, length(spectra(object)))
        rt_adj[object@sampleDataLinks[["spectra"]][, 2L]] <-
            unlist(res, use.names = FALSE)
        object@spectra$rtime_adjusted <- rt_adj
        if (!is(object, "XcmsExperiment"))
            object <- as(object, "XcmsExperiment")
        if (hasChromPeaks(object)) {
            fidx <- as.factor(fromFile(object))
            object@chromPeaks <- .applyRtAdjToChromPeaks(
                .chromPeaks(object),
                rtraw = split(rtime(object, adjusted = FALSE), fidx),
                rtadj = split(rt_adj, fidx))
        }
        ph <- XProcessHistory(param = param,
                              type. = .PROCSTEP.RTIME.CORRECTION,
                              fileIndex. = seq_along(object),
                              msLevel = msLevel)
        object <- addProcessHistory(object, ph)
        validObject(object)
        object
})

#' @rdname adjustRtime
setMethod(
    "adjustRtime", signature(object = "MsExperiment",
                             param = "PeakGroupsParam"),
    function(object, param, msLevel = 1L, ...) {
        if (!inherits(object, "XcmsExperiment"))
            object <- as(object, "XcmsExperiment")
        if (hasAdjustedRtime(object))
            stop("Alignment results already present. Please either remove ",
                 "them with 'dropAdjustedRtime' in order to perform an ",
                 "alternative, new, alignment, or use 'applyAdjustedRtime'",
                 " prior 'adjustRtime' to perform a second round of ",
                 "alignment.")
        if (any(msLevel != 1L))
            stop("Alignment is currently only supported for MS level 1")
        if (!nrow(peakGroupsMatrix(param))) {
            if (!hasFeatures(object))
                stop("No feature definitions present in 'object'. Please ",
                     "perform first a correspondence analysis using ",
                     "'groupChromPeaks'")
            peakGroupsMatrix(param) <- adjustRtimePeakGroups(
                object, param = param)
        }
        fidx <- as.factor(fromFile(object))
        rt_raw <- split(rtime(object), fidx)
        rt_adj <- .adjustRtime_peakGroupsMatrix(
            rt_raw, peakGroupsMatrix(param), smooth = smooth(param),
            span = span(param), family = family(param),
            subset = subset(param), subsetAdjust = subsetAdjust(param))
        pt <- vapply(object@processHistory, processType, character(1))
        idx_pg <- .match_last(.PROCSTEP.PEAK.GROUPING, pt, nomatch = -1L)
        if (idx_pg > 0)
            ph <- object@processHistory[idx_pg]
        else ph <- list()
        object <- dropFeatureDefinitions(object)
        object@spectra$rtime_adjusted <- unlist(rt_adj, use.names = FALSE)
        object@chromPeaks <- .applyRtAdjToChromPeaks(
            .chromPeaks(object), rtraw = rt_raw, rtadj = rt_adj)
        xph <- XProcessHistory(
            param = param, type. = .PROCSTEP.RTIME.CORRECTION,
            fileIndex = seq_along(object), msLevel = msLevel)
        object@processHistory <- c(object@processHistory, ph, list(xph))
        validObject(object)
        object
})

#'@rdname LamaParama
setMethod(
    "adjustRtime", signature(object = "XcmsExperiment", param = "LamaParama"),
    function(object, param, BPPARAM = bpparam(), ...) {
        if (!hasChromPeaks(object))
            stop("'object' needs to have detected chromPeaks. ",
                 "Run 'findChromPeaks()' first")
        if (hasAdjustedRtime(object))
            stop("Alignment results already present. Please either remove ",
                 "them with 'dropAdjustedRtime' in order to perform an ",
                 "alternative, new, alignment, or use 'applyAdjustedRtime'",
                 " prior 'adjustRtime' to perform a second round of ",
                 "alignment.")
        fidx <- as.factor(fromFile(object))
        rt_raw <- split(rtime(object), fidx)
        idx <- seq_along(object)

        # Check if user as ran matching lama vs chrompeaks beforehand
        if (length(param@rtMap) == 0)
            param <- matchLamasChromPeaks(object, param)
        rtMap <- param@rtMap
        if (length(rtMap) != length(object))
            stop("Mismatch between the number of files matched to lamas: ",
                 length(rtMap), " and files in the object: ", length(object))

        # Make model and adjust retention for each file
        rt_adj <- bpmapply(rtMap, rt_raw, idx, FUN = function(x, y, i, param) {
            if (nrow(x) >= 10) { # too strict ? Gam always throws error when less than that and loess does not work that well either.
                .adjust_rt_model(y, method = param@method,
                                 rt_map = x, span = param@span,
                                 resid_ratio = param@outlierTolerance,
                                 zero_weight = param@zeroWeight,
                                 bs = param@bs)
            } else {
                warning("Too few chrompeaks could be assigned to external",
                        " reference peaks (lamas) for sample ", i,
                        ". Skipping alignment for this sample.")
                y
            }
        }, SIMPLIFY = FALSE, BPPARAM = BPPARAM, MoreArgs = list(param = param))

        # post processing housekeeping steps
        pt <- vapply(object@processHistory, processType, character(1))
        idx_pg <- .match_last(.PROCSTEP.PEAK.GROUPING, pt,
                                     nomatch = -1L)
        if (idx_pg > 0)
            ph <- object@processHistory[idx_pg]
        else ph <- list()
        object <- dropFeatureDefinitions(object)
        object@spectra$rtime_adjusted <- unlist(rt_adj, use.names = FALSE)
        object@chromPeaks <-.applyRtAdjToChromPeaks(
            .chromPeaks(object), rtraw = rt_raw, rtadj = rt_adj)
        xph <- XProcessHistory(
            param = param, type. = .PROCSTEP.RTIME.CORRECTION,
            fileIndex = seq_along(object))
        object@processHistory <- c(object@processHistory, ph, list(xph))
        validObject(object)
        object
    })

#' @rdname XcmsExperiment
setMethod("dropAdjustedRtime", "XcmsExperiment", function(object) {
    if (!hasAdjustedRtime(object))
        return(object)
    ptype <- vapply(object@processHistory, processType, character(1))
    nom <- length(ptype) + 1L
    idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, ptype, nomatch = nom)
    idx_co <- .match_last(.PROCSTEP.PEAK.GROUPING, ptype, nomatch = nom)
    if (hasChromPeaks(object)) {
        fidx <- as.factor(fromFile(object))
        object@chromPeaks <- .applyRtAdjToChromPeaks(
            object@chromPeaks,
            rtraw = split(rtime(object, adjusted = TRUE), fidx),
            rtadj = split(rtime(object, adjusted = FALSE), fidx))
    }
    svs <- unique(c(spectraVariables(object@spectra), "mz", "intensity"))
    object@spectra <- selectSpectraVariables(
        object@spectra, svs[svs != "rtime_adjusted"])
    object@processHistory <- dropProcessHistoriesList(
        object@processHistory, type = .PROCSTEP.RTIME.CORRECTION, num = 1L)
    if (hasFeatures(object) && idx_co > idx_al) {
        warning("Had to remove feature definitions along with the adjusted ",
                "retention times because of the dependency between them.")
        object <- dropFeatureDefinitions(object)
    }
    object
})

#' @rdname XcmsExperiment
setMethod("hasAdjustedRtime", "MsExperiment", function(object) {
    if (length(spectra(object)))
        any(spectraVariables(spectra(object)) == "rtime_adjusted")
    else FALSE
})

#' @rdname XcmsExperiment
setMethod(
    "rtime", "XcmsExperiment",
    function(object, adjusted = hasAdjustedRtime(object)) {
        if (adjusted && hasAdjustedRtime(object))
            spectra(object)$rtime_adjusted
        else rtime(spectra(object))
    })

#' @rdname XcmsExperiment
setMethod("adjustedRtime", "XcmsExperiment", function(object) {
    rtime(object, adjusted = TRUE)
})

################################################################################
## correspondence
################################################################################

#' @rdname groupChromPeaks
setMethod(
    "groupChromPeaks",
    signature(object = "XcmsExperiment", param = "Param"),
    function(object, param, msLevel = 1L, add = FALSE) {
        msLevel <- unique(msLevel)
        if (length(msLevel) != 1)
            stop("Can only perform the correspondence analysis on one MS",
                 " level at a time. Please repeat for other MS levels ",
                 "with parameter `add = TRUE`.")
        if (!hasChromPeaks(object, msLevel))
            stop("No chromatographic peak for MS level ", msLevel,
                 " present. Please perform first a peak detection ",
                 "using the 'findChromPeaks' method.", call. = FALSE)
        if (hasFeatures(object) && !add)
            object <- dropFeatureDefinitions(object)
        cps <- chromPeaks(object, msLevel = msLevel)
        res <- .xmse_group_cpeaks(
            cps, param = param,
            index = match(rownames(cps), rownames(.chromPeaks(object))))
        if (!nrow(res))
            return(object)
        res$ms_level <- as.integer(msLevel)
        if (add && hasFeatures(object)) {
            sf <- max(
                as.integer(sub("FT", "", rownames(object@featureDefinitions))))
            rownames(res) <- .featureIDs(nrow(res), from = (sf + 1L))
            object@featureDefinitions <- rbindFill(
                object@featureDefinitions, res)
        } else {
            rownames(res) <- .featureIDs(nrow(res))
            object@featureDefinitions <- res
        }
        xph <- XProcessHistory(param = param, type. = .PROCSTEP.PEAK.GROUPING,
                               fileIndex = seq_along(object), msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

#' @rdname XcmsExperiment
setMethod(
    "hasFeatures", "XcmsExperiment",
    function(object, msLevel = integer()) {
        if (length(msLevel))
            any(object@featureDefinitions$ms_level %in% msLevel)
        else as.logical(nrow(object@featureDefinitions))
    })

#' @rdname XcmsExperiment
setReplaceMethod("featureDefinitions", "XcmsExperiment",
                 function(object, value) {
                     object@featureDefinitions <- value
                     object
                 })

#' @rdname XcmsExperiment
setMethod(
    "featureDefinitions", "XcmsExperiment",
    function(object, mz = numeric(), rt = numeric(), ppm = 0,
             type = c("any", "within", "apex_within"), msLevel = integer()) {
        if (length(msLevel)) {
            fdef <- object@featureDefinitions[
                               object@featureDefinitions$ms_level %in% msLevel,
                             , drop = FALSE]
        } else fdef <- object@featureDefinitions
        type <- match.arg(type)
        .subset_feature_definitions(fdef, mz = mz, rt = rt,
                                    ppm = ppm, type = type)
    })

#' @rdname XcmsExperiment
setMethod(
    "dropFeatureDefinitions", "XcmsExperiment",
    function(object, keepAdjustedRtime = FALSE) {
        if (!hasFeatures(object))
            return(object)
        ptype <- vapply(object@processHistory, processType, character(1))
        nom <- length(ptype) + 1L
        idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, ptype, nomatch = nom)
        idx_co <- .match_last(.PROCSTEP.PEAK.GROUPING, ptype, nomatch = nom)
        object@processHistory <- dropProcessHistoriesList(
            object@processHistory, type = .PROCSTEP.PEAK.GROUPING, num = 1L)
        object@featureDefinitions <- .empty_feature_definitions()
        if (.hasFilledPeaks(object)) {
            object <- .filter_chrom_peaks(
                object, which(!.chromPeakData(object)$is_filled))
            object@processHistory <- dropProcessHistoriesList(
                object@processHistory, type = .PROCSTEP.PEAK.FILLING)
        }
        if (!keepAdjustedRtime && hasAdjustedRtime(object) && idx_al > idx_co) {
            object <- dropAdjustedRtime(object)
        }
        object
    })

#' @rdname manualChromPeaks
setMethod(
    "manualFeatures", "XcmsExperiment",
    function(object, peakIdx = list(), msLevel = 1L) {
        if (!length(peakIdx))
            return(object)
        if (length(msLevel) > 1L)
            stop("Can only define features for one MS level at a time")
        if (!hasChromPeaks(object))
            stop("No chromatographic peaks present. ",
                 "Please run 'findChromPeaks' first.")
        res <- .manual_feature_definitions(.chromPeaks(object), peakIdx)
        res$ms_level <- as.integer(msLevel)
        if (hasFeatures(object)) {
            maxi <- max(as.integer(
                sub("FT", "", rownames(featureDefinitions(object)))))
            rownames(res) <- .featureIDs(nrow(res), from = maxi + 1)
            object@featureDefinitions <- rbindFill(
                object@featureDefinitions, res)
        } else {
            rownames(res) <- .featureIDs(nrow(res))
            object@featureDefinitions <- res
        }
        object
    })

#' @rdname chromPeakChromatograms
setMethod(
    "chromPeakChromatograms", "XcmsExperiment",
    function(object, expandRt = 0, expandMz = 0, aggregationFun = "max",
             peaks = character(),
             return.type = c("XChromatograms", "MChromatograms"),
             ..., progressbar = TRUE) {
        return.type <- match.arg(return.type)
        if (!hasChromPeaks(object))
            stop("'object' does not have any 'chromPeaks'. Please run ",
                 "'findChromPeaks' first.")
        pks <- .chromPeaks(object)
        pkd <- .chromPeakData(object)
        if (length(peaks)) {
            msg <- paste0("If provided, 'peaks' is expected to be a character ",
                          "vector with the IDs (row names) of the ",
                          "chromatographic peaks.")
            if (!is.character(peaks))
                stop(msg)
            if (!all(peaks %in% rownames(.chromPeaks(object))))
                stop("'peaks' don't match row names of 'chromPeaks'. ", msg)
            pks <- pks[peaks, , drop = FALSE]
            pkd <- pkd[peaks, ]
        }
        pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                               "total (:percent) in ",
                                               ":elapsed"),
                               total = (length(object) + 1), clear = FALSE)
        pb$tick(0)
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        ph <- object@processHistory
        object <- as(object, "MsExperiment")
        res <- lapply(seq_along(object), function(z) {
            idx <- which(pks[, "sample"] == z)
            if (length(idx)) {
                mzr <- pks[idx, c("mzmin", "mzmax"), drop = FALSE]
                rtr <- pks[idx, c("rtmin", "rtmax"), drop = FALSE]
                if (expandMz != 0) {
                    mzr[, 1] <- mzr[, 1] - expandMz
                    mzr[, 2] <- mzr[, 2] + expandMz
                }
                if (expandRt != 0) {
                    rtr[, 1] <- rtr[, 1] - expandRt
                    rtr[, 2] <- rtr[, 2] + expandRt
                }
                chrs <- .mse_chromatogram(
                    object[z], rt = rtr, mz = mzr,
                    aggregationFun = aggregationFun,
                    msLevel = pkd$ms_level[idx],
                    isolationWindow = pkd$isolationWindow[idx],
                    chunkSize = 1L, progressbar = FALSE,
                    BPPARAM = SerialParam())
                fData(chrs)$sample_index <- z # report sample index
                pb$tick()
                rownames(chrs) <- rownames(pks)[idx]
                rownames(fData(chrs)) <- rownames(chrs)
                chrs
            } else {
                pb$tick()
                NULL
            }
        })
        res <- as(do.call(c, res[lengths(res) > 0]), return.type)
        pData(res)[,] <- NA             # it's not from a single file.
        ## re-order the result - if needed.
        if (any(rownames(res) != rownames(pks)))
            res <- res[match(rownames(pks), rownames(res)), 1L]
        if (return.type == "XChromatograms") {
            idx <- seq_along(res)
            pks <- split.data.frame(pks, idx)
            pkd <- split.data.frame(pkd, idx)
            for (i in seq_along(res)) {
                tmp <- res@.Data[i, 1L][[1L]]
                slot(tmp, "chromPeaks", check = FALSE) <- pks[[i]]
                slot(tmp, "chromPeakData", check = FALSE) <-
                    as(pkd[[i]], "DataFrame")
                res@.Data[i, 1L][[1L]] <- tmp
            }
            res@.processHistory <- ph
        }
        pb$tick()
        res
    })

#' @rdname featureChromatograms
setMethod(
    "featureChromatograms", "XcmsExperiment",
    function(object, expandRt = 0, expandMz = 0, aggregationFun = "max",
             features = character(), return.type = "XChromatograms",
             chunkSize = 2L, mzmin = min, mzmax = max, rtmin = min,
             rtmax = max, ..., progressbar = TRUE, BPPARAM = bpparam()) {
        return.type <- match.arg(return.type)
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        area <- featureArea(object, mzmin = mzmin, mzmax = mzmax, rtmin = rtmin,
                            rtmax = rtmax, features = features)
        if (expandRt != 0) {
            area[, "rtmin"] <- area[, "rtmin"] - expandRt
            area[, "rtmax"] <- area[, "rtmax"] + expandRt
        }
        if (expandMz != 0) {
            area[, "mzmin"] <- area[, "mzmin"] - expandMz
            area[, "mzmax"] <- area[, "mzmax"] + expandMz
        }
        fts <- featureDefinitions(object)[rownames(area), ]
        chrs <- as(.mse_chromatogram(
            as(object, "MsExperiment"),
            rt = area[, c("rtmin", "rtmax"), drop = FALSE],
            mz = area[, c("mzmin", "mzmax"), drop = FALSE],
            aggregationFun = aggregationFun, msLevel = fts$ms_level,
            chunkSize = chunkSize, progressbar = progressbar,
            BPPARAM = BPPARAM), "XChromatograms")
        ## Populate with chrom peaks.
        nf <- nrow(fts)
        js <- seq_len(ncol(chrs))
        pks_empty <- .chromPeaks(object)[integer(), ]
        pkd_empty <- as(.chromPeakData(object)[integer(), ], "DataFrame")
        tmp <- chrs@.Data
        for (i in seq_len(nf)) {
            idx <- fts$peakidx[[i]]
            smpl <- .chromPeaks(object)[idx, "sample"]
            for (j in js) {
                keep <- smpl == j
                tmp_i <- tmp[i, j][[1L]]
                if (any(keep)) {
                    slot(tmp_i, "chromPeaks", check = FALSE) <-
                        .chromPeaks(object)[idx[keep], , drop = FALSE]
                    slot(tmp_i, "chromPeakData",
                         check = FALSE) <- as(
                        object@chromPeakData[idx[keep], ], "DataFrame")
                } else {
                    slot(tmp_i, "chromPeaks", check = FALSE) <- pks_empty
                    slot(tmp_i, "chromPeakData", check = FALSE) <- pkd_empty
                }
                tmp[i, j][[1L]] <- tmp_i
            }
        }
        chrs@.Data <- tmp
        ## Update peakidx in feature definitions.
        fts$row <- seq_len(nf)
        pkid_all <- rownames(.chromPeaks(object))
        pkid <- chromPeaks(chrs)[, c("row", "column")]
        pkid <- cbind(pkid, index = seq_len(nrow(pkid)))
        pkidl <- split.data.frame(pkid, pkid[, "row"])
        fts$peakidx <- lapply(fts$row, function(z) {
            unname(pkidl[[z]][pkid_all[fts$peakidx[[z]]], "index"])
        })
        colnames(chrs) <- basename(fileNames(object))
        rownames(chrs@phenoData) <- colnames(chrs)
        chrs@featureDefinitions <- DataFrame(fts)
        chrs@.processHistory <- object@processHistory
        chrs
    })

#' @rdname XcmsExperiment
setMethod(
    "filterFeatureDefinitions", "XcmsExperiment",
    function(object, features = integer()) {
        if (!length(features))
            return(object)
        if (!hasFeatures(object))
            stop("No feature definitions present! Please run ",
                 "'groupChromPeaks' first.")
        idx <- .i2index(features, ids = rownames(object@featureDefinitions),
                        name = "features")
        object@featureDefinitions <- object@featureDefinitions[idx, ]
        validObject(object)
        object
    })

#' @rdname featureSpectra
setMethod(
    "featureSpectra", "XcmsExperiment",
    function(object, msLevel = 2L, expandRt = 0, expandMz = 0, ppm = 0,
             skipFilled = FALSE, return.type = c("Spectra", "List"),
             features = character(), ...) {
        return.type <- match.arg(return.type)
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run ",
                 "'groupChromPeaks' first.")
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        features_all <- rownames(featureDefinitions(object))
        if (!length(features))
            features <- features_all
        findex <- .i2index(features, ids = features_all, name = "features")
        features <- features_all[findex]
        findex <- unique(findex)
        ufeatures <- features_all[findex]
        pindex <- unlist(featureDefinitions(object)$peakidx[findex],
                         use.names = FALSE)
        sps <- .mse_spectra_for_peaks(
            object, msLevel = msLevel, expandRt = expandRt,
            expandMz = expandMz, ppm = ppm, skipFilled = skipFilled,
            peaks = unique(pindex), ...)
        mtch <- as.matrix(
            findMatches(sps$peak_id, rownames(.chromPeaks(object))[pindex]))
        sps <- sps[mtch[, 1L]]
        fid <- rep(
            ufeatures, lengths(featureDefinitions(object)$peakidx[findex]))
        sps$feature_id <- fid[mtch[, 2L]]
        if (return.type == "List") {
            sps <- List(split(sps, f = factor(sps$feature_id,
                                              levels = ufeatures)))
            sps[features]
        } else sps
    })

################################################################################
## gap filling
################################################################################

#' @rdname XcmsExperiment
setMethod("hasFilledChromPeaks", "XcmsExperiment", function(object) {
    any(.chromPeakData(object)$is_filled)
})

#' @rdname fillChromPeaks
setMethod(
    "fillChromPeaks",
    signature(object = "XcmsExperiment", param = "ChromPeakAreaParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        if (length(msLevel) != 1)
            stop("Can only perform peak filling for one MS level at a time.")
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No feature definitions for MS level ", msLevel, " present.")
        ## Define region to integrate from for each file
        feature_ids <- rownames(featureDefinitions(object, msLevel = msLevel))
        fr <- .features_ms_region(object, mzmin = param@mzmin,
                                  mzmax = param@mzmax, rtmin = param@rtmin,
                                  rtmax = param@rtmax, features = feature_ids)
        fr <- cbind(
            fr, mzmed = featureDefinitions(object, msLevel = msLevel)$mzmed)
        fvals <- featureValues(object, value = "index", msLevel = msLevel)
        ## For each sample, keep features with some missing values.
        pal <- lapply(seq_len(ncol(fvals)), function(i) {
            fr[is.na(fvals[, i]), , drop = FALSE]
        })
        names(pal) <- seq_along(pal)
        ## Get integration function and other info.
        ph <- .xmse_process_history(object, .PROCSTEP.PEAK.DETECTION,
                                    msLevel = msLevel)
        fill_fun <- .history2fill_fun(ph)
        mzf <- "wMean"
        if (length(ph) && inherits(ph[[1L]], "XProcessHistory")) {
            prm <- ph[[1L]]@param
            if (any(slotNames(prm) == "mzCenterFun"))
                mzf <- prm@mzCenterFun
        } else
            prm <- MatchedFilterParam()
        mzf <- paste0("mzCenter.", gsub("mzCenter.", "", mzf, fixed = TRUE))
        ## Manual chunk processing because we have to split `object` and `pal`
        idx <- seq_along(object)
        chunks <- split(idx, ceiling(idx / chunkSize))
        pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                               "total (:percent) in ",
                                               ":elapsed"),
                               total = length(chunks) + 1L, clear = FALSE)
        pb$tick(0)
        res <- lapply(chunks, function(z, ...) {
            pb$tick()
            .xmse_integrate_chrom_peaks(
                .subset_xcms_experiment(object, i = z, keepAdjustedRtime = TRUE,
                                        ignoreHistory = TRUE),
                pal = pal[z], intFun = fill_fun, mzCenterFun = mzf,
                param = prm, BPPARAM = BPPARAM)
        })
        res <- do.call(rbind, res)
        ## Update feature definitions
        i_res <- seq((nrow(.chromPeaks(object)) + 1L), length.out = nrow(res))
        i_res <- split(i_res, rownames(res))
        i_ft <- match(names(i_res), rownames(featureDefinitions(object)))
        for (i in seq_along(i_res))
            object@featureDefinitions$peakidx[[i_ft[i]]] <-
                sort(c(object@featureDefinitions$peakidx[[i_ft[i]]],i_res[[i]]))
        ## Add results
        nr <- nrow(res)
        maxi <- max(as.integer(sub("CP", "", rownames(.chromPeaks(object)))))
        rownames(res) <- .featureIDs(nr, "CP", maxi + 1)
        cpd <- data.frame(ms_level = rep(as.integer(msLevel), nr),
                          is_filled = rep(TRUE, nr))
        rownames(cpd) <- rownames(res)
        object@chromPeaks <- rbind(object@chromPeaks, res)
        object@chromPeakData <- rbindFill(object@chromPeakData, cpd)
        pb$tick()
        ## Need to update the index in the featureDefinitions
        ph <- XProcessHistory(param = param,
                              date. = date(),
                              type. = .PROCSTEP.PEAK.FILLING,
                              fileIndex = seq_along(object),
                              msLevel = msLevel)
        object <- addProcessHistory(object, ph)
        validObject(object)
        object
    })

#' @rdname XcmsExperiment
setMethod("dropFilledChromPeaks", "XcmsExperiment", function(object) {
    if (!.hasFilledPeaks(object))
        return(object)
    keep_pks <- which(!.chromPeakData(object)$is_filled)
    object <- .filter_chrom_peaks(object, keep_pks)
    object@processHistory <- dropProcessHistoriesList(
        object@processHistory, type = .PROCSTEP.PEAK.FILLING)
                type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                         .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                         .PROCSTEP.PEAK.REFINEMENT)
    validObject(object)
    object
})

################################################################################
## results
################################################################################

#' @rdname XcmsExperiment
setMethod("quantify", "XcmsExperiment", function(object, ...) {
    if (!hasFeatures(object))
        stop("No correspondence results present. Pease run ",
             "'groupChromPeaks' first.")
    fd <- featureDefinitions(object)
    SummarizedExperiment(
        assays = list(raw = featureValues(object, ...)),
        rowData = fd[, colnames(fd) != "peakidx"],
        colData = sampleData(object),
        metadata = processHistory(object))
})

#' @rdname XcmsExperiment
setMethod(
    "featureValues", "XcmsExperiment",
    function(object, method = c("medret", "maxint", "sum"), value = "into",
             intensity = "into", filled = TRUE, missing = NA_real_,
             msLevel = integer()) {
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No feature definitions for MS level(s) ", msLevel," present.")
        method <- match.arg(method)
        if (method == "sum" && !(value %in% c("into", "maxo")))
            stop("method 'sum' is only allowed if value is set to 'into'",
                 " or 'maxo'")
        if (is.character(missing) && !(missing %in% c("rowmin_half")))
            stop("if 'missing' is not 'NA' or a numeric it should",
                 " be one of: \"rowmin_half\".")
        fNames <- basename(fileNames(object))
        pks <- .chromPeaks(object)
        ## issue #157: replace all values for filled-in peaks with NA
        if (!filled)
            pks[.chromPeakData(object)$is_filled, ] <- NA_real_
        .feature_values(
            pks = pks, fts = featureDefinitions(object, msLevel = msLevel),
            method = method, value = value, intensity = intensity,
            colnames = fNames, missing = missing)
    })

################################################################################
## utility and unsorted methods
################################################################################

#' @rdname XcmsExperiment
setMethod(
    "chromatogram", "XcmsExperiment",
    function(object, rt = matrix(nrow = 0, ncol = 2),
             mz = matrix(nrow = 0, ncol = 2), aggregationFun = "sum",
             msLevel = 1L, chunkSize = 2L, isolationWindowTargetMz = NULL,
             return.type = c("XChromatograms", "MChromatograms"),
             include = character(),
             chromPeaks = c("apex_within", "any", "none"),
             BPPARAM = bpparam()) {
        if (!is.matrix(rt)) rt <- matrix(rt, ncol = 2L)
        if (!is.matrix(mz)) mz <- matrix(mz, ncol = 2L)
        if (length(include)) {
            warning("Parameter 'include' is deprecated, please use ",
                    "'chromPeaks' instead")
            chromPeaks <- include
        }
        if (nrow(mz) && !nrow(rt))
            rt <- cbind(rep(-Inf, nrow(mz)), rep(Inf, nrow(mz)))
        if (nrow(rt) && !nrow(mz))
            mz <- cbind(rep(-Inf, nrow(rt)), rep(Inf, nrow(rt)))
        return.type <- match.arg(return.type)
        chromPeaks <- match.arg(chromPeaks)
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        .xmse_extract_chromatograms_old(
            object, rt = rt, mz = mz, aggregationFun = aggregationFun,
            msLevel = msLevel, isolationWindow = isolationWindowTargetMz,
            chunkSize = chunkSize, chromPeaks = chromPeaks,
            return.type = return.type, BPPARAM = BPPARAM)
    })

#' @rdname XcmsExperiment
setMethod("processHistory", "XcmsExperiment", function(object, type) {
    ph <- object@processHistory
    if (length(ph) && !missing(type))
        ph <- ph[vapply(ph, function(z) processType(z) %in% type, logical(1))]
    ph
})

setMethod("addProcessHistory", "XcmsExperiment", function(object, ph) {
    if (!inherits(ph, "ProcessHistory"))
        stop("Argument 'ph' has to be of type 'ProcessHistory' or a class ",
             "extending it!")
    object@processHistory[[(length(object@processHistory) + 1)]] <- ph
    object
})

#' @rdname XcmsExperiment
setMethod(
    "filterFile", "XcmsExperiment",
    function(object, file, keepAdjustedRtime = hasAdjustedRtime(object),
             keepFeatures = FALSE, ...) {
        if (missing(file)) return(object)
        object[i = sort(unique(file)), keepAdjustedRtime = keepAdjustedRtime,
               keepFeatures = keepFeatures, ...]
    })
