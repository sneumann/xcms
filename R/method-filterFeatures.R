#' @title Filtering of features based on conventional quality assessment
#'
#' @description
#'
#' When dealing with metabolomics results, it is often necessary to filter
#' features based on certain criteria. These criteria are typically derived
#' from statistical formulas applied to full rows of data, where each row
#' represents a feature and its abundance of signal in each samples.
#' The `filterFeatures` function filters features based on these conventional
#' quality assessment criteria. Multiple types of filtering are implemented and
#' can be defined by the `filter` argument.
#'
#' Supported `filter` arguments are:
#'
#' - [`RsdFilter`]: Calculates the relative standard deviation
#'  (i.e. coefficient of variation) in abundance for each feature in QC
#'  (Quality Control) samples and filters them in the input object according to
#'  a provided threshold.
#'
#' - [`DratioFilter`]: Computes the D-ratio or *dispersion ratio*, defined as
#'  the standard deviation in abundance for QC samples divided by the standard
#'  deviation for biological test samples, for each feature and filters them
#'  according to a provided threshold
#'
#' - [`PercentMissingFilter`]: Determines the percentage of missing values for
#'  each feature in the various sample groups and filters them according to a
#'  provided threshold.
#'
#' - [`BlankFlag`]: Identifies features where the mean abundance in test samples
#'  is lower than a specified multiple of the mean abundance of blank samples.
#'  This can be used to flag features that result from contamination in the
#'  solvent of the samples. A new column `possible_contaminants` is added to the
#'  `featureDefinitions` (`XcmsExperiment` object) or `rowData`
#'  (`SummarizedExperiment` object) reflecting this.
#'
#' For specific examples, see the help pages of the individual parameter classes
#' listed above.
#'
#' @param object `XcmsExperiment` or `SummarizedExperiment`. For an
#' `XcmsExperiment` object, the `featureValues(object)` will be evaluated, and
#' for `Summarizedesxperiment` the `assay(object, assay)`. The object will be
#' filtered.
#'
#' @param filter The parameter object selecting and configuring the type of
#' filtering. It can be one of the following classes: [`RsdFilter`],
#' [`DratioFilter`], [`PercentMissingFilter`] or [`BlankFlag`].
#'
#' @param assay For filtering of `SummarizedExperiment` objects only. Indicates
#' which assay the filtering will be based on. Note that the features for the
#' entire object will be removed, but the computations are performed on a single
#' assay. Default is 1, which means the first assay of the `object` will
#' be evaluated.
#'
#' @param ... Optional parameters. For `object` being an `XcmsExperiment`:
#' parameters for the [featureValues()] call.
#'
#' @name filterFeatures
#'
#' @author Philippine Louail
#'
#' @importFrom MetaboCoreUtils rowRsd rowDratio rowPercentMissing rowBlank
#'
#' @importMethodsFrom ProtGenerics filterFeatures
#'
#' @references
#'
#' Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR,
#' Dunn WB. Guidelines and considerations for the use of system suitability
#' and quality control samples in mass spectrometry assays applied in
#' untargeted clinical metabolomic studies. Metabolomics. 2018;14(6):72.
#' doi: 10.1007/s11306-018-1367-3. Epub 2018 May 18. PMID: 29805336;
#' PMCID: PMC5960010.
#'
#' @examples
#' ## See the vignettes for more detailed examples
#' library(MsExperiment)
#'
#' ## Load a test data set with features defined.
#' test_xcms <- loadXcmsData()
#' ## Set up parameter to filter based on coefficient of variation. By setting
#' ## the filter such as below, features that have a coefficient of variation
#' ## superior to 0.3 in QC samples will be removed from the object `test_xcms`
#' ## when calling the `filterFeatures` function.
#'
#' rsd_filter <- RsdFilter(threshold = 0.3,
#'                         qcIndex = sampleData(test_xcms)$sample_type == "QC")
#'
#' filtered_data_rsd <- filterFeatures(object = test_xcms, filter = rsd_filter)
#'
#' ## Set up parameter to filter based on D-ratio. By setting the filter such
#' ## as below, features that have a D-ratio computed based on their abundance
#' ## between QC and study samples superior to 0.5 will be removed from the
#' ## object `test_xcms`.
#'
#' dratio_filter <- DratioFilter(threshold = 0.5,
#'                  qcIndex = sampleData(test_xcms)$sample_type == "QC",
#'                  studyIndex = sampleData(test_xcms)$sample_type == "study")
#'
#' filtered_data_dratio <- filterFeatures(object = test_xcms,
#'                                        filter = dratio_filter)
#'
#' ## Set up parameter to filter based on the percent of missing data.
#' ## Parameter f should represent the sample group of samples, for which the
#' ## percentage of missing values will be evaluated. As the setting is defined
#' ## bellow, if a feature as less (or equal) to 30% missing values in one
#' ## sample group, it will be kept in the `test_xcms` object.
#'
#' missing_data_filter <- PercentMissingFilter(threshold = 30,
#'                                        f = sampleData(test_xcms)$sample_type)
#'
#' filtered_data_missing <- filterFeatures(object = test_xcms,
#'                                         filter = missing_data_filter)
#'
#' ## Set up parameter to flag possible contaminants based on blank samples'
#' ## abundance. By setting the filter such as below, features that have mean
#' ## abundance ratio between blank(here use study as an example) and QC
#' ## samples less than 2 will be marked as `TRUE` in an extra column named
#' ## `possible_contaminants` in the `featureDefinitions` table of the object
#' ## `test_xcms`.
#'
#' filter <- BlankFlag(threshold = 2,
#'                     qcIndex = sampleData(test_xcms)$sample_type == "QC",
#'                     blankIndex = sampleData(test_xcms)$sample_type == "study")
#' filtered_xmse <- filterFeatures(test_xcms, filter)
#'
#' @md
NULL

#' @title Filter features based on their coefficient of variation
#'
#' @name RsdFilter
#'
#' @export
#'
#' @family Filter features in xcms
#'
#' @description
#' The `RsdFilter` class and methods enable users to filter features from an
#' `XcmsExperiment` or `SummarizedExperiment` object based on their relative
#' standard deviation (coefficient of variation) for a specified threshold.
#'
#' This `filter` is part of the possible dispatch of the generic function
#' `filterFeatures`. Features *above* (`>`) the user-input threshold will be
#' removed from the entire dataset.
#'
#' @param threshold `numeric` value representing the threshold. Features with a
#' coefficient of variation *strictly higher* (`>`) than this will be removed
#' from the entire dataset.
#'
#' @param qcIndex `integer` (or `logical`) vector corresponding to the indices
#' of QC samples.
#'
#' @param na.rm `logical` indicates whether missing values (`NA`) should be
#' removed prior to the calculations.
#'
#' @param mad `logical` indicates whether the *Median Absolute Deviation* (MAD)
#' should be used instead of the standard deviation. This is suggested for
#' non-gaussian distributed data.
#'
#' @inheritParams filterFeatures
#'
#' @return For `RsdFilter`: a `RsdFilter` class. `filterFeatures` return the
#' input object minus the features that did not met the user input threshold.
#'
#' @note
#' It is assumed that the abundance values are in natural scale. Abundances in
#' log scale should be first transformed to natural scale before calculating
#' the RSD.
#'
#' @author Philippine Louail
#'
#' @importFrom MetaboCoreUtils rowRsd
#'
NULL

#' @noRd
setClass("RsdFilter",
         slots = c(threshold = "numeric",
                   qcIndex = "integer",
                   na.rm = "logical",
                   mad = "logical"),
         contains = "Param",
         prototype = prototype(
             threshold = numeric(),
             qcIndex = integer(),
             na.rm = logical(),
             mad = logical()),
         validity = function(object) {
             msg <- NULL
             if (length(object@threshold) != 1)
                 msg <- c("'threshold' has to be a numeric of length 1")
             msg
         })

#' @rdname RsdFilter
#'
#' @export
RsdFilter <- function(threshold = 0.3, qcIndex = integer(),
                      na.rm = TRUE, mad = FALSE) {
    if (is.logical(qcIndex))
        qcIndex <- which(qcIndex)
    if (is.numeric(qcIndex))
        qcIndex <- as.integer(qcIndex)
    new("RsdFilter", threshold = threshold, qcIndex = qcIndex, na.rm = na.rm,
        mad = mad)
}

#' @rdname RsdFilter
setMethod("filterFeatures",
          signature(object = "XcmsResult",
                    filter = "RsdFilter"),
          function(object, filter, ...){
            .check_index_range(filter@qcIndex, length(object), name = "qcIndex")
            vals <- featureValues(object, ...)[, filter@qcIndex]
            vals <- rowRsd(x = vals,  na.rm = filter@na.rm, mad = filter@mad)
            fts_idx <- which(vals <= filter@threshold)
            message(length(vals) - length(fts_idx), " features were removed")
            ph <- XProcessHistory(param = filter,
                                  date. = date(),
                                  type. = .PROCSTEP.FEATURE.FILTERING,
                                  fileIndex = seq_along(object))
            object <- addProcessHistory(object, ph)
            filterFeatureDefinitions(object, features = fts_idx)
          }
)

#' @rdname RsdFilter
setMethod("filterFeatures",
          signature(object = "SummarizedExperiment",
                    filter = "RsdFilter"),
          function(object, filter, assay = 1){
              .check_index_range(filter@qcIndex, ncol(object), name = "qcIndex")
              vals <- assay(object, assay)[, filter@qcIndex]
              vals <- rowRsd(vals, na.rm = filter@na.rm, mad = filter@mad)
              fts_idx <- which(vals <= filter@threshold)
              message(length(vals) - length(fts_idx), " features were removed")
              object[fts_idx]
          }
)

#' @title Filter features based on the dispersion ratio
#'
#' @name DratioFilter
#'
#' @export
#'
#' @family Filter features in xcms
#'
#' @description
#' The `DratioFilter` class and method enable users to filter features from an
#' `XcmsExperiment` or `SummarizedExperiment` object based on the D-ratio or
#' *dispersion ratio*. This is defined as the standard deviation for QC
#' samples divided by the standard deviation for biological test samples, for
#' each feature of the object (Broadhurst et al.).
#'
#' This `filter` is part of the possible dispatch of the generic function
#' `filterFeatures`.  Features *above* (`>`) the user-input threshold will be
#' removed from the entire dataset.
#'
#' @param threshold `numeric` value representing the threshold. Features with a
#' D-ratio *strictly higher* (`>`) than this will be removed from the entire
#' dataset.
#'
#' @param qcIndex `integer` (or `logical`) vector corresponding to the indices
#' of QC samples.
#'
#' @param studyIndex `integer` (or `logical`) vector corresponding of the
#' indices of study samples.
#'
#' @param na.rm `logical` Indicates whether missing values (`NA`) should be
#' removed prior to the calculations.
#'
#' @param mad `logical` Indicates whether the *Median Absolute Deviation*
#' (MAD) should be used instead of the standard deviation. This is suggested
#' for non-gaussian distributed data.
#'
#' @inheritParams filterFeatures
#'
#' @return For `DratioFilter`: a `DratioFilter` class. `filterFeatures` return
#' the input object minus the features that did not met the user input threshold
#'
#' @author Philippine Louail
#'
#' @importFrom MetaboCoreUtils rowDratio
#'
#' @references
#'
#' Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR,
#' Dunn WB. Guidelines and considerations for the use of system suitability
#' and quality control samples in mass spectrometry assays applied in
#' untargeted clinical metabolomic studies. Metabolomics. 2018;14(6):72.
#' doi: 10.1007/s11306-018-1367-3. Epub 2018 May 18. PMID: 29805336;
#' PMCID: PMC5960010.
#'
NULL

#' @noRd
setClass("DratioFilter",
         slots = c(threshold = "numeric",
                   qcIndex = "integer",
                   studyIndex = "integer",
                   na.rm = "logical",
                   mad = "logical"),
         contains = "Param",
         prototype = prototype(
             threshold = numeric(),
             qcIndex = integer(),
             studyIndex = integer(),
             na.rm = logical(),
             mad = logical()),
         validity = function(object) {
             msg <- NULL
             if (length(object@threshold) != 1)
                 msg <- c("'threshold' has to be a numeric of length 1")
             msg
         })

#' @rdname DratioFilter
#'
#' @export
DratioFilter <- function(threshold = 0.5,
                         qcIndex = integer(),
                         studyIndex = integer(),
                         na.rm = TRUE, mad = FALSE) {
    if (is.logical(qcIndex))
        qcIndex <- which(qcIndex)
    if(is.numeric(qcIndex))
        qcIndex <- as.integer(qcIndex)
    if (is.logical(studyIndex))
        studyIndex <- which(studyIndex)
    if(is.numeric(studyIndex))
        studyIndex <- as.integer(studyIndex)
    new("DratioFilter", threshold = threshold, qcIndex = qcIndex,
        studyIndex = studyIndex, na.rm = na.rm, mad = mad)
}

#' @rdname DratioFilter
setMethod("filterFeatures",
          signature(object = "XcmsResult",
                    filter = "DratioFilter"),
          function(object, filter, ...){
              .check_index_range(filter@studyIndex, length(object),
                                 name = "studyIndex")
              .check_index_range(filter@qcIndex, length(object), name = "qcIndex")
              x <- featureValues(object, ...)[, filter@studyIndex]
              y <- featureValues(object, ...)[, filter@qcIndex]
              vals <- rowDratio(x = x, y = y,
                                na.rm = filter@na.rm,
                                mad = filter@mad)
              fts_idx <- which(vals <= filter@threshold)
              message(length(vals) - length(fts_idx), " features were removed")
              ph <- XProcessHistory(param = filter,
                                    date. = date(),
                                    type. = .PROCSTEP.FEATURE.FILTERING,
                                    fileIndex = seq_along(object))
              object <- addProcessHistory(object, ph)
              filterFeatureDefinitions(object, features = fts_idx)
          }
)

#' @rdname DratioFilter
setMethod("filterFeatures",
          signature(object = "SummarizedExperiment",
                    filter = "DratioFilter"),
          function(object, filter, assay = 1){
              .check_index_range(filter@studyIndex, ncol(object),
                                 name = "studyIndex")
              .check_index_range(filter@qcIndex, ncol(object), name = "qcIndex")
              x <- assay(object, assay)[, filter@studyIndex]
              y <- assay(object, assay)[, filter@qcIndex]
              vals <- rowDratio(x = x, y = y,
                                na.rm = filter@na.rm,
                                mad = filter@mad)
              fts_idx <- which(vals <= filter@threshold)
              message(length(vals) - length(fts_idx), " features were removed")
              object[fts_idx]
          }
)

#' @title Filter features based on the percentage of missing data
#'
#' @name PercentMissingFilter
#'
#' @export
#'
#' @family Filter features in xcms
#'
#' @description
#' The `PercentMissingFilter` class and method enable users to filter features
#' from an `XcmsExperiment` or `SummarizedExperiment` object based on the
#' percentage (values from 1 to 100) of missing values for each features in
#' different sample groups and filters them according to a provided threshold.
#'
#' This `filter` is part of the possible dispatch of the generic function
#' `filterFeatures`. Features with a percentage of missing values *higher* (`>`)
#' than the user input threshold in all sample groups will be removed (i.e.
#' features for which the proportion of missing values is below (`<=`) the
#' threshold in at least one sample group will be retained).
#'
#' @param f `vector` of the same length as the `object`, specifying the sample
#' type for each sample in the dataset. The percentage of missing values per
#' feature will be computed within each of these sample groups. Parameter `f`,
#' if not already a `factor`, will be converted to one using the factor function.
#' Samples with an `NA` as their value in `f` will be excluded from calculation.
#'
#' @param threshold `numeric` percentage (between 0 and 100) of accepted missing
#' values for a feature in one sample group.
#'
#' @inheritParams filterFeatures
#'
#' @return For `PercentMissingFilter`: a `PercentMissingFilter` class.
#' `filterFeatures` return the input object minus the features that did not met
#' the user input threshold
#'
#' @author Philippine Louail
#'
#' @importFrom MetaboCoreUtils rowPercentMissing
#'
NULL

#' @noRd
setClass("PercentMissingFilter",
         slots = c(threshold = "numeric",
                   f = "factor"),
         contains = "Param",
         prototype = prototype(
             threshold = numeric(),
             f = factor()),
         validity = function(object) {
             msg <- NULL
             if (length(object@threshold) != 1)
                 msg <- c("'threshold' has to be a numeric of length 1")
             msg
         })

#' @rdname PercentMissingFilter
#'
#' @export
PercentMissingFilter <- function(threshold = 30, f = factor()) {
    if (!is.factor(f))
        f <- factor(f)
    new("PercentMissingFilter", threshold = threshold, f = f)
}

#' @rdname PercentMissingFilter
setMethod("filterFeatures",
          signature(object = "XcmsResult",
                    filter = "PercentMissingFilter"),
          function(object, filter, ...){
              if (length(filter@f) != length(object))
                  stop("'f' must be same lenght as object")
              fts_idx <- c()
              for (i in levels(filter@f)){
                  spl_idx <- which(filter@f == i)
                  vals <- rowPercentMissing(featureValues(object, ...)[, spl_idx])
                  fts_idx <- c(fts_idx, which(vals <= filter@threshold))
              }
              fts_idx <- order(unique(fts_idx))
              message(length(vals) - length(fts_idx), " features were removed")
              ph <- XProcessHistory(param = filter,
                                    date. = date(),
                                    type. = .PROCSTEP.FEATURE.FILTERING,
                                    fileIndex = seq_along(object))
              object <- addProcessHistory(object, ph)
              filterFeatureDefinitions(object, features = fts_idx)
          }
)

#' @rdname PercentMissingFilter
setMethod("filterFeatures",
          signature(object = "SummarizedExperiment",
                    filter = "PercentMissingFilter"),
          function(object, filter, assay = 1){
              if (length(filter@f) != ncol(object))
                  stop("'f' must be same lenght as object")
              fts_idx <- c()
              for (i in levels(filter@f)){
                  spl_idx <- which(filter@f == i)
                  vals <- rowPercentMissing(assay(object, assay)[, spl_idx])
                  fts_idx <- c(fts_idx, which(vals <= filter@threshold))
              }
              fts_idx <- order(unique(fts_idx))
              message(length(vals) - length(fts_idx), " features were removed")
              object[fts_idx]
          }
)

#' @title Flag features based on the intensity in blank samples
#'
#' @name BlankFlag
#'
#' @export
#'
#' @family Filter features in xcms
#'
#' @description
#' The `BlankFlag` class and method enable users to flag features of an
#' `XcmsExperiment` or `SummarizedExperiment` object based on the relationship
#' between the intensity of a feature in blanks compared to the intensity in the
#' samples.
#'
#' This class and method are part of the possible dispatch of the
#' generic function `filterFeatures`. Features *below* (`<`) the user-input
#' threshold will be flagged by calling the `filterFeatures` function. This
#' means that an extra column will be created in `featureDefinitions` or
#' `rowData` called `possible_contaminants` with a logical value for each
#' feature.
#'
#' @param threshold `numeric` indicates the minimum difference
#' required between the mean abundance of a feature in samples compared to the
#' mean abundance of the same feature in blanks for it to not be considered a
#' possible contaminant. For example, the default threshold of 2 signifies that
#' the mean abundance of the features in samples has to be at least twice the
#' mean abundance in blanks for it to not be flagged as a possible contaminant.
#'
#' @param blankIndex `integer` (or `logical`) vector corresponding to the
#' indices of blank samples.
#'
#' @param qcIndex `integer` (or `logical`) vector corresponding to the
#' indices of quality control (QC) samples.
#'
#' @param na.rm `logical` indicates whether missing values (`NA`) should be
#' removed prior to the calculations.
#'
#' @inheritParams filterFeatures
#'
#' @return For `BlankFlag`: a `BlankFlag` class. `filterFeatures` returns
#' the input object with an added column in the features metadata called
#' `possible_contaminants` with a logical value for each feature. This is added
#' to `featureDefinitions` for `XcmsExperiment` objects and `rowData` for
#' `SummarizedExperiment` objects.
#'
#' @author Philippine Louail
#'
#' @importFrom MetaboCoreUtils rowBlank
#'
NULL

#' @noRd
setClass("BlankFlag",
         slots = c(threshold = "numeric", blankIndex = "integer",
                   qcIndex = "integer", na.rm = "logical"),
         contains = "Param",
         prototype = prototype(
             threshold = numeric(),
             blankIndex = integer(),
             qcIndex = integer(),
             na.rm = logical()),
         validity = function(object) {
             msg <- NULL
             if (length(object@threshold) != 1)
                 msg <- c("'threshold' has to be a numeric of length 1")
             msg
         })

#' @rdname BlankFlag
#'
#' @export
BlankFlag <- function(threshold = 2, blankIndex = integer(),
                         qcIndex = integer(), na.rm = TRUE) {
    if (is.logical(blankIndex))
        blankIndex <- which(blankIndex)
    if(is.numeric(blankIndex))
        blankIndex <- as.integer(blankIndex)
    if (is.logical(qcIndex))
        qcIndex <- which(qcIndex)
    if(is.numeric(qcIndex))
        qcIndex <- as.integer(qcIndex)
    new("BlankFlag", threshold = threshold, blankIndex = blankIndex,
        qcIndex = qcIndex, na.rm = na.rm)
}

#' @rdname BlankFlag
setMethod("filterFeatures",
          signature(object = "XcmsResult",
                    filter = "BlankFlag"),
          function(object, filter, ...){
              .check_index_range(filter@blankIndex, length(object), name = "blankIndex")
              .check_index_range(filter@qcIndex, length(object), name = "qcIndex")
              x <- featureValues(object, ...)[, filter@qcIndex]
              y <- featureValues(object, ...)[, filter@blankIndex]
              vals <- rowBlank(x = x, y = y,
                               na.rm = filter@na.rm,
                               threshold = filter@threshold)
              message(sum(vals, na.rm = TRUE) - nrow(featureValues(object)),
                            " features were flagged")
              featureDefinitions(object)$possible_contaminants <- vals
              ph <- XProcessHistory(param = filter,
                                    date. = date(),
                                    type. = .PROCSTEP.FEATURE.FILTERING,
                                    fileIndex = seq_along(object))
              object <- addProcessHistory(object, ph)
              validObject(object)
              object
          }
)

#' @rdname BlankFlag
setMethod("filterFeatures",
          signature(object = "SummarizedExperiment",
                    filter = "BlankFlag"),
          function(object, filter, assay = 1){
              .check_index_range(filter@blankIndex, ncol(object),
                                 name = "blankIndex")
              .check_index_range(filter@qcIndex, ncol(object), name = "qcIndex")
              x <- assay(object, assay)[, filter@qcIndex]
              y <- assay(object, assay)[, filter@blankIndex]
              vals <- rowBlank(x = x, y = y,
                               na.rm = filter@na.rm,
                               threshold = filter@threshold)
              message(sum(vals, na.rm = TRUE) - length(object),
                      " features were flagged")
              rowData(object)$possible_contaminants <- vals
              object
          }
)

#' @noRd
.check_index_range <- function(x, l, name = "") {
    if (!all(x %in% seq_len(l)) | length(x) == 0)
        stop(name, " should be between 1 and ", l)
}
