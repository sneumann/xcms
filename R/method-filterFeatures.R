#' @title Filtering of features based on conventional quality assessment
#'
#' @description
#'
#' When dealing with metabolomics results, it is often necessary to filter
#' features based on certain criteria. These criteria are typically derived
#' from statistical formulas applied to full rows of data, where each row
#' represents a feature.
#' The `filterFeatures` function filters features based on these conventional
#' quality assessment criteria. Multiple type of filtering are implemented and\
#' can be defined by the `filter` argument.
#'
#' Supported `filter` arguments are:
#'
#' - [`RsdFilter`]: Calculates the relative standard deviation
#'  (i.e. coefficient of variation) of each features in QC (Quality Control)
#'  samples and filter them in the input object according to a provided
#'  threshold.
#'
#' - [`DratioFilter`]: Computes the D-ratio or *dispersion ratio*, defined as
#'  the standard deviation for QC samples divided by the standard deviation for
#'  biological test samples, for each features and filter them according to a
#'  provided threshold
#'
#' - [`PercentMissingFilter`]: Determine the percentage of missing values for
#'  each features in  and filters them according to a provided threshold.
#'
#' - [`BlankFlag`]: Identifies features where the mean of test samples
#'  is lower than a specified multiple of the mean of blank samples. This can
#'  be used to flag features that result from contamination in the solvent of
#'  the samples. A new column `possible_contaminant` is added to the
#'  `featureDefinitions` reflecting this.
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
#' assay. Default is 1, which mean the first assay of the `object` will
#' be evaluated.
#'
#' @param ... Optional parameters.
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
#'
#' ## Load a test data set with features
#' test_xcms <- loadXcmsData()
#'
#' ## Set up parameter to filter based on coefficient of variation
#' rsd_filter <- RsdFilter(qcIndex = sampleData(test_xcms)$sample_type == "QC")
#'
#' filtered_data_rsd <- filterFeatures(object = test_xcms, filter = rsd_filter)
#'
#' ## Set up parameter to filter based on D-ratio and filter
#' dratio_filter <- DratioFilter(qcIndex = sampleData(test_xcms)$sample_type == "QC",
#'                               studyIndex = sampleData(test_xcms)$sample_type == "study")
#'
#' filtered_data_dratio <- filterFeatures(object = test_xcms,
#' filter = dratio_filter)
#'
#' ## Set up parameter to filter based on the percent of missing data
#' missing_data_filter <- PercentMissingFilter(f = sampleData(test_xcms)$sample_type)
#'
#' filtered_data_missing <- filterFeatures(object = test_xcms,
#' filter = missing_data_filter)
#'
#' ## Set up parameter to flag possible contaminant based on blank samples
#' filter <- BlankFlag(
#'     qcIndex = sampleData(test_xcms)$sample_type == "QC",
#'     blankIndex = sampleData(test_xcms)$sample_type == "study")
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
#' `XcmsExperiment` or `SummarizedExperiment` object based on their coefficent
#' of variation for a specified threshold.
#'
#' This `filter` is part of the possible dispatch of the generic function
#' `filterFeatures`. Features *above* the user-input threshold will be
#' removed from the entire dataset.
#'
#' @param threshold `numeric` value representing the threshold. Features with a
#' coefficient of variation higher than this will be removed from the entire
#' dataset.
#'
#' @param qcIndex `integer` (or logical) vector corresponding to the indexes of
#' QC samples.
#'
#' @param na.rm `logical` indicates whether missing values (`NA`) should be
#' removed prior to the calculations.
#'
#' @param mad `logical` indicates whether the *Median Absolute Deviation*
#' (MAD) should be used instead of the standard deviation. This is suggested
#' for non-gaussian distributed data.
#'
#' @inheritParams filterFeatures
#'
#' @return For `RsdFilter`: a `RsdFilter` class. `filterFeatures` return the
#' input object minus the features that did not met the user input threshold.
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
          function(object, filter){
              if (length(filter@qcIndex) < 1 |
                  length(filter@qcIndex) > length(object))
                  stop("'qcIndex' must be within object length range")
            vals <- featureValues(object)[, filter@qcIndex]
            vals <- rowRsd(x = vals,  na.rm = filter@na.rm, mad = filter@mad)
            fts_idx <- which(vals <= filter@threshold)
            message(paste(length(vals) - length(fts_idx),
                          "features were removed"))
            ph <- XProcessHistory(param = filter,
                                  date. = date(),
                                  type. = .PROCSTEP.FEATURE.FILTERING,
                                  fileIndex = seq_along(object))
            object <- addProcessHistory(object, ph)
            object <- filterFeatureDefinitions(object, features = fts_idx)
          }
)

#' @rdname RsdFilter
setMethod("filterFeatures",
          signature(object = "SummarizedExperiment",
                    filter = "RsdFilter"),
          function(object, filter, assay = 1){
              if (length(filter@qcIndex) < 1 | length(filter@qcIndex) > ncol(object))
                  stop("'qcIndex' must be within object length range")
              vals <- assay(object, assay)[, filter@qcIndex]
              vals <- rowRsd(vals,  na.rm = filter@na.rm, mad = filter@mad)
              fts_idx <- which(vals <= filter@threshold)
              message(paste(length(vals) - length(fts_idx), "features were removed"))
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
#' each feature of the object.
#'
#' This `filter` is part of the possible dispatch of the
#' generic function `filterFeatures`.  Features *above* the user-input
#' threshold will be removed from the entire dataset.
#'
#' @param threshold `numeric` value representing the threshold. Features with a
#' coefficient of variation higher than this will be removed from the entire
#' dataset.
#'
#' @param qcIndex `integer` (or `logical`) Vector corresponding to the indexes of
#' QC samples.
#'
#' @param studyIndex `integer` (or `logical`) Vector corresponding of the
#' indexes of study samples.
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
          function(object, filter){
              if (length(filter@studyIndex) < 1 |
                  length(filter@studyIndex) > length(object))
                  stop("'studyIndex' must be within object length range")
              if (length(filter@qcIndex) < 1 |
                  length(filter@qcIndex) > length(object))
                  stop("'qcIndex' must be within object length range")
              x <- featureValues(object)[, filter@studyIndex]
              y <- featureValues(object)[, filter@qcIndex]
              vals <- rowDratio(x = x, y = y,
                                na.rm = filter@na.rm,
                                mad = filter@mad)
              fts_idx <- which(vals <= filter@threshold)
              message(paste(length(vals) - length(fts_idx),
                            "features were removed"))
              ph <- XProcessHistory(param = filter,
                                    date. = date(),
                                    type. = .PROCSTEP.FEATURE.FILTERING,
                                    fileIndex = seq_along(object))
              object <- addProcessHistory(object, ph)
              object <- filterFeatureDefinitions(object, features = fts_idx)
          }
)

#' @rdname DratioFilter
setMethod("filterFeatures",
          signature(object = "SummarizedExperiment",
                    filter = "DratioFilter"),
          function(object, filter, assay = 1){
              if (length(filter@studyIndex) < 1 |
                  length(filter@studyIndex) > ncol(object))
                  stop("'studyIndex' must be within object length range")
              if (length(filter@qcIndex) < 1 |
                  length(filter@qcIndex) > ncol(object))
                  stop("'qcIndex' must be within object length range")
              x <-  assay(object, assay)[, filter@studyIndex]
              y <- assay(object, assay)[, filter@qcIndex]
              vals <- rowDratio(x = x, y = y,
                                na.rm = filter@na.rm,
                                mad = filter@mad)
              fts_idx <- which(vals <= filter@threshold)
              message(paste(length(vals) - length(fts_idx),
                            "features were removed"))
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
#' percentage of missing values for each features and filters them according to
#' a provided threshold.
#'
#' This `filter` is part of the possible dispatch of the generic function
#' `filterFeatures`. The features *above* the user input threshold will be
#' removed by calling the `filterFeatures` function.
#'
#' @param f `vector` of the same length as the `object`, specifying the sample
#' type for each sample in the dataset. The percentage of missing values will
#' be computed for each of the sample types.
#'
#' @param threshold `numeric`  percent of accepted missing data for one feature.
#' If a feature has a percentage of missing data *above* the threshold, it
#' will be removed.
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
#' @examples
#'
NULL

#' @noRd
setClass("PercentMissingFilter",
         slots = c(threshold = "numeric",
                   f = "character"),
         contains = "Param",
         prototype = prototype(
             threshold = numeric(),
             f = character()),
         validity = function(object) {
             msg <- NULL
             if (length(object@threshold) != 1)
                 msg <- c("'threshold' has to be a numeric of length 1")
             msg
         })

#' @rdname PercentMissingFilter
#'
#' @export
PercentMissingFilter <- function(threshold = 30, f = character()) {
    new("PercentMissingFilter", threshold = threshold, f = f)
}

#' @rdname PercentMissingFilter
setMethod("filterFeatures",
          signature(object = "XcmsResult",
                    filter = "PercentMissingFilter"),
          function(object, filter){
              if (length(filter@f) != length(object))
                  stop("'f' must be same lenght as object")
              f <- factor(filter@f)
              fts_idx <- c()
              for (i in levels(f)){
                  spl_idx <- which(f == i)
                  vals <- rowPercentMissing(featureValues(object)[, spl_idx])
                  fts_idx <- c(fts_idx, which(vals <= filter@threshold))
              }
              fts_idx <- order(unique(fts_idx))
              message(paste(length(vals) - length(fts_idx),
                            "features were removed"))
              ph <- XProcessHistory(param = filter,
                                    date. = date(),
                                    type. = .PROCSTEP.FEATURE.FILTERING,
                                    fileIndex = seq_along(object))
              object <- addProcessHistory(object, ph)
              object <- filterFeatureDefinitions(object, features = fts_idx)
          }
)

#' @rdname PercentMissingFilter
setMethod("filterFeatures",
          signature(object = "SummarizedExperiment",
                    filter = "PercentMissingFilter"),
          function(object, filter, assay = 1){
              if (length(filter@f) != ncol(object))
                  stop("'f' must be same lenght as object")
              f <- factor(filter@f)
              fts_idx <- c()
              for (i in levels(f)){
                  spl_idx <- which(f == i)
                  vals <- rowPercentMissing(assay(object, assay)[, spl_idx])
                  fts_idx <- c(fts_idx, which(vals <= filter@threshold))
              }
              fts_idx <- order(unique(fts_idx))
              message(paste(length(vals) - length(fts_idx),
                            "features were removed"))
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
#' generic function `filterFeatures`. Features *above* the user-input
#' threshold will be flagged by calling the `filterFeatures` function. This
#' means that an extra column will be created in `featureDefinitions` or
#' `rowData` called `possible_contaminants` with a logical value for each
#' feature.
#'
#' @param threshold `numeric` indicates the minimum difference
#' required between the mean of a feature in samples compared to the mean of
#' the same feature in blanks for it to not be considered a possible
#' contaminant. For example, the default threshold of 2 signifies that the mean
#' of the features in samples has to be at least twice the mean in blanks for
#' it not to be flagged as a possible contaminant.
#'
#' @param blankIndex `integer` (or `logical`) vector corresponding to the
#' indexes of blank samples.
#'
#' @param qcIndex `integer` (or `logical`) vector corresponding to the
#' indexes of quality control (QC) samples.
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
          function(object, filter){
              if (length(filter@blankIndex) < 1 |
                  length(filter@blankIndex) > length(object))
                  stop("'blankIndex' must be within object length range")
              if (length(filter@qcIndex) < 1 |
                  length(filter@qcIndex) > length(object))
                  stop("'qcIndex' must be within object length range")
              x <- featureValues(object)[, filter@qcIndex]
              y <- featureValues(object)[, filter@blankIndex]
              vals <- rowBlank(x = x, y = y,
                               na.rm = filter@na.rm,
                               threshold = filter@threshold)
              message(paste(length(featureValues(object)) - sum(vals),
                            "features were flagged"))
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
              if (length(filter@blankIndex) < 1 |
                  length(filter@blankIndex) > ncol(object))
                  stop("'blankIndex' must be within object length range")
              if (length(filter@qcIndex) < 1 |
                  length(filter@qcIndex) > ncol(object))
                  stop("'qcIndex' must be within object length range")
              x <- assay(object, assay)[, filter@qcIndex]
              y <- assay(object, assay)[, filter@blankIndex]
              vals <- rowBlank(x = x, y = y,
                               na.rm = filter@na.rm,
                               threshold = filter@threshold)
              message(paste(length(object) - sum(vals, na.rm = TRUE),
                            "features were flagged"))
              rowData(object)$possible_contaminants <- vals
              object
          }
)

