#' @title Compounding of LC-MS features
#'
#' @name feature-grouping
#'
#' @description
#'
#' Feature *compounding* aims at identifying and grouping LC-MS features
#' representing different ions or adducts (including isotopes) of the same
#' originating compound.
#' The [MsFeatures](https://bioconductor.org/packages/MsFeatures) package
#' provides a general framework and functionality to group features based on
#' different properties. The `groupFeatures` methods for [XCMSnExp-class]
#' objects implemented in `xcms` extend these to enable the *compounding* of
#' LC-MS data. Note that these functions simply define feature groups but don't
#' actually *aggregate* or combine the features.
#'
#' See [MsFeatures::groupFeatures()] for an overview on the general feature
#' grouping concept as well as details on the individual settings and
#' parameters.
#'
#' The available options for `groupFeatures` on `xcms` preprocessing results
#' (i.e. on `XCMSnExp` objects after correspondence analysis with
#' [groupChromPeaks()]) are:
#'
#' - Grouping by similar retention times: [groupFeatures-similar-rtime()].
#'
#' - Grouping by similar feature values across samples:
#'   [AbundanceSimilarityParam()].
#'
#' - Grouping by similar peak shape of extracted ion chromatograms:
#'   [EicSimilarityParam()].
#'
#' An ideal workflow grouping features should sequentially perform the above
#' methods (in the listed order).
#'
#' Compounded feature groups can be accessed with the `featureGroups` function.
#'
#' @param object an [XCMSnExp()] object.
#'
#' @param value for `featureGroups<-`: replacement for the feature groups in
#'     `object`. Has to be of length 1 or length equal to the number of features
#'     in `object`.
#'
#' @author Johannes Rainer, Mar Garcia-Aloy, Vinicius Veri Hernandes
#'
#' @seealso [plotFeatureGroups()] for visualization of grouped features.
#'
#' @md
NULL

#' @rdname feature-grouping
#'
#' @export
setMethod("featureGroups", "XCMSnExp", function(object) {
    if (!hasFeatures(object))
        stop("No feature definitions present. Please run 'groupChromPeak'",
             call. = FALSE)
    if (any(colnames(featureDefinitions(object)) == "feature_group"))
        as.character(featureDefinitions(object)$feature_group)
    else rep(NA_character_, nrow(featureDefinitions(object)))
})

#' @rdname feature-grouping
#'
#' @export
setReplaceMethod("featureGroups", "XCMSnExp", function(object, value) {
    if (!hasFeatures(object))
        stop("No feature definitions present. Please run 'groupChromPeak'",
             call. = FALSE)
    if (!(length(value) == 1 |
          length(value) == nrow(featureDefinitions(object))))
        stop("'value' has to be either of length 1 or equal to the number ",
             "of features in object")
    featureDefinitions(object)$feature_group <- as.character(value)
    object
})


#' @title Compounding/feature grouping based on similar retention times
#'
#' @name groupFeatures-similar-rtime
#'
#' @description
#'
#' Group features based on similar retention time. This method is supposed to be
#' used as an initial *crude* grouping of features based on the median retention
#' time of all their chromatographic peaks. All features with a difference in
#' their retention time which is `<=` parameter `diffRt` of the parameter object
#' are grouped together. If a column `"feature_group"` is found in
#' [xcms::featureDefinitions()] this is further sub-grouped by this method.
#'
#' See [MsFeatures::SimilarRtimeParam()] in `MsFeatures` for more details.
#'
#' @param msLevel `integer(1)` defining the MS level on which the features
#'     should be grouped.
#'
#' @param object [XCMSnExp()] object containing also correspondence results.
#'
#' @param param `SimilarRtimeParam` object with the settings for the method. See
#'     [MsFeatures::SimilarRtimeParam()] for details and options.
#'
#' @param ... passed to the `groupFeatures` function on numeric values.
#'
#' @return input `XCMSnExp` with feature groups added (i.e. in column
#'     `"feature_group"` of its `featureDefinitions` data frame.
#'
#' @family feature grouping methods
#'
#' @rdname groupFeatures-similar-rtime
#'
#' @importClassesFrom ProtGenerics Param
#'
#' @importFrom MsFeatures SimilarRtimeParam AbundanceSimilarityParam
#'
#' @importClassesFrom MsFeatures SimilarRtimeParam AbundanceSimilarityParam
#'
#' @importMethodsFrom MsFeatures groupFeatures featureGroups featureGroups<-
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' library(MsFeatures)
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Group chromatographic peaks across samples
#' xodg <- groupChromPeaks(faahko_sub, param = PeakDensityParam(sampleGroups = rep(1, 3)))
#'
#' ## Group features based on similar retention time (i.e. difference <= 2 seconds)
#' xodg_grp <- groupFeatures(xodg, param = SimilarRtimeParam(diffRt = 2))
#'
#' ## Feature grouping get added to the featureDefinitions in column "feature_group"
#' head(featureDefinitions(xodg_grp)$feature_group)
#'
#' table(featureDefinitions(xodg_grp)$feature_group)
#' length(unique(featureDefinitions(xodg_grp)$feature_group))
#'
#' ## Using an alternative groupiing method that creates larger groups
#' xodg_grp <- groupFeatures(xodg,
#'     param = SimilarRtimeParam(diffRt = 2, groupFun = MsCoreUtils::group))
#'
#' length(unique(featureDefinitions(xodg_grp)$feature_group))
NULL

#' @rdname groupFeatures-similar-rtime
#'
#' @importMethodsFrom xcms hasFeatures featureDefinitions featureDefinitions<-
#'
#' @importFrom MsCoreUtils group
#'
#' @exportMethod groupFeatures
#'
#' @importClassesFrom xcms XCMSnExp XProcessHistory
setMethod(
    "groupFeatures",
    signature(object = "XCMSnExp", param = "SimilarRtimeParam"),
    function(object, param, msLevel = 1L, ...) {
        fgs <- featureGroups(object)
        if (length(msLevel) > 1)
            stop("Currently only grouping of features from a single MS level",
                 " is supported.", call. = FALSE)
        if (any(colnames(featureDefinitions(object)) == "ms_level"))
            is_msLevel <- featureDefinitions(object)$ms_level == msLevel
        else is_msLevel <- rep(TRUE, nrow(featureDefinitions(object)))
        if (all(is.na(fgs)))
            fgs <- rep("FG", length(fgs))
        fgs[!is_msLevel] <- NA
        nas <- is.na(fgs)
        fgs <- factor(fgs, levels = unique(fgs))
        rtl <- split(featureDefinitions(object)$rtmed, fgs)
        res <- lapply(
            rtl, function(z, param)
                MsFeatures:::.format_id(groupFeatures(z, param = param, ...)),
            param = param, ...)
        res <- paste(fgs, unsplit(res, f = fgs), sep = ".")
        if (any(nas))
            res[nas] <- NA_character_
        featureGroups(object) <- res
        xph <- new("XProcessHistory", param = param, date = date(),
                   type = .PROCSTEP.FEATURE.GROUPING,
                   fileIndex = 1:length(fileNames(object)),
                   msLevel = as.integer(msLevel))
        object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
        validObject(object)
        object
    })

#' @title Compounding/feature grouping based on similarity of abundances across samples
#'
#' @name groupFeatures-abundance-correlation
#'
#' @description
#'
#' Features from the same originating compound are expected to have similar
#' intensities across samples. This method this groups features based on
#' similarity of abundances (i.e. *feature values*) across samples.
#' See also [AbundanceSimilarityParam()] for additional information and details.
#'
#' This help page lists parameters specific for `xcms` result objects (i.e. the
#' [XCMSnExp()] object). Documentation of the parameters for the similarity
#' calculation is available in the [AbundanceSimilarityParam()] help page in
#' the `MsFeatures` package.
#'
#' @param filled `logical(1)` whether filled-in values should be included in
#'     the correlation analysis. Defaults to `filled = TRUE`.
#'
#' @param intensity `character(1)` passed to the `featureValues` call. See
#'     [featureValues()] for details. Defaults to `intensity = "into"`.
#'
#' @param method `character(1)` passed to the `featureValues` call. See
#'     [featureValues()] for details. Defaults to `method = "medret"`.
#'
#' @param msLevel `integer(1)` defining the MS level on which the features
#'     should be grouped.
#'
#' @param object [XCMSnExp()] object containing also correspondence results.
#'
#' @param param `AbudanceSimilarityParam` object with the settings for the
#'     method. See [AbundanceSimilarityParam()] for details on the grouping
#'     method and its parameters.
#'
#' @param value `character(1)` passed to the `featureValues` call. See
#'     [featureValues()] for details. Defaults to `value = "into"`.
#'
#' @param ... additional parameters passed to the `groupFeatures` method for
#'     `matrix`.
#'
#' @return input `XCMSnExp` with feature group definitions added to a column
#'     `"feature_group"` in its `featureDefinitions` data frame.
#'
#' @family feature grouping methods
#'
#' @rdname groupFeatures-abundance-correlation
#'
#' @importClassesFrom MsFeatures AbundanceSimilarityParam
#'
#' @importFrom MsFeatures AbundanceSimilarityParam
#'
#' @author Johannes Rainer
#'
#' @seealso feature-grouping for a general overview.
#'
#' @md
#'
#' @examples
#'
#' library(MsFeatures)
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Group chromatographic peaks across samples
#' xodg <- groupChromPeaks(faahko_sub, param = PeakDensityParam(sampleGroups = rep(1, 3)))
#'
#' ## Group features based on correlation of feature values (integrated
#' ## peak area) across samples. Note that there are many missing values
#' ## in the feature value which influence grouping of features in the present
#' ## data set.
#' xodg_grp <- groupFeatures(xodg,
#'     param = AbundanceSimilarityParam(threshold = 0.8))
#' table(featureDefinitions(xodg_grp)$feature_group)
#'
#' ## Group based on the maximal peak intensity per feature
#' xodg_grp <- groupFeatures(xodg,
#'     param = AbundanceSimilarityParam(threshold = 0.8, value = "maxo"))
#' table(featureDefinitions(xodg_grp)$feature_group)
NULL

#' @rdname groupFeatures-abundance-correlation
#'
setMethod(
    "groupFeatures",
    signature(object = "XCMSnExp", param = "AbundanceSimilarityParam"),
    function(object, param, msLevel = 1L, method = c("medret", "maxint", "sum"),
             value = "into", intensity = "into", filled = TRUE, ...) {
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run ",
                 "first 'groupChromPeaks'")
        if (length(msLevel) > 1)
            stop("Currently only grouping of features from a single MS level",
                 " is supported.")
        fgs <- featureGroups(object)
        fgs_orig <- fgs
        if (any(colnames(featureDefinitions(object)) == "ms_level"))
            is_msLevel <- featureDefinitions(object)$ms_level == msLevel
        else is_msLevel <- rep(TRUE, nrow(featureDefinitions(object)))
        if (all(is.na(fgs)))
            fgs <- rep("FG", length(fgs))
        fgs[!is_msLevel] <- NA
        nas <- is.na(fgs)
        fgs <- factor(fgs, levels = unique(fgs))
        l <- split.data.frame(
            featureValues(object, method = method, value = value,
                          intensity = intensity, filled = filled), fgs)
        res <- lapply(
            l, function(z, param)
                MsFeatures:::.format_id(groupFeatures(z, param = param, ...)),
            param = param, ...)
        res <- paste(fgs, unsplit(res, f = fgs), sep = ".")
        if (any(nas))
            res[nas] <- fgs_orig[nas]
        featureGroups(object) <- res
        xph <- new("XProcessHistory", param = param, date = date(),
                   type = .PROCSTEP.FEATURE.GROUPING,
                   fileIndex = 1:length(fileNames(object)),
                   msLevel = as.integer(msLevel))
        object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
        validObject(object)
        object
    })


#' @title Plot feature groups in the m/z-retention time space
#'
#' @description
#'
#' `plotFeatureGroups` visualizes defined feature groups in the m/z by
#' retention time space. Features are indicated by points with features from
#' the same feature group being connected by a line. See [featureGroups()]
#' for details on and options for feature grouping.
#'
#' @param x [XCMSnExp()] object with grouped features (i.e. after calling
#'     [groupFeatures()].
#'
#' @param xlim `numeric(2)` with the lower and upper limit for the x-axis.
#'
#' @param ylim `numeric(2)` with the lower and upper limit for the y-axis.
#'
#' @param xlab `character(1)` with the label for the x-axis.
#'
#' @param ylab `character(1)` with the label for the y-axis.
#'
#' @param pch the plotting character. Defaults to `pch = 4` i.e. plotting
#'     features as crosses. See [par()] for more information.
#'
#' @param col color to be used to draw the features. At present only a single
#'     color is supported.
#'
#' @param type plotting type (see [par()]). Defaults to `type = "o"` which
#'     draws each feature as a point and connecting the features of the same
#'     feature group with a line.
#'
#' @param main `character(1)` with the title of the plot.
#'
#' @param featureGroups optional `character` of feature group IDs to draw only
#'     specified feature group(s). If not provided, all feature groups are
#'     drawn.
#'
#' @importFrom graphics lines
#'
#' @md
#'
#' @export
#'
#' @author Johannes Rainer
plotFeatureGroups <- function(x, xlim = numeric(), ylim = numeric(),
                              xlab = "retention time", ylab = "m/z",
                              pch = 4, col = "#00000060", type = "o",
                              main = "Feature groups",
                              featureGroups = character()) {
    if (!inherits(x, "XCMSnExp"))
        stop("'x' is supposed to be an xcms result object ('XCMSnExp')")
    if (!length(featureGroups(x)))
        stop("No feature groups present. Please run 'groupFeatures' first")
    fts <- factor(featureGroups(x))
    if (!length(featureGroups))
        featureGroups <- levels(fts)
    fts <- fts[fts %in% featureGroups]
    fts <- droplevels(fts)
    if (!length(fts))
        stop("None of the specified feature groups found")
    fdef <- featureDefinitions(x)[featureGroups(x) %in% fts, ]
    rts <- split(fdef$rtmed, fts)
    mzs <- split(fdef$mzmed, fts)
    xy <- cbind(
        x = unlist(lapply(rts, function(z) c(z, NA)), use.names = FALSE),
        y = unlist(lapply(mzs, function(z) c(z, NA)), use.names = FALSE))
    if (length(xlim) != 2)
        xlim <- range(unlist(rts, use.names = FALSE))
    if (length(ylim) != 2)
        ylim <- range(unlist(mzs, use.names = FALSE))
    plot(3, 3, pch = NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
    lines(xy, type = type, col = col, pch = pch)
}

## #' @title Extract spectra for feature groups
## #'
## #' @description
## #'
## #' `featureGroupSpectra` allows to extract a `Spectrum` object for each feature
## #' group in `x`. Based on the specified function `FUN` different *types* of
## #' spectra can be returned:
## #'
## #' - `featureGroupPseudoSpectrum` creates a *pseudo* spectrum based on the
## #'   feature values (defined by `value`) of all features within a feature group
## #'   (i.e. each feature is represented as a mass peak in the resulting
## #'   spectrum). The reported m/z values will be the `"mzmed"` of the respective
## #'   feature from the [featureDefinitions()] data frame. The associated
## #'   intensity is calculated from the values of the features from the feature
## #'   group: by default, for each feature, the median intensity across all
## #'   samples part of `subset` is reported. Parameters `value` and `filled` are
## #'   passed to the internal call to [featureValues()] that returns the features'
## #'   values which are used in these calculations. Parameter `n` allows to
## #'   further restrict the samples being considered in the calculations: for each
## #'   feature group samples are first ordered by the sum of signal of the
## #'   features of the group and then only the *top n* samples are used in the
## #'   calculations.
## #'
## #'   Examples:
## #'   To report the mean intensity of each feature in the 10 samples with the
## #'   highest signal for the feature group use `n = 10` and
## #'   `intensityFun = mean`. The m/z values reported in the `Spectrum` object
## #'   of a feature group will be the `"mzmed"` of the features, the intensity
## #'   values the mean intensity (`value = "maxo"`) across the 10 samples with
## #'   the highest signal for that feature group.
## #'
## #'   To report the maximal intensity (`value = "maxo"` of each feature in
## #'   samples 1, 4, 8 and 10 use `subset = c(1, 4, 8, 10)` and
## #'   `intensityFun = max`. More examples in the examples section.
## #'
## #' - `featureGroupFullScan`: reports the full MS1 spectrum (full scan) in the
## #'   sample with the highest total signal (defined by `value`) for the feature
## #'   group at the retention time closest to the median `"rtmed"` across all
## #'   features of the feature group.
## #'
## #' @param x [XCMSnExp()] object with available `featureGroups()`.
## #'
## #' @param featureGroup `character` with the IDs of the feature group(s) for
## #'     which the spectra should be returned. Defaults to all feature groups
## #'     defined in `x`. Only `featureGroupSpectra` supports
## #'     `length(featureGroup)` to be of length > 1.
## #'
## #' @param filled for `featureGroupPseudoSpectra`: `logical(1)` whether
## #'     filled-in values should also be considered. See [featureValues()] for
## #'     details.
## #'
## #' @param FUN `function` to be used to define the spectrum for each feature
## #'     group. Can be `featureGroupPseudoSpectrum`, `featureGroupFullScan` or
## #'     any function taking parameters `featureGroup`, `x`, `fvals`.
## #'
## #' @param fvals for `featureGroupPseudoSpectra` and `featureGroupFullScan`:
## #'     `matrix` with feature values (rows being features and columns samples)
## #'     such as returned by [featureValues()].
## #'
## #' @param intensityFun for `featureGroupPseudoSpectra`: `function` that should
## #'     be applied across samples (defined by `subset`) of the feature value
## #'     matrix to calculate the intensity for each mass peak of the returned
## #'     pseudo spectrum. By default (`intensityFun = median`) the median
## #'     intensity of a feature across all samples (defined by `subset` and `n`)
## #'     is used. See description section for examples.
## #'
## #' @param n for `featureGroupPseudoSpectra`: `integer(1)` defining the *top n*
## #'     samples (in `subset`) on which spectra should be defined. Samples are
## #'     ordered based on the sum of signal (defined by parameter `value`) from
## #'     the features of each feature group. See description section for more
## #'     details.
## #'
## #' @param subset `integer` with indices of specific samples if spectra should
## #'     only be defined on a subset of samples. See description section for
## #'     details.
## #'
## #' @param value `character(1)` specifying the column in `chromPeaks` matrix to
## #'     be used as *feature values* for each feature. This parameter is passed
## #'     to the [featureValues()] call.
## #'
## #' @param ... additional parameters passed down to the function specifyed with
## #'     `FUN`.
## #'
## #' @return for `featureGroupSpectra`: `MSpectra` object of length equal to the
## #'     number of feature groups in `x` and each element being one spectrum.
## #'     For all other functions: a `Spectrum` object.
## #'
## #' @author Johannes Rainer
## #'
## #' @importMethodsFrom xcms filterFile hasAdjustedRtime hasFeatures rtime
## #'
## #' @importFrom xcms applyAdjustedRtime
## #'
## #' @importFrom S4Vectors DataFrame
## #'
## #' @importFrom IRanges CharacterList
## #'
## #' @importFrom MSnbase MSpectra
## #'
## #' @importFrom stats median
## #'
## #' @export
## #' @md
## #' @examples
## #'
## #' ## Load test data set from xcms
## #' library(xcms)
## #' data(faahko_sub)
## #' ## Update the path to the files for the local system
## #' dirname(faahko_sub) <- system.file("cdf/KO/", package = "faahKO")
## #'
## #' ## Group chromatographic peaks across samples
## #' xodg <- groupChromPeaks(faahko_sub, param = PeakDensityParam(sampleGroups = rep(1, 3)))
## #' ## Perform correspondence analysis
## #' xdata <- groupChromPeaks(faahko_sub,
## #'     param = PeakDensityParam(sampleGroup = rep(1, 3)))
## #'
## #' ## Group features
## #' xdata <- groupFeatures(xdata, param = SimilarRtimeParam(4))
## #' xdata <- groupFeatures(xdata,
## #'     param = AbundanceSimilarityParam(threshold = 0.3))
## #'
## #' sort(table(featureGroups(xdata)))
## #'
## #' ################
## #' ## featureGroupSpectra
## #' ##
## #'
## #' ## Get a pseudo spectrum for each feature group
## #' res <- featureGroupSpectra(xdata)
## #' res
## #'
## #' ## Get a full scan spectrum for a subset of the feature groups
## #' ## considering only the subset of the last two samples
## #' res <- featureGroupSpectra(xdata,
## #'     featureGroup = unique(featureGroups(xdata))[1:4],
## #'     FUN = featureGroupFullScan, subset = 2:3)
## #' res
## #'
## #' ################
## #' ## Pseudo Spectrum
## #' ##
## #'
## #' ## Get the pseudo spectrum for one feature group reporting the per-feature
## #' ## maximal "maxo" value across samples as the spectrum's intensities
## #' res <- featureGroupPseudoSpectrum(featureGroup = "FG.010.001", xdata,
## #'     fvals = featureValues(xdata, value = "maxo"), intensityFun = max)
## #'
## #' intensity(res)
## #' mz(res)
## #'
## #' ## Get the pseudo spectrum using the values in the one sample with the
## #' ## highest total sum of signal ("maxo") for the feature group.
## #' res <- featureGroupPseudoSpectrum(featureGroup = "FG.010.001", xdata,
## #'     fvals = featureValues(xdata, value = "maxo"), n = 1)
## #'
## #' intensity(res)
## #' mz(res)
## #'
## #'
## #' ################
## #' ## Full Scan Spectrum
## #' ##
## #'
## #' ## Get the full MS1 spectrum from the sample with the highest total signal
## #' ## of one specific feature group
## #' res <- featureGroupFullScan(featureGroup = "FG.010.001", xdata,
## #'     fvals = featureValues(xdata, value = "maxo"))
## #'
## #' plot(mz(res), intensity(res), type = "h", xlab = "m/z", ylab = "intensity")
## #' ## Highlight the peaks for the features of the group.
## #' idx <- which(featureGroups(xdata) == "FG.001.001")
## #' points(x = featureDefinitions(xdata)$mzmed[idx],
## #'     y = rep(0, length(idx)), pch = 4, col = "red")
## featureGroupSpectra <- function(x, featureGroup = featureGroups(x),
##                                 FUN = featureGroupPseudoSpectrum,
##                                 value = "maxo", filled = TRUE,
##                                 subset = seq_along(fileNames(x)),
##                                 ...) {
##     if (!all(subset %in% seq_along(fileNames(x))))
##         stop("'subset' is expected to be an integer vector with values ",
##              "between 1 and ", length(fileNames(x)))
##     subset <- unique(subset)
##     if (!hasFeatures(x))
##         stop("No feature definitions present. Please run 'groupChromPeaks' first")
##     featureGroup <- unique(featureGroup)
##     featureGroup <- featureGroup[!is.na(featureGroup)]
##     if (!length(featureGroup))
##         stop("No feature groups present. Please run 'groupFeatures' first")
##     if (!all(featureGroup %in% featureGroups(x)))
##         stop("Not all feature groups defined with parameter 'featureGroup' ",
##              "found in 'featureGroups(x)'")
##     if (length(subset) < length(fileNames(x)))
##         x <- filterFile(x, subset, keepFeatures = TRUE)
##     fvals <- featureValues(x, method = "maxint", intensity = value,
##                            value = value, filled = filled)
##     res <- lapply(featureGroup, FUN, x = x, fvals = fvals, ...)
##     fids <- split(rownames(featureDefinitions(x)), featureGroups(x))
##     MSnbase::MSpectra(res, elementMetadata = DataFrame(
##                               feature_group = featureGroup,
##                               feature_id = CharacterList(fids[featureGroup],
##                                                          compress = FALSE)))
## }

## #' @rdname featureGroupSpectra
## #'
## #' @importClassesFrom MSnbase Spectrum Spectrum1 Spectrum2 MSpectra
## #'
## #' @importMethodsFrom MSnbase polarity
## #'
## #' @export
## featureGroupPseudoSpectrum <- function(featureGroup = character(), x,
##                                        fvals = featureValues(x),
##                                        n = ncol(fvals),
##                                        intensityFun = median, ...) {
##     if (n < 1 || n > ncol(fvals))
##         stop("'n' has to be an integer between 1 and ", ncol(fvals))
##     ft_idx <- which(featureGroups(x) == featureGroup)
##     ft_fvals <- fvals[ft_idx, , drop = FALSE]
##     ordr <- order(colSums(ft_fvals, na.rm = TRUE), decreasing = TRUE)
##     ft_fvals <- ft_fvals[, ordr, drop = FALSE][, 1:n, drop = FALSE]
##     ft_fdef <- extractROWS(featureDefinitions(x), ft_idx)
##     if (any(colnames(ft_fdef) == "ms_level") && all(ft_fdef$ms_level == 1))
##         cls <- "Spectrum1"
##     else cls <- "Spectrum2"
##     sp <- new(cls)
##     sp@rt <- median(ft_fdef$rtmed)
##     sp@mz <- ft_fdef$mzmed
##     sp@intensity <- apply(ft_fvals, MARGIN = 1, FUN = intensityFun, na.rm = TRUE)
##     sp@peaksCount <- length(ft_idx)
##     sp@centroided <- TRUE
##     sp@polarity <- polarity(x)[1]
##     sp
## }

## #' @rdname featureGroupSpectra
## #'
## #' @importFrom S4Vectors extractROWS
## #'
## #' @export
## featureGroupFullScan <- function(featureGroup = character(), x,
##                                  fvals = featureValues(x), ...) {
##     ft_idx <- which(featureGroups(x) == featureGroup)
##     ft_fvals <- fvals[ft_idx, , drop = FALSE]
##     samp_idx <- which.max(colSums(ft_fvals, na.rm = TRUE))
##     ft_fdef <- extractROWS(featureDefinitions(x), ft_idx)
##     if (hasAdjustedRtime(x))
##         x <- applyAdjustedRtime(x)
##     x <- filterFile(as(x, "OnDiskMSnExp"), samp_idx)
##     rtmed <- median(ft_fdef$rtmed)
##     sp <- x[[which.min(abs(rtime(x) - rtmed))]]
##     sp@fromFile <- samp_idx
##     sp
## }

#' @title Compounding/feature grouping based on similarity of extracted ion chromatograms
#'
#' @aliases EicSimilarityParam-class
#'
#' @name groupFeatures-eic-similarity
#'
#' @description
#'
#' Features from the same originating compound are expected to share their
#' elution pattern (i.e. chromatographic peak shape) with it.
#' Thus, this methods allows to group features based on similarity of their
#' extracted ion chromatograms (EICs). The similarity calculation is performed
#' separately for each sample with the similarity score being aggregated across
#' samples for the final generation of the similarity matrix on which the
#' grouping (considering parameter `threshold`) will be performed.
#'
#' The [compareChromatograms()] function is used for similarity calculation
#' which by default calculates the Pearson's correlation coefficient. The
#' settings for `compareChromatograms` can be specified with parameters
#' `ALIGNFUN`, `ALIGNFUNARGS`, `FUN` and `FUNARGS`. `ALIGNFUN` defaults to
#' [alignRt()] and is the function used to *align* the chromatograms before
#' comparison. `ALIGNFUNARGS` allows to specify additional arguments for the
#' `ALIGNFUN` function. It defaults to
#' `ALIGNFUNARGS = list(tolerance = 0, method = "closest")` which ensures that
#' data points from the same spectrum (scan, i.e. with the same retention time)
#' are compared between the EICs from the same sample. Parameter `FUN` defines
#' the function to calculate the similarity score and defaults to `FUN = cor`
#' and `FUNARGS` allows to pass additional arguments to this function (defaults
#' to `FUNARGS = list(use = "pairwise.complete.obs")`. See also
#' [compareChromatograms()] for more information.
#'
#' The grouping of features based on the EIC similarity matrix is performed
#' with the function specified with parameter `groupFun` which defaults to
#' `groupFun = groupSimilarityMatrix` which groups all rows (features) in the
#' similarity matrix with a similarity score larger than `threshold` into the
#' same cluster. This creates clusters of features in which **all** features
#' have a similarity score `>= threshold` with **any** other feature in that
#' cluster. See [groupSimilarityMatrix()] for details. Additional parameters to
#' that function can be passed with the `...` argument.
#'
#' This feature grouping should be called **after** an initial feature
#' grouping by retention time (see [SimilarRtimeParam()]). The feature groups
#' defined in columns `"feature_group"` of `featureDefinitions(object)` (for
#' features matching `msLevel`) will be used and refined by this method.
#' Features with a value of `NA` in `featureDefinitions(object)$feature_group`
#' will be skipped/not considered for feature grouping.
#'
#' @note
#'
#' While being possible to be performed on the full data set without prior
#' feature grouping, this is not suggested for the following reasons: I) the
#' selection of the top `n` samples with the highest signal for the
#' *feature group* will be biased by very abundant compounds as this is
#' performed on the full data set (i.e. the samples with the highest overall
#' intensities are used for correlation of all features) and II) it is
#' computationally much more expensive because a pairwise correlation between
#' all features has to be performed.
#'
#' It is also suggested to perform the correlation on a subset of samples
#' per feature with the highest intensities of the peaks (for that feature)
#' although it would also be possible to run the correlation on all samples by
#' setting `n` equal to the total number of samples in the data set. EIC
#' correlation should however be performed ideally on samples in which the
#' original compound is highly abundant to avoid correlation of missing values
#' or noisy peak shapes as much as possible.
#'
#' By default also the signal which is outside identified chromatographic peaks
#' is excluded from the correlation.
#'
#' @param ALIGNFUN `function` defining the function to be used to *align*
#'     chromatograms prior similarity calculation. Defaults to
#'     `ALIGNFUN = alignRt`. See [alignRt()] and [compareChromatograms()] for
#'     more information.
#'
#' @param ALIGNFUNARGS **named** `list` with arguments for `ALIGNFUN`.
#'     Defaults to `ALIGNFUNARGS = list(tolerance = 0, method = "closest")`.
#'
#' @param FUN `function` defining the function to be used to calculate a
#'     similarity between (aligned) chromatograms. Defaults to `FUN = cor`.
#'     See [cor()] and [compareChromatograms()] for more information.
#'
#' @param FUNARGS **named** `list` with arguments for `FUN`. Defaults to
#'     `FUN = list(use = "pairwise.complete.obs")`.
#'
#' @param groupFun `function` defining the function to be used to group rows
#'     based on a pairwise similarity matrix. Defaults to
#'     [groupSimilarityMatrix()].
#'
#' @param msLevel `integer(1)` defining the MS level on which the features
#'     should be grouped.
#'
#' @param n `numeric(1)` defining the total number of samples per feature group
#'     on which this similarity calculation should be performed. This value is
#'     rounded up to the next larger integer value.
#'
#' @param object [XCMSnExp()] object containing also correspondence results.
#'
#' @param onlyPeak `logical(1)` whether the correlation should be performed only
#'     on the signals within the identified chromatographic peaks
#'     (`onlyPeak = TRUE`, default) or all the signal from the extracted ion
#'     chromatogram.
#'
#' @param param `EicSimilarityParam` object with the settings for the method.
#'
#' @param threshold `numeric(1)` with the minimal required similarity score to
#'     group featues. This is passed to the `groupFun` function.
#'
#' @param value `character(1)` defining whether samples should be grouped based
#'     on the sum of the maximal peak intensity (`value = "maxo"`, the default)
#'     or the integrated peak area (`value = "into"`) for a feature.
#'
#' @param ... for `EicSimilarityParam`: additional arguments to be passed to
#'     `groupFun` and `featureChromatograms` (such as `expandRt` to expand the
#'     retention time range of each feature).
#'
#' @return input `XCMSnExp` with feature groups added (i.e. in column
#'     `"feature_group"` of its `featureDefinitions` data frame.
#'
#' @family feature grouping methods
#'
#' @seealso feature-grouping for a general overview.
#'
#' @rdname groupFeatures-eic-similarity
#'
#' @exportClass EicSimilarityParam
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' library(MsFeatures)
#' ## Load a test data set with detected peaks
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Group chromatographic peaks across samples
#' xodg <- groupChromPeaks(faahko_sub, param = PeakDensityParam(sampleGroups = rep(1, 3)))
#'
#' ## Performing a feature grouping based on EIC similarities on a single
#' ## sample
#' xodg_grp <- groupFeatures(xodg, param = EicSimilarityParam(n = 1))
#'
#' table(featureDefinitions(xodg_grp)$feature_group)
#'
#' ## Usually it is better to perform this correlation on pre-grouped features
#' ## e.g. based on similar retention time.
#' xodg_grp <- groupFeatures(xodg, param = SimilarRtimeParam(diffRt = 4))
#' xodg_grp <- groupFeatures(xodg_grp, param = EicSimilarityParam(n = 1))
#'
#' table(featureDefinitions(xodg_grp)$feature_group)
NULL

setClass(
    "EicSimilarityParam",
    slots = c(threshold = "numeric",
              n = "numeric",
              onlyPeak = "logical",
              value = "character",
              groupFun = "function",
              ALIGNFUN = "function",
              FUN = "function",
              ALIGNFUNARGS = "list",
              FUNARGS = "list",
              dots = "list"),
    contains = "Param",
    prototype = prototype(
        threshold = 0.9,
        n = 1,
        onlyPeak = TRUE,
        value = "maxo",
        groupFun = groupSimilarityMatrix,
        ALIGNFUN = alignRt,
        ALIGNFUNARGS = list(tolerance = 0, method = "closest"),
        FUN = cor,
        FUNARGS = list(use = "pairwise.complete.obs"),
        dots = list()
    ),
    validity = function(object) {
        msg <- NULL
        if (length(object@threshold) != 1 || object@threshold < 0)
            msg <- "'threshold' has to be a positive numeric of length 1"
        if (length(object@n) != 1 || object@n < 0)
            msg <- c(msg, "'n' has to be a positive numeric of length 1")
        if (length(object@onlyPeak) != 1)
            msg <- c(msg, "'onlyPeak' has to a logical of length 1")
        if (length(object@value) != 1 && !(object@value %in%
                                           c("maxo", "into")))
            msg <- c(msg, "'value' has to be either \"maxo\" or \"into\"")
        msg
    })

#' @rdname groupFeatures-eic-similarity
#'
#' @export
EicSimilarityParam <- function(threshold = 0.9, n = 1, onlyPeak = TRUE,
                               value = c("maxo", "into"),
                               groupFun = groupSimilarityMatrix,
                               ALIGNFUN = alignRt,
                               ALIGNFUNARGS = list(tolerance = 0,
                                                   method = "closest"),
                               FUN = cor,
                               FUNARGS = list(use = "pairwise.complete.obs"),
                               ...) {
    value <- match.arg(value)
    groupFun <- match.fun(groupFun)
    ALIGNFUN <- match.fun(ALIGNFUN)
    FUN <- match.fun(FUN)
    new("EicSimilarityParam", threshold = threshold, n = ceiling(n),
        onlyPeak = onlyPeak, value = value, groupFun = groupFun,
        ALIGNFUN = ALIGNFUN, ALIGNFUNARGS = ALIGNFUNARGS, FUN = FUN,
        FUNARGS = FUNARGS, dots = list(...))
}

#' @rdname groupFeatures-eic-similarity
setMethod(
    "groupFeatures",
    signature(object = "XCMSnExp", param = "EicSimilarityParam"),
    function(object, param, msLevel = 1L) {
        if (!requireNamespace("progress", quietly = TRUE))
            stop("Package 'progress' is required. Please install with ",
                 "'BiocManager::install(\"progress\")'")
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run ",
                 "first 'groupChromPeaks'")
        if (length(msLevel) > 1)
            stop("Currently only grouping of features from a single MS level",
                 " is supported.")
        n <- ceiling(param@n)
        nc <- length(fileNames(object))
        if (n > nc)
            stop("'n' should be smaller than or equal to the number of ",
                 "samples (", nc, ")")
        if (any(colnames(featureDefinitions(object)) == "ms_level"))
            is_msLevel <- featureDefinitions(object)$ms_level == msLevel
        else is_msLevel <- rep(TRUE, nrow(featureDefinitions(object)))
        if (any(colnames(featureDefinitions(object)) == "feature_group")) {
            f <- featureDefinitions(object)$feature_group
            f_new <- as.character(f)
        } else {
            f <- rep("FG", nrow(featureDefinitions(object)))
            f_new <- rep(NA_character_, length(f))
        }
        f[!is_msLevel] <- NA
        if (is.factor(f)) {
            f <- droplevels(f)
            fgroups <- levels(f)
        } else {
            fgroups <- unique(f)
            f <- factor(f, levels = fgroups)
        }
        fvals <- featureValues(object, method = "maxint", value = param@value)
        ffun <- function(z, na.rm = TRUE)
            quantile(z, probs = 0.75, na.rm = na.rm)
        pb <- progress::progress_bar$new(format = paste0("[:bar] :current/:",
                                                         "total (:percent) in ",
                                                         ":elapsed"),
                                         total = length(fgroups),
                                         clear = FALSE, force = TRUE)
        pb$tick(0)
        for (fg in fgroups) {
            idx <- which(f == fg)
            idxl <- length(idx)
            if (idxl > 1) {
                vals <- apply(fvals[idx, ], MARGIN = 2, sum, na.rm = TRUE)
                sample_idx <- order(vals, decreasing = TRUE)[seq_len(n)]
                obj_sub <- .filter_file_XCMSnExp(object, sample_idx,
                                                 keepFeatures = TRUE)
                ## Can happen that some of the features are not present in the
                ## subsetted object. Will put them into their own individual
                ## groups later.
                idx_miss <- which(!rownames(fvals)[idx] %in%
                                  rownames(featureDefinitions(obj_sub)))
                if (length(idx_miss)) {
                    tmp <- idx[idx_miss]
                    idx <- idx[-idx_miss]
                    idx_miss <- tmp
                }
                if (length(idx) > 1) {
                    eics <- do.call(
                        featureChromatograms,
                        args = c(list(obj_sub, features = rownames(fvals)[idx],
                                      filled = TRUE), param@dots))
                    if (param@onlyPeak)
                        eics <- removeIntensity(eics,
                                                which = "outside_chromPeak")
                    res <- do.call(
                        .group_eic_similarity,
                        args = c(list(as(eics, "MChromatograms"),
                                      aggregationFun = ffun,
                                      threshold = param@threshold,
                                      ALIGNFUN = param@ALIGNFUN,
                                      ALIGNFUNARGS = param@ALIGNFUNARGS,
                                      FUN = param@FUN,
                                      FUNARGS = param@FUNARGS,
                                      groupFun = param@groupFun), param@dots))
                } else res <- factor(1)
                f_new[idx] <- paste0(fg, ".", MsFeatures:::.format_id(res))
                if (length(idx_miss))
                    f_new[idx_miss] <- paste0(
                        fg, ".", MsFeatures:::.format_id(
                                                  seq((length(levels(res)) + 1),
                                                      length.out = length(
                                                          idx_miss))))
            } else
                f_new[idx] <- paste0(fg, ".001")
            pb$tick(1)
        }
        featureDefinitions(object)$feature_group <- f_new
        xph <- new("XProcessHistory", param = param, date = date(),
                   type = .PROCSTEP.FEATURE.GROUPING,
                   fileIndex = seq_along(fileNames(object)),
                   msLevel = as.integer(msLevel))
        object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
        validObject(object)
        object
    })

#' @title Group EICs based on their correlation
#'
#' @description
#'
#' `groupEicCorrelation` groups (extracted ion) chromatograms (EICs) based on
#' their similarity with each other. If this correlation is `>=` than the
#' provided `threshold` they are grouped.
#'
#' If `x` is a [MChromatograms()] object with more than one column (sample),
#' pairwise similarity is calculated between EICs first for each column
#' (sample) of `x` separately and subsequently aggregated across samples using
#' `aggregationFun`. If `x` is a `MChromatograms` with 4 rows (EICs) and 3
#' columns (samples), pairwise correlations are first calculated between all
#' 4 EICs in each of the 3 columns resulting in 3 correlation matrices (of
#' dimension 4x4). These correlation matrices are combined into a single matrix
#' by combining the 3 correlation values per comparison with
#' `aggregationFun`. By default the mean of the correlation value between e.g.
#' EIC 1 and EIC 2 in each of the 3 columns is used as the final correlation
#' value. Similar to the one-column case EICs are grouped if their (aggregated)
#' correlation coefficient is larger than `threshold`.
#'
#' @param x [MChromatograms()] object.
#'
#' @param aggregationFun `function` to combine the correlation values between
#'     pairs of EICs across samples (columns). See description for details.
#'
#' @inheritParams groupFeatures-eic-similarity
#'
#' @return `factor` same length as `nrow(x)` (if `x` is a `MChromatograms`
#'     object) or `length(x)` (if `x` is a `list`) with the group each EIC
#'     is assigned to.
#'
#' @family grouping operations
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' library(MSnbase)
#' set.seed(123)
#' chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
#' chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
#'     intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
#' chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
#'     intensity = c(53, 80, 130, 15, 5, 3, 2))
#' chrs <- MChromatograms(list(chr1, chr2, chr3))
#'
#' groupEicCorrelation(chrs)
#'
#' ## With a MChromatograms with two columns, use the maximal correlation
#' ## coefficient found in each of the columns
#' chrs <- MChromatograms(list(chr1, chr2, chr3, chr1, chr2, chr3), ncol = 2)
#' groupEicCorrelation(chrs, aggregationFun = max)
.group_eic_similarity <- function(x, aggregationFun = mean,
                                  threshold = 0.8, ALIGNFUN = alignRt,
                                  ALIGNFUNARGS = list(tolerance = 0), FUN = cor,
                                  FUNARGS = list(use = "pairwise.complete.obs"),
                                  groupFun = groupSimilarityMatrix,
                                  ...) {
    nr <- nrow(x)
    nc <- ncol(x)
    res <- array(NA_real_, dim = c(nr, nr, nc))
    for (i in seq_len(nc))
        res[, , i] <- compareChromatograms(x[, i],
                                           ALIGNFUN = ALIGNFUN,
                                           ALIGNFUNARGS = ALIGNFUNARGS,
                                           FUN = FUN, FUNARGS = FUNARGS)
    suppressWarnings(
        res <- apply(res, c(1, 2), aggregationFun, na.rm = TRUE)
    )
    ## Ensure diagonal is always TRUE to not drop any features!
    res[cbind(1:nr, 1:nr)] <- 1
    as.factor(do.call(groupFun,
                      args = c(list(res, threshold = threshold),
                               list(...))))
}
