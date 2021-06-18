#' @title Grouping of LC-MS features
#'
#' @name feature-grouping
#'
#' @description
#'
#' After correspondence analysis ([xcms::groupChromPeaks()]) the identified
#' LC-MS features can be further grouped based on different criteria into
#' *feature groups* which contain ideally features representing
#' ions (adducts, isotopes) of the same compound/metabolite.
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
#'   [EicCorrelationParam()].
#'
#' An ideal workflow grouping features should sequentially perform the above
#' methods (in the listed order).
#'
#' Feature groups can be accessed with the `featureGroups` function.
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


#' @title Group features based on similar retention times
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

#' @title Group features based on similarity of abundances across samples
#'
#' @name groupFeatures-abundance-correlation
#'
#' @description
#'
#' This method groups features based on similarity of abundances (i.e.
#' *feature values*) across samples. See also [AbundanceSimilarityParam()] for
#' additional information and details.
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
## #'
## #' @examples
## #'
## #' ## Load test data set from xcms
## #' library(xcms)
## #' data(faahko_sub)
## #' ## Update the path to the files for the local system
## #' dirname(faahko_sub) <- system.file("cdf/KO/", package = "faahKO")
## #'
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

## #' @title Group features based on correlation of extracted ion chromatograms
## #'
## #' @name groupFeatures-eic-correlation
## #'
## #' @description
## #'
## #' Group features based on correlation of their extracted ion chromatograms.
## #' This correlation is performed separately for each sample with the correlation
## #' coefficients being aggregated across samples for the final comparison with
## #' parameter `threshold` (the 75% quantile of the per-sample correlation values
## #' is used for the comparison with `threshold`).
## #'
## #' This feature grouping should be called **after** an initial feature
## #' grouping by retention time (see [SimilarRtimeParam()]). The feature groups
## #' defined in columns `"feature_group"` of `featureDefinitions(object)` (for
## #' features matching `msLevel`) will be used and refined by this method.
## #' Features with a value of `NA` in `featureDefinitions(object)$feature_group`
## #' will be skipped/not considered for feature grouping.
## #'
## #' While being possible to be performed on the full data set without prior
## #' feature grouping , this is not suggested for the following reasons: I) the
## #' selection of the top `n` samples with the highest signal for the
## #' *feature group* will be biased by very abundant compounds as this is
## #' performed on the full data set (i.e. the samples with the highest overall
## #' intensities are used for correlation of all features) and II) it is
## #' computationally much more expensive because a pairwise correlation between
## #' all features has to be performed.
## #'
## #' It is also suggested to perform the correlation on a subset of samples
## #' per feature with the highest intensities of the peaks (for that feature)
## #' although it would also be possible to run the correlation on all samples by
## #' setting `n` equal to the total number of samples in the data set. EIC
## #' correlation should however be performed ideally on samples in which the
## #' original compound is highly abundant to avoid correlation of missing values
## #' or noisy peak shapes as much as possible.
## #'
## #' By default also the signal which is outside identified chromatographic peaks
## #' is excluded from the correlation  (parameter `clean`).
## #'
## #' @param clean `logical(1)` whether the correlation should be performed only
## #'     on the signals within the identified chromatographic peaks
## #'     (`clean = TRUE`, default) or all the signal from the extracted ion
## #'     chromatogram.
## #'
## #' @param inclusive `logical(1)` which grouping algorithm should be used: one that
## #'     creates small groups of highly correlated features (`inclusive = FALSE`, the
## #'     default) or whether features should be grouped that have at least one
## #'     correlation with any other member of the group in common
## #'     (`inclusive = TRUE`). See [groupByCorrelation()] for more information.
## #'
## #' @param msLevel `integer(1)` defining the MS level on which the features
## #'     should be grouped.
## #'
## #' @param n `numeric(1)` defining the total number of samples per feature group
## #'     on which this correlation should be performed. This value is rounded up
## #'     to the next larger integer value.
## #'
## #' @param object [XCMSnExp()] object containing also correspondence results.
## #'
## #' @param param `EicCorrelationParam` object with the settings for the method.
## #'
## #' @param threshold `numeric(1)` with the minimal required correlation
## #'     coefficient to group featues.
## #'
## #' @param value `character(1)` defining whether samples should be grouped based
## #'     on the sum of the maximal peak intensity (`value = "maxo"`, the default)
## #'     or the integrated peak area (`value = "into"`) for a feature.
## #'
## #' @return input `XCMSnExp` with feature groups added (i.e. in column
## #'     `"feature_group"` of its `featureDefinitions` data frame.
## #'
## #' @family feature grouping methods
## #'
## #' @seealso feature-grouping for a general overview.
## #'
## #' @rdname groupFeatures-eic-correlation
## #'
## #' @exportClass EicCorrelationParam
## #'
## #' @author Johannes Rainer
## #'
## #' @examples
## #'
## #' library(MsFeatures)
## #' ## Load a test data set with detected peaks
## #' data(faahko_sub)
## #' ## Update the path to the files for the local system
## #' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
## #'
## #' ## Disable parallel processing for this example
## #' register(SerialParam())
## #'
## #' ## Performing a feature grouping based on EIC correlation on a single
## #' ## sample
## #' xodg_grp <- groupFeatures(faahko_sub, param = EicCorrelationParam(n = 1))
## #'
## #' table(featureDefinitions(xodg_grp)$feature_group)
## #'
## #' ## Usually it is better to perform this correlation on pre-grouped features
## #' ## e.g. based on similar retention time.
## #' xodg_grp <- groupFeatures(faahko_sub, param = SimilarRtimeParam(diffRt = 4))
## #' xodg_grp <- groupFeatures(xodg_grp, param = EicCorrelationParam(n = 1))
## #'
## #' table(featureDefinitions(xodg_grp)$feature_group)
## NULL

## setClass("EicCorrelationParam",
##          slots = c(threshold = "numeric",
##                    n = "numeric",
##                    clean = "logical",
##                    value = "character",
##                    inclusive = "logical"),
##          contains = "Param",
##          prototype = prototype(
##              threshold = 0.9,
##              n = 1,
##              clean = TRUE,
##              value = "maxo",
##              inclusive = FALSE
##          ),
##          validity = function(object) {
##              msg <- NULL
##              if (length(object@threshold) != 1 || object@threshold < 0)
##                  msg <- "'threshold' has to be a positive numeric of length 1"
##              if (length(object@n) != 1 || object@n < 0)
##                  msg <- c(msg, "'n' has to be a positive numeric of length 1")
##              if (length(object@clean) != 1)
##                  msg <- c(msg, "'clean' has to a logical of length 1")
##              if (length(object@value) != 1 && !(object@value %in%
##                                                 c("maxo", "into")))
##                  msg <- c(msg, "'value' has to be either \"maxo\" or \"into\"")
##              msg
##          })

## #' @rdname groupFeatures-eic-correlation
## #'
## #' @export
## EicCorrelationParam <- function(threshold = 0.9, n = 1, clean = TRUE,
##                                 value = c("maxo", "into"), inclusive = FALSE) {
##     value = match.arg(value)
##     new("EicCorrelationParam", threshold = threshold, n = n, clean = clean,
##         value = value, inclusive = FALSE)
## }

## #' @rdname groupFeatures-eic-correlation
## #'
## #' @importFrom utils txtProgressBar setTxtProgressBar
## #'
## #' @importMethodsFrom MSnbase fileNames
## #'
## #' @importMethodsFrom xcms filterFile featureValues removeIntensity
## #'
## #' @importFrom xcms featureChromatograms
## setMethod(
##     "groupFeatures",
##     signature(object = "XCMSnExp", param = "EicCorrelationParam"),
##     function(object, param, msLevel = 1L) {
##         if (!hasFeatures(object))
##             stop("No feature definitions present. Please run ",
##                  "first 'groupChromPeaks'")
##         if (length(msLevel) > 1)
##             stop("Currently only grouping of features from a single MS level",
##                  " is supported.")
##         n <- ceiling(param@n)
##         nc <- length(fileNames(object))
##         if (n > nc)
##             stop("'n' should be smaller or equal than the number of ",
##                  "samples (", nc, ")")
##         if (any(colnames(featureDefinitions(object)) == "ms_level"))
##             is_msLevel <- featureDefinitions(object)$ms_level == msLevel
##         else is_msLevel <- rep(TRUE, nrow(featureDefinitions(object)))
##         if (any(colnames(featureDefinitions(object)) == "feature_group")) {
##             f <- featureDefinitions(object)$feature_group
##             f_new <- as.character(f)
##         } else {
##             f <- rep("FG", nrow(featureDefinitions(object)))
##             f_new <- rep(NA_character_, length(f))
##         }
##         f[!is_msLevel] <- NA
##         if (is.factor(f))
##             f <- droplevels(f)
##         else
##             f <- factor(f, levels = unique(f))
##         fvals <- featureValues(object, method = "maxint", value = param@value)
##         ffun <- function(z, na.rm = TRUE)
##             quantile(z, probs = 0.75, na.rm = na.rm)
##         pb <- txtProgressBar(min = 0, max  = nrow(fvals), style = 3)
##         setTxtProgressBar(pb, 0)
##         counter <- 0
##         for (fg in levels(f)) {
##             idx <- which(f == fg)
##             idxl <- length(idx)
##             counter <- counter + idxl
##             if (idxl > 1) {
##                 vals <- apply(fvals[idx, ], MARGIN = 2, sum, na.rm = TRUE)
##                 sample_idx <- order(vals, decreasing = TRUE)[seq_len(n)]
##                 obj_sub <- filterFile(object, sample_idx, keepFeatures = TRUE)
##                 ## Can happen that some of the features are not present in the
##                 ## subsetted object. Will put them into their own individual
##                 ## groups later.
##                 idx_miss <- which(!rownames(fvals)[idx] %in%
##                                   rownames(featureDefinitions(obj_sub)))
##                 if (length(idx_miss)) {
##                     tmp <- idx[idx_miss]
##                     idx <- idx[-idx_miss]
##                     idx_miss <- tmp
##                 }
##                 if (length(idx) > 1) {
##                     eics <- featureChromatograms(
##                         obj_sub, features = rownames(fvals)[idx], filled = TRUE)
##                     if (param@clean)
##                         eics <- removeIntensity(eics,
##                                                 which = "outside_chromPeak")
##                     res <- groupEicCorrelation(
##                         as(eics, "MChromatograms"), aggregationFun = ffun,
##                         threshold = param@threshold, tolerance = 0,
##                         inclusive = param@inclusive)
##                 } else res <- factor(1)
##                 f_new[idx] <- paste0(fg, ".", MsFeatures:::.format_id(res))
##                 if (length(idx_miss))
##                     f_new[idx_miss] <- paste0(
##                         fg, ".", MsFeatures:::.format_id(
##                                                   seq((length(levels(res)) + 1),
##                                                       length.out = length(
##                                                           idx_miss))))
##             } else
##                 f_new[idx] <- paste0(fg, ".001")
##             setTxtProgressBar(pb, counter)
##         }
##         setTxtProgressBar(pb, nrow(fvals))
##         close(pb)
##         featureDefinitions(object)$feature_group <- f_new
##         xph <- new("XProcessHistory", param = param, date = date(),
##                    type = .PROCSTEP.FEATURE.GROUPING,
##                    fileIndex = 1:length(fileNames(object)),
##                    msLevel = as.integer(msLevel))
##         object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
##         validObject(object)
##         object
##     })
