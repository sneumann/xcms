## Alignment functions.
#' @include functions-Params.R

##' @title Align features across samples based on time dimension feature densities
##'
##' @description The \code{do_groupFeatures_density} function performs feature
##' alignment based on the density (distribution) of features, found in different
##' samples, along the retention time axis in slices of overlapping mz ranges.
##'
##' @details For overlapping slices along the mz dimension, the function
##' calculates the density distribution of identified features along the
##' retention time axis and groups features from the same or different samples
##' that are close to each other. See [Smith 2006] for more details.
##'
##' @param features A \code{matrix} or \code{data.frame} with the mz values and
##' retention times of the identified features in all samples of an experiment.
##' Required columns are \code{"mz"}, \code{"rt"}.
##'
##' @param sampleGroups A vector of the same length than samples defining the
##' sample group assignments.
##'
##' @param bw numeric(1) defining the bandwidth (standard deviation ot the
##' smoothing kernel) to be used. This argument is passed to the
##' \code{\link{density}} method.
##'
##' @param minFraction numeric(1) defining the minimum fraction of samples in at
##' least one sample group in which the features have to be present to be
##' considered as a feature group.
##'
##' @param minSamples numeric(1) with the minimum number of samples in at least
##' one sample group in which the features have to be detected to be considered
##' as a feature group.
##'
##' @param binSize numeric(1) defining the size of the overlapping slices in mz
##' dimension.
##'
##' @param maxFeatures numeric(1) with the maximum number of feature groups to
##' be identified in a single mz slice.
##'
##' @return A \code{list} with elements \code{"featureGroups"} and
##' \code{"featureIndex"}. \code{"featureGroups"} is a \code{matrix}, each row
##' representing an aligned feature group and with columns:
##' \describe{
##' \item{"mzmed"}{median of the features' apex mz values.}
##' \item{"mzmin"}{smallest mz value of all features' apex within the feature
##' group.}
##' \item{"mzmax"}{largest mz value of all features' apex within the feature
##' group.}
##' \item{"rtmed"}{the median of the features' retention times.}
##' \item{"rtmin"}{the smallest retention time of the features in the group.}
##' \item{"rtmax"}{the largest retention time of the features in the group.}
##' \item{"npeaks"}{the total number of features assigned to the feature group.
##' Note that this number can be larger than the total number of samples, since
##' multiple features from the same sample could be assigned to a group.}
##' }
##' \code{"featureIndex"} is a \code{list} with the indices of all features in a
##' feature group in the \code{features} input matrix.
##'
##' @family core feature alignment algorithms
##'
##' @references
##' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
##' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
##' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
##' \emph{Anal. Chem.} 2006, 78:779-787.
##'
##' @author Colin Smith, Johannes Rainer
##'
##' @examples
##' ## Load the test data set
##' library(faahKO)
##' data(faahko)
##'
##' ## Extract the matrix with the identified features from the xcmsSet:
##' fts <- peaks(faahko)
##'
##' ## Perform the feature alignment with default settings:
##' res <- do_groupFeatures_density(fts, sampleGroups = sampclass(faahko))
##'
##' ## The feature groups:
##' head(res$featureGroups)
##'
##' ## The assignment of features from the input matrix to the feature groups
##' head(res$featureIndex)
do_groupFeatures_density <- function(features, sampleGroups,
                                     bw = 30, minFraction = 0.5, minSamples = 1,
                                     binSize = 0.25, maxFeatures = 50) {
    if (missing(sampleGroups))
        stop("Parameter 'sampleGroups' is missing! This should be a vector of ",
             "length equal to the number of samples specifying the group ",
             "assignment of the samples.")
    if (missing(features))
        stop("Parameter 'peaks' is missing!")
    if (!is.matrix(features) | is.data.frame(features))
        stop("Peaks has to be a 'matrix' or a 'data.frame'!")
    ## Check that we've got all required columns
    .reqCols <- c("mz", "rt", "sample")
    if (!all(.reqCols %in% colnames(features)))
        stop("Required columns ",
             paste0("'", .reqCols[!.reqCols %in% colnames(features)],"'",
                    collapse = ", "), " not found in 'features' parameter")

    ## samples <- sampnames(object)
    ## classlabel <- sampclass(object)
    ## classnames <- as.character(unique(sampclass(object)))
    ## classlabel <- as.vector(unclass(classlabel))
    ## classnum <- table(classlabel)

    sampleGroups <- as.character(sampleGroups)
    sampleGroupNames <- unique(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)

    ## peakmat <- peaks(object)
    ## porder <- order(peakmat[,"mz"])
    ## peakmat <- peakmat[porder,, drop=FALSE]
    ## rownames(peakmat) <- NULL
    ## retrange <- range(peakmat[,"rt"])

    ## Order features matrix by mz
    featureOrder <- order(features[, "mz"])
    features <- features[featureOrder, .reqCols, drop = FALSE]
    rownames(features) <- NULL
    rtRange <- range(features[, "rt"])

    ## Define the mass slices and the index in the features matrix with an mz
    ## value >= mass[i].
    mass <- seq(features[1, "mz"], features[nrow(features), "mz"] + binSize,
                by = binSize / 2)
    masspos <- findEqualGreaterM(features[,"mz"], mass)

    groupmat <- matrix(nrow = 512, ncol = 7 + nSampleGroups)
    groupindex <- vector("list", 512)

    endIdx <- 0
    num <- 0
    gcount <- integer(nSampleGroups)
    for (i in seq_len(length(mass)-2)) {
        ## That's identifying overlapping mz slices.
        startIdx <- masspos[i]
        endIdx <- masspos[i + 2] - 1
        if (endIdx - startIdx < 0)
            next
        curMat <- features[startIdx:endIdx, , drop = FALSE]
        den <- density(curMat[,"rt"], bw, from = rtRange[1] - 3 * bw,
                       to = rtRange[2] + 3 * bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2 + gcount)
        snum <- 0
        ## What's that 20 there?
        while (deny[maxy <- which.max(deny)] > maxden / 20 && snum < maxFeatures) {
            grange <- descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(curMat[,"rt"] >= den$x[grange[1]] &
                          curMat[,"rt"] <= den$x[grange[2]])
            ## Determine the sample group of the samples in which the features
            ## were detected and check if they correspond to the required limits.
            tt <- table(sampleGroups[unique(curMat[gidx, "sample"])])
            if (!any(tt / sampleGroupTable[names(tt)] >= minFraction &
                     tt >= minSamples))
                next
            snum <- snum + 1
            num <- num + 1
            ## Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat,
                                  matrix(nrow = nrow(groupmat),
                                         ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            gcount <- rep(0, length(sampleGroupNames))
            names(gcount) <- sampleGroupNames
            gcount[names(tt)] <- as.numeric(tt)
            groupmat[num, 1] <- median(curMat[gidx, "mz"])
            groupmat[num, 2:3] <- range(curMat[gidx, "mz"])
            groupmat[num, 4] <- median(curMat[gidx, "rt"])
            groupmat[num, 5:6] <- range(curMat[gidx, "rt"])
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7 + seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(featureOrder[(startIdx:endIdx)[gidx]])
        }
    }

    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", sampleGroupNames)

    groupmat <- groupmat[seq_len(num), , drop = FALSE]
    groupindex <- groupindex[seq_len(num)]

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[, (match("npeaks",
                                         colnames(groupmat))+1):ncol(groupmat),
                                drop = FALSE])
    uorder <- order(-numsamp, groupmat[, "npeaks"])

    uindex <- rectUnique(groupmat[, c("mzmin","mzmax","rtmin","rtmax"),
                                  drop = FALSE],
                         uorder)

    return(list(featureGroups = groupmat[uindex, , drop = FALSE],
                featureIndex = groupindex[uindex]))
}

