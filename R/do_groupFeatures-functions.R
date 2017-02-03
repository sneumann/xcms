## Alignment functions.
#' @include functions-Params.R

##' @title Core API function for feature density based feature alignment
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
##' @note The default settings might not be appropriate for all LC/GC-MS setups,
##' especially the \code{bw} and \code{binSize} parameter should be adjusted
##' accordingly.
##' 
##' @param features A \code{matrix} or \code{data.frame} with the mz values and
##' retention times of the identified features in all samples of an experiment.
##' Required columns are \code{"mz"}, \code{"rt"} and \code{"sample"}. The latter
##' should contain \code{numeric} values representing the index of the sample in
##' which the feature was found.
##'
##' @inheritParams groupFeatures-density
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

    sampleGroups <- as.character(sampleGroups)
    sampleGroupNames <- unique(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)

    ## Check that sample groups matches with sample column.
    if (max(features[, "sample"]) > length(sampleGroups))
        stop("Sample indices in 'features' are larger than there are sample",
             " groups specified with 'sampleGroups'!")
    
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

    densFrom <- rtRange[1] - 3 * bw
    densTo <- rtRange[2] + 3 * bw
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
        den <- density(curMat[, "rt"], bw = bw, from = densFrom, to = densTo)
        maxden <- max(den$y)
        deny <- den$y
        ## gmat <- matrix(nrow = 5, ncol = 2 + gcount)
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

## Just to check if we could squeeze a little bit more out using parallel
## processing...
do_groupFeatures_density_par <- function(features, sampleGroups,
                                         bw = 30, minFraction = 0.5,
                                         minSamples = 1, binSize = 0.25,
                                         maxFeatures = 50) {
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

    sampleGroups <- as.character(sampleGroups)
    sampleGroupNames <- unique(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)

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

    ## Create the list of feature data subsets.
    ftsL <- vector("list", length(mass))
    for (i in seq_len(length(mass) - 2)) {
        startIdx <- masspos[i]
        endIdx <- masspos[i + 2] - 1
        ftsL[[i]] <- cbind(features[startIdx:endIdx, , drop = FALSE],
                         idx = startIdx:endIdx)
    }
    ftsL <- ftsL[lengths(ftsL) > 0]
    ## Here we can run bplapply:
    res <- bplapply(ftsL, function(z, rtr, bw, maxF, sampleGrps,
                                   sampleGroupTbl, minFr, minSmpls,
                                   sampleGroupNms, featureOrdr) {
        den <- density(z[, "rt"], bw = bw, from = rtr[1] - 3 * bw,
                       to = rtr[2] + 3 * bw)
        maxden <- max(den$y)
        deny <- den$y
        snum <- 0
        tmpL <- vector("list", maxF)
        tmpL2 <- tmpL
        while (deny[maxy <- which.max(deny)] > maxden / 20 && snum < maxF) {
            grange <- xcms:::descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(z[,"rt"] >= den$x[grange[1]] &
                          z[,"rt"] <= den$x[grange[2]])
            ## Determine the sample group of the samples in which the features
            ## were detected and check if they correspond to the required limits.
            tt <- table(sampleGrps[unique(z[gidx, "sample"])])
            if (!any(tt / sampleGroupTbl[names(tt)] >= minFr &
                     tt >= minSmpls))
                next
            snum <- snum + 1
            gcount <- rep(0, length(sampleGroupNms))
            names(gcount) <- sampleGroupNms
            gcount[names(tt)] <- as.numeric(tt)

            tmpL[[snum]] <- c(median(z[gidx, "mz"]),
                              range(z[gidx, "mz"]),
                              median(z[gidx, "rt"]),
                              range(z[gidx, "rt"]),
                              length(gidx),
                              gcount)
            tmpL2[[snum]] <- sort(featureOrdr[z[, "idx"][gidx]])
        }
        tmpL <- tmpL[lengths(tmpL) > 0]
        tmpL2 <- tmpL2[lengths(tmpL2) > 0]
        if (length(tmpL))
            return(list(grps = do.call(rbind, tmpL), idx = tmpL2))
    }, rtr = rtRange, bw = bw, maxF = maxFeatures, sampleGrps = sampleGroups,
    sampleGroupTbl = sampleGroupTable, minFr = minFraction,
    minSmpls = minSamples, sampleGroupNms = sampleGroupNames,
    featureOrdr = featureOrder)

    res <- res[lengths(res) > 0]
    ## Now we have to process that list of results.
    groupmat <- do.call(rbind, lapply(res, function(z) z[["grps"]]))
    groupidx <- unlist(lapply(res, function(z) z[["idx"]]), recursive = FALSE)
    
    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", sampleGroupNames)

    ## groupmat <- groupmat[seq_len(num), , drop = FALSE]
    ## groupindex <- groupindex[seq_len(num)]

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[, (match("npeaks",
                                         colnames(groupmat))+1):ncol(groupmat),
                                drop = FALSE])
    uorder <- order(-numsamp, groupmat[, "npeaks"])

    uindex <- rectUnique(groupmat[, c("mzmin","mzmax","rtmin","rtmax"),
                                  drop = FALSE],
                         uorder)

    return(list(featureGroups = groupmat[uindex, , drop = FALSE],
                featureIndex = groupidx[uindex]))
}

##' @title Core API function for feature alignment using mzClust
##'
##' @description The \code{do_groupFeatures_mzClust} function performs high
##' resolution alignment on single spectra samples.
##' 
##' @inheritParams groupFeatures-density
##' @inheritParams do_groupFeatures_density
##' @inheritParams groupFeatures-mzClust
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
##' \item{"rtmed"}{always \code{-1}.}
##' \item{"rtmin"}{always \code{-1}.}
##' \item{"rtmax"}{always \code{-1}.}
##' \item{"npeaks"}{the total number of features assigned to the feature group.
##' Note that this number can be larger than the total number of samples, since
##' multiple features from the same sample could be assigned to a group.}
##' }
##' \code{"featureIndex"} is a \code{list} with the indices of all features in a
##' feature group in the \code{features} input matrix.
##'
##' @family core feature alignment algorithms
##'
##' @references Saira A. Kazmi, Samiran Ghosh, Dong-Guk Shin, Dennis W. Hill
##' and David F. Grant\cr \emph{Alignment of high resolution mass spectra:
##' development of a heuristic approach for metabolomics}.\cr Metabolomics,
##' Vol. 2, No. 2, 75-83 (2006)
do_groupFeatures_mzClust <- function(features, sampleGroups, ppm = 20,
                                     absMz = 0, minFraction = 0.5,
                                     minSamples = 1) {
    if (missing(sampleGroups))
        stop("Parameter 'sampleGroups' is missing! This should be a vector of ",
             "length equal to the number of samples specifying the group ",
             "assignment of the samples.")
    if (missing(features))
        stop("Parameter 'peaks' is missing!")
    if (!is.matrix(features) | is.data.frame(features))
        stop("Peaks has to be a 'matrix' or a 'data.frame'!")
    ## Check that we've got all required columns
    .reqCols <- c("mz", "sample")
    if (!all(.reqCols %in% colnames(features)))
        stop("Required columns ",
             paste0("'", .reqCols[!.reqCols %in% colnames(features)],"'",
                    collapse = ", "), " not found in 'features' parameter")
    if (!is.factor(sampleGroups))
        sampleGroups <- factor(sampleGroups, levels = unique(sampleGroups))
    sampleGroupNames <- levels(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)
    ##sampleGroups <- as.numeric(sampleGroups)
    
    ## Check that sample groups matches with sample column.
    if (max(features[, "sample"]) > length(sampleGroups))
        stop("Sample indices in 'features' are larger than there are sample",
             " groups specified with 'sampleGroups'!")

    ##features <- features[, .reqCols, drop = FALSE]
    grps <- mzClustGeneric(features[, .reqCols, drop = FALSE],
                           sampclass = sampleGroups,
                           mzppm = ppm,
                           mzabs = absMz,
                           minsamp = minSamples,
                           minfrac = minFraction)
    grpmat <- grps$mat
    if (is.null(nrow(grpmat))) {
        matColNames <- names(grpmat)
        grpmat <- matrix(grpmat, ncol = length(grpmat), byrow = FALSE)
        colnames(grpmat) <- matColNames
    }
    rts <- rep(-1, nrow(grpmat))
    cns <- colnames(grpmat)
    grpmat <- cbind(grpmat[, 1:3, drop = FALSE], rts, rts, rts,
                    grpmat[, 4:ncol(grpmat), drop = FALSE])
    colnames(grpmat) <- c(cns[1:3], c("rtmed", "rtmin", "rtmax"),
                          cns[4:length(cns)])
    return(list(featureGroups = grpmat, featureIndex = grps$idx))    
}

##' @title Core API function for feature alignment using a nearest neighbor approach
##'
##' @description The \code{do_groupFeatures_nearest} function groups features
##' across samples by creating a master feature list and assigning corresponding
##' features from all samples to each feature group. The method is inspired by
##' the alignment algorithm of mzMine [Katajamaa 2006].
##' 
##' @inheritParams do_groupFeatures_density
##' @inheritParams groupFeatures-nearest
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
##' \item{"npeaks"}{the total number of features assigned to the feature group.}
##' }
##' \code{"featureIndex"} is a \code{list} with the indices of all features in a
##' feature group in the \code{features} input matrix.
##'
##' @family core feature alignment algorithms
##'
##' @references Katajamaa M, Miettinen J, Oresic M: MZmine: Toolbox for
##' processing and visualization of mass spectrometry based molecular profile
##' data. \emph{Bioinformatics} 2006, 22:634-636. 
do_groupFeatures_nearest <- function(features, sampleGroups, mzVsRtBalance = 10,
                                     absMz = 0.2, absRt = 15, kNN = 10) {
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
    if (!is.factor(sampleGroups))
        sampleGroups <- factor(sampleGroups, levels = unique(sampleGroups))
    sampleGroupNames <- levels(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)

    ## sampleGroups == classlabel
    ## nSampleGroups == gcount
    ## features == peakmat

    features <- features[, .reqCols, drop = FALSE]
    
    parameters <- list(mzVsRTBalance = mzVsRtBalance, mzcheck = absMz,
                       rtcheck = absRt, knn = kNN)

    ptable <- table(features[,"sample"])
    pord <- ptable[order(ptable, decreasing = TRUE)]
    sid <- as.numeric(names(pord))
    pn <- as.numeric(pord)

    ## environment - probably not a good idea for parallel processing - we
    ## would like to have data copying there (or better just provide the data
    ## chunk it needs to process).
    mplenv <- new.env(parent = .GlobalEnv)
    mplenv$mplist <- matrix(0, pn[1], length(sid))
    mplenv$mplist[, sid[1]] <- which(features[,"sample"] == sid[1])
    mplenv$mplistmean <- data.frame(features[which(features[,"sample"] == sid[1]),
                                            c("mz", "rt")])
    mplenv$peakmat <- features
    assign("peakmat", features, envir = mplenv)  ## double assignment?

    sapply(sid[2:length(sid)], function(sample, mplenv){
        message("Processing sample number ", sample, " ... ", appendLF = FALSE)
        ## require(parallel)
        ## cl <- makeCluster(getOption("cl.cores", nSlaves))
        ## clusterEvalQ(cl, library(RANN))
        ## parSapply(cl, 2:length(samples), function(sample,mplenv, object){
        ## might slightly improve on this for loop.
        ## Calculating for each row (peak) the mean mz or rt for features
        ## assigned yet to this feature group.
        for (mml in seq(mplenv$mplist[,1])) {
            mplenv$mplistmean[mml, "mz"] <-
                mean(mplenv$peakmat[mplenv$mplist[mml, ], "mz"])
            mplenv$mplistmean[mml, "rt"] <-
                mean(mplenv$peakmat[mplenv$mplist[mml, ], "rt"])
        }

        mplenv$peakIdxList <- data.frame(
            peakidx = which(mplenv$peakmat[, "sample"] == sample),
            isJoinedPeak = FALSE
        )
        if (length(mplenv$peakIdxList$peakidx) == 0)
            message("Warning: No peaks in sample number ", sample)
        
        ## this really doesn't take a long time not worth parallel version here.
        ## but make an apply loop now faster even with rearranging the data :D : PB
        scoreList <- sapply(mplenv$peakIdxList$peakidx,
                            function(currPeak, para, mplenv){
                                xcms:::patternVsRowScore(currPeak, para, mplenv)
                            }, parameters, mplenv, simplify = FALSE)
        scoreList <- do.call(rbind, scoreList)

        ## Browse scores in order of descending goodness-of-fit
        scoreListcurr <- scoreList[order(scoreList[, "score"]), ]
        if (nrow(scoreListcurr) > 0) {
            for (scoreIter in 1:nrow(scoreListcurr)) {

                iterPeak <- scoreListcurr[scoreIter, "peak"]
                iterRow <- scoreListcurr[scoreIter, "mpListRow"]

                ## Check if master list row is already assigned with peak
                if (scoreListcurr[scoreIter, "isJoinedRow"] == TRUE)
                    next

                ## Check if peak is already assigned to some master list row
                if (scoreListcurr[scoreIter, "isJoinedPeak"] == TRUE)
                    next

                ##  Check if score good enough
                ## Assign peak to master peak list row
                mplenv$mplist[iterRow, sample] <- iterPeak

                ## Mark peak as joined
                setTrue <- which(scoreListcurr[, "mpListRow"] == iterRow)
                scoreListcurr[setTrue, "isJoinedRow"] <- TRUE
                setTrue <- which(scoreListcurr[, "peak"] == iterPeak)
                scoreListcurr[setTrue, "isJoinedPeak"] <- TRUE
                mplenv$peakIdxList[which(mplenv$peakIdxList$peakidx == iterPeak),
                                   "isJoinedPeak"] <- TRUE
            }
        }
        notJoinedPeaks <- mplenv$peakIdxList[which(mplenv$peakIdxList$isJoinedPeak == FALSE), "peakidx"]
        
        for (notJoinedPeak in notJoinedPeaks) {
            mplenv$mplist <- rbind(mplenv$mplist,
                                   matrix(0, 1, dim(mplenv$mplist)[2]))
            mplenv$mplist[length(mplenv$mplist[,1]), sample] <- notJoinedPeak
        }

        ## Clear "Joined" information from all master peaklist rows
        rm(list = "peakIdxList", envir = mplenv)
        message("OK")
    }, mplenv)

    groupmat <- matrix( 0, nrow(mplenv$mplist),  7 + nSampleGroups)
    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", sampleGroupNames)
    groupindex <- vector("list", nrow(mplenv$mplist))
    ## Variable to count samples for a feature
    sampCounts <- rep_len(0, nSampleGroups)
    names(sampCounts) <- sampleGroupNames
    ## gcount <- integer(nSampleGroups)
    ## Can we vectorize that below somehow?
    for (i in 1:nrow(mplenv$mplist)) {
        groupmat[i, "mzmed"] <- median(features[mplenv$mplist[i, ], "mz"])
        groupmat[i, c("mzmin", "mzmax")] <- range(features[mplenv$mplist[i, ], "mz"])
        groupmat[i, "rtmed"] <- median(features[mplenv$mplist[i, ], "rt"])
        groupmat[i, c("rtmin", "rtmax")] <- range(features[mplenv$mplist[i, ], "rt"])

        groupmat[i, "npeaks"] <- length(which(features[mplenv$mplist[i, ]] > 0))

        ## Now summarizing the number of samples in which the peak was identified
        sampCounts[] <- 0
        tbl <- table(sampleGroups[features[mplenv$mplist[i, ], "sample"]])
        sampCounts[names(tbl)] <- as.numeric(tbl)
        groupmat[i, 7 + seq_len(nSampleGroups)] <- sampCounts
        ## gnum <- sampleGroups[unique(features[mplenv$mplist[i, ], "sample"])]
        ## for (j in seq(along = gcount))
        ##     gcount[j] <- sum(gnum == j)
        ## groupmat[i, 7 + seq(along = gcount)] <- gcount
        groupindex[[i]] <- mplenv$mplist[i, (which(mplenv$mplist[i,]>0))]
    }
    
    return(list(featureGroups = groupmat, featureIndex = groupindex))    
}
