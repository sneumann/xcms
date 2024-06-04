## Correspondence functions.
#' @include functions-Params.R

#' @title Core API function for peak density based chromatographic peak
#' grouping
#'
#' @description
#'
#' The `do_groupChromPeaks_density` function performs chromatographic peak
#' grouping based on the density (distribution) of peaks, found in different
#' samples, along the retention time axis in slices of overlapping m/z ranges.
#' By default (with parameter `ppm = 0`) these m/z ranges have all the same
#' (constant) size (depending on parameter `binSize`). For values of `ppm`
#' larger than 0 the m/z bins (ranges or slices) will have increasing sizes
#' depending on the m/z value. This better models the m/z-dependent
#' measurement error/precision seen on some MS instruments.
#'
#' @details
#'
#' For overlapping slices along the mz dimension, the function
#' calculates the density distribution of identified peaks along the
#' retention time axis and groups peaks from the same or different samples
#' that are close to each other. See (Smith 2006) for more details.
#'
#' @note
#'
#' The default settings might not be appropriate for all LC/GC-MS setups,
#' especially the `bw` and `binSize` parameter should be adjusted
#' accordingly.
#'
#' @param peaks A `matrix` or `data.frame` with the mz values and
#'     retention times of the identified chromatographic peaks in all samples
#'     of an experiment. Required columns are `"mz"`, `"rt"` and
#'     `"sample"`. The latter should contain `numeric` values representing
#'     the index of the sample in which the peak was found.
#'
#' @param index An optional `integer` providing the indices of the peaks in the
#'     original peak matrix.
#'
#' @inheritParams groupChromPeaks
#'
#' @param sleep `numeric(1)` defining the time to *sleep* between
#'     iterations and plot the result from the current iteration.
#'
#' @return
#'
#' A `data.frame`, each row representing a (mz-rt) feature (i.e. a peak group)
#' with columns:
#'
#' - `"mzmed"`: median of the peaks' apex mz values.
#' - `"mzmin"`: smallest mz value of all peaks' apex within the feature.
#' - `"mzmax"`:largest mz value of all peaks' apex within the feature.
#' - `"rtmed"`: the median of the peaks' retention times.
#' - `"rtmin"`: the smallest retention time of the peaks in the group.
#' - `"rtmax"`: the largest retention time of the peaks in the group.
#' - `"npeaks"`: the total number of peaks assigned to the feature.
#' - `"peakidx"`: a `list` with the indices of all peaks in a feature in the
#'   `peaks` input matrix.
#'
#' Note that this number can be larger than the total number of samples, since
#' multiple peaks from the same sample could be assigned to a feature.
#'
#' @references
#'
#' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#' Anal. Chem. 2006, 78:779-787.
#'
#' @author Colin Smith, Johannes Rainer
#'
#' @family core peak grouping algorithms
#'
#' @md
#'
#' @examples
#' ## Load the test file
#' library(xcms)
#' library(MsExperiment)
#' faahko_sub <- loadXcmsData("faahko_sub2")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## Extract the matrix with the identified peaks from the xcmsSet:
#' pks <- chromPeaks(faahko_sub)
#'
#' ## Perform the peak grouping with default settings:
#' res <- do_groupChromPeaks_density(pks, sampleGroups = rep(1, 3))
#'
#' ## The feature definitions:
#' head(res)
do_groupChromPeaks_density <- function(peaks, sampleGroups,
                                       bw = 30, minFraction = 0.5,
                                       minSamples = 1, binSize = 0.25,
                                       maxFeatures = 50, sleep = 0,
                                       index = seq_len(nrow(peaks)),
                                       ppm = 0) {
    if (missing(sampleGroups))
        stop("Parameter 'sampleGroups' is missing! This should be a vector of ",
             "length equal to the number of samples specifying the group ",
             "assignment of the samples.")
    if (missing(peaks))
        stop("Parameter 'peaks' is missing!")
    if (!(is.matrix(peaks) | is.data.frame(peaks)))
        stop("'peaks' has to be a 'matrix' or a 'data.frame'!")
    ## Check that we've got all required columns
    .reqCols <- c("mz", "rt", "sample")
    if (sleep > 0)
        .reqCols <- c(.reqCols, "into")
    if (!all(.reqCols %in% colnames(peaks)))
        stop("Required columns ",
             paste0("'", .reqCols[!.reqCols %in% colnames(peaks)],"'",
                    collapse = ", "), " not found in 'peaks' parameter")

    ## With a `factor` we also support excluding samples/groups, i.e. samples
    ## with an NA are not considered in the feature definition.
    if (!is.factor(sampleGroups))
        sampleGroups <- factor(sampleGroups)
    sampleGroupNames <- levels(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)

    ## Check that sample groups matches with sample column.
    if (max(peaks[, "sample"]) > length(sampleGroups))
        stop("Sample indices in 'peaks' are larger than there are sample",
             " groups specified with 'sampleGroups'!")

    peaks <- cbind(peaks[, .reqCols, drop = FALSE],
                   index = index)

    ## Order peaks matrix by mz
    rownames(peaks) <- NULL
    peaks <- peaks[order(peaks[, "mz"]), , drop = FALSE]
    rtRange <- range(peaks[, "rt"])

    ## Define the mass slices and the index in the peaks matrix with an mz
    ## value >= mass[i]. If ppm != 0 the size of the individual bins will
    ## be dependend on the m/z value.
    mass <- breaks_ppm(peaks[1, "mz"], peaks[nrow(peaks), "mz"] + binSize,
                       by = binSize / 2, ppm = ppm / 2)
    masspos <- findEqualGreaterM(peaks[, "mz"], mass)

    densFrom <- rtRange[1] - 3 * bw
    densTo <- rtRange[2] + 3 * bw
    ## Increase the number of sampling points for the density distribution.
    densN <- max(512, 2 * 2^(ceiling(log2(diff(rtRange) / (bw / 2)))))
    endIdx <- 0
    niter <- length(mass) - 2L
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = 100,
                           clear = FALSE)
    pb_perc <- floor(seq(1, niter, length.out = 100))
    pb$tick(0)
    resL <- vector("list", niter)
    for (i in seq_len(niter)) {
        ## That's identifying overlapping mz slices.
        startIdx <- masspos[i]
        endIdx <- masspos[i + 2] - 1
        if (any(i == pb_perc))
            pb$tick()
        if (endIdx - startIdx < 0)
            next
        resL[[i]] <- .group_peaks_density(
            peaks[startIdx:endIdx, , drop = FALSE], bw = bw,
            densFrom = densFrom, densTo = densTo, densN = densN,
            sampleGroups = sampleGroups, sampleGroupTable = sampleGroupTable,
            minFraction = minFraction, minSamples = minSamples,
            maxFeatures = maxFeatures, sleep = sleep)
    }
    res <- do.call(rbind, resL)
    if (nrow(res)) {
        ## Remove groups that overlap with more "well-behaved" groups
        numsamp <- rowSums(
            as.matrix(res[, (match("npeaks", colnames(res)) +1):(ncol(res) -1),
                          drop = FALSE]))
        uorder <- order(-numsamp, res[, "npeaks"])

        uindex <- rectUnique(
            as.matrix(res[, c("mzmin", "mzmax", "rtmin", "rtmax"),
                          drop = FALSE]), uorder)
        res <- res[uindex, , drop = FALSE]
        rownames(res) <- NULL
    }
    res
}

#' @title Core API function for peak grouping using mzClust
#'
#' @description
#'
#' The `do_groupPeaks_mzClust` function performs high resolution
#' correspondence on single spectra samples.
#'
#' @inheritParams groupChromPeaks
#'
#' @inheritParams do_groupChromPeaks_density
#'
#' @inheritParams groupChromPeaks
#'
#' @return A `list` with elements `"featureDefinitions"` and
#' `"peakIndex"`. `"featureDefinitions"` is a `matrix`, each row
#' representing an (mz-rt) feature (i.e. peak group) with columns:
#' - `"mzmed"`: median of the peaks' apex mz values.
#' - `"mzmin"`: smallest mz value of all peaks' apex within the feature.
#' - `"mzmax"`: largest mz value of all peaks' apex within the feature.
#' - `"rtmed"`: always `-1`.
#' - `"rtmin"`: always `-1`.
#' - `"rtmax"`: always `-1`.
#' - `"npeaks"`: the total number of peaks assigned to the feature. Note that
#'   this number can be larger than the total number of samples, since
#'   multiple peaks from the same sample could be assigned to a group.
#'
#' `"peakIndex"` is a `list` with the indices of all peaks in a peak group in
#' the `peaks` input matrix.
#'
#' @md
#'
#' @family core peak grouping algorithms
#'
#' @references Saira A. Kazmi, Samiran Ghosh, Dong-Guk Shin, Dennis W. Hill
#' and David F. Grant\cr \emph{Alignment of high resolution mass spectra:
#' development of a heuristic approach for metabolomics}.\cr Metabolomics,
#' Vol. 2, No. 2, 75-83 (2006)
do_groupPeaks_mzClust <- function(peaks, sampleGroups, ppm = 20,
                                  absMz = 0, minFraction = 0.5,
                                  minSamples = 1) {
    if (missing(sampleGroups))
        stop("Parameter 'sampleGroups' is missing! This should be a vector of ",
             "length equal to the number of samples specifying the group ",
             "assignment of the samples.")
    if (missing(peaks))
        stop("Parameter 'peaks' is missing!")
    if (!(is.matrix(peaks) | is.data.frame(peaks)))
        stop("Peaks has to be a 'matrix' or a 'data.frame'!")
    ## Check that we've got all required columns
    .reqCols <- c("mz", "sample")
    if (!all(.reqCols %in% colnames(peaks)))
        stop("Required columns ",
             paste0("'", .reqCols[!.reqCols %in% colnames(peaks)],"'",
                    collapse = ", "), " not found in 'peaks' parameter")
    if (!is.factor(sampleGroups))
        sampleGroups <- factor(sampleGroups, levels = unique(sampleGroups))
    sampleGroupNames <- levels(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)
    ##sampleGroups <- as.numeric(sampleGroups)

    ## Check that sample groups matches with sample column.
    if (max(peaks[, "sample"]) > length(sampleGroups))
        stop("Sample indices in 'peaks' are larger than there are sample",
             " groups specified with 'sampleGroups'!")
    peaks <- .fix_mz_clust_peaks(peaks)
    grps <- mzClustGeneric(peaks[, .reqCols, drop = FALSE],
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
    return(list(featureDefinitions = grpmat, peakIndex = grps$idx))
}

#' @title Core API function for chromatic peak grouping using a nearest
#' neighbor approach
#'
#' @description
#'
#' The `do_groupChromPeaks_nearest` function groups peaks across samples by
#' creating a master peak list and assigning corresponding peaks from all
#' samples to each peak group (i.e. feature). The method is inspired by the
#' correspondence algorithm of mzMine (Katajamaa 2006).
#'
#' @inheritParams do_groupChromPeaks_density
#'
#' @inheritParams groupChromPeaks
#'
#' @return A `list` with elements `"featureDefinitions"` and
#' `"peakIndex"`. `"featureDefinitions"` is a `matrix`, each row
#' representing an (mz-rt) feature (i.e. peak group) with columns:
#'
#' - `"mzmed"`: median of the peaks' apex mz values.
#' - `"mzmin"`: smallest mz value of all peaks' apex within the feature.
#' - `"mzmax"`:largest mz value of all peaks' apex within the feature.
#' - `"rtmed"`: the median of the peaks' retention times.
#' - `"rtmin"`: the smallest retention time of the peaks in the feature.
#' - `"rtmax"`: the largest retention time of the peaks in the feature.
#' - `"npeaks"`: the total number of peaks assigned to the feature.
#'
#' `"peakIndex"` is a `list` with the indices of all peaks in a feature in the
#' `peaks` input matrix.
#'
#' @family core peak grouping algorithms
#'
#' @md
#'
#' @references Katajamaa M, Miettinen J, Oresic M: MZmine: Toolbox for
#' processing and visualization of mass spectrometry based molecular profile
#' data. Bioinformatics 2006, 22:634-636.
do_groupChromPeaks_nearest <- function(peaks, sampleGroups, mzVsRtBalance = 10,
                                       absMz = 0.2, absRt = 15, kNN = 10) {
    if (missing(sampleGroups))
        stop("Parameter 'sampleGroups' is missing! This should be a vector of ",
             "length equal to the number of samples specifying the group ",
             "assignment of the samples.")
    if (missing(peaks))
        stop("Parameter 'peaks' is missing!")
    if (!(is.matrix(peaks) | is.data.frame(peaks)))
        stop("Peaks has to be a 'matrix' or a 'data.frame'!")
    ## Check that we've got all required columns
    .reqCols <- c("mz", "rt", "sample")
    if (!all(.reqCols %in% colnames(peaks)))
        stop("Required columns ",
             paste0("'", .reqCols[!.reqCols %in% colnames(peaks)],"'",
                    collapse = ", "), " not found in 'peaks' parameter")
    if (!is.factor(sampleGroups))
        sampleGroups <- factor(sampleGroups, levels = unique(sampleGroups))
    sampleGroupNames <- levels(sampleGroups)
    sampleGroupTable <- table(sampleGroups)
    nSampleGroups <- length(sampleGroupTable)

    ## sampleGroups == classlabel
    ## nSampleGroups == gcount
    ## peaks == peakmat

    peaks <- peaks[, .reqCols, drop = FALSE]

    parameters <- list(mzVsRTBalance = mzVsRtBalance, mzcheck = absMz,
                       rtcheck = absRt, knn = kNN)

    ptable <- table(peaks[,"sample"])
    pord <- ptable[order(ptable, decreasing = TRUE)]
    sid <- as.numeric(names(pord))
    pn <- as.numeric(pord)

    ## environment - probably not a good idea for parallel processing - we
    ## would like to have data copying there (or better just provide the data
    ## chunk it needs to process).
    mplenv <- new.env(parent = .GlobalEnv)
    mplenv$mplist <- matrix(0, pn[1], length(sid))
    mplenv$mplist[, sid[1]] <- which(peaks[,"sample"] == sid[1])
    mplenv$mplistmean <- data.frame(peaks[which(peaks[,"sample"] == sid[1]),
                                            c("mz", "rt")])
    mplenv$peakmat <- peaks
    assign("peakmat", peaks, envir = mplenv)  ## double assignment?

    sapply(sid[2:length(sid)], function(sample, mplenv){
        message("Processing sample number ", sample, " ... ", appendLF = FALSE)
        ## require(parallel)
        ## cl <- makeCluster(getOption("cl.cores", nSlaves))
        ## clusterEvalQ(cl, library(RANN))
        ## parSapply(cl, 2:length(samples), function(sample,mplenv, object){
        ## might slightly improve on this for loop.
        ## Calculating for each row (peak) the mean mz or rt for peaks
        ## assigned yet to this peak group.
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
    ## Variable to count samples for a peak
    sampCounts <- rep_len(0, nSampleGroups)
    names(sampCounts) <- sampleGroupNames
    ## gcount <- integer(nSampleGroups)
    ## Can we vectorize that below somehow?
    for (i in 1:nrow(mplenv$mplist)) {
        groupmat[i, "mzmed"] <- median(peaks[mplenv$mplist[i, ], "mz"])
        groupmat[i, c("mzmin", "mzmax")] <- range(peaks[mplenv$mplist[i, ], "mz"])
        groupmat[i, "rtmed"] <- median(peaks[mplenv$mplist[i, ], "rt"])
        groupmat[i, c("rtmin", "rtmax")] <- range(peaks[mplenv$mplist[i, ], "rt"])

        groupmat[i, "npeaks"] <- length(which(peaks[mplenv$mplist[i, ]] > 0))

        ## Now summarizing the number of samples in which the peak was identified
        sampCounts[] <- 0
        tbl <- table(sampleGroups[peaks[mplenv$mplist[i, ], "sample"]])
        sampCounts[names(tbl)] <- as.numeric(tbl)
        groupmat[i, 7 + seq_len(nSampleGroups)] <- sampCounts
        ## gnum <- sampleGroups[unique(peaks[mplenv$mplist[i, ], "sample"])]
        ## for (j in seq(along = gcount))
        ##     gcount[j] <- sum(gnum == j)
        ## groupmat[i, 7 + seq(along = gcount)] <- gcount
        groupindex[[i]] <- mplenv$mplist[i, (which(mplenv$mplist[i,]>0))]
    }

    return(list(featureDefinitions = groupmat, peakIndex = groupindex))
}

#' Low level function to group chromatographic peaks within a m/z slice.
#'
#' @param x `matrix` such as the one returned by `chromPeaks,XCMSnExp`, just
#'     with the peaks within one m/z slice. Note that we require in addition
#'     a column `"index"` with the index of the peak within the full peak table.
#'
#' @param return `data.frame`
#'
#' @author Johannes Rainer
#'
#' @noRd
.group_peaks_density <- function(x, bw, densFrom, densTo, densN, sampleGroups,
                                 sampleGroupTable, minFraction,
                                 minSamples, maxFeatures, sleep = 0) {
    den <- density(x[, "rt"], bw = bw, from = densFrom, to = densTo,
                   n = densN)
    maxden <- max(den$y)
    deny <- den$y
    sampleGroupNames <- names(sampleGroupTable)
    nSampleGroups <- length(sampleGroupNames)
    col_nms <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                 "npeaks", sampleGroupNames)
    res_mat <- matrix(nrow = 0, ncol = length(col_nms),
                      dimnames = list(character(), col_nms))
    res_idx <- list()
    while (deny[maxy <- which.max(deny)] > maxden / 20 && nrow(res_mat) <
           maxFeatures) {
        grange <- descendMin(deny, maxy)
        deny[grange[1]:grange[2]] <- 0
        gidx <- which(x[,"rt"] >= den$x[grange[1]] &
                      x[,"rt"] <= den$x[grange[2]])
        ## Determine the sample group of the samples in which the peaks
        ## were detected and check if they correspond to the required limits.
        tt <- table(sampleGroups[unique(x[gidx, "sample"])])
        if (!any(tt / sampleGroupTable[names(tt)] >= minFraction &
                 tt >= minSamples))
            next
        gcount <- rep(0, length(sampleGroupNames))
        names(gcount) <- sampleGroupNames
        gcount[names(tt)] <- as.numeric(tt)
        res_mat <- rbind(res_mat,
                         c(median(x[gidx, "mz"]),
                           range(x[gidx, "mz"]),
                           median(x[gidx, "rt"]),
                           range(x[gidx, "rt"]),
                           length(gidx),
                           gcount)
                         )
        res_idx <- c(res_idx, list(unname(sort(x[gidx, "index"]))))
    }
    if (sleep > 0) {
        ## Plot the density
        plot(den, main = paste(round(min(x[,"mz"]), 2), "-",
                               round(max(x[,"mz"]), 2)))
        ## Highlight peaks per sample group.
        for (j in seq_len(nSampleGroups)) {
            ## Which peaks belong to this sample group.
            cur_group_samples <- which(sampleGroups == sampleGroupNames[j])
            idx <- x[, "sample"] %in% cur_group_samples
            points(x[idx, "rt"], x[idx, "into"] /
                                 max(x[, "into"]) * maxden,
                   col = j, pch=20)
        }
        for (j in seq_len(nrow(res_mat)))
            abline(v = res_mat[j, 5:6], lty = "dashed", col = j)
        Sys.sleep(sleep)
    }
    res <- as.data.frame(res_mat)
    res$peakidx <- res_idx
    res
}

#' @description
#'
#' Check the input peaks table eventually replacing `NA` values in column `"mz"`
#' with the mean of columns `"mzmin"` and `"mzmax"` (if present).
#' This fixes issue #416.
#'
#' @param x peaks `matrix`.
#'
#' @return peaks `matrix`
#'
#' @noRd
#'
#' @md
.fix_mz_clust_peaks <- function(x) {
    ## Issue #416: fix for peaks with an m/z of NA.
    nas <- is.na(x[, "mz"])
    if (any(nas)) {
        ## if we have mzmin and mzmax use mean of them.
        if (all(c("mzmin", "mzmax") %in% colnames(x)) &&
            !any(is.na(x[nas, c("mzmin", "mzmax")]))) {
            warning("Got ", sum(nas), " peaks with missing values in column ",
                    "'mz'. Replaced them with the mean of values in columns ",
                    "'mzmin' and 'mzmax' values.")
            x[nas, "mz"] <- rowMeans(x[nas, c("mzmin", "mzmax")])
        } else {
            stop("Got ", sum(nas), " peaks with missing values in column 'mz'.")
        }
    }
    x
}
