## Functions for xcmsSet objects.
#' @include DataClasses.R do_findChromPeaks-functions.R

## The "constructor"
## The "new" xcmsSet method using BiocParallel.
xcmsSet <- function(files = NULL, snames = NULL, sclass = NULL,
                    phenoData = NULL, profmethod = "bin",
                    profparam = list(), polarity = NULL,
                    lockMassFreq=FALSE, mslevel=NULL, nSlaves=0,
                    progressCallback=NULL, scanrange=NULL,
                    BPPARAM=bpparam(), stopOnError = TRUE, ...) {

    if (nSlaves != 0) {
        message("Use of argument 'nSlaves' is deprecated,",
                " please use 'BPPARAM' instead.")
        options(mc.cores = nSlaves)
    }
    if (!is.logical(stopOnError))
        stop("'stopOnError' has to be a logical.")
    ## Overwriting the stop.on.error in BPPARAM:
    ## bpstopOnError(BPPARAM) <- stopOnError
    orig <- bpstopOnError(BPPARAM)
    on.exit(BPPARAM@.xData$stop.on.error <- orig)
    BPPARAM@.xData$stop.on.error <- stopOnError

    object <- new("xcmsSet")

    ## initialise progress information
    xcms.options <- getOption("BioC")$xcms
    xcms.methods <- c(paste("group", xcms.options$group.methods,sep="."),
                      paste("findPeaks", xcms.options$findPeaks.methods,sep="."),
                      paste("retcor", xcms.options$retcor.methods,sep="."),
                      paste("fillPeaks", xcms.options$fillPeaks.methods,sep="."))
    eval(parse(text=paste("object@progressInfo <- list(",paste(xcms.methods,"=0",
                                                               sep="",collapse=","),")") ))

    if (is.function(progressCallback))
        object@progressCallback <- progressCallback

    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
        files <- getwd()
    info <- file.info(files)
    if (any(info$isdir)) {
        message("Scanning files in directory ", files[info$isdir], " ... ",
                appendLF = FALSE)
        listed <- list.files(files[info$isdir], pattern = filepattern,
                             recursive = TRUE, full.names = TRUE)
        message("found ", length(listed), " files")
        files <- c(files[!info$isdir], listed)
    }
    ## try making paths absolute
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]
    if (length(files) == 0 | all(is.na(files)))
        stop("No NetCDF/mzXML/mzML files were found.\n")

    if(lockMassFreq==TRUE){
        ## remove the 02 files if there here
        lockMass.files<-grep("02.CDF", files)
        if(length(lockMass.files) > 0){
            files<-files[-lockMass.files]
        }
    }
    filepaths(object) <- files

    ## determine experimental design
    fromPaths <- phenoDataFromPaths(files)
    if (is.null(snames)) {
        snames <- rownames(fromPaths)
    } else {
        rownames(fromPaths) <- snames
    }
    pdata <- phenoData
    if (is.null(pdata)) {
        pdata <- sclass
        if (is.null(pdata))
            pdata <- fromPaths
    }else{
        if(class(pdata)=="AnnotatedDataFrame")
            pdata <- as(pdata, "data.frame")
        if(class(pdata)!="data.frame")
            stop("phenoData has to be a data.frame or AnnotatedDataFrame!")
    }
    phenoData(object) <- pdata
    if (is.null(phenoData))
        rownames(phenoData(object)) <- snames
    rtlist <- list(raw = vector("list", length(snames)),
                   corrected = vector("list", length(snames)))

    if ("step" %in% names(profparam)) {
        if ("step" %in% names(list(...)) && profparam$step != list(...)$step) {
            stop("different step values defined in profparam and step arguments")
        }
        profstep <- profparam$step
        profparam <- profparam[names(profparam) != "step"]
    } else if ("step" %in% names(list(...))) {
        profstep <- list(...)$step
    } else {
        profstep <- 0.1
    }

    if ("method" %in% names(profparam)) {
        if (profparam$method != profmethod) {
            stop("different method values defined in profparam and profmethod arguments")
        }
        profmethod <- profparam$method
        profparam <- profparam[names(profparam) != "method"]
    }

    profinfo(object) <- c(list(method = profmethod, step = profstep), profparam)

    object@polarity <- as.character(polarity)
    includeMSn=FALSE

    ## implicitely TRUE if selecting MSn
    includeMSn <- !is.null(mslevel) &&  mslevel>1

    ## implicitely TRUE if MS1 parent peak picking
    xcmsSetArgs <- as.list(match.call())
    if (!is.null(xcmsSetArgs$method)) {
        if (xcmsSetArgs$method=="MS1") {
            includeMSn=TRUE
        }
    }

    params <- list(...)
    params$profmethod <- profmethod
    params$profparam <- profparam
    params$includeMSn <- includeMSn
    params$scanrange <- scanrange

    params$mslevel <- mslevel ## Actually, this is
    params$lockMassFreq <- lockMassFreq

    ft <- cbind(file = files,id = 1:length(files))
    argList <- apply(ft ,1 ,function(x) list(file = x["file"],
                                           id = as.numeric(x["id"]),
                                           params = params))
    ## Use BiocParallel: bplapply along with bptry
    res <- bptry(bplapply(argList, findPeaksPar, BPPARAM = BPPARAM))
    ## Catch the error.
    if (stopOnError) {
        if (inherits(res, "bperror"))
            stop(res)
    }

    ## Processing the results
    ## Feature detection in <file>: xy peaks. -> info
    ## Error identifying features in <>: Error:... -> info + error
    isOK <- bpok(res)
    if (all(!isOK))
        stop("Chromatographic peak detection failed for all files!",
             " The first error was: ", res[[1]])
    if (any(!isOK)) {
        ## Use scantime from a working file one of the failing ones.
        scnt <- res[[which(isOK)[1]]]$scantime
    }
    nres <- length(res)
    peaklist <- vector("list", length = nres)
    proclist <- vector("list", length = nres)
    scntlist <- vector("list", length = nres)
    ## Loop trough the results and gather information.
    for (i in 1:nres) {
        if (isOK[i]) {
            scntlist[[i]] <- res[[i]]$scantime
            pks <- res[[i]]$peaks
            peaklist[[i]] <- pks
            if (is.null(pks))
                warning("No peaks found in sample ", snames[i], ".")
            else if (nrow(pks) == 0)
                warning("No peaks found in sample ", snames[i], ".")
            else if (nrow(pks) == 1)
                warning("Only 1 peak found in sample ", snames[i], ".")
            else if (nrow(pks) < 5)
                warning("Only ", nrow(pks), " found in sample ", snames[i], ".")
            proclist[[i]] <- ProcessHistory(info. = paste0("Peak detection in '",
                                                          basename(files[i]),
                                                          "': ", nrow(pks),
                                                          " peaks identified."),
                                            date. = res[[i]]$date,
                                            type. = .PROCSTEP.PEAK.DETECTION,
                                            fileIndex. = i)
        } else {
            scntlist[[i]] <- scnt
            peaklist[[i]] <- NULL
            proclist[[i]] <- ProcessHistory(info. = paste0("Error identifying",
                                                          " peaks in '",
                                                          basename(files[i]),
                                                          "': ", res[[i]]),
                                            error. = res[[i]],
                                            type. = .PROCSTEP.PEAK.DETECTION,
                                            fileIndex. = i)
            warning("Peak detection failed in '", files[i], "':", res[[i]])
        }
    }
    ## peaklist <- lapply(res, function(x) x$peaks)
    ## rtlist$raw <-  rtlist$corrected <- lapply(res, function(x) x$scantime)
    if(lockMassFreq){
        object@dataCorrection[1:length(files)] <- 1
    }

    rtlist$raw <- rtlist$corrected <- scntlist
    peaks(object) <- do.call(rbind, peaklist)
    object@rt <- rtlist
    object@.processHistory <- proclist

    mslevel(object) <- as.numeric(mslevel)
    scanrange(object) <- as.numeric(scanrange)
    OK <- .validProcessHistory(object)
    if (!is.logical(OK))
       stop(OK)
    object
}

############################################################
## c
c.xcmsSet <- function(...) {
    lcsets <- list(...)
    object <- new("xcmsSet")

    peaklist <- vector("list", length(lcsets))
    namelist <- vector("list", length(lcsets))
    if (any(duplicated(unlist(namelist)))) {
        stop("Duplicated sample names\n")
    }

    classlist <- vector("list", length(lcsets))
    cdflist <- vector("list", length(lcsets))
    rtraw <- vector("list", 0)
    rtcor <- vector("list", 0)
    procHist <- vector("list", 0)
    startIdx <- 1
    nsamp <- 0
    for (i in seq(along = lcsets)) {
        peaklist[[i]] <- peaks(lcsets[[i]])
        namelist[[i]] <- sampnames(lcsets[[i]])
        classlist[[i]] <- sampclass(lcsets[[i]])
        classlist[[i]] <- levels(classlist[[i]])[classlist[[i]]]
        cdflist[[i]] <- filepaths(lcsets[[i]])
        rtraw <- c(rtraw, lcsets[[i]]@rt$raw)
        rtcor <- c(rtcor, lcsets[[i]]@rt$corrected)

        ## Update samples only if we've got any peaks. Issue #133
        if (nrow(peaks(lcsets[[i]]))) {
            sampidx <- seq(along = namelist[[i]]) + nsamp
            peaklist[[i]][,"sample"] <- sampidx[peaklist[[i]][,"sample"]]
            ## Don't increment if we don't have any peaks
            nsamp <- nsamp + length(namelist[[i]])
        }
        if (.hasSlot(lcsets[[i]], ".processHistory")) {
            ph <- .getProcessHistory(lcsets[[i]])
            if (length(ph) > 0) {
                num_files <- length(namelist[[i]])
                ph <- lapply(ph, updateFileIndex, old = 1:num_files,
                             new = startIdx:(startIdx+num_files-1))
            } else {
                ph <- list()
            }
            procHist <- c(procHist, ph)
        } else {
            procHist <- c(procHist, list())
        }
        startIdx <- startIdx + length(namelist[[i]])
    }

    peaks(object) <- do.call(rbind, peaklist)
    sampnames(object) <- unlist(namelist)
    classlist <- unlist(classlist)
    sampclass(object) <- factor(classlist, unique(classlist))
    filepaths(object) <- unlist(cdflist)
    profinfo(object) <- profinfo(lcsets[[1]])
    object@rt <- list(raw = rtraw, corrected = rtcor)

    object@.processHistory <- procHist
    OK <- .validProcessHistory(object)
    if (!is.logical(OK))
        stop(OK)
    invisible(object)
}

############################################################
## split
split.xcmsSet <- function(x, f, drop = TRUE, ...) {
    ## Update the object
    x <- updateObject(x)
    if (!is.factor(f))
        f <- factor(f)
    sampidx <- unclass(f)
    peakmat <- peaks(x)
    samples <- sampnames(x)
    classlabel <- sampclass(x)
    cdffiles <- filepaths(x)
    prof <- profinfo(x)
    rtraw <- x@rt$raw
    rtcor <- x@rt$corrected

    lcsets <- vector("list", length(levels(f)))
    names(lcsets) <- levels(f)

    ## get the phenoData and all other parameters we want to pass
    ## down to the splitted objects
    pd <- phenoData(x)
    dataCor <- x@dataCorrection
    suppressWarnings(
        msL <- mslevel(x)
        )
    suppressWarnings(
        scanR <- scanrange(x)
        )
    pol <- x@polarity
    procHist <- x@.processHistory
    for (i in unique(sampidx)) {
        samptrans = which(sampidx == i)
        samptrans = samptrans[samptrans <= nrow(x@phenoData)]

        if (length(samptrans) < 1) next

        lcsets[[i]] <- new("xcmsSet")

        cpeaks = peakmat[peakmat[,"sample"] %in% samptrans, ,drop=F]
        cpeaks[,"sample"] <- as.numeric(factor(cpeaks[,"sample"]))
        peaks(lcsets[[i]]) <- cpeaks

        ## don't need these, since I'm going to use the phenoData instead.
        ## sampnames(lcsets[[i]]) <- samples[samptrans]
        ## sampclass(lcsets[[i]]) <- classlabel[samptrans, drop = TRUE]
        phenoData(lcsets[[i]]) <- droplevels(pd[samptrans, , drop=FALSE])
        ## set also all other settings.
        if(length(dataCor) > 1){
            lcsets[[i]]@dataCorrection <- dataCor[samptrans]
        }
        ## suppressWarnings, as "old" xcmsSet objects don't have these
        ## slots.
        suppressWarnings(
            mslevel(lcsets[[i]]) <- as.numeric(msL)
            )
        suppressWarnings(
            scanrange(lcsets[[i]]) <- as.numeric(scanR)
            )
        lcsets[[i]]@polarity <- pol
        filepaths(lcsets[[i]]) <- cdffiles[samptrans]
        profinfo(lcsets[[i]]) <- prof
        lcsets[[i]]@rt$raw <- rtraw[samptrans]
        lcsets[[i]]@rt$corrected <- rtcor[samptrans]
        ## .processHistory
        procHist <- .getProcessHistory(x, fileIndex = samptrans)
        procHist <- lapply(procHist, updateFileIndex, old = samptrans,
                           new = 1:length(samptrans))
        lcsets[[i]]@.processHistory <- procHist
        OK <- .validProcessHistory(lcsets[[i]])
        if (!is.logical(OK))
            stop(OK)
    }

    if (drop)
        lcsets <- lcsets[!sapply(lcsets, is.null)]

    lcsets
}

############################################################
## phenoDataFromPaths
## derive experimental design from set of file paths
#' @title Derive experimental design from file paths
#'
#' @description The `phenoDataFromPaths` function builds a `data.frame`
#'     representing the experimental design from the folder structure in which
#'     the files of the experiment are located.
#'
#' @note This function is used by the *old* `xcmsSet` function to guess
#'     the experimental design (i.e. group assignment of the files) from the
#'     folders in which the files of the experiment can be found.
#'
#' @param paths `character` representing the file names (including the full
#'     path) of the experiment's files.
#'
#' @md
#'
#' @examples
#' ## List the files available in the faahKO package
#' base_dir <- system.file("cdf", package = "faahKO")
#' cdf_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)
phenoDataFromPaths <- function(paths) {
    ## create factors from filesystem hierarchy
    sclass <- gsub("^\\.$", "sample", dirname(paths))
    lev <- strsplit(sclass, "/")
    levlen <- sapply(lev, length)
    if(length(lev) > 1 && !all(levlen[1] == levlen))
        stop("Directory tree must be level")
    pdata <- as.data.frame(matrix(unlist(lev), nrow=length(lev), byrow=TRUE))
    redundant <- apply(pdata, 2, function(col) length(unique(col)) == 1)
    if (!any(!redundant)) {
        redundant[length(redundant)] <- FALSE
    }
    pdata <- pdata[,!redundant,drop=FALSE]
    if (ncol(pdata) == 1) { ## if not multiple factors, behave as before
        ## Make the default group names less redundant
        scomp <- strsplit(substr(sclass, 1, min(nchar(sclass))), "")
        scomp <- matrix(c(scomp, recursive = TRUE), ncol = length(scomp))
        i <- 1
        while(all(scomp[i,1] == scomp[i,-1]) && i < nrow(scomp))
            i <- i + 1
        i <- min(i, tail(c(0, which(scomp[1:i,1] == .Platform$file.sep)), n = 1) + 1)
        if (i > 1 && i <= nrow(scomp))
            sclass <- substr(sclass, i, max(nchar(sclass)))
        pdata <- data.frame(factor(sclass))
        colnames(pdata) <- "class"
    }
    rownames(pdata) <- gsub("\\.[^.]*$", "", basename(paths))
    pdata
}

############################################################
## patternVsRowScore
patternVsRowScore <- function(currPeak, parameters, mplenv)
{
    mplistmeanCurr <- mplenv$mplistmean[, c("mz", "rt")]
    mplistmeanCurr[, "mz"] <- mplistmeanCurr[, "mz"] * parameters$mzVsRTBalance
    peakmatCurr <- mplenv$peakmat[currPeak, c("mz", "rt"), drop = FALSE]
    peakmatCurr[, "mz"] <- peakmatCurr[, "mz"] * parameters$mzVsRTBalance

    nnDist <- nn2(mplistmeanCurr, peakmatCurr[, c("mz", "rt"), drop = FALSE],
                  k = min(length(mplistmeanCurr[, 1]), parameters$knn))

    scoreListcurr <- data.frame(score = numeric(0),
                                peak = integer(0),
                                mpListRow = integer(0),
                                isJoinedPeak = logical(0),
                                isJoinedRow = logical(0))

    rtTolerance = parameters$rtcheck

    for (mplRow in 1:length(nnDist$nn.idx)) {
        mplistMZ <- mplenv$mplistmean[nnDist$nn.idx[mplRow], "mz"]
        mplistRT <- mplenv$mplistmean[nnDist$nn.idx[mplRow], "rt"]

        ## Calculate differences between M/Z and RT values of current peak and
        ## median of the row
        diffMZ = abs(mplistMZ - mplenv$peakmat[[currPeak, "mz"]])
        diffRT = abs(mplistRT - mplenv$peakmat[[currPeak, "rt"]])

        ## Calculate if differences within tolerancdiffRT < rtTolerance)es
        if ( (diffMZ < parameters$mzcheck) & (diffRT < rtTolerance) ) {
            scoreListcurr <- rbind(scoreListcurr,
                                   data.frame(score = nnDist$nn.dists[mplRow],
                                              peak = currPeak,
                                              mpListRow = nnDist$nn.idx[mplRow],
                                              isJoinedPeak = FALSE,
                                              isJoinedRow = FALSE))
            ## goodEnough = true
            return(scoreListcurr)
        }
    }

    return(scoreListcurr) ## empty
}

############################################################
## getSpecWindow
getSpecWindow <- function(xs, gidxs, borderwidth=1){
    groupidx <- groupidx(xs)
    if (length(groupidx)==0) {
        stop("no groups found in xcmsSet object.")
    }
    minmaxs <- matrix(ncol=2, nrow=length(gidxs))
    cat("Processing data from sample: ")
    for (a in 1:length(gidxs)){
        ## 1st step: getting boundaries
        mzmin <- min(peaks(xs)[groupidx[[gidxs[a]]],"mzmin"]) ## lower bound
        mzmax <- max(peaks(xs)[groupidx[[gidxs[a]]],"mzmax"]) ## upper bound
        mzw <- mzmax-mzmin
        minmaxs[a,1] <- mzmin - borderwidth*mzw
        minmaxs[a,2] <- mzmax + borderwidth*mzw ## a peakwidth left and right
    }
    mzlistlist=list()
    for (s in 1:length(gidxs)) mzlistlist[[s]] <- list()
    mzlistlist$minmax <- minmaxs
    for (s in 1:length(sampnames(xs))){
        ## 2nd Step: getting each sample
        cat(" ",s)
        xr <- xcmsRaw(xs@filepaths[s])
        for (a in 1:length(gidxs)){
            pmin <- min(which(abs(xr@env$mz - minmaxs[a,1]) == min(abs(xr@env$mz - minmaxs[a,1]))))
            pmax <- min(which(abs(xr@env$mz - minmaxs[a,2]) == min(abs(xr@env$mz - minmaxs[a,2]))))
            mzlistlist[[a]][[s]] <- matrix(ncol=2, nrow=(pmax-pmin+1),
                                           data=c(xr@env$mz[pmin:pmax],xr@env$intensity[pmin:pmax]))
        }
    }
    cat("\n")
    invisible(mzlistlist)
}

############################################################
## plotSpecWindow
plotSpecWindow <- function(xs, gidxs, borderwidth=1){
    if (length(groupidx(xs))==0) {
        stop("no groups found in xcmsSet object.")
    }

    mzll <- getSpecWindow(xs, gidxs, borderwidth)
    minmax <- mzll$minmax
    groupidx <- groupidx(xs)
    pcolors <- c("black","darkred", "green","blue")
    ecolors <- c("grey","red","lightgreen","lightblue")

    cat("\ngroup: ")
    for (a in 1:length(gidxs)){
        cat(" ",gidxs[a])
        nsa <- length(levels(sampclass(xs)))
        pcol <- pcolors[which(levels(sampclass(xs)) == sampclass(xs)[1])]
        ecol <- ecolors[which(levels(sampclass(xs)) == sampclass(xs)[1])]
        mzmin <- minmax[a,1]
        mzmax <- minmax[a,2]
        maxints <- NA

        for (s in 1:length(mzll[[a]])) {
            maxints[s] <- max(mzll[[a]][[s]][,2])
        }

        maxo <- max(maxints)
        for (n in 1:length(mzll[[a]])){
            par(lty=1)
            pcol <- pcolors[which(levels(sampclass(xs)) == sampclass(xs)[n])]
            ecol <- ecolors[which(levels(sampclass(xs)) == sampclass(xs)[n])]
            if(n==1) {
                plot(mzll[[a]][[n]][,1],
                     mzll[[a]][[n]][,2],
                     xlim=c(mzmin,mzmax), ylim=c(0,(maxo+10000)),
                     type='l', col=ecol, xlab="", ylab="") ## complete raw-range
            } else {
                lines(mzll[[a]][[n]][,1],
                      mzll[[a]][[n]][,2],
                      xlim=c(mzmin,mzmax), ylim=c(0,(maxo+10000)),
                      type='l', col=ecol, xlab="", ylab="") ## complete raw-range}
            }
            if (length(which(peaks(xs)[groupidx[[gidxs[a]]],"sample"] == n))>0){
                ## peak entry of the first sample
                apeak <- groupidx[[gidxs[a]]][min(which(peaks(xs)[groupidx[[gidxs[a]]],"sample"] == n))]
                if (apeak %in% xs@filled) {
                    par(lty=2)
                }else{
                    par(lty=1)
                }
                ppmin <- min(which(abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmin"]) ==
                                   min(      abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmin"]))))
                ppmax <- min(which(abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmax"]) ==
                                   min(      abs(mzll[[a]][[n]][,1] - peaks(xs)[apeak,"mzmax"]))))
                lines(mzll[[a]][[n]][ppmin:ppmax,1],mzll[[a]][[n]][ppmin:ppmax,2],
                      xlim=c(mzmin,mzmax), col=pcol, type='l')
            }
        }

        title(main=paste("m/z vs. intensity for group",gidxs[a]), xlab="m/z", ylab="intensity")
        legend("topright", as.vector(levels(sampclass(xs))),col=pcolors[1:nsa],lty=rep(1,nsa))
    }
    cat("\n")
}

##
## before starting conversion to metaboAnalyst format:
## grouping, opt. retcoring (+grouping) and peak filling the xcmsObject
##
.write.metaboanalyst <- function(object, filename, phenoDataColumn=NULL, value="into", ...) {

    if (! "value" %in% names(list(...))) {
        p <-groupval(object, value="into", ... )
    } else {
        p <-groupval(object, value=value, ... )
    }

    if (missing(phenoDataColumn)) {
        labels <- as.character(sampclass(object))
    } else {
        labels <- as.character(phenoData(object)[, phenoDataColumn])
    }

    if(any(table(labels)<3))
        stop(paste("The classes", paste(names(which(table(labels)<3)), collapse=", "), "have less than 3 samples"))

    p <- rbind(Sample=sampnames(object),
               Label=labels,
               p)

    write.table(p, file = filename,
                dec=".", sep=",", qmethod="double",
                col.names=F, row.names = T)

}

############################################################
## xcmsBoxPlot
xcmsBoxPlot<-function(values, className, dirpath, pic, width=640, height=480){

    if (pic == "png"){
        png(filename = file.path(dirpath, "%003d.png"), width = width,
            height = height, units = "px")
    } else{
        pdf(file.path(dirpath, "%003d.pdf"), width = width/72, height = height/72, onefile = FALSE)
    }

    ind<-which(colnames(values) != "metlin")
    for (i in 1:nrow(values)){
        boxplot(as.numeric(values[i,ind]) ~ className, col="blue",
                outline=FALSE, main=paste("Feature ", row.names(values)[i] ))
    }
    if (length(values) > 0) {
        dev.off()
    }
}

############################################################
## retexp
retexp <- function(peakrange, width = 200) {

    retmean <- rowMeans(peakrange[,c("rtmin", "rtmax"),drop=FALSE])
    peakrange[,"rtmin"] <- retmean-width/2
    peakrange[,"rtmax"] <- retmean+width/2

    peakrange
}

############################################################
## na.flatfill
## Fill in NA values with the minimum and maximum value
na.flatfill <- function(x) {

    realloc <- which(!is.na(x))
    if (realloc[1] > 1)
        x[1:(realloc[1]-1)] <- x[realloc[1]]
    if (realloc[length(realloc)] < length(x))
        x[(realloc[length(realloc)]+1):length(x)] <- x[realloc[length(realloc)]]
    x
}

############################################################
## pval
pval <- function(X, classlabel, teststat) {

    n1 <- rowSums(!is.na(X[,classlabel == 0]))
    n2 <- rowSums(!is.na(X[,classlabel == 1]))
    A <- apply(X[,classlabel == 0], 1, sd, na.rm=TRUE)^2/n1 ## sd(t(X[,classlabel == 0]), na.rm = TRUE)^2/n1
    B <- apply(X[,classlabel == 1], 1, sd, na.rm=TRUE)^2/n2 ## sd(t(X[,classlabel == 1]), na.rm = TRUE)^2/n2
    df <- (A+B)^2/(A^2/(n1-1)+B^2/(n2-1))

    pvalue <- 2 * (1 - pt(abs(teststat), df))
    invisible(pvalue)
}

############################################################
## panel.cor
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex)
}

############################################################
## defineColAndTy
## defines the colors and ltys to be used in the plot functions...
## the code combines original code taken from the plot functions and
## adds the possibility to sumbmit colors.
## col: colors to be used for the samples
## ty: lty to be used for the samples
## classlabels: factor with class labels (group definitions), length
##              being equal to the number of samples.
defineColAndTy <- function(col=NULL, ty=NULL, classlabel){
    if(missing(classlabel))
        stop("classlabel is required")
    n <- length(classlabel)
    if(!is.factor(classlabel))
        classlabel <- factor(classlabel)
    ## want to transform the class labels to integer values, i.e. skip the levels
    classlabel <- as.numeric(classlabel)
    if(is.null(col)) {
        col <- integer(n)
        for (i in 1:max(classlabel))
            col[classlabel == i] <- 1:sum(classlabel == i)
    }else{
        ## check if we do have the same number of colors than samples
        if(length(col) != n){
            warning("Less colors than samples! Using the first color for all samples.")
            col <- rep(col[1], n)
        }
    }
    if(is.null(ty)) {
        ## allow col being not just integers...
        col.int <- as.numeric(factor(col))
        ty <- integer(n)
        for (i in 1:max(col.int))
            ty[col.int == i] <- 1:sum(col.int == i)
    }else{
        if(length(ty) != n){
            warning("Less line types than samples! Using the first type for all samples.")
            ty <- rep(ty[1], n)
        }
    }
    ## if col is a character vector (e.g. colors defined by RColorBrewer)
    if(!is.numeric(col)){
        ## define the mypal... that's the color vector as used below.
        mypal <- col
        ## col is now a numeric vector
        col <- 1:length(mypal)
    }else{
        if (length(palette()) < max(col))
            mypal <- rainbow(max(col), end = 0.85)
        else
            mypal <- palette()[1:max(col)]
    }
    return(list(col=col, ty=ty, mypal=mypal ))
}

############################################################
## filtfft
filtfft <- function(y, filt) {

    yfilt <- numeric(length(filt))
    yfilt[1:length(y)] <- y
    yfilt <- fft(fft(yfilt, inverse = TRUE) * filt)

    Re(yfilt[1:length(y)])
}

############################################################
## getProcessErrors
## get ProcessHistory objects with an error.
.getProcessErrors <- function(x, PROCSTEP = .PROCSTEPS) {
    if (!missing(PROCSTEP)) {
        PROCSTEP <- PROCSTEP %in% .PROCSTEPS
        if (length(PROCSTEP) == 0)
            stop("PROCSTEP not OK.")
    }
    if (.hasSlot(x, ".processHistory")) {
        errs <- lapply(x@.processHistory, function(z) {
            if (!is.null(z@error)) {
                return(z)
            }
        })
        errs <- errs[lengths(errs) > 0]
        return(errs)
    }
    return(list())
}

############################################################
## .getProcessHistory
## Get ProcessHistory objects allowing to retrieve either ProcessHistory
## objects for selected fileIndex, or for a specific PROCSTEP
.getProcessHistory <- function(x, fileIndex, PROCSTEP = .PROCSTEPS) {
    if (!missing(PROCSTEP)) {
        PROCSTEP <- PROCSTEP %in% .PROCSTEPS
        if (length(PROCSTEP) == 0)
            stop("PROCSTEP not OK.")
    }
    if (missing(fileIndex)) {
        fileIndex <- 1:length(filepaths(x))
    } else {
        if (!all(fileIndex %in% 1:length(filepaths(x))))
            stop("'fileIndex' does not match the number of files.")
    }
    if (.hasSlot(x, ".processHistory")) {
        phs <- lapply(x@.processHistory, function(z) {
            if (any(z@fileIndex %in% fileIndex) & z@type %in% PROCSTEP) {
                return(z)
            }
        })
        phs <- phs[lengths(phs) > 0]
        return(phs)
    }
    return(list())
}

############################################################
## .validProcessHistory
## Check the validity of the .processHistory slot.
.validProcessHistory <- function(x) {
    msg <- character()
    if (.hasSlot(x, ".processHistory")) {
        if (length(x@.processHistory) > 0) {
            ## All elements have to inherit from ProcessHistory
            if (!all(unlist(lapply(x@.processHistory, function(z) {
                return(inherits(z, "ProcessHistory"))
            }))))
                msg <- c(msg, paste0("All objects in slot .processHistory",
                                     " have to be 'ProcessHistory' objects!"))
            ## Each element has to be valid
            vals <- lapply(x@.processHistory, validObject)
            for (i in seq_along(vals)) {
                if (!is.logical(vals[[i]]))
                    msg <- c(msg, vals[[i]])
            }
            ## The fileIndex has to be within 1:length(filepaths(x))
            fidx <- 1:length(filepaths(x))
            for (z in x@.processHistory) {
                if (length(z@fileIndex) == 0 |
                    !(all(z@fileIndex %in% fidx)))
                    msg <- c(msg, paste0("Value of 'fileIndex' slot of some",
                                         " ProcessHistory objects does not",
                                         " match the number of available",
                                         " files!"))
            }
        }
    }
    if (length(msg)) msg
    else TRUE
}
