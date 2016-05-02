test.getEICxraw <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    step <- 0.1
    xraw <- xcmsRaw(file, profstep=step)
    e <- getEIC(xraw, rtrange=cbind(3000,3500), mzrange=cbind(200,201))
    ## calculate the EIC manually...
    mass <- seq(floor(xraw@mzrange[1]/step)*step,
                ceiling(xraw@mzrange[2]/step)*step, by = step)
    rtIdx <- which(xraw@scantime >= 3000 & xraw@scantime <= 3500)
    mzIdx <- which(mass >= 200 & mass <= 201)
    ecalc <- apply(xraw@env$profile[mzIdx, rtIdx, drop=FALSE], MARGIN=2, max)
    checkEqualsNumeric(ecalc, e@eic[[1]][[1]][, 2])
    ## extract the /full/ EIC, i.e. the base peak chromatogram (BPC).
    ## here we have to use the "new" method...
    BioC <- getOption("BioC")
    BioC$xcms$getEIC.method <- "getEICNew"
    options(BioC=BioC)
    e <- getEIC(xraw, mzrange=matrix(xraw@mzrange, nrow=1),
                rtrange=matrix(range(xraw@scantime), nrow=1))
    rtIdx <- which(xraw@scantime >= min(xraw@scantime) &
                       xraw@scantime <= max(xraw@scantime))
    mzIdx <- which(mass >= xraw@mzrange[1] & mass <= xraw@mzrange[2])
    ecalc <- apply(xraw@env$profile[mzIdx, rtIdx, drop=FALSE], MARGIN=2, max)
    checkEqualsNumeric(ecalc, e@eic[[1]][[1]][, 2])
    ## for two ranges...
    mzrange <- rbind(c(200, 201), c(300, 310), c(300, 402))
    rtrange <- rbind(c(3000, 3500), c(4000, 4300), c(2600, 3000))
    e <- getEIC(xraw, mzrange=mzrange, rtrange=rtrange)
    ## manually calculate...
    for(i in 1:nrow(mzrange)){
        rtIdx <- which(xraw@scantime >= rtrange[i, 1] & xraw@scantime <= rtrange[i, 2])
        mzIdx <- which(mass >= mzrange[i, 1] & mass <= mzrange[i, 2])
        ecalc <- apply(xraw@env$profile[mzIdx, rtIdx, drop=FALSE], MARGIN=2, max)
        checkEqualsNumeric(ecalc, e@eic[[1]][[i]][, 2])
    }
    ## restoring the setting...
    BioC <- getOption("BioC")
    BioC$xcms$getEIC.method <- "getEICOld"
    options(BioC=BioC)
}


## This test case checks for the problem reported by Alan Smith in issue #7
## on github.
## Apparently (and actually), the length of the retention time corrected and
## original retention times differ, e.g. when scanrange is used. Using the
## getXcmsRaw method fixes that, since within that method we ensure that all
## vectors, times etc are aligned (i.e. match in length and ordering).
test.issue7 <- function(){
    library(xcms)
    library(faahKO)

    cdfpath <- system.file("cdf", package = "faahKO")
    list.files(cdfpath, recursive = TRUE)
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)

    xset<-xcmsSet(cdffiles,profmethod = "bin", method="centWave",
                  ppm=30, peakwidth=c(5,50), snthresh=10, prefilter=c(7,1000),
                  integrate=1, mzdiff= -0.001, fitgauss=FALSE,
                  scanrange= c(0,500),mzCenterFun="wMean",noise=300,
                  sleep=0, verbose.columns=T); # parallel disabled: ,nSlaves=3);
    xset2<-retcor(xset, method = "obiwarp",distFunc="cor_opt",
                  profStep=1,plot="deviation",gapInit=0.8, gapExtend=3.3)
    gxset2<-group(xset2,method="density", bw=8, minfrac = .5,
                  minsamp = 2, mzwid = .01, max = 50, sleep = 0)
    xset3<-fillPeaks(gxset2)

    peaktab <-data.frame(ID=groupnames(xset3),peakTable(xset3))
    rownames(peaktab) <- peaktab$ID

    pkid <- c("M205T2790","M392T2927")
    a <- getEIC(xset3,groupidx=pkid)
    ## Get the retention times and eic for the first file.
    rts <- a@eic[[1]][[1]]
    ## rt has to be increasing!
    rtidx <- order(rts[, "rt"])
    checkEquals(rtidx, 1:length(rtidx))

    ## Reproduce what's happening within the getEIC for xcmsSet:
    ## That's how the xcmsRaw is loaded.
    braw <- xcmsRaw(cdffiles[1], profmethod=profMethod(xset3), profstep=0)
    ## Note that the vector lengths are different.
    checkTrue(length(braw@scantime) != length(xset3@rt$corrected[[1]]))
    ## That's the problem! That's how the corrected scantime is applied.
    braw@scantime <- xset3@rt$corrected[[1]]

    ## Now the lengths are different and we get the error.
    rtrange <- as.matrix(peaktab[pkid, c("rtmin", "rtmax")])
    rtrange[,1] <- rtrange[,1] - 200
    rtrange[,2] <- rtrange[,2] + 200
    beic <- getEIC(braw, mzrange=as.matrix(peaktab[pkid, c("mzmin", "mzmax")]),
                   rtrange=rtrange)
    plot(beic)
    rts <- beic@eic[[1]][[1]]
    ## rt has to be increasing!
    rtidx <- order(rts[, "rt"])
    checkTrue(all(rtidx != 1:length(rtidx)))

    ## Applying also the scanrange (as done in getXcmsRaw) fix it.
    craw <- getXcmsRaw(xset3)
    ceic <- getEIC(craw, mzrange=as.matrix(peaktab[pkid, c("mzmin", "mzmax")]),
                   rtrange=rtrange)
    plot(ceic)
    rts <- ceic@eic[[1]][[1]]
    ## rt has to be increasing!
    rtidx <- order(rts[, "rt"])
    checkEquals(rtidx, 1:length(rtidx))

    ## Interestingly, using getEICNew also helps:
    BioC <- getOption("BioC")
    BioCBefore <- getOption("BioC")
    BioC$xcms$getEIC.method <- "getEICNew"
    options(BioC=BioC)
    zic <- getEIC(xset3, groupidx=pkid[1])
    plot(zic)
    rts <- zic@eic[[1]][[1]]
    ## rt has to be increasing!
    rtidx <- order(rts[, "rt"])
    checkEquals(rtidx, 1:length(rtidx))

    ## Reset options.
    options(BioC=BioCBefore)
}

test.getEICxset <- function() {
    xset <- fillPeaks(group(faahko))
    e <- getEIC(xset, sampleidx=c(1,2), groupidx=c(1,2), rtrange=200)
    plot(e)
}

test.getEICretcor <- function() {
    xset <- fillPeaks(group(retcor(group(faahko))))
    opt.warn <- options("warn")$warn
    options("warn" = 2) ## turns warning into errors
    e <- getEIC(xset, sampleidx=c(1,2), groupidx=c(1,2),
                rt="corrected", rtrange=200)
    options("warn" = opt.warn)
    plot(e)
}

test.plotEIC <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)
    e <- getEIC(xraw, rtrange=cbind(3000,3500), mzrange=cbind(200,201))
    plot(e)
}
