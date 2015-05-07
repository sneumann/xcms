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
