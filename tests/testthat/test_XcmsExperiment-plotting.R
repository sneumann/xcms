library(MsExperiment)
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))
mse <- readMsExperiment(spectraFiles = fls, sampleData = df)
p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
xmse <- findChromPeaks(mse, param = p)
pdp <- PeakDensityParam(sampleGroups = rep(1, 3))
xmseg <- groupChromPeaks(xmse, param = pdp, add = FALSE)
xmsegr <- adjustRtime(xmseg, param = PeakGroupsParam(span = 0.4))

test_that(".plot_adjusted_rtime works, .plot_peak_groups works", {
    .plot_adjusted_rtime(rtime(xmsegr, adjusted = FALSE), rtime(xmsegr),
                         from_file = fromFile(xmsegr))
    ph <- processHistory(xmsegr, type = xcms:::.PROCSTEP.RTIME.CORRECTION)[[1L]]

    rt <- split(rtime(xmsegr, adjusted = FALSE), fromFile(xmsegr))
    rtadj <- split(rtime(xmsegr), fromFile(xmsegr))
    .plot_peak_groups(rt, rtadj, peakGroupsMatrix(ph@param))

    ## simulating a subset:
    idx <- c(1, 3)
    .plot_peak_groups(rt[idx], rtadj[idx],
                      peakGroupsMatrix(ph@param)[, idx, drop = FALSE],
                      col = "red")
})

test_that("plotAdjustedRtime,XcmsExperiment works", {
    expect_warning(plotAdjustedRtime(xmseg), "results present")
    expect_error(plotAdjustedRtime(mse), "XcmsExperiment")
    plotAdjustedRtime(xmsegr, col = "red")
    plotAdjustedRtime(xmsegr, adjustedRtime = FALSE)
})
