library(MsExperiment)
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))
mse <- readMsExperiment(spectraFiles = fls, sampleData = df)
mse_ms2 <- readMsExperiment(msdata::proteomics(full.names=TRUE)[4])
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

test_that("plotChromPeaks,XcmsExperiment works", {
    expect_true(plotChromPeaks(xmse, 2))
})

test_that("plotChromPeakImage,XcmsExperiment works", {
    expect_true(plotChromPeakImage(xmse))
})

test_that("plot,XcmsExperiment and .xmse_plot_xic works", {
    tmp <- filterRt(mse, c(3000, 3100))
    plot(tmp, cex = 0.1)
    plot(tmp, msLevel = 2L)
    tmp <- filterRt(xmse, c(3000, 3100))
    .xmse_plot_xic(tmp, lwd = 3, col = NA, cex = 0.2)

    tmp <- filterMz(filterRt(xmse, rt = c(2550, 2800)), mz = c(342.5, 344.5))
    plot(tmp)
})

test_that("plotPrecursorIons works", {
    expect_error(plotPrecursorIons(3), "MsExperiment")
    fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
                      package = "msdata")
    a <- readMsExperiment(fl)
    plotPrecursorIons(a, main = "SWATH")

    fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
                      package = "msdata")
    a <- readMsExperiment(fl)
    plotPrecursorIons(a)
})

test_that(".xmse_plot_xic works with ms2 data", {
  tmp <-  filterMz(filterRt(mse_ms2, rt= c(2160, 2190)), mz = c(990,1000))
  plot(tmp)
})
