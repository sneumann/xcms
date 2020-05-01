test_that("groupval returns expected output on object converted from XCMSnExp", {

    cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
        recursive = TRUE)[1:3]
    xset <- xcmsSet(cdfs, method = 'centWave', ppm = 25, peakwidth = c(20, 80),
        snthresh = 10, prefilter = c(3,100), integrate = 1, mzdiff = -0.001,
        verbose.columns = FALSE, fitgauss = FALSE, noise = 5000)
    xset <- xcms::group(xset, method = "density", bw = 30, minfrac = 0.5,
        minsamp = 1)
    
    raw_data <- readMSData(files = cdfs,
        pdata = new("NAnnotatedDataFrame"), mode = "onDisk")
    cwp <- xcms::CentWaveParam(peakwidth = c(20, 80), noise = 5000)
    xdata <- xcms::findChromPeaks(raw_data, param = cwp)
    xdata <- xcms::groupChromPeaks(xdata,
        xcms::PeakDensityParam(sampleGroups = rep("KO", 3)))
    
    expect_equivalent(xcms::peaks(xset), xcms::chromPeaks(xdata))
    expect_equivalent(xcms::groups(xset),
        as.matrix(xcms::featureDefinitions(xdata)[, 1:8]))
    
    expect_true(is(xdata, "XCMSnExp"))
    
    # Convert XCMSnExp to xcmsSet
    xdata <- as(xdata, "xcmsSet")
    expect_true(is(xdata, "xcmsSet"))
    
    expect_equivalent(xcms::peaks(xset), xcms::peaks(xdata))
    expect_equivalent(xcms::groups(xset), xcms::groups(xdata)[, 1:8])
    
    expect_equivalent(groupval(xset), groupval(xdata))
})