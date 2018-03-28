test_that(".getPeaks_xxx functions works", {
    ## Compare the old and new getPeaks implementations.
    xs_m <- xcmsSet(faahko_3_files[1])

    pks_range <- peaks(xs_m)[1:200, ]
    ## Extend the range
    pks_range[, "mzmin"] <- pks_range[, "mzmin"] - 0.05
    pks_range[, "mzmax"] <- pks_range[, "mzmax"] + 0.05
    expect_warning(pks_o <- .getPeaks_orig(faahko_xr_1, peakrange = pks_range))
    pks_n <- .getPeaks_new(faahko_xr_1, peakrange = pks_range)
    expect_equal(pks_o, pks_n)

    pks_tmp <- pks_o
    ## Force it to use different step.
    expect_warning(pks_o <- .getPeaks_orig(faahko_xr_1, peakrange = pks_range,
                                           step = 0.3))
    pks_n <- .getPeaks_new(faahko_xr_1, peakrange = pks_range, step = 0.3)
    expect_equal(pks_o, pks_n)
    expect_true(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
    
    ## Change profile generation settings.
    tmp <- deepCopy(faahko_xr_1)
    tmp@profmethod <- "binlin"
    expect_warning(pks_o <- .getPeaks_orig(tmp, peakrange = pks_range,
                                           step = 0.2))
    pks_n <- .getPeaks_new(tmp, peakrange = pks_range, step = 0.2)
    ## Can not expect identical values because of differences in binlin
    ## See issues #46 and #49.
    expect_true(cor(pks_o[, "into"], pks_n[, "into"]) > 0.999)
    expect_true(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
    pks_tmp <- pks_o
    
    ## Change profile generation settings.
    tmp@profmethod <- "binlinbase"
    expect_warning(pks_o <- .getPeaks_orig(tmp, peakrange = pks_range,
                                           step = 0.2))
    pks_n <- .getPeaks_new(tmp, peakrange = pks_range, step = 0.2)
    expect_equal(pks_o, pks_n)
    expect_true(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
    pks_tmp <- pks_o

    tmp@profmethod <- "intlin"
    pks_o <- .getPeaks_orig(tmp, peakrange = pks_range, step = 0.2)
    pks_n <- .getPeaks_new(tmp, peakrange = pks_range, step = 0.2)
    expect_equal(pks_o, pks_n)
    expect_true(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
})

test_that("xcmsSet can handle MS2 data", {
    filename <- system.file('iontrap/extracted.mzData', package = "msdata")
    expect_warning(xs2 <- xcmsSet(filename, snthresh = 4, mslevel = 2))
})

test_that("xcmsSet works with MS2... again", {
    filename <- system.file('iontrap/extracted.mzData', package = "msdata")
    expect_warning(xs2 <- xcmsSet(filename, method="centWave", mslevel = 2))
})

test_that("phenoDataFromPaths and others don't fail", {
    files <- system.file(c("cdf/KO/ko15.CDF", "cdf/KO/ko16.CDF",
                           "cdf/KO/ko18.CDF", "cdf/KO/ko19.CDF",
                           "cdf/KO/ko21.CDF", "cdf/KO/ko22.CDF",
                           "cdf/WT/wt15.CDF", "cdf/WT/wt16.CDF",
                           "cdf/WT/wt18.CDF", "cdf/WT/wt19.CDF",
                           "cdf/WT/wt21.CDF", "cdf/WT/wt22.CDF"),
                         package = "faahKO")
    pd <- phenoDataFromPaths(files)
    xs <- faahko

    ##xcms::phenoData(xs) <- pd
    ## https://stat.ethz.ch/pipermail/r-devel/2008-April/049184.html
    xs <- xcms::`phenoData<-`(xs, pd)
    
    xsg <- group(xs)

    pd <- phenoDataFromPaths(files)
    xs <- faahko

    ##xcms::phenoData(xs) <- pd
    ## https://stat.ethz.ch/pipermail/r-devel/2008-April/049184.html
    xs <- xcms::`phenoData<-`(xs, pd)
    xs <- group(xs)
    ## Setting the filepaths again; otherwise we will have problem finding these
    ## files ... obviously.
    filepaths(xs) <- files
    xs <- fillPeaks(xs)
    dr <- diffreport(xs, class1="KO", class2="WT")
    expect_true(nrow(dr) > 0)
})

test_that("showError works", {
    data(xs)
    errs <- .getProcessErrors(xs)
    expect_equal(length(errs), 0)
    ph <- .getProcessHistory(xs)
    expect_equal(length(ph), 0)
    xs <- updateObject(xs)
    expect_true(.hasSlot(xs, ".processHistory"))

    errs <- .getProcessErrors(xs)
    expect_equal(length(errs), 0)

    errs <- showError(xs)
    expect_equal(length(errs), 0)

    ph <- .getProcessHistory(xs, fileIndex = 3)
    expect_equal(length(ph), 0)
})

test_that("split.xcmsSet works", {
    xsl <- split(faahko, sampclass(faahko))
    expect_equal(length(xsl), 2)
    xsl <- split(faahko, sampnames(faahko))
    expect_equal(length(xsl), length(sampnames(faahko)))
    xsl <- split(faahko, c(1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
    expect_equal(length(xsl), 2)
    expect_equal(length(sampnames(xsl[[1]])), 1)
    expect_equal(length(sampnames(xsl[[2]])), 11)
    expect_equal(sampnames(xsl[[1]]), sampnames(faahko)[1])
    xsl <- split(faahko, c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1))
    expect_equal(length(xsl), 2)
    expect_equal(length(sampnames(xsl[[1]])), 1)
    expect_equal(length(sampnames(xsl[[2]])), 11)
    expect_equal(sampnames(xsl[[1]]), sampnames(faahko)[12])
    xsl <- split(faahko, c(2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2))
    expect_equal(length(xsl), 2)
    expect_equal(sampnames(xsl[[1]]), sampnames(faahko)[6])
    xsl <- split(faahko,c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))
    expect_equal(length(xsl), 1)
    xsl <- split(faahko, sampclass(faahko))
    expect_equal(length(sampnames(c(xsl[[1]], xsl[[2]]))), 12)
    ## short factor:
    f = c(1, 2, 2)
    xsl <- split(faahko, f)
    for (x in unique(f)) {
        num.samps = sum(x==f)
        expect_equal(length(xsl[[x]]@filepaths), num.samps)
        expect_equal(nrow(xsl[[x]]@phenoData), num.samps)
        expect_equal(length(xsl[[x]]@rt$raw), num.samps)
    }
    ## long factor
    f = c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3)
    actual.f = f[1:nrow(faahko@phenoData)]
    xsl <- split(faahko, f)
    for (x in unique(actual.f)) {
        num.samps = sum(x==actual.f)
        expect_equal(length(xsl[[x]]@filepaths), num.samps)
        expect_equal(nrow(xsl[[x]]@phenoData), num.samps)
        expect_equal(length(xsl[[x]]@rt$raw), num.samps)
    }
    ## check phenoData
    xset <- faahko
    ## make a dummy phenoDate data.frame
    pd <- data.frame(class = xset$class, exp = rep(1:6, 2))
    rownames(pd) <- sampnames(xset)
    phenoData(xset) <- pd
    ## split by sampleclass
    xsetList <- split(xset, sampclass(xset))
    expect_equal(xsetList[[1]]$class, droplevels(pd[pd$class=="KO", "class"]))
    expect_equal(xsetList[[2]]$exp, pd[pd$class=="WT", "exp"])
    ## now split by the experiment
    xsetList <- split(xset, xset$exp)
    ## just check some stuff...
    expect_equal(length(xsetList[[1]]$class), 2)
    expect_equal(nrow(phenoData(xsetList[[4]])), 2)
    ## have to make sure we also do have the mslevel, profinfo,
    ## scanrange, dataCorrection and polarity set.
    expect_equal(mslevel(xset), mslevel(xsetList[[1]]))
    expect_equal(scanrange(xset), scanrange(xsetList[[1]]))
    expect_equal(xset@polarity, xsetList[[1]]@polarity)
})

test_that("split.xcmsSet preserves ProcessHistory", {
    spl <- split(faahko, sampclass(faahko))
    expect_true(length(spl[[1]]@.processHistory) == 0)

    spl <- split(xs_1, f = c(1, 2, 1, 2))
    ph <- xs_1@.processHistory[c(1, 3)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 2L
    expect_equal(spl[[1]]@.processHistory, ph)

    ph <- xs_1@.processHistory[c(2, 4)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 2L
    expect_equal(spl[[2]]@.processHistory, ph)

    ## Add fake ProcessHistory steps.
    ph <- xs_1@.processHistory
    ph <- c(list(ProcessHistory(fileIndex = 1:4)), ph,
            list(ProcessHistory()))
    xs_2 <- xs_1
    xs_2@.processHistory <- ph
    ##
    spl <- split(xs_2, f = c(2, 1, 1, 1))
    ph <- xs_2@.processHistory[c(1, 3, 4, 5)]
    ph[[1]]@fileIndex <- 1L:3L
    ph[[2]]@fileIndex <- 1L
    ph[[3]]@fileIndex <- 2L
    ph[[4]]@fileIndex <- 3L
    expect_equal(spl[[1]]@.processHistory, ph)

    ph <- xs_2@.processHistory[c(1, 2)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 1L
    expect_equal(spl[[2]]@.processHistory, ph)
})

test_that("c.xcmsSet preserves ProcessHistory", {
    spl <- split(faahko, sampclass(faahko))
    conc <- do.call("c", spl)
    expect_equal(length(conc@.processHistory), 0)

    spl <- split(xs_1, c(1, 1, 2, 2))
    conc <- c(spl[[1]], spl[[2]])
    expect_equal(xs_1@.processHistory, conc@.processHistory)

    ## Different ordering
    xs_1.1 <- xs_1[, c(1, 3)]
    xs_1.2 <- xs_1[, 2]
    xs_1.3 <- xs_1[, 4]
    ## Add a fake processing for the second one.
    ph <- xs_1.2@.processHistory
    ph <- c(list(ProcessHistory(fileIndex = 1)), ph)
    xs_1.2@.processHistory <- ph
    expect_true(.validProcessHistory(xs_1.2))
    ## Combine them.
    conc <- c(xs_1.1, xs_1.2, xs_1.3)
    ## 1st
    expect_equal(conc@.processHistory[[1]], xs_1@.processHistory[[1]])
    ## 2nd
    ph <- xs_1@.processHistory[[3]]
    ph@fileIndex <- 2L
    expect_equal(conc@.processHistory[[2]], ph)
    ## 3rd
    ph <- xs_1.2@.processHistory[[1]]
    ph@fileIndex <- 3L
    expect_equal(conc@.processHistory[[3]], ph)
    ## 4th
    ph <- xs_1@.processHistory[[2]]
    ph@fileIndex <- 3L
    expect_equal(conc@.processHistory[[4]], ph)
    ## 5th
    ph <- xs_1@.processHistory[[4]]
    ph@fileIndex <- 4L
    expect_equal(conc@.processHistory[[5]], ph)
    ## empty
    library(msdata)
    suppressWarnings(
        xs <- xcmsSet(system.file("microtofq/MM8.mzML", package="msdata"),
                      method="centWave", ppm=25, peakwidth=c(20, 50))
    )
    xs2 <- xcmsSet(system.file("microtofq/MM14.mzML", package="msdata"),
                   method="centWave", ppm=25, peakwidth=c(20, 50))
    comb <- c(xs, xs2)
    expect_true(nrow(peaks(comb)) == 0)
})

test_that(".getProcessHistory works", {
    expect_equal(.getProcessHistory(xs_1, fileIndex = 2:3),
                 xs_1@.processHistory[2:3])
})

test_that("xcmsSet centWave works", {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xset1 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      profparam = list(step = 0))
    xset2 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      profparam = list(step = 0), scanrange=c(1,112))
    xset3 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      profparam = list(step = 0), scanrange=c(1,80))
    expect_true(nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))
    expect_true(nrow((peaks(xset1)@.Data))  > nrow((peaks(xset3)@.Data)))
})

test_that("xcmsSet matchedFilter works", {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xset1 <- xcmsSet(files=file, method="matchedFilter", fwhm=10)
    xset2 <- xcmsSet(files=file, method="matchedFilter", fwhm=10,
                     scanrange=c(1,112))
    xset3 <- xcmsSet(files=file, method="matchedFilter", fwhm=10,
                     scanrange=c(1,80))
    expect_true(nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))
    expect_true(nrow((peaks(xset1)@.Data))  > nrow((peaks(xset3)@.Data)))
})

test_that("xcmsSet parallel works", {
    file <- system.file('microtofq/MM14.mzdata', package = "msdata")
    xset1 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      scanrange=c(1,80), profparam = list(step = 0))
    xset2 <-  xcmsSet(files=file, method="centWave", peakwidth=c(5,12),
                      scanrange=c(1,80), profparam = list(step = 0))
    ## parallel disabled: , nSlaves=2)
    expect_true(nrow((peaks(xset1)@.Data)) == nrow((peaks(xset2)@.Data)))
})

test_that("phenoDataFromPaths works", {
    base_dir <- system.file("cdf", package = "faahKO")
    cdf_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)
    pd <- phenoDataFromPaths(cdf_files)
    expect_true(colnames(pd) == "class")
    expect_equal(levels(pd$class), c("KO", "WT"))
})
