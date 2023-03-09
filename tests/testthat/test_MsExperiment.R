library(MsExperiment)
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))
mse <- readMsExperiment(spectraFiles = fls, sampleData = df)

test_that("filterRt,MsExperiment works", {
    res <- filterRt(mse, rt = c(2700, 2900))
    expect_true(all(rtime(spectra(res)) > 2700 & rtime(spectra(res)) < 2900))
    b <- spectra(mse[2])
    B <- spectra(res[2])
    expect_equal(rtime(B), rtime(filterRt(b, rt = c(2700, 2900))))

    res <- filterRt(mse, rt = c(2700, 2900), msLevel = 2L)
    expect_equal(rtime(spectra(res)), rtime(spectra(mse)))
})

test_that("filterFile,MsExperiment works", {
    res <- filterFile(mse)
    expect_s4_class(res, "MsExperiment")
    expect_true(length(res) == 0)

    res <- filterFile(mse, 2)
    expect_equal(res, mse[2])
    res <- filterFile(mse, c(3, 1))
    expect_equal(res, mse[c(1, 3)])
})

test_that("rtime,MsExperiment works", {
    res <- rtime(MsExperiment())
    expect_identical(res, numeric())

    res <- rtime(mse)
    expect_true(is.numeric(res))
    expect_identical(res, rtime(spectra(mse)))
})

test_that("fromFile,MsExperiment works", {
    expect_identical(fromFile(MsExperiment()), integer())

    res <- fromFile(mse)
    expect_equal(res, mse@sampleDataLinks[["spectra"]][, 1L])
})

test_that("fileNames,MsExperiment works", {
    expect_identical(fileNames(MsExperiment()), character())

    res <- fileNames(mse)
    expect_equal(res, unique(dataOrigin(spectra(mse))))
})

test_that("chromatogram,MsExperiment works", {
    res <- chromatogram(mse)
    expect_s4_class(res, "MChromatograms")
    expect_equal(length(mse), ncol(res))
    expect_equal(1L, nrow(res))
    ref <- chromatogram(faahko_od)
    expect_equal(rtime(res[1, 1]), unname(rtime(ref[1, 1])))
    expect_equal(rtime(res[1, 2]), unname(rtime(ref[1, 2])))
    expect_equal(rtime(res[1, 3]), unname(rtime(ref[1, 3])))
    expect_equal(intensity(res[1, 1]), unname(intensity(ref[1, 1])))
    expect_equal(intensity(res[1, 2]), unname(intensity(ref[1, 2])))
    expect_equal(intensity(res[1, 3]), unname(intensity(ref[1, 3])))

    ## Subset.
    res <- chromatogram(mse, rt = c(10, 3000))
    expect_s4_class(res, "MChromatograms")
    expect_equal(length(mse), ncol(res))
    expect_equal(1L, nrow(res))
    expect_true(all(rtime(res[1, 1]) <= 3000))

    res <- chromatogram(mse, msLevel = 2L)
    expect_s4_class(res, "MChromatograms")
    expect_equal(length(mse), ncol(res))
    expect_equal(1L, nrow(res))
    expect_equal(intensity(res[1, 1]), numeric())
    expect_equal(intensity(res[1, 2]), numeric())
    expect_equal(intensity(res[1, 2]), numeric())
})

test_that("uniqueMsLevels,MsExperiment works", {
    expect_equal(uniqueMsLevels(mse), 1L)
})

test_that("filterMzRange,MsExperiment works", {
    res <- filterMzRange(mse[1L], mz = c(100, 500))
    mzs <- unlist(mz(res))
    expect_true(all(mzs <= 500))
    expect_true(all(mzs >= 100))
    expect_warning(
        res <- filterMzRange(mse[1L], mz = c(100, 500), msLevel. = 2L),
        "not available")
    mzs <- unlist(mz(res))
    expect_true(any(mzs > 500))
})
