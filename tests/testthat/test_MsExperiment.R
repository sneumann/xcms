library(MsExperiment)
mse <- MsExperiment()
fls <- normalizePath(faahko_3_files)
df <- data.frame(mzML_file = basename(fls),
                 dataOrigin = fls,
                 sample = c("ko15", "ko16", "ko18"))

spectra(mse) <- Spectra::Spectra(fls)
sampleData(mse) <- DataFrame(df)
## Link samples to spectra.
mse <- linkSampleData(mse, with = "sampleData.dataOrigin = spectra.dataOrigin")

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
