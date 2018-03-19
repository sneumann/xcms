test_that("do_findPeaks_MSW works", {
    first_file <- filterFile(fticr, file = 1)
    spctr <- spectra(first_file)
    expect_true(length(spctr) == 1)
    mzs <- unname(mz(spctr[[1]]))
    ints <- unname(intensity(spctr[[1]]))
    feats1 <- do_findPeaks_MSW(mz = mzs[10000:20000],
                               int = ints[10000:20000],
                               snthresh = 100)
    feats2 <- do_findPeaks_MSW(mz = mzs[10000:20000],
                               int = ints[10000:20000],
                               snthresh = 50)
    expect_true(nrow(feats2) > nrow(feats1))
})
