context("Spectra-methods")

sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
           precursorMz = 2, rt = 1.232446)
sp3 <- new("Spectrum1", mz = c(1, 2, 3, 5, 6), intensity = c(6:10),
           rt = 1.232445)
spl_ <- Spectra(sp1, sp2, sp3)

test_that("mz, intensity, rtime work", {
    spl <- spl_
    expect_true(length(mz(spl)) == 3)
    expect_equal(mz(spl[[2]]), mz(spl)[[2]])
    expect_equal(mz(spl)[[3]], mz(sp3))
    
    expect_true(length(intensity(spl)) == 3)
    expect_equal(intensity(spl[[2]]), intensity(spl)[[2]])
    expect_equal(intensity(spl)[[3]], intensity(sp3))
    
    expect_true(length(rtime(spl)) == 3)
    expect_equal(rtime(spl[[2]]), rtime(spl)[[2]])
    expect_equal(rtime(spl)[[3]], rtime(sp3))
    
    spl <- c(spl, Spectra(new("Spectrum2")))
    expect_true(lengths(mz(spl))[4] == 0)
    expect_true(lengths(intensity(spl))[4] == 0)
    expect_equal(rtime(spl)[4], NA_real_)

    ## Put names on it.
    names(spl) <- c("a", "b", "c", "d")
    expect_equal(names(rtime(spl)), c("a", "b", "c", "d"))
    expect_equal(names(mz(spl)), c("a", "b", "c", "d"))
    expect_equal(names(intensity(spl)), c("a", "b", "c", "d"))

    ## Empty spectra
    spl <- Spectra(new("Spectrum1"), new("Spectrum2"))
    expect_equal(rtime(spl), c(NA_real_, NA_real_))
    expect_true(length(mz(spl)) == 2)
    expect_true(all(lengths(mz(spl)) == 0))
    expect_true(length(intensity(spl)) == 2)
    expect_true(all(lengths(intensity(spl)) == 0))
})

test_that("precursor* work", {
    sp1 <- new("Spectrum2", precursorMz = 123.3, precursorCharge = 1L,
               precursorIntensity = 1234.4)
    sp2 <- new("Spectrum2", precursorMz = NA_real_, precursorCharge = integer(),
               precursorIntensity = NA_real_, precScanNum = 34L)
    spl <- Spectra(sp1, sp2, sp3)
    expect_equal(precursorMz(spl), c(123.3, NA, NA))
    expect_equal(precursorCharge(spl), c(1L, NA, NA))
    expect_true(is.integer(precursorCharge(spl)))
    expect_equal(precursorIntensity(spl), c(1234.4, NA, NA))
    expect_equal(precScanNum(spl), c(NA_integer_, 34L, NA_integer_))
    expect_true(is.integer(precScanNum(spl)))
    
    expect_equal(precursorMz(spl_), c(NA_real_, 2, NA_real_))
    expect_equal(precursorCharge(spl_), rep(NA_integer_, length(spl_)))
    expect_equal(precursorIntensity(spl_), rep(NA_real_, length(spl_)))
    expect_true(is.integer(precScanNum(spl_)))
})

test_that("acquisitionNum and scanIndex work", {
    sp1 <- new("Spectrum2", acquisitionNum = 2L, scanIndex = 1L)
    sp2 <- new("Spectrum2", acquisitionNum = 4L)
    spl <- Spectra(sp1, sp2)
    expect_identical(acquisitionNum(spl), c(2L, 4L))
    expect_identical(scanIndex(spl), c(1L, NA_integer_))
        
    expect_equal(acquisitionNum(spl_), rep(NA_integer_, length(spl_)))
    expect_equal(scanIndex(spl_), rep(NA_integer_, length(spl_)))
    expect_true(is.integer(acquisitionNum(spl_)))
    expect_true(is.integer(scanIndex(spl_)))
})

test_that("peaksCount, msLevel, tic and ionCount work", {
    sp1 <- new("Spectrum2", msLevel = 3L, tic = 5)
    sp2 <- new("Spectrum2")
    spl <- Spectra(sp1, sp2)

    expect_true(is.integer(peaksCount(spl)))
    expect_equal(peaksCount(spl), c(0, 0))
    expect_true(is.integer(msLevel(spl)))
    expect_equal(msLevel(spl), c(3, 2))
    expect_true(is.numeric(tic(spl)))
    expect_equal(tic(spl), c(5, 0))
    expect_true(is.numeric(ionCount(spl)))
    expect_equal(ionCount(spl), c(0, 0))

    expect_equal(peaksCount(spl_), c(3, 4, 5))
    expect_equal(msLevel(spl_), c(2, 2, 1))
    expect_equal(tic(spl_), c(11, 15, 40))
    expect_equal(ionCount(spl_), unlist(lapply(intensity(spl_), sum)))
})

test_that("collisionEnergy works", {
    sp1 <- new("Spectrum2")
    sp2 <- new("Spectrum2", collisionEnergy = 23.3)
    spl <- Spectra(sp1, sp2)

    expect_true(is.numeric(collisionEnergy(spl)))
    expect_equal(collisionEnergy(spl), c(NA, 23.3))

    expect_true(is.numeric(collisionEnergy(spl_)))
    expect_equal(collisionEnergy(spl_), c(NA_real_, NA_real_, NA_real_))
})

test_that("fromFile and polarity work", {
    sp1 <- new("Spectrum2", polarity = 1L, fromFile = 5L)
    sp2 <- new("Spectrum2", fromFile = 3L)
    spl <- Spectra(sp1, sp2)

    expect_true(is.integer(fromFile(spl)))
    expect_equal(fromFile(spl), c(5, 3))
    expect_true(is.integer(polarity(spl)))
    expect_equal(polarity(spl), c(1, NA))

    expect_true(is.integer(fromFile(spl_)))
    expect_true(all(is.na(fromFile(spl_))))
    expect_true(is.integer(polarity(spl_)))
    expect_equal(polarity(spl_), rep(NA_integer_, 3))
})

test_that("smoothed, isEmpty, centroided and isCentroided work", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(4, 2, 4, 5))
    sp2 <- new("Spectrum2", centroided = TRUE, smoothed = TRUE)
    spl <- Spectra(sp1, sp2)

    expect_true(is.logical(smoothed(spl)))
    expect_equal(smoothed(spl), c(NA, TRUE))
    expect_true(is.logical(isEmpty(spl)))
    expect_equal(isEmpty(spl), c(FALSE, TRUE))
    expect_true(is.logical(centroided(spl)))
    expect_equal(centroided(spl), c(NA, TRUE))
    expect_true(is.logical(isCentroided(spl)))
    expect_equal(isCentroided(spl), c(NA, NA))

    expect_true(is.logical(smoothed(spl_)))
    expect_equal(smoothed(spl_), rep(NA, length(spl_)))
    expect_true(is.logical(isEmpty(spl_)))
    expect_true(all(!isEmpty(spl_)))
    expect_true(is.logical(centroided(spl_)))
    expect_true(all(is.na(centroided(spl_))))
    expect_true(is.logical(isCentroided(spl_)))
    expect_equal(isCentroided(spl_), c(NA, NA, NA))
})

test_that("writeMgfData,Spectra works", {
    tmpf <- tempfile()

    writeMgfData(spl_, tmpf)
    res <- readLines(tmpf)
    ## No additional fields.
    expect_equal(res[7], "1 4")
    expect_equal(res[17], "1 5")
    expect_equal(res[27], "1 6")

    expect_error(writeMgfData(spl_, tmpf))

    spl <- spl_
    mcols(spl) <- DataFrame(index = 1:3, some_id = c("sp_1", "sp_2", "sp_3"))
    file.remove(tmpf)

    writeMgfData(spl, tmpf)
    res_2 <- readLines(tmpf)
    expect_equal(res_2[7], "INDEX=1")
    expect_equal(res_2[8], "SOME_ID=sp_1")
    expect_equal(res_2[19], "INDEX=2")
    expect_equal(res_2[20], "SOME_ID=sp_2")
    expect_equal(res_2[31], "INDEX=3")
    expect_equal(res_2[32], "SOME_ID=sp_3")
    expect_equal(res[-1], res_2[-c(1, 7, 8, 19, 20, 31, 32)])
})
