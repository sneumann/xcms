test_that(".which_mz_in_range works", {
    skip_on_os(os = "windows", arch = "i386")

    mz <- 3.4
    lowerMz <- c(1, 3, 8, 12)
    upperMz <- c(3, 7, 11, 15)
    res <- .which_mz_in_range(mz, lowerMz, upperMz)
    expect_equal(res, 2L)
    res <- .which_mz_in_range(c(3, 3.4, 9), lowerMz, upperMz)
    expect_equal(res, list(c(1L, 2L), 2L, 3L))
})

test_that(".which_chrom_peak_overlap_rt works", {
    skip_on_os(os = "windows", arch = "i386")

    pks <- cbind(rtmin = c(1, 2, 3, 4, 5), rtmax = c(2, 3, 4, 5, 6))
    res <- .which_chrom_peak_overlap_rt(c(rtmin = 1.1, rtmax = 2.2), pks)
    expect_equal(res, c(1L, 2L))
    res <- .which_chrom_peak_overlap_rt(c(rtmin = 3.1, rtmax = 3.2), pks)
    expect_equal(res, 3L)
})

test_that(".which_chrom_peak_diff_rt works", {
    skip_on_os(os = "windows", arch = "i386")

    pks <- cbind(rt = c(1, 2, 3, 4, 5), rtmax = c(2, 3, 4, 5, 6))
    res <- .which_chrom_peak_diff_rt(c(rt = 1.1), pks, diffRt = 2)
    expect_equal(res, c(1L, 2L, 3L))
    res <- .which_chrom_peak_diff_rt(c(rt = 3.1), pks, diffRt = 0.101)
    expect_equal(res, 3L)
})

test_that(".reconstruct_ms2_for_chrom_peak works", {
    skip_on_os(os = "windows", arch = "i386")

    ## Fluopicolide, exact mass = 381.965430576, [M+H]+ = 382.972706
    pk <- chromPeaks(pest_swth, mz = 382.972706, ppm = 10)
    res <- .reconstruct_ms2_for_chrom_peak(pk, pest_swth, fromFile = 7L,
                                           expandRt = 3, diffRt = 2,
                                           minCor = 0.8)
    expect_true(is(res, "MSpectra"))
    expect_equal(unname(fromFile(res)), 7L)
    expect_true(!is.na(rtime(res)))
    expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 1)
    expect_true(all(mcols(res)$ms2_peak_cor[[1]] > 0.8))
    expect_equal(length(intensity(res[[1]])), 1)

    res <- .reconstruct_ms2_for_chrom_peak(pk, pest_swth, fromFile = 7L,
                                           expandRt = 3, diffRt = 2,
                                           minCor = 1)
    expect_true(is(res, "MSpectra"))
    expect_equal(unname(fromFile(res)), 7L)
    expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 0)
    expect_true(all(isEmpty(res)))

    ## Some random other peak.
    pk <- chromPeaks(pest_swth, msLevel = 1)[1, ]
    res <- .reconstruct_ms2_for_chrom_peak(pk, pest_swth)
    expect_equal(length(intensity(res[[1]])), 0)
    expect_true(is(mcols(res)$ms2_peak_id, "CharacterList"))
    expect_true(is(mcols(res)$ms2_peak_cor, "NumericList"))
    expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 0)
    expect_equal(length(mcols(res)$ms2_peak_cor[[1]]), 0)

    pk <- chromPeaks(pest_swth, msLevel = 1)[3, ]
    res <- .reconstruct_ms2_for_chrom_peak(pk, pest_swth)
    expect_equal(length(intensity(res[[1]])), 0)
    expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 0)
    expect_equal(length(mcols(res)$ms2_peak_cor[[1]]), 0)
})

test_that(".reconstruct_ms2_for_peaks_file works", {
    skip_on_os(os = "windows", arch = "i386")

    ## No MS2 level peaks: expect empty.
    tmp <- filterRt(filterFile(xod_x, 1), rt = c(2500, 3000))
    res <- xcms:::.reconstruct_ms2_for_peaks_file(tmp)
    expect_true(all(isEmpty(res)))
    expect_true(is(mcols(res)$ms2_peak_id, "CharacterList"))
    expect_true(is(mcols(res)$ms2_peak_cor, "NumericList"))
    expect_identical(rownames(chromPeaks(tmp)),
                     mcols(res)$peak_id)

    res_2 <- .reconstruct_ms2_for_peaks_file(pest_swth, fromFile = 2L)
    expect_true(isEmpty(res_2)[1])
    expect_true(length(intensity(res_2[[7]])) == 5)
    expect_true(all(fromFile(res_2) == 2L))

    res_3 <- .reconstruct_ms2_for_peaks_file(
                        pest_swth, fromFile = 2L,
                        peakId = c("CP03", "CP04", "CP07"))
    expect_identical(intensity(res_3)[1L], intensity(res_2[3]))
    expect_identical(intensity(res_3)[2L], intensity(res_2[4]))
    expect_identical(intensity(res_3)[[3L]], intensity(res_2[7])[[1L]])
    expect_identical(mcols(res_3)$ms2_peak_id[[1L]],
                     mcols(res_2)$ms2_peak_id[[3L]])
    expect_identical(mcols(res_3)$ms2_peak_id[[2L]],
                     mcols(res_2)$ms2_peak_id[[4L]])

    ## Same with return.type = "Spectra"
    res_s <- .reconstruct_ms2_for_peaks_file(tmp, return.type = "Spectra")
    expect_true(is(res_s, "Spectra"))
    expect_true(all(msLevel(res_s) == 2L))
    expect_equal(res_s$peak_id, mcols(res)$peak_id)
    expect_equal(rtime(res_s), unname(rtime(res)))
    expect_equal(precursorMz(res_s), unname(precursorMz(res)))
    expect_true(all(isEmpty(res_s)))

    res_s2 <- .reconstruct_ms2_for_peaks_file(pest_swth, fromFile = 2L,
                                                     return.type = "Spectra")
    expect_true(isEmpty(res_s2)[1])
    expect_true(length(intensity(res_s2)[[7L]]) == 5)
    expect_true(all(res_s2$fromFile == 2L))

    res_s3 <- .reconstruct_ms2_for_peaks_file(
                        pest_swth, fromFile = 2L,
                        peakId = c("CP03", "CP04", "CP07"),
                        return.type = "Spectra")
    expect_identical(intensity(res_s3)[1L], intensity(res_s2[3]))
    expect_identical(intensity(res_s3)[2L], intensity(res_s2[4]))
    expect_identical(intensity(res_s3)[[3L]], intensity(res_s2[7])[[1L]])
    expect_identical(res_s3$ms2_peak_id[[1L]],
                     mcols(res_2)$ms2_peak_id[[3L]])
    expect_identical(res_s3$ms2_peak_id[[2L]],
                     mcols(res_2)$ms2_peak_id[[4L]])
})

test_that(".data2spectra works", {
    skip_on_os(os = "windows", arch = "i386")

    ## empty one, no mz, no intensity, just rt fromFile etc.
    res <- .data2spectra(rt = 3.2, polarity = 1L, precursorMz = 12.3,
                         precursorIntensity = 1234, return.type = "MSpectra")
    expect_true(is(res, "MSpectra"))
    expect_true(validObject(res))
    expect_equal(unname(rtime(res)), 3.2)
    expect_equal(mcols(res)$ms2_peak_cor[[1L]], numeric())

    res <- .data2spectra(rt = 3.2, polarity = 1L, precursorMz = 12.3,
                         precursorIntensity = 1234, return.type = "MSpectra",
                         mz = c(1, 2, 3, 4), intensity = c(12, 13, 14, 15),
                         peak_id = c("a", "b", "c", "d"))
    expect_true(validObject(res))
    expect_equal(mz(res)[[1L]], c(1, 2, 3, 4))
    expect_equal(mcols(res)$ms2_peak_id[[1L]], c("a", "b", "c", "d"))
    expect_true(centroided(res))
    expect_equal(unname(polarity(res)), 1L)

    ## DataFrame
    res <- .data2spectra(rt = 3.2, polarity = 1L, precursorMz = 12.3,
                         precursorIntensity = 1234, return.type = "Spectra",
                         mz = c(1, 2, 3, 4), intensity = c(12, 13, 14, 15),
                         peak_id = c("a", "b", "c", "d"))
    expect_true(is(res, "DataFrame"))
    expect_equal(res$mz, NumericList(c(1, 2, 3, 4), compress = FALSE))

    res <- .data2spectra(rt = 3.2, polarity = 1L, precursorMz = 12.3,
                         precursorIntensity = 1234, return.type = "Spectra")
    expect_true(is(res, "DataFrame"))
    expect_equal(res$mz, NumericList(c(), compress = FALSE))
})
