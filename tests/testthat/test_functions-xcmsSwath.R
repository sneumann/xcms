test_that(".which_mz_in_range works", {
    mz <- 3.4
    lowerMz <- c(1, 3, 8, 12)
    upperMz <- c(3, 7, 11, 15)
    res <- .which_mz_in_range(mz, lowerMz, upperMz)
    expect_equal(res, 2L)
    res <- .which_mz_in_range(c(3, 3.4, 9), lowerMz, upperMz)
    expect_equal(res, list(c(1L, 2L), 2L, 3L))
})

test_that(".which_chrom_peak_overlap_rt works", {
    pks <- cbind(rtmin = c(1, 2, 3, 4, 5), rtmax = c(2, 3, 4, 5, 6))
    res <- .which_chrom_peak_overlap_rt(c(rtmin = 1.1, rtmax = 2.2), pks)
    expect_equal(res, c(1L, 2L))
    res <- .which_chrom_peak_overlap_rt(c(rtmin = 3.1, rtmax = 3.2), pks)
    expect_equal(res, 3L)
})

test_that(".which_chrom_peak_diff_rt works", {
    pks <- cbind(rt = c(1, 2, 3, 4, 5), rtmax = c(2, 3, 4, 5, 6))
    res <- .which_chrom_peak_diff_rt(c(rt = 1.1), pks, diffRt = 2)
    expect_equal(res, c(1L, 2L, 3L))
    res <- .which_chrom_peak_diff_rt(c(rt = 3.1), pks, diffRt = 0.101)
    expect_equal(res, 3L)
})

test_that(".reconstruct_ms2_for_chrom_peak works", {
    ## Fluopicolide, exact mass = 381.965430576, [M+H]+ = 382.972706
    pk <- chromPeaks(pest_swth, mz = 382.972706, ppm = 10)
    res <- xcms:::.reconstruct_ms2_for_chrom_peak(pk, pest_swth, fromFile = 7L,
                                           expandRt = 3, diffRt = 2,
                                           minCor = 0.8)
    expect_true(is(res, "MSpectra"))
    expect_equal(unname(fromFile(res)), 7L)
    expect_true(!is.na(rtime(res)))
    expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 14)
    expect_true(all(mcols(res)$ms2_peak_cor[[1]] > 0.8))
    expect_equal(length(intensity(res[[1]])), 14)

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
    expect_equal(length(intensity(res[[1]])), 8)
    expect_equal(length(mcols(res)$ms2_peak_id[[1]]), 8)
    expect_equal(length(mcols(res)$ms2_peak_cor[[1]]), 8)
})

test_that(".reconstruct_ms2_for_peaks_file works", {
    ## No MS2 level peaks: expect empty.
    res <- .reconstruct_ms2_for_peaks_file(filterFile(xod_x, 1))
    expect_true(all(isEmpty(res)))
    expect_true(is(mcols(res)$ms2_peak_id, "CharacterList"))
    expect_true(is(mcols(res)$ms2_peak_cor, "NumericList"))
    expect_identical(rownames(chromPeaks(filterFile(xod_x, 1))),
                     mcols(res)$peak_id)

    res <- .reconstruct_ms2_for_peaks_file(pest_swth, fromFile = 2L)
    expect_true(isEmpty(res)[1])
    expect_true(length(intensity(res[[3]])) == 8)
    expect_true(all(fromFile(res) == 2L))

    res_3 <- .reconstruct_ms2_for_peaks_file(pest_swth, fromFile = 2L,
                                                    peakId = "CP03")
    expect_identical(intensity(res_3), intensity(res[3]))
    expect_identical(mcols(res_3)$ms2_peak_id[[1]], mcols(res)$ms2_peak_id[[3]])
})
