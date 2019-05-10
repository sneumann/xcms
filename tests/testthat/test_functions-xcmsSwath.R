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
})
