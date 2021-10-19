test_that("matchpeaks doesn't fail", {
    skip_on_os(os = "windows", arch = "i386")

    faahko_file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")

    faahko_xs <- xcmsSet(faahko_file, profparam = list(step = 0),
                         method = "centWave", noise = 10000, snthresh = 40,
                         prefilter = c(6, 10000))
    pks <- peaks(faahko_xs)

    calibs <- pks[c(3, 5, 7, 13, 17, 29), "mz"]

    res <- xcms:::matchpeaks(pks, calibs)
    res_2 <- xcms:::matchpeaks(pks[order(pks[, "mz"]), ], calibs)
})
