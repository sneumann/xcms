
test_that(".correlate_chromatogram works", {
    set.seed(112)
    ## create
    rtime1 <- seq(5,25, by = 1)
    intensity1 <- dnorm(rtime1, mean=14, sd=1.0)*200

    rtime2 <- seq(5.1, 25.1, by = 1)
    intensity2 <- dnorm(rtime2, mean=14, sd=1.0)*500

    ## bogus chromatograms
    ch1 <- new("Chromatogram",
               rtime = rtime1,
               intensity = intensity1)

    ch2 <- new("Chromatogram",
               rtime = rtime2,
               intensity = intensity2)

    ## check that correlation with NA values fails
    expect_equal(xcms:::.correlate_chromatogram(ch1, ch2),
                 cor(intensity1, intensity2))
    expect_equal(xcms:::.correlate_chromatogram(ch1, ch2),
                 xcms:::.correlate_chromatogram(ch2, ch1))

    expect_equal(
        xcms:::.correlate_chromatogram(ch1, ch2, interpolate = TRUE),
        xcms:::.correlate_chromatogram(ch2, ch1, interpolate = TRUE),
        tolerance = 0.0001
    )

    ch2 <- filterRt(ch2, rt = c(5.1, 23))
    res <- xcms:::.correlate_chromatogram(ch1, ch2)
    expect_equal(res, cor(intensity(ch2), intensity(ch1)[1:length(ch2)]))

    res <- xcms:::.correlate_chromatogram(ch1, ch2, interpolate = TRUE,
                                          use = "everything")
    expect_equal(res, NA_real_)
})
