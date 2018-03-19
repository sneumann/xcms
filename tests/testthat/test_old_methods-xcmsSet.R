test_that("diffreport etc works", {
    g <- group(faahko)
    f <- fillPeaks(g)
    d <- diffreport(f)
    ## Fake xcmsSet with 1 sample in class 1
    sampclass(f) <- c(1,2,2,3,3,3,3,3,3,3,3,3)
    d <- diffreport(f, class1=1, class2=2)

    ## Test anova
    d <- diffreport(f)
    sampclass(f) <- rep(1:4,3)
    d <- diffreport(f)
    d$anova
})

test_that("findPeaks.MSW works", {
    ## We do expect an error if we have multiple spectra (issue #237)
    expect_error(findPeaks.MSW(faahko_xr_1))
})
