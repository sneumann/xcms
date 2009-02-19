testDiffreport <- function() {
    g <- group(faahko)
    f <- fillPeaks(g)
    d <- diffreport(f)
}

testAnova <- function() {
    g <- group(faahko)
    f <- fillPeaks(g)
    sampclass(f) <- rep(1:4,3)
    d <- diffreport(f)
    d$anova
}
