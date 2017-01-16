g <- group(faahko)
f <- fillPeaks(g)

testDiffreport <- function() {
    ## g <- group(faahko)
    ## f <- fillPeaks(g)
    d <- diffreport(f)
}

testDiffreportMinSamples <- function() {
    ## Test the case where one of the classes has just one sample
    ## Should skip the t-test, and proceed w/o error
    ## Reported by  Tomas Mikula
    ## g <- group(faahko)
    ## f <- fillPeaks(g)

    ## Fake xcmsSet with 1 sample in class 1
    sampclass(f) <- c(1,2,2,3,3,3,3,3,3,3,3,3)
    d <- diffreport(f, class1=1, class2=2)

    ## Test anova
    d <- diffreport(f)
}

testAnova <- function() {
    ## g <- group(faahko)
    ## f <- fillPeaks(g)
    sampclass(f) <- rep(1:4,3)
    d <- diffreport(f)
    d$anova
}
