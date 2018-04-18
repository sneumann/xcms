test_that("present absent sum works", {
    xsg <- group(faahko)
    a <- length(which(absent(xsg, class="WT", minfrac=0)))
    expect_equal(a, 407)
    p <- length(which(absent(xsg, class="KO", minfrac=0)))
    expect_equal(p, 407)
    ## Because minfrac = 0.5 translates to >= 3 samples,
    ## those peak groups with *exactly* 3 Na values
    ## are counted as both absent and present
    p <- length(which(present(xsg, class="WT", minfrac=0.5)))
    a <- length(which(absent(xsg , class="WT", minfrac=0.5)))
    expect_equal(a+p, 479)
    p <- length(which(present(xsg, class="KO", minfrac=0.333333333)))
    a <- length(which(absent(xsg,  class="KO", minfrac=0.666666666)))
    expect_true(a + p > nrow(groups(xsg)))
    ## minfrac=0.8,0.2, translate to 1.2 and 4.8 samples
    ##
    p <- length(which(present(xsg, class="WT", minfrac=0.8)))
    a <- length(which(absent(xsg, class="WT", minfrac=0.2)))
    expect_equal(a+p, 407)
    p <- length(which(present(xsg, class="WT", minfrac=0.2)))
    a <- length(which(absent(xsg, class="WT", minfrac=0.8)))
    expect_equal(a+p, 407)
    p <- length(which(present(xsg, class="KO", minfrac=0.75)))
    a <- length(which(absent(xsg, class="KO", minfrac=0.25)))
    expect_equal(a+p, 407)
    p <- length(which(present(xsg, class="KO", minfrac=0.25)))
    a <- length(which(absent(xsg, class="KO", minfrac=0.75)))
    expect_equal(a+p, 407)
})

##
## Same as above, with fillPeaks()
##
test_that("PresentAbsentSumAfterFillPeaks works", {
    xsg <- faahko_grouped_filled
    ## xsg <- faahko_grouped_filled
    a <- length(which(absent(xsg, class="WT", minfrac=0)))
    expect_equal(a, 407)
    p <- length(which(absent(xsg, class="KO", minfrac=0)))
    expect_equal(p, 407)
    ## Because minfrac = 0.5 translates to >= 3 samples,
    ## those peak groups with *exactly* 3 Na values
    ## are counted as both absent and present
    p <- length(which(present(xsg, class="WT", minfrac=0.5)))
    a <- length(which(absent(xsg , class="WT", minfrac=0.5)))
    expect_equal(a+p, 479)
    p <- length(which(present(xsg, class="KO", minfrac=0.333333333)))
    a <- length(which(absent(xsg,  class="KO", minfrac=0.666666666)))
    expect_true(a + p > nrow(groups(xsg)))
    ## minfrac=0.8,0.2, translate to 1.2 and 4.8 samples
    ##
    p <- length(which(present(xsg, class="WT", minfrac=0.8)))
    a <- length(which(absent(xsg, class="WT", minfrac=0.2)))
    expect_equal(a+p, 407)
    p <- length(which(present(xsg, class="WT", minfrac=0.2)))
    a <- length(which(absent(xsg, class="WT", minfrac=0.8)))
    expect_equal(a+p, 407)
    p <- length(which(present(xsg, class="KO", minfrac=0.75)))
    a <- length(which(absent(xsg, class="KO", minfrac=0.25)))
    expect_equal(a+p, 407)
    p <- length(which(present(xsg, class="KO", minfrac=0.25)))
    a <- length(which(absent(xsg, class="KO", minfrac=0.75)))
    expect_equal(a+p, 407)
})
