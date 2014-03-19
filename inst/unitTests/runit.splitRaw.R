testSplitRawEven <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)

    xrl <- split(xraw, f=xraw@scanindex%%2)
    checkEqualsNumeric(length(xrl), 2)

    checkTrue(length(xraw@scanindex) == length(xrl[[1]]@scanindex) + length(xrl[[2]]@scanindex))
    checkTrue(length(xraw@env$mz) == length(xrl[[1]]@env$mz) + length(xrl[[2]]@env$mz))
}

testSplitRawOdd <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)

    xrl <- split(xraw, f=(xraw@scanindex+1)%%2)
    checkEqualsNumeric(length(xrl), 2)

    checkTrue(length(xraw@scanindex) == length(xrl[[1]]@scanindex) + length(xrl[[2]]@scanindex))
    checkTrue(length(xraw@env$mz) == length(xrl[[1]]@env$mz) + length(xrl[[2]]@env$mz))
}


testSplitRawFirst <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)

    xrl <- split(xraw, f=c(1,rep(2,length(xraw@scanindex)-1)))

    checkEqualsNumeric(length(xrl), 2)

}

testSplitRawLast <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)

    xrl <- split(xraw, f=c(rep(1,length(xraw@scanindex)-1), 2))

    checkEqualsNumeric(length(xrl), 2)
}

testSplitRawNone <- function() {

    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xraw <- xcmsRaw(file)

    xrl <- split(xraw, f=rep(1,length(xraw@scanindex)))

    checkEqualsNumeric(length(xrl), 1)

    xraw2 <- xrl[[1]]
    checkTrue(length(xraw@scanindex) == length(xrl[[1]]@scanindex))
    checkTrue(all(xraw@scanindex == xraw2@scanindex))

}
