testAssignNames <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(files=file, snames=c("A"))
}

testAssignClass <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(files=file, sclass=c("A"))
}
