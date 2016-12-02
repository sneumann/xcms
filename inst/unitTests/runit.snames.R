testAssignNames <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(files=file, snames=c("A"), method = "centWave", noise = 10000)
}

testAssignClass <- function() {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(files=file, sclass=c("A"), method = "centWave", noise = 10000)
}
