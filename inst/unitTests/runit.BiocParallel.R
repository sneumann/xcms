## Testing BiocParallel against the old parallel processing setup with
## manual setup of the cluster etc.
## NOTE: This test file might be removed after the old parallel processing
## code has been removed.

test_BiocParallel <- function() {

    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))

    ## Default.
    library(BiocParallel)
    bpp <- bpparam()
    bpprogressbar(bpp) <- TRUE
    xs2 <- xcmsSet(fs, BPPARAM=bpp)

    bpp <- SerialParam()
    xs3 <- xcmsSet(fs, BPPARAM=bpp)
    checkIdentical(xs2, xs3)
}
