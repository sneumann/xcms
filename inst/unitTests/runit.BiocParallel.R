## Testing BiocParallel against the old parallel processing setup with
## manual setup of the cluster etc.
## NOTE: This test file might be removed after the old parallel processing
## code has been removed.

test_BiocParallel <- function() {

    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))

    xs <- xcmsSet(fs, nSlaves = 2)

    ## Default.
    library(BiocParallel)
    bpp <- bpparam()
    bpprogressbar(bpp) <- TRUE
    xs2 <- xcms:::xcmsSet2(fs, BPPARAM=bpp)
    checkIdentical(xs, xs2)

    bpp <- SerialParam()
    bpprogressbar(bpp) <- TRUE
    xs2 <- xcms:::xcmsSet2(fs, BPPARAM=bpp)
    checkIdentical(xs, xs2)

    bpp <- MulticoreParam()
    bpprogressbar(bpp) <- TRUE
    xs2 <- xcms:::xcmsSet2(fs, BPPARAM=bpp)
    checkIdentical(xs, xs2)

    ## With fill peaks.
    xs <- group(xs)
    xsf <- fillPeaks.chrom(xs, nSlaves = 2)
    xsf2 <- xcms:::fillPeaks.chrom2(xs)
    checkIdentical(xsf, xsf2)

    ## Serial processing
    bpp <- SerialParam()
    bpprogressbar(bpp) <- TRUE
    xsf2 <- xcms:::fillPeaks.chrom2(xs, BPPARAM=bpp)
    checkIdentical(xsf, xsf2)

}
