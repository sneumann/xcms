testFilledFlag <- function() {
    xsg <- group(faahko)
    xsgf <- fillPeaks(xsg, method="chrom")

    checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled), nrow(peaks(xsgf)))
}

dontrun_testFillPeaksPar <- function() {
    xsg <- group(faahko)

    xsgfSerial <- fillPeaks(xsg, method="chrom")
    xsgfParallel <- fillPeaks(xsg, method="chrom") # parallel disabled: , nSlaves=2)

    checkEqualsNumeric(nrow(peaks(xsgfSerial)),
                       nrow(peaks(xsgfParallel)))
}



test.fillPeaksColumns <- function() {
    xsg <- group(faahko)
    xsg <- group(faahko_xs)
    peaks(xsg) <- cbind(peaks(xsg), anotherColumn=4711)

    oldCnames <- colnames(peaks(xsg))
    xsgf <- fillPeaks(xsg) # parallel disabled: , nSlaves=2)

    newCnames <- colnames(peaks(xsgf))
    checkEquals(oldCnames, newCnames)

    ## Check dims if nothing to do
    oldDims <- dim(peaks(xsgf))
    xsgf2 <- fillPeaks(xsgf) # parallel disabled: , nSlaves=2)
    newDims <- dim(peaks(xsgf2))
    checkEquals(oldDims, newDims)

    ## Case where only some samples have NA values
    xsg <- group(faahko_xs, minfrac=1)
    xsgf <- fillPeaks(xsg) # parallel disabled: , nSlaves=2)
    sampclass(xsgf) <- c(rep("KO", 1), rep("WT", 2))
    xsgf <- group(xsgf, minfrac=1)
    xsgf <- fillPeaks(xsgf) # parallel disabled: , nSlaves=2)
}

test.getPeaks_implementation <- function() {
    ## Compare the old and new getPeaks implementations.
    xs_m <- xcmsSet(faahko_3_files[1])

    pks_range <- peaks(xs_m)[1:200, ]
    ## Extend the range
    pks_range[, "mzmin"] <- pks_range[, "mzmin"] - 0.05
    pks_range[, "mzmax"] <- pks_range[, "mzmax"] + 0.05
    suppressWarnings(
        pks_o <- xcms:::.getPeaks_orig(faahko_xr_1, peakrange = pks_range)
    )
    pks_n <- xcms:::.getPeaks_new(faahko_xr_1, peakrange = pks_range)
    checkEquals(pks_o, pks_n)

    pks_tmp <- pks_o
    ## Force it to use different step.
    suppressWarnings(
        pks_o <- xcms:::.getPeaks_orig(faahko_xr_1, peakrange = pks_range, step = 0.3)
    )
    pks_n <- xcms:::.getPeaks_new(faahko_xr_1, peakrange = pks_range, step = 0.3)
    checkEquals(pks_o, pks_n)
    checkTrue(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
    
    ## Change profile generation settings.
    tmp <- deepCopy(faahko_xr_1)
    tmp@profmethod <- "binlin"
    suppressWarnings(
        pks_o <- xcms:::.getPeaks_orig(tmp, peakrange = pks_range, step = 0.2)
    )
    pks_n <- xcms:::.getPeaks_new(tmp, peakrange = pks_range, step = 0.2)
    ## Can not expect identical values because of differences in binlin
    ## See issues #46 and #49.
    checkTrue(cor(pks_o[, "into"], pks_n[, "into"]) > 0.999)
    checkTrue(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
    pks_tmp <- pks_o
    
    ## Change profile generation settings.
    tmp@profmethod <- "binlinbase"
    suppressWarnings(
        pks_o <- xcms:::.getPeaks_orig(tmp, peakrange = pks_range, step = 0.2)
    )
    pks_n <- xcms:::.getPeaks_new(tmp, peakrange = pks_range, step = 0.2)
    checkEquals(pks_o, pks_n)
    checkTrue(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
    pks_tmp <- pks_o

    tmp@profmethod <- "intlin"
    suppressWarnings(
        pks_o <- xcms:::.getPeaks_orig(tmp, peakrange = pks_range, step = 0.2)
    )
    pks_n <- xcms:::.getPeaks_new(tmp, peakrange = pks_range, step = 0.2)
    checkEquals(pks_o, pks_n)
    checkTrue(sum(pks_o[, "into"] != pks_tmp[, "into"]) > 0)
 }

## Compare the results we get when running the old and new fillPeaks.
dontrun_test.fillPeaks_old_vs_new <- function() {
    xsg <- group(faahko, minfrac = 1)

    register(SerialParam())
    res_n <- fillPeaks(xsg)
    useOriginalCode(TRUE)
    res_o <- fillPeaks(xsg)
    useOriginalCode(FALSE)
    pks_n <- peaks(res_n)[res_n@filled, ]
    pks_o <- peaks(res_o)[res_o@filled, ]
    checkTrue(cor(pks_o[pks_o[, "sample"] == 7, "into"],
                  pks_n[pks_n[, "sample"] == 7, "into"]) > 0.999)
    ## plot(pks_n[, "into"], pks_o[, "into"])
    
    profinfo(xsg) <- list(method = "binlin", step = 0.2)
    res_n <- fillPeaks(xsg)
    useOriginalCode(TRUE)
    res_o <- fillPeaks(xsg)
    useOriginalCode(FALSE)
    pks_n <- peaks(res_n)[res_n@filled, ]
    pks_o <- peaks(res_o)[res_o@filled, ]
    checkTrue(cor(pks_o[pks_o[, "sample"] == 7, "into"],
                  pks_n[pks_n[, "sample"] == 7, "into"]) > 0.999)
}


## testFilledFlagMSW <- function() {

##   xsg <- group(ham)
##   xsgf <- fillPeaks(xsg)

##   checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled, nrow(peaks(xsgf))) )
## }
