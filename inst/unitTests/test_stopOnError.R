############################################################
## Unit tests related to stopIfNot parameters (translated to
## into the stop.if.not BPPARAM argument.

dontrun_test_xcmsSet_stopIfNot <- function() {

    ## Run the feature detection using xcmsSet in parallel mode
    ## and evaluate the stopOnError parameter
    ## We need some problematic test files here; we're using the ones
    ## from issue #55

    mzfs <- paste0("../../local_data/mzML-files/",
                   c("123.mzML", "261.mzML", "263.mzML"))

    ## Default: stops on the first encountered error.
    checkException(xcmsSet(mzfs, method = "centWave", ppm = 2.5,
                           peakwidth = c(2.5, 9), mzdiff = -0.001,
                           snthresh = 10, stopOnError = TRUE))
    suppressWarnings(
        res <- xcmsSet(mzfs, method = "centWave", ppm = 2.5,
                       peakwidth = c(2.5, 9), mzdiff = -0.001,
                       snthresh = 10, stopOnError = FALSE)
    )
    checkTrue(all(res@peaks[, "sample"] == 1))
}


