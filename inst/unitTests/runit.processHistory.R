############################################################
## Tests related to the processing history.
library(xcms)
library(RUnit)

############################################################
## Check constructor, errors etc.
test_ProcessHistory_class <- function() {
    ph <- xcms:::ProcessHistory()
    checkTrue(inherits(ph, "ProcessHistory"))

    ph@type <- "FAILURE"
    checkException(validObject(ph))
    checkException(xcms:::ProcessHistory(type = "any"))

    ## Accessor methods.
    checkException(xcms:::processType(ph) <- "OOOO")
    xcms:::processType(ph) <- xcms:::.PROCSTEP.UNKNOWN
    checkEquals(processType(ph), xcms:::.PROCSTEP.UNKNOWN)
    checkException(xcms:::processDate(ph) <- c("a", "b"))
    xcms:::processDate(ph) <- "A"
    checkEquals(processDate(ph), "A")
    checkException(xcms:::processInfo(ph) <- c("a", "b"))
    xcms:::processInfo(ph) <- "B"
    checkEquals(processInfo(ph), "B")
    xcms:::fileIndex(ph) <- 1:3
    checkEquals(fileIndex(ph), 1:3)
}

############################################################
## That's testing on real failing files (from issue #55)
dontrun_test_ProcessHistory_failing_files <- function() {
    mzfs <- paste0("../../local_data/mzML-files/",
                   c("123.mzML", "261.mzML", "263.mzML"))
    suppressWarnings(
        res <- xcmsSet(mzfs, method = "centWave", ppm = 2.5,
                       peakwidth = c(2.5, 9), mzdiff = -0.001,
                       snthresh = 10, stopOnError = FALSE)
    )
    ## We expect that the last two files generated an error.
    ## Get the ProcessHistory objects with an error.
    errs <- xcms:::.getProcessErrors(res)
    checkEquals(length(errs), 2)
    ## evaluate the showError method
    errs <- showError(res)
    checkEquals(length(errs), 2)
    checkTrue(all(is.character(unlist(errs))))
}

############################################################
## Check showError method
test_showError <- function() {

    library(msdata)
    data(xs)

    errs <- xcms:::.getProcessErrors(xs)
    checkEquals(length(errs), 0)
    ph <- xcms:::.getProcessHistory(xs)
    checkEquals(length(ph), 0)
    xs <- updateObject(xs)
    checkTrue(.hasSlot(xs, ".processHistory"))

    errs <- xcms:::.getProcessErrors(xs)
    checkEquals(length(errs), 0)

    errs <- showError(xs)
    checkEquals(length(errs), 0)

    ph <- xcms:::.getProcessHistory(xs, fileIndex = 3)
    checkEquals(length(ph), 0)

}

############################################################
## Check concatenating xcmsSet objects.
test_ProcessHistory_c_xcmsSet <- function() {
}


############################################################
## Test XProcessHistory
test_XProcessHistory_class <- function() {
    ph <- xcms:::XProcessHistory()
    checkTrue(is(ph, "XProcessHistory"))
    checkTrue(inherits(ph, "ProcessHistory"))

    ph <- xcms:::XProcessHistory(info = "some info",
                                 type = xcms:::.PROCSTEP.PEAK.DETECTION)
    checkEquals(ph@info, "some info")
    checkEquals(ph@type, xcms:::.PROCSTEP.PEAK.DETECTION)

    ph@type <- "other"
    checkException(validObject(ph))

    ph <- xcms:::XProcessHistory(info = "some info",
                                 type = xcms:::.PROCSTEP.PEAK.DETECTION,
                                 param = CentWaveParam())

    checkTrue(is(ph@param, "CentWaveParam"))
    checkTrue(is(processParam(ph), "CentWaveParam"))
}
