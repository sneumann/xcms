#' @rdname XcmsExperiment
setMethod("filterRt", "MsExperiment",
          function(object, rt = numeric(), ...) {
              message("Filter spectra")
              object <- .mse_filter_spectra(object, filterRt, rt = rt, ...)
              object
          })

#' @rdname XcmsExperiment
setMethod("filterMzRange", "MsExperiment",
          function(object, mz = numeric(), msLevel. = uniqueMsLevels(object)) {
              message("Filter spectra")
              object@spectra <- filterMzRange(object@spectra, mz, msLevel.)
              object
          })

#' @rdname XcmsExperiment
setMethod("filterMz", "MsExperiment",
          function(object, mz = numeric(), msLevel. = uniqueMsLevels(object)) {
              filterMzRange(object, mz, msLevel.)
          })

#' @rdname XcmsExperiment
setMethod("uniqueMsLevels", "MsExperiment", function(object) {
    uniqueMsLevels(spectra(object))
})

#' @rdname XcmsExperiment
setMethod("filterFile", "MsExperiment", function(object,
                                                 file = integer(), ...) {
    object[i = sort(unique(file)), ...]
})

#' @rdname profMat-xcmsSet
setMethod("profMat", "MsExperiment", function(object,
                                              method = "bin",
                                              step = 0.1,
                                              baselevel = NULL,
                                              basespace = NULL,
                                              mzrange. = NULL,
                                              fileIndex = seq_along(object),
                                              chunkSize = 1L, msLevel = 1L,
                                              BPPARAM = bpparam(), ...) {
    .mse_profmat_chunks(object, msLevel = msLevel, method = method, step = step,
                        baselevel = baselevel, basespace = basespace,
                        mzrange. = mzrange., fileIndex = fileIndex,
                        chunkSize = chunkSize, BPPARAM = BPPARAM, ...)
})

################################################################################
## These functions below are needed to re-use code from the xcms package
## developed for OnDiskMSnExp/XCMSnExp objects for MsExperiment objects. They
## are NOT indended to go to the MsExperiment package as they do not make full
## use of the new data structure.

#' @rdname XcmsExperiment
setMethod("rtime", "MsExperiment", function(object) {
    if (length(spectra(object)))
        rtime(spectra(object))
    else numeric()
})

#' @rdname XcmsExperiment
setMethod("fromFile", "MsExperiment", function(object) {
    if (length(spectra(object))) {
        .mse_check_spectra_sample_mapping(object)
        object@sampleDataLinks[["spectra"]][, 1L]
    } else integer()
})

#' @rdname XcmsExperiment
setMethod("fileNames", "MsExperiment", function(object) {
    if (length(spectra(object)))
        unique(dataOrigin(spectra(object)))
    else character()
})

#' @rdname XcmsExperiment
setMethod("polarity", "MsExperiment", function(object) {
    if (length(spectra(object)))
        polarity(spectra(object))
    else integer()
})

#' @rdname XcmsExperiment
setMethod(
    "filterIsolationWindow", "MsExperiment", function(object, mz = numeric()) {
        .mse_filter_spectra(object, filterIsolationWindow, mz = mz)
    })

#' @rdname XcmsExperiment
setMethod(
    "chromatogram", "MsExperiment",
    function(object, rt = matrix(nrow = 0, ncol = 2),
             mz = matrix(nrow = 0, ncol = 2), aggregationFun = "sum",
             msLevel = 1L, isolationWindowTargetMz = NULL, chunkSize = 2L,
             return.type = "MChromatograms", BPPARAM = bpparam()) {
        if (!is.matrix(rt))
            rt <- matrix(rt, ncol = 2L)
        if (!is.matrix(mz))
            mz <- matrix(mz, ncol = 2L)
        .mse_chromatogram(
            object, rt = rt, mz = mz, aggregationFun = aggregationFun,
            msLevel = msLevel, isolationWindow = isolationWindowTargetMz,
            chunkSize = chunkSize, BPPARAM = BPPARAM)
    })
