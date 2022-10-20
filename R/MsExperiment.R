#' @rdname XcmsExperiment
setMethod("filterRt", "MsExperiment",
          function(object, rt = numeric(), ...) {
              message("Filter spectra")
              object <- .mse_filter_spectra(object, filterRt, rt = rt, ...)
              object
          })

#' @rdname XcmsExperiment
setMethod("filterFile", "MsExperiment", function(object, file = integer()) {
    object[i = sort(unique(file))]
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

## .mse_profmat_chunks <- function(x, method = "bin", step = 0.1,
##                                 baseleve = NULL, basespace = NULL,
##                                 mzrange. = NULL, fileIndex = seq_along(x),
##                                 chunkSize = 1L, ...) {
##     if (!length(spectra(x)))
##         stop("No spectra available.")
##     if (!any(names(x@sampleDataLinks) == "spectra"))
##         stop("No link between samples and spectra found. Please use ",
##              "'linkSampleData' to define which spectra belong to which ",
##              "samples.")
##     if (!all(fileIndex %in% seq_along(x)))
##         stop("fileIndex out of bounds", call. = FALSE)
##     ## Need to do it chunk-wise, again.
##     chunks <- split(fileIndex, ceiling(fileIndex / chunkSize))
##     keep <- x@sampleDataLinks[["spectra"]][, 1L] %in% fileIndex
##     sps <- spectra(x)[x@sampleDataLinks[["spectra"]][keep, 2L]]
##     sps$.SAMPLE_IDX <- x@sampleDataLinks[["spectra"]][keep, 1L]
##     res <- unlist(lapply(chunks, function(z) {
##         .mse_find_chrom_peaks_chunk(
##             sps[sps$.SAMPLE_IDX %in% z],
##             msLevel = msLevel, param = param,
##             BPPARAM = BPPARAM)
##     }), recursive = FALSE, use.names = FALSE)
##     sidx <- vapply(res, function(z) if (is.matrix(z)) nrow(z) else length(z),
##                    integer(1), USE.NAMES = FALSE)
##     res <- cbind(do.call(rbind, res), sample = rep(seq_along(x), sidx))
##     if (!length(res))
##         .empty_chrom_peaks()
##     else res
## }

## .mse_profmat_chunk <- function(x, method = "bin", step = 0.1,
##                                baseleve = NULL, basespace = NULL,
##                                mzrange. = NULL, fileIndex = seq_along(object),
##                                chunkSize = 1L) {
##     sidx <- unique(x@.SAMPLE_IDX)
##     x <- filterMsLevel(x, msLevel = 1L)
##     if (length(x))
##         f <- factor(x$.SAMPLE_IDX, levels = sidx)
##     else f <- factor(integer(), levels = sidx)
##     bplapply(
##         split(Spectra::peaksData(x, columns = c("mz", "intensity"),
##                                  BPPARAM = SerialParam()), f),
##         function(z) {
##             pk_count <- vapply(z, nrow, integer(1), USE.NAMES = FALSE)
##             empty_spectra <- which(!pk_count)
##             ## check profMat method.
##         }
##     )
## }
