
setClass("XcmsExperimentHdf5",
         contains = "XcmsExperiment",
         slots = c(hdf5_file = "character",
                   hdf5_mod_count = "integer",
                   sample_id = "integer"))

.h5_subset_xcms_experiment <- function(x, i = integer(),
                                       keepChromPeaks = TRUE,
                                       keepAdjustedRtime = FALSE,
                                       keepFeatures = FALSE,
                                       ignoreHistory = FALSE,
                                       keepSampleIndex = FALSE,
                                       ...) {
    i <- i2index(i, length(x))
    if (any(i < 0)) {
        if (all(i < 0))
            i <- seq_along(x)[i]
        else stop("Mixing positive and negative indices is not supported.")
    }
    stop("Needs to be implemented")
}

#' chromPeaks needs to iterate through all - will need also a chunkSize for
#' that.
