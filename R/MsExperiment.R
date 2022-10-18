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
