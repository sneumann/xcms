#' @rdname XChromatogram
setMethod("show", "XChromatogram", function(object) {
    callNextMethod()
    cat("Identified chromatographic peaks:\n")
    cat(" rt\trtmin\trtmax\tinto\tmaxo\tsn\n")
    for (i in seq_len(nrow(object@chromPeaks)))
        cat(" ", paste(object@chromPeaks[i, ], collapse = "\t"), sep = "")
    cat("\n")
})
