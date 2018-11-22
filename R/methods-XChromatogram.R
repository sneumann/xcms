#' @rdname XChromatogram
setMethod("show", "XChromatogram", function(object) {
    callNextMethod()
    cat("Identified chromatographic peaks:\n")
    cat(" rt\trtmin\trtmax\tinto\tmaxo\tsn\n")
    for (i in seq_len(nrow(object@chromPeaks)))
        cat(" ", paste(object@chromPeaks[i, ], collapse = "\t"), sep = "")
    cat("\n")
})

#' @rdname XChromatogram
#'
#' @description
#'
#' - `chromPeaks`, `chromPeaks<-`: extract or set the matrix with the
#'   chromatographic peak definitions. Parameter `rt` allows to specify a
#'   retention time range for which peaks should be returned along with
#'   parameter `type` that defines how *overlapping* is defined (parameter
#'   description for details).
#'
#' @param rt For `chromPeaks`: `numeric(2)` defining the retention time range
#'     for which chromatographic peaks should be returned.
#'
#' @param type For `chromPeaks`: `character(1)` defining which peaks to return
#'     if `rt` is provided: `"any"` (default) return all peaks that are even
#'     partially overlapping with `rt`, `"within"` return peaks that are
#'     completely within `rt` and `"apex_within"` return peaks which apex
#'     is within `rt`.
#'
#' @param value For `chromPeaks<-`: a numeric `matrix` with required columns
#'     `"rt"`, `"rtmin"`, `"rtmax"`, `"into"` and `"maxo"`.
#'
#' @return
#'
#' For `chromPeaks`: a `matrix` with columns:
#' - `"rt"`: the retention time of the peak apex.
#' - `"rtmin"`: the lower peak boundary.
#' - `"rtmax"`: the upper peak boundary.
#' - `"into"`: the integrated peak signal (area of the peak).
#' - `"maxo"`: the maximum intensity of the peak.
#'
#' ## Extract the chromatographic peaks
#' chromPeaks(xchr)
setMethod("chromPeaks", "XChromatogram", function(object, rt = numeric(),
                                                  type = c("any", "within",
                                                           "apex_within")) {
    type <- match.arg(type)
    pks <- object@chromPeaks
    if (length(rt) && nrow(pks)) {
        rt <- range(rt)
        pks <- switch(type,
                      any = pks[which(pks[, "rtmin"] <= rt[2] &
                                      pks[, "rtmax"] >= rt[1]), , drop = FALSE],
                      within = pks[which(pks[, "rtmin"] >= rt[1] &
                                         pks[, "rtmax"] <= rt[2]), ,
                                   drop = FALSE],
                      apex_within = pks[which(pks[, "rt"] >= rt[1] &
                                              pks[, "rt"] <= rt[2]), ,
                                        drop = FALSE]
                      )
    }
    pks
})

#' @rdname XChromatogram
setReplaceMethod("chromPeaks", "XChromatogram", function(object, value) {
    if (!is.matrix(value))
        stop("'value' should be a numeric matrix")
    object@chromPeaks <- value
    if (validObject(object)) object
})
