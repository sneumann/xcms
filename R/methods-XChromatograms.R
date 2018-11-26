setAs("Chromatograms", "XChromatograms", function(from) {
    res <- new("XChromatograms")
    res@.Data <- matrix(lapply(from, function(z) {
        if (is(z, "Chromatogram"))
            as(z, "XChromatogram")
        else z
    }), nrow = nrow(from), ncol = ncol(from), dimnames = dimnames(from))
    res@phenoData <- from@phenoData
    res@featureData <- from@featureData
    if (validObject(resetClass)) res
})

#' @rdname XChromatogram
setMethod("show", "XChromatograms", function(object) {
    nr <- nrow(object)
    nc <- ncol(object)
    cat(class(object), " with ",
        nr, ifelse(nr == 1, " row and ", " rows and "),
        nc, ifelse(nc == 1, " column\n", " columns\n"),
        sep = "")
    sumFun <- function(z) {
        paste0("peaks: ", nrow(z[[1]]@chromPeaks))
    }
    if (nr > 0 && nc > 0) {
        if (nr <= 4) {
            out <- apply(object, MARGIN = c(1, 2), sumFun)
            rownames(out) <- paste0("[", 1:nrow(out), ",]")
        }
        else {
            out <- rbind(
                apply(object[c(1, 2), , drop = FALSE], MARGIN = c(1, 2), sumFun),
                rep(" ... ", ncol(object)),
                apply(object[nrow(object) - c(1, 0), , drop = FALSE],
                      MARGIN = c(1, 2), sumFun)
            )
            rownames(out) <- c("[1,]", "[2,]", "...",
                               paste0("[", c(nrow(object) - c(1, 0)), ",]"))
        }
        rn <- rownames(out)
        out <- rbind(rep("<XChromatogram>", ncol(out)), out)
        rownames(out) <- c("", rn)
        print(out, quote = FALSE, right = TRUE)
    }
    cat("phenoData with", length(varLabels(object@phenoData)), "variables\n")
    cat("featureData with", length(fvarLabels(object)), "variables\n")
})

#' @rdname XChromatogram
setMethod("hasChromPeaks", "XChromatograms", function(object) {
    matrix(vapply(object, hasChromPeaks, logical(1)), ncol = ncol(object),
           dimnames = dimnames(object))
})
