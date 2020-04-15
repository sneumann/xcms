## This was copied from MsCoreUtils - remove after the package was added to
## Bioconductor!

rbindFill <- function(...) {
    l <- list(...)

    if (length(l) == 1L)
        l <- l[[1L]]

    cnms <- c("matrix", "data.frame", "DataFrame", "DFrame")
    cls <- vapply(l, inherits, integer(length(cnms)), what = cnms, which = TRUE)
    rownames(cls) <- cnms

    if (any(!as.logical(colSums(cls))))
        stop("'rbindFill' just works for ", paste(cls, collapse = ", "))

    ## convert matrix to data.frame for easier and equal subsetting and class
    ## determination
    isMatrix <- as.logical(cls["matrix",])
    l[isMatrix] <- lapply(l[isMatrix], as.data.frame)

    allcl <- unlist(
        lapply(l, function(ll) {
            vapply(ll, function(lll)class(lll)[1L], character(1),
                   USE.NAMES = TRUE)
        })
    )
    allnms <- unique(names(allcl))
    allcl <- allcl[allnms]

    for (i in seq(along=l)) {
        diffcn <- setdiff(allnms, names(l[[i]]))
        if (length(diffcn))
            l[[i]][, diffcn] <- lapply(allcl[diffcn], as, object = NA)
    }
    r <- do.call(rbind, l)

    ## if we had just matrices as input we need to convert our temporary
    ## data.frame back to a matrix
    if (all(isMatrix))
        r <- as.matrix(r)
    r
}

## helper function to allow lapply(..., as, ...) in rbindFill
setAs("logical", "factor", function(from, to) factor(from))
