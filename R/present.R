setGeneric("present", function(object, class, minfrac) standardGeneric("present"))

setMethod("present", "xcmsSet", function(object, class, minfrac) {
    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]

    sampidx <- which(classlabel %in% class)

    if (length(sampidx) == 0) {
        stop("Class ", class, "not found")
    }

    classnum <- length(sampidx)
    minpresent <- classnum * minfrac

    filled <- rep(FALSE, nrow(peaks(object)))

    ## exists(object@filled) always returns FALSE ??
    sloti <- try(slot(object, "filled"), silent = TRUE)
    if (class(sloti) != "try-error") {
        filled[object@filled] <- TRUE
    }

    apply (groupval(object), 1, function(x) {
        length(which(  (!(is.na(x[sampidx]) | is.nan(x[sampidx])))
                     & !filled[x[sampidx]])) >= minpresent
    })
})

setGeneric("absent", function(object, class, minfrac) standardGeneric("absent"))
setMethod("absent", "xcmsSet", function(object, class, minfrac) {
    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]

    sampidx <- which(classlabel %in% class)

    if (length(sampidx) == 0) {
        stop("Class ", class, "not found")
    }

    classnum <- length(sampidx)
    minabsent <- classnum * minfrac

    filled <- rep(FALSE, nrow(peaks(object)))

    ## exists(object@filled) always returns FALSE ??
    sloti <- try(slot(object, "filled"), silent = TRUE)
    if (class(sloti) != "try-error") {
        filled[object@filled] <- TRUE
    }

    apply (groupval(object), 1, function(x) {
        length(which(is.na(x[sampidx]) | is.nan(x[sampidx]) | filled[x[sampidx]])) >= minabsent
    })
})
