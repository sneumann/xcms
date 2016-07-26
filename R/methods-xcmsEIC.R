## Methods for xcmsEIC
#' @include functions-xcmsEIC.R

############################################################
## show
setMethod("show", "xcmsEIC", function(object) {

    cat("An \"xcmsEIC\" object with", length(object@eic), "samples\n\n")

    cat("Time range: ", paste(round(range(object@rtrange), 1), collapse = "-"),
        " seconds (", paste(round(range(object@rtrange)/60, 1), collapse = "-"),
        " minutes)\n", sep = "")
    cat("Mass range:", paste(round(range(object@mzrange, na.rm = TRUE), 4), collapse = "-"),
        "m/z\n")
    cat("EICs per sample:", nrow(object@mzrange), "\n")
    cat("Retention times:", object@rt, "\n")
    cat("xcmsSet group names: ")
    if (length(object@groupnames))
        cat("present")
    else
        cat("absent")
    cat("\n\n")
    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
})

############################################################
## groupnames
setMethod("groupnames", "xcmsEIC", function(object) {

    object@groupnames
})

############################################################
## sampnames
setMethod("sampnames", "xcmsEIC", function(object) {

    names(object@eic)
})

############################################################
## mzrange
setMethod("mzrange", "xcmsEIC", function(object) {

    object@mzrange
})

############################################################
## rtrange
setMethod("rtrange", "xcmsEIC", function(object) {

    object@rtrange
})
