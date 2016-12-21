## Methods for the XCMSnExp object representing untargeted metabolomics
## results

setMethod("initialize", "XCMSnExp", function(.Object, ...) {
    classVersion(.Object)["XCMSnExp"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

##' @rdname XCMSnExp-class
##' @noRd
setMethod("show", "XCMSnExp", function(object) {
    cat("Object of class: ", class(object), "\n")
    callNextMethod()
})
