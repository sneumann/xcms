## Methods for xcmsPeaks
#' @include DataClasses.R

#' @rdname hidden_aliases
setMethod("show", "xcmsPeaks", function(object) {
    cat("A matrix of", nrow(object), "peaks\n")
    cat("Column names:\n")
    print(colnames(object))
})
