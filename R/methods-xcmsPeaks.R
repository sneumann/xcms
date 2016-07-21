## Methods for xcmsPeaks
#' @include DataClasses.R

setMethod("show", "xcmsPeaks", function(object) {
    cat("A matrix of", nrow(object), "peaks\n")
    cat("Column names:\n")
    print(colnames(object))
})
