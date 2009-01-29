testGroup <- function() {
    classes <- sampclass(faahko)
    levels(classes) <- c(levels(classes), "Leftover")

    sampclass(faahko) <- classes
    xsg <- group(faahko)

    checkEquals(sampclass(faahko), sampclass(xsg))    
}
