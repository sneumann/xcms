testGroup <- function() {
    classes <- sampclass(faahko)
    levels(classes) <- c(levels(classes), "Leftover")

    xs <- faahko
    sampclass(xs) <- classes
    xsg <- group(xs)

    checkEquals(sampclass(xs), sampclass(xsg))
}
