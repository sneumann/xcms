testGroupDensityOnePeak <- function() {
    xs <- faahko
    p <- peaks(xs)
    peaks(xs) <- p[1,,drop=FALSE]
    g <- group(xs, minsamp=1, minfrac=0.001, method="density")
}

## testGroupNearestOnePeak <- function() {
##   xs <- faahko
##   p <- peaks(xs)
##   peaks(xs) <- p[1,,drop=FALSE]
##   g <- group(xs, method="nearest")
## }

testGroupNearest<- function() {
  xs <- faahko
  p <- peaks(xs)
  g <- group(xs, method="nearest")
  checkEquals(range(unlist(g@groupidx))[1],  1)
}

