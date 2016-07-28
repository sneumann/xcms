## Test functions to evaluate binning of values.

.notrun_yet <- function() {
    x <- abs(rnorm(10000, mean = 200, sd = 300))

    ## R binning.
    binR <- function(x, binSize = 1L) {
        nVals <- length(x)
        br <- seq(floor(min(x)), ceiling(max(x)), by = binSize)
        idx <- findInterval(x, br)
        out <- double(length(br))
        ## split does already do the sorting for us.
        out[unique(idx)] <- unlist(lapply(split(x, idx), mean), use.names = FALSE)
        return(cbind(bin = brMidpoint(br), x = out))
    }

    brMidpoint <- function(x) {
        return(c((x[-length(x)]+x[-1L])/2L, x[length(x)]))
    }

    brMidpoint2 <- function(x, binS) {
        if(missing(binS))
            binS <- unique(diff(x))
        return((x+binS/2)[-length(x)])
    }
}
