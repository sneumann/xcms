## Test updateObject.

test_updateObject_xcmsSet <- function() {
    library(msdata)
    data(xs)

    newXs <- updateObject(xs)
    checkIdentical(xs@peaks, newXs@peaks)
    ## xs might be an old object.
    if (!.hasSlot(xs, "scanrange"))
        checkIdentical(newXs@scanrange, numeric(0))
    checkIdentical(xs@groups, newXs@groups)
}
