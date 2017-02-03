## Tests related to profile matrix generation and update. Tests in this file
## cover
## o .createProfileMatrix internal function.
## o profMat method for xcmsRaw.
## o profStep<- method for xcmsRaw.
## o profMethod<- method for xcmsRaw.
## o profMethod for XCMSnExp and OnDiskMSnExp

## library(faahKO)
fs <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
## xr <- deepCopy(faahko_xr_1)
## xr <- xcmsRaw(fs, profstep = 0)

## .createProfileMatrix
test_createProfileMatrix <- function() {
    xr <- deepCopy(faahko_xr_1)
    mz <- xr@env$mz
    int <- xr@env$intensity
    numPerSc <- diff(c(xr@scanindex, length(xr@env$mz)))
    ## Testing all properties.
    ## o bin
    pm <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                      valsPerSpect = numPerSc,
                                      method = "bin", step = 2)
    ## o binlin
    pm_2 <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                        valsPerSpect = numPerSc,
                                        method = "binlin", step = 2)
    checkEquals(dim(pm), dim(pm_2))
    ## o binlinbase
    pm_3 <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                        valsPerSpect = numPerSc,
                                        method = "binlinbase", step = 2)
    checkEquals(dim(pm), dim(pm_3))
    checkEquals(sum(pm == 0), sum(pm_3 == 35))
    ##   setting parameter: baselevel
    pm_3_2 <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                        valsPerSpect = numPerSc,
                                        method = "binlinbase", step = 2,
                                        baselevel = 666666)
    checkEquals(sum(pm_3_2 == 666666), sum(pm == 0))
    ##   setting parameter: basespace
    pm_3_3 <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                        valsPerSpect = numPerSc,
                                        method = "binlinbase", step = 2,
                                        basespace = 0.5)
    checkEquals(pm_3_3, pm_3)
    pm_3_3 <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                        valsPerSpect = numPerSc,
                                        method = "binlinbase", step = 2,
                                        basespace = 300)
    checkTrue(!all(pm_3_3 == pm_3))
    ## o intlin
    pm_4 <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                        valsPerSpect = numPerSc,
                                        method = "intlin", step = 2)
    checkEquals(dim(pm), dim(pm_4))
}

## profMat
test_profMat <- function() {
    ## Create a new profile matrix:

    xr <- deepCopy(faahko_xr_1)
    checkException(pm <- profMat(xr))

    xr_2 <- xcmsRaw(fs, profstep = 2)
    checkEquals(profMat(xr_2), profMat(xr, step = 2))
    xr_3 <- xcmsRaw(fs, profstep = 2, profmethod = "binlinbase",
                    profparam = list(baselevel = 666666))
    checkEquals(xr_3@env$profile, profMat(xr, step = 2, method = "binlinbase",
                                          baselevel = 666666))
    checkEquals(xr_3@env$profile, profMat(xr_3))
}

test_profMat_OnDiskMSnExp <- function() {
    ## Get it from all 3 files in one go.
    res <- profMat(faahko_od, step = 2)
    res_2 <- profMat(xcmsRaw(faahko_3_files[2], profstep = 0), step = 2)
    checkEquals(res_2, res[[2]])
    res_2 <- profMat(xcmsRaw(faahko_3_files[3], profstep = 0), step = 2)
    checkEquals(res_2, res[[3]])

    res_2 <- profMat(faahko_xod, step = 2)
    checkEquals(res, res_2)

    res <- profMat(faahko_od, step = 2, method = "binlin", fileIndex = 2)
    res_2 <- profMat(xcmsRaw(faahko_3_files[2], profstep = 0), step = 2,
                     method = "binlin")
    checkEquals(res_2, res[[1]])
}


## profStep<-
test_profStepReplace <- function() {
    ## Profile matrix will be generated/replaced if the step parameter is > 0
    ## and differs from the one within the object.
    ## xr <- xcmsRaw(fs, profstep = 0)
    xr <- deepCopy(faahko_xr_1)
    xr_2 <- xr
    checkTrue(length(xr_2@env$profile) == 0)
    profStep(xr_2) <- 2
    checkTrue(length(xr_2@env$profile) > 0)
    checkEquals(profMat(xr_2), profMat(xr, step = 2))
    checkEquals(profMat(xr_2, step = 2), xr_2@env$profile)
}

## profMethod<-
test_profMethodReplace <- function() {
    ## Profile matrix will be generated/replaced if profMethod is changed.
    xr <- deepCopy(faahko_xr_1)
    xr_2 <- xr
    checkTrue(length(xr_2@env$profile) == 0)
    ## Just setting profMethod doesn't help here
    profMethod(xr_2) <- "binlin"
    checkEquals(profMethod(xr_2), "binlin")
    checkTrue(length(xr_2@env$profile) == 0)
    profStep(xr_2) <- 2
    xr_3 <- xcmsRaw(fs, profstep = 2, profmethod = "binlin")
    checkEquals(profMat(xr_3), profMat(xr_2))
    ## binlinbase
    xr_3@profparam <- list(baselevel = 666666)
    profMethod(xr_3) <- "binlinbase"
    checkEquals(xr_3@env$profile, profMat(xr_3, method = "binlinbase",
                                          baselevel = 666666))
    xr_4 <- xcmsRaw(fs, profstep = 2, profmethod = "binlinbase",
                    profparam = list(baselevel = 666666))
    checkEquals(xr_4@env$profile, xr_3@env$profile)
}

## This test compares the new methods to the old deprecated one.
dontrun_createProfileMatrix <- function() {
    xr <- deepCopy(faahko_xr_1)
    mz <- xr@env$mz
    int <- xr@env$intensity
    numPerSc <- diff(c(xr@scanindex, length(xr@env$mz)))
    ## profile matrix:
    pm <- xcms:::.createProfileMatrix(mz = mz, int = int,
                                      valsPerSpect = numPerSc,
                                      method = "bin", step = 0.1)
    ## Do the same with profBin.
    minmass <- round(min(xr@env$mz) / 0.1) * 0.1
    maxmass <- round(max(xr@env$mz) / 0.1) * 0.1
    num <- (maxmass - minmass) / 0.1 + 1
    pm_2 <- xcms:::profBinM(xr@env$mz, xr@env$intensity, xr@scanindex, num)
    checkEquals(pm, pm_2)

    xr_2 <- xcmsRaw(fs, profstep = 2, profmethod = "binlinbase",
                    profparam = list(basespace = 0.5))
    profinfo(xr_2)
}

