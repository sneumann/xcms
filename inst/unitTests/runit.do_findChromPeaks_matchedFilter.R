## Testo findChromPeaks matchedFilter
## library(xcms)
## library(RUnit)

## library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))
## library(msdata)
## f <- msdata::proteomics(full.names = TRUE, pattern = "TMT_Erwinia")
## mzf <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
##          system.file("microtofq/MM8.mzML", package = "msdata"))


test_do_findChromPeaks_matchedFilter <- function() {
    ## xr <- xcmsRaw(fs[1], profstep = 0)
    xr <- deepCopy(faahko_xr_1)
    ## We expect that changing a parameter has an influence on the result.
    mzVals <- xr@env$mz
    intVals <- xr@env$intensity
    ## Define the values per spectrum:
    valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
    res1 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = xr@scantime,
                                            valsPerSpect,
                                            binSize = 10)
    res2 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = xr@scantime,
                                            valsPerSpect,
                                            binSize = 10,
                                            snthresh = 100)
    checkTrue(nrow(res1) > nrow(res2))
    res2 <- do_findChromPeaks_matchedFilter(mz = mzVals,
                                            int = intVals,
                                            scantime = xr@scantime,
                                            valsPerSpect,
                                            binSize = 20)
    checkTrue(nrow(res1) > nrow(res2))
}

## Evaluate the peak detection method using matchedFilter on MSnExp and
## OnDiskMSnExp objects. For now we can't read CDF files, so we have to restrict
## to provided mzML files.
test_findChromPeaks_matchedFilter <- function() {
    library(MSnbase)
    mfp <- MatchedFilterParam(binSize = 20, impute = "lin")
    res <- xcmsSet(fs[1], method = "matchedFilter", profmethod = "binlin",
                   step = binSize(mfp))
    ## onDisk
    ## onDisk <- readMSData(fs[1], mode = "onDisk")
    onDisk <- filterFile(faahko_od, file = 1)
    res_o <- findChromPeaks(onDisk, param = mfp, return.type = "xcmsSet")
    checkEquals(peaks(res_o), peaks(res))
    checkEquals(res_o@rt$raw, res@rt$raw, checkNames = FALSE)

    checkException(findChromPeaks(onDisk, param = mfp, msLevel = 2))
    ## inMem
    ## inMem <- readMSData(mzf, msLevel. = 1)
    ## res_i <- findChromPeaks(inMem, param = mfp, return.type = "xcmsSet")
    ## checkEquals(peaks(res_i), peaks(res))
    ## checkEquals(res_i@rt$raw, res@rt$raw, checkNames = FALSE)

    ## xs <- xcmsSet(fs, , method = "matchedFilter", profmethod = "binlin",
    ##               step = binSize(mfp))
    ## onDisk <- readMSData(fs, mode = "onDisk")
    ## res <- findChromPeaks(onDisk, param = mfp)
    ## checkTrue(hasChromPeaks(res))
    ## checkTrue(!hasAdjustedRtime(res))
    ## checkTrue(!hasFeatures(res))
    ## checkEquals(peaks(xs)@.Data, chromPeaks(res))
    ## checkEquals(processParam(processHistory(res)[[1]]), mfp)
}

## Some benchmarks
dontrun_benchmark_detecfFeatures_matchedFilter <- function() {
    library(microbenchmark)
    library(MSnbase)
    mfp <- MatchedFilterParam(binSize = 0.2, impute = "lin")
    onDisk <- readMSData(mzf, mode = "onDisk")
    inMem <- readMSData(mzf, msLevel. = 1)
    microbenchmark(xcmsSet(mzf, method = "matchedFilter", profmethod = "binlin",
                           step = binSize(mfp)),
                   findChromPeaks(onDisk, param = mfp, return.type = "xcmsSet"),
                   findChromPeaks(inMem, param = mfp, return.type = "xcmsSet"),
                   times = 3)
    ## netCDF.
    onDisk <- readMSData(fs, mode = "onDisk")
    inMem <- readMSData(fs, msLevel. = 1)
    microbenchmark(xcmsSet(fs, method = "matchedFilter", profmethod = "binlin",
                           step = binSize(mfp)),
                   findChromPeaks(onDisk, param = mfp, return.type = "xcmsSet"),
                   findChromPeaks(inMem, param = mfp, return.type = "xcmsSet"),
                   times = 3)
}

############################################################
## Compare each individual function to the original one changing
## settings.
## Comparing each of the functions to the original one:
## A: do_findChromPeaks_matchedFilter (original code)
## B: .matchedFilter_binYonX_iter
## C: .matchedFilter_no_iter
## D: .matchedFilter_binYonX_no_iter
## This is also discussed on issue #47 on github:
## https://github.com/sneumann/xcms/issues/47
## A description of the results is provided in section "Implementation and
## comparison for matchedFilter" section of "new_functionality.org".
dontrun_test_do_findChromPeaks_matchedFilter_impl <- function() {

    library(xcms)
    library(RUnit)

    library(faahKO)
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))
    i <- 1

    cat("Comparison of results from different implementations:\n")
    cat("- orig: the original findPeaks.matchedFilter method.\n")
    cat("- A: do_findChromPeaks_matchedFilter (containing original code).\n")
    cat(paste0("- B: .matchedFilter_binYonX_iter: new function using binYonX",
               " for binning and imputeLinInterpol for interpolation. Uses",
               " iterative buffering like the original code."))
    cat(paste0("- C: .matchedFilter_no_iter: original code but without",
               " iterative buffering."))
    cat(paste0("- D: .matchedFilter_binYonX_no_iter: new code without",
               " iterative buffering.\n\n"))

    for (i in 1:length(fs)) {
        cat("============================================================\n")
        cat("|   file", i, ":", fs[i], "\n")
        cat("------------------------------------------------------------")

        xr <- xcmsRaw(fs[i])
        mz <- xr@env$mz
        int <- xr@env$intensity
        scantime <- xr@scantime
        scanindex <- xr@scanindex
        valsPerSpect <- diff(c(scanindex, length(mz)))

        cat("\n------------------------------------------------------------\n")
        cat("|    Impute: none\n\n")
        impute <- "none"
        fwhm <- 30
        sigma <- fwhm / 2.3548
        max <- 5
        snthresh <- 10
        step <- 0.1
        steps <- 2
        mzdiff <- 0.8 - step * steps
        index <- FALSE

        ## bin
        res <- .compare_em(xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step,
                           steps = steps, mzdiff = mzdiff,
                           index = index)    # OK
        steps <- 4
        mzdiff <- 0.8 - step * steps
        res <- .compare_em(xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step,
                           steps = steps, mzdiff = mzdiff,
                           index = index)    # OK
        steps <- 5
        mzdiff <- 0.8 - step * steps
        res <- .compare_em(xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step,
                           steps = steps, mzdiff = mzdiff,
                           index = index)    # OK
        ## THIS WILL FAIL!!! REASONE: rounding error in bin definition.
        step <- 0.2
        steps <- 2
        mzdiff <- 0.8 - step * steps
        res <- .compare_em(xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step,
                           steps = steps, mzdiff = mzdiff,
                           index = index,
                           stopOnError = FALSE)    # FAIL
        step <- 0.1999
        steps <- 2
        mzdiff <- 0.8 - step * steps
        res <- .compare_em(xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step,
                           steps = steps, mzdiff = mzdiff,
                           index = index)    # OK
        snthresh <- 4
        res <- .compare_em(xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step,
                           steps = steps, mzdiff = mzdiff,
                           index = index)    # OK
        fwhm <- 15
        sigma <- fwhm / 2.3548
        res <- .compare_em(xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index)    # OK

        ## ######
        ## binLin EXPECT DIFFERENCES HERE:
        ##  - imputation on the full matrix could be different to subsets.
        ##  - imputation of binLin is wrong (first and last bin).
        cat("\n------------------------------------------------------------\n")
        cat("|    Impute: lin\n\n")
        impute <- "lin"
        fwhm <- 30
        sigma <- fwhm / 2.3548
        max <- 5
        snthresh <- 10
        step <- 0.1
        steps <- 2
        mzdiff <- 0.8 - step * steps
        index <- FALSE

        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)
        ## Differences: B != orig, C != orig, C != orig, B != C, B != D, C != D
        ## But results are comparable.
        steps <- 4
        mzdiff <- 0.8 - step * steps
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)
        ## Differences: B != orig, C != orig, C != orig, B != C, B != D, C != D
        step <- 0.1
        mzdiff <- 0.8 - step * steps
        snthresh <- 100
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)
        ## Differences: B != orig, C != orig, C != orig, B != C, B != D
        fwhm <- 15
        sigma <- fwhm / 2.3548
        snthresh <- 40
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)
        step <- 10 ## Smaller step (0.01 etc) results in larger differences
        fwhm <- 30
        sigma <- fwhm / 2.3548
        snthresh <- 20
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)

        ## ######
        ## binLinBase
        cat("\n------------------------------------------------------------\n")
        cat("|    Impute: linbase\n\n")
        impute <- "linbase"
        fwhm <- 30
        sigma <- fwhm / 2.3548
        max <- 5
        snthresh <- 10
        step <- 0.1
        steps <- 2
        mzdiff <- 0.8 - step * steps
        index <- FALSE

        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = TRUE)  # OK; but there was no interpolation
        ## Changing baseValue.
        baseV <- 233
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           baseValue = baseV,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = TRUE)  ## OK; but also no interpolation
        ## Changing distance:
        distance <- 1L
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           baseValue = baseV, distance = distance,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)  ## A vs orig and C vs D OK
        ## and not passing baseVal
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           distance = distance,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)  ## A vs orig and C vs D OK
        distance <- 5L
        ## and not passing baseVal
        res <- .compare_em(xr = xr, mz = mz, int = int, scantime = scantime,
                           valsPerSpect = valsPerSpect, impute = impute,
                           distance = distance,
                           fwhm = fwhm, sigma = sigma, max = max,
                           snthresh = snthresh, binSize = step, steps = steps,
                           mzdiff = mzdiff, index = index,
                           stopOnError = FALSE)  ## A vs orig and C vs D OK
    }
}


.compare_em <- function(xr, mz, int, scantime, valsPerSpect, impute, baseValue,
                        distance, fwhm, sigma, max, snthresh, binSize, steps, mzdiff,
                        index, stopOnError = TRUE) {
    if (impute == "none") {
        profMethod(xr) <- "bin"
    }
    if (impute == "lin") {
        profMethod(xr) <- "binlin"
    }
    if (impute == "linbase") {
        profMethod(xr) <- "binlinbase"
    }
    if (impute == "intlin") {
        profMethod(xr) <- "intlin"
    }
    profparam <- list()
    if (!missing(baseValue))
        profparam$baselevel <- baseValue
    if (!missing(distance)) {
        mass <- seq(floor(min(mz) / binSize) * binSize,
                    ceiling(max(mz) / binSize) * binSize, by = binSize)
        bin_size <- (max(mass) - min(mass)) / (length(mass) - 1)
        profparam$basespace <- distance * bin_size
    }
    xr@profparam <- profparam
    ## The reference is the old code.
    ## Have to use the _orig method here, since the "official" one uses
    ## already do_findChromPeaks...
    orig <- xcms:::findPeaks.matchedFilter_orig(xr,
                                                fwhm = fwhm,
                                                sigma = sigma,
                                                max = max,
                                                snthresh = snthresh,
                                                step = binSize,
                                                steps = steps,
                                                mzdiff = mzdiff,
                                                index = index)@.Data

    ## A
    A <- xcms:::.matchedFilter_orig(mz = mz, int = int,
                                    scantime = scantime,
                                    valsPerSpect = valsPerSpect,
                                    binSize = binSize,
                                    impute = impute,
                                    baseValue,
                                    distance,
                                    fwhm = fwhm,
                                    sigma = sigma,
                                    max = max,
                                    snthresh = snthresh,
                                    steps = steps,
                                    mzdiff = mzdiff,
                                    index = index)
    res <- all.equal(orig, A)
    if (is.character(res)) {
        msg <- paste0("A vs original FAILED! ", res, "\n")
        cat(msg)
        .compare_peaks(A, orig)
        if (stopOnError)
            stop(msg)
    } else {
        cat("A vs original OK\n")
    }
    ## B
    B <- xcms:::.matchedFilter_binYonX_iter(mz = mz,
                                            int = int,
                                            scantime = scantime,
                                            valsPerSpect = valsPerSpect,
                                            binSize = binSize,
                                            impute = impute,
                                            baseValue,
                                            distance,
                                            fwhm = fwhm,
                                            sigma = sigma,
                                            max = max,
                                            snthresh = snthresh,
                                            steps = steps,
                                            mzdiff = mzdiff,
                                            index = index
                                            )
    res <- all.equal(orig, B)
    if (is.character(res)) {
        msg <- paste0("B vs original FAILED! ", res, "\n")
        cat(msg)
        .compare_peaks(B, orig)
        if (stopOnError)
            stop(msg)
    } else {
        cat("B vs original OK\n")
    }
    ## C
    C <- xcms:::.matchedFilter_no_iter(mz = mz, int = int,
                                       scantime = scantime,
                                       valsPerSpect = valsPerSpect,
                                       binSize = binSize,
                                       impute = impute,
                                       baseValue,
                                       distance,
                                       fwhm = fwhm,
                                       sigma = sigma,
                                       max = max,
                                       snthresh = snthresh,
                                       steps = steps,
                                       mzdiff = mzdiff,
                                       index = index
                                       )
    res <- all.equal(orig, C)
    if (is.character(res)) {
        msg <- paste0("C vs original FAILED! ", res, "\n")
        cat(msg)
        .compare_peaks(C, orig)
        if (stopOnError)
            stop(msg)
    } else {
        cat("C vs original OK\n")
    }
    ## D
    D <- xcms:::.matchedFilter_binYonX_no_iter(mz = mz,
                                               int = int,
                                               scantime = scantime,
                                               valsPerSpect = valsPerSpect,
                                               binSize = binSize,
                                               impute = impute,
                                               baseValue,
                                               distance,
                                               fwhm = fwhm,
                                               sigma = sigma,
                                               max = max,
                                               snthresh = snthresh,
                                               steps = steps,
                                               mzdiff = mzdiff,
                                               index = index
                                               )
    res <- all.equal(orig, D)
    if (is.character(res)) {
        msg <- paste0("D vs original FAILED! ", res, "\n")
        cat(msg)
        .compare_peaks(D, orig)
        if (stopOnError)
            stop(msg)
    } else {
        cat("D vs original OK\n")
    }
    ## Now compare between:
    res <- all.equal(B, C)
    if (is.character(res)) {
        msg <- paste0("B vs C FAILED! ", res, "\n")
        cat(msg)
        .compare_peaks(B, C)
    } else {
        cat("B vs C OK\n")
    }
    res <- all.equal(B, D)
    if (is.character(res)) {
        msg <- paste0("B vs D FAILED! ", res, "\n")
        cat(msg)
        .compare_peaks(B, D)
    } else {
        cat("B vs D OK\n")
    }
    res <- all.equal(C, D)
    if (is.character(res)) {
        msg <- paste0("C vs D FAILED! ", res, "\n")
        cat(msg)
        .compare_peaks(C, D)
    } else {
        cat("C vs D OK\n")
    }
    return(list(A = A, B = B, C = C, D = D, orig = orig))
}

.compare_peaks <- function(a, b, cols = c("into", "intf", "maxo", "maxf",
                                          "i", "sn")) {
    ## find peaks with same mz and rt:
    a_num <- nrow(a)
    b_num <- nrow(b)
    rownames(a) <- paste(a[, "mz"], a[, "rt"], sep = ":")
    rownames(b) <- paste(b[, "mz"], b[, "rt"], sep = ":")
    common <- intersect(rownames(a), rownames(b))
    cat("-----------------------------\n")
    cat("| Peaks: a: ", a_num, " b: ", b_num, " common: ",
        length(common), "\n", sep = "")
    ## Comparing peaks.
    if (length(common) > 0) {
        a <- a[common, , drop = FALSE]
        b <- b[common, , drop = FALSE]
        for (theCol in cols) {
            cat("| '", theCol, "' comparison: ", sep = "")
            if (is.character(all.equal(a[, theCol], b[, theCol]))) {
                same <- mapply(a[, theCol], b[, theCol], FUN = identical)
                cat(sum(same), " equal (", sum(same)/nrow(a)*100,"%)\n",
                    sep ="")
            } else{
                cat("OK\n")
            }
        }
    }
    cat("-----------------------------\n")
}


############################################################
## Lengthy explanation:
## Point is that the division taking place inside the ProfBinLinBase C-function
## to determine how many neighboring bins should be included in the linear
## interpolation can become unstable if step and basespace are set to the
## same value; this does seem to take place specifically if iterative
## buffer creation is used.
notrun_odd_matchedFilter_behaviour <- function() {

    library(xcms)
    library(RUnit)
    library(faahKO)
    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"))
    xr <- xcmsRaw(fs)
    profparam <- list()
    ## That's stable
    profparam$basespace <- 0.10000001
    step <- 0.1
    xr@profparam <- profparam
    profMethod(xr) <- "binlinbase"
    res_stable <- findPeaks.matchedFilter(xr, step = step)

    ## Do the same with 0.1
    profparam$basespace <- 0.1
    xr@profparam <- profparam
    res_unstable <- findPeaks.matchedFilter(xr, step = step)

    checkTrue(nrow(res_stable) != nrow(res_unstable))

    ## Checking the number of peaks we would get if we're switching on or off
    ## the interpolation.
    profparam$basespace <- 0.05
    xr@profparam <- profparam
    res_no_inter <- findPeaks.matchedFilter(xr, step = step)
    profparam$basespace <- 0.13
    xr@profparam <- profparam
    res_with_inter <- findPeaks.matchedFilter(xr, step = step)

    ## Check what we've got:
    checkEquals(res_with_inter, res_stable)
    checkTrue(nrow(res_no_inter) > nrow(res_with_inter))
}

############################################################
## Check what happens with the max parameter if we provide a file with
## presumably large number of peaks.
dontrun_check_matchedFilter_max_param <- function() {

    f <- "/Volumes/Ext64/data/2016/2016-06/PILOT_POS/140616_POOL_IntraP_PE_POS_2.mzML"
    xr <- xcmsRaw(f, step = 0)

    res1 <- xcms:::.matchedFilter_orig(mz = xr@env$mz,
                                       int = xr@env$intensity,
                                       scantime = xr@scantime,
                                       valsPerSpect = diff(c(xr@scanindex,
                                                             length(xr@env$mz))),
                                       binSize = 0.1,
                                       max = 5)
    res2 <- xcms:::.matchedFilter_binYonX_no_iter(mz = xr@env$mz,
                                                  int = xr@env$intensity,
                                                  scantime = xr@scantime,
                                                  valsPerSpect = diff(c(xr@scanindex,
                                                                        length(xr@env$mz))),
                                                  binSize = 0.1,
                                                  max = 5)
    checkEquals(res1, res2)
    ## max 20
    res1 <- xcms:::.matchedFilter_orig(mz = xr@env$mz,
                                       int = xr@env$intensity,
                                       scantime = xr@scantime,
                                       valsPerSpect = diff(c(xr@scanindex,
                                                             length(xr@env$mz))),
                                       binSize = 0.2,
                                       max = 20)
    res2 <- xcms:::.matchedFilter_binYonX_no_iter(mz = xr@env$mz,
                                                  int = xr@env$intensity,
                                                  scantime = xr@scantime,
                                                  valsPerSpect = diff(c(xr@scanindex,
                                                                        length(xr@env$mz))),
                                                  binSize = 0.2,
                                                  max = 20)
    checkEquals(res1, res2)
    ## Just a performance comparison:
    library(microbenchmark)
    microbenchmark(xcms:::.matchedFilter_orig(mz = xr@env$mz,
                                              int = xr@env$intensity,
                                              scantime = xr@scantime,
                                              valsPerSpect = diff(c(xr@scanindex,
                                                                    length(xr@env$mz))),
                                              binSize = 0.1,
                                              max = 5),
                   xcms:::.matchedFilter_binYonX_no_iter(mz = xr@env$mz,
                                                         int = xr@env$intensity,
                                                         scantime = xr@scantime,
                                                         valsPerSpect = diff(c(xr@scanindex,
                                                                               length(xr@env$mz))),
                                                         binSize = 0.1,
                                                         max = 5),
                   times = 10
                   )
    ## Result: original code: 18.8 sec, mine: 5.2 secs.
}
