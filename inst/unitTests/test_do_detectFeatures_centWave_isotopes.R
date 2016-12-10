## Tests related to the findPeaks.centWaveWithAddIsotopeROIs.
library(faahKO)
fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
        system.file('cdf/KO/ko18.CDF', package = "faahKO"),
        system.file('cdf/KO/ko19.CDF', package = "faahKO"))
xr <- xcmsRaw(fs[1], profstep = 0)
mzVals <- xr@env$mz
intVals <- xr@env$intensity
## Define the values per spectrum:
valsPerSpect <- diff(c(xr@scanindex, length(mzVals)))
res1 <- do_detectFeatures_centWave(mz = mzVals,
                                   int = intVals,
                                   scantime = xr@scantime,
                                   valsPerSpect,
                                   ##snthresh = 40,
                                   ##noise = 10000,
                                   verboseColumns = TRUE)

## Test during implementation of do_define_isotopes to check whether results
## fit with the ones from the original code.
dontrun_test_do_define_isotopes <- function() {
    ## Subset to columns we want, transform to data.frame.
    req_cols <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "intb", "scale")
    roiL <- split(res1, f = 1:nrow(res1))  ## is splitting a matrix faster?
    cns <- colnames(res1)
    roiL <- lapply(roiL, function(z) {
        names(z) <- cns
        return(z)
    })
    ## Original, returns a list which is appended to the original one.
    orig_ <- xcms:::do_define_isotopes_orig(roiList = roiL)
    new_ <- xcms:::do_define_isotopes(res1)
    checkIdentical(new_, orig_)
}

## Test during implementation of do_define_adducts to check whether results
## fit with the ones from the original code.
dontrun_test_do_define_adducts <- function() {
    ## Subset to columns we want, transform to data.frame.
    req_cols <- c("mz", "mzmin", "mzmax", "scmin", "scmax", "intb", "scale")
    roiL <- split(res1, f = 1:nrow(res1))  ## is splitting a matrix faster?
    cns <- colnames(res1)
    roiL <- lapply(roiL, function(z) {
        names(z) <- cns
        return(z)
    })
    ## Original, returns a list which is appended to the original one.
    orig_ <- xcms:::do_define_adducts_orig(roiList = roiL)
    new_ <- xcms:::do_define_adducts(res1)
    checkIdentical(new_, orig_)
}


dontrun_test_impl <- function() {

    ## Needed variables
    maxCharge <- 3
    maxIso <- 5
    mzIntervalExtension <- TRUE
    ppm <- 25

    mz <- mzVals
    int <- intVals
    scantime <- xr@scantime

    if (!is.double(mz))
        mz <- as.double(mz)
    if (!is.double(int))
        int <- as.double(int)

    feats_1 <- res1

    addNewIsotopeROIs <- TRUE
    addNewAdductROIs <- FALSE
    ## Extend the mzmin and mzmax if needed.
    tittle <- feats_1[, "mz"] * (ppm / 2) / 1E6
    expand_mz <- (feats_1[, "mzmax"] - feats_1[, "mzmin"]) < (tittle * 2)
    if (any(expand_mz)) {
        feats_1[expand_mz, "mzmin"] <- feats_1[expand_mz, "mz"] - tittle[expand_mz]
        feats_1[expand_mz, "mzmax"] <- feats_1[expand_mz, "mz"] + tittle[expand_mz]
    }
    ## Add predicted ROIs
    if (addNewIsotopeROIs) {
        iso_ROIs <- xcms:::do_define_isotopes(features. = feats_1,
                                       maxcharge = maxCharge,
                                       maxiso = maxIso,
                                       mzIntervalExtension = mzIntervalExtension)
    } else {
        iso_ROIs <- matrix(nrow = 0, ncol = 8)
        colnames(iso_ROIs) <- c("mz", "mzmin", "mzmax", "scmin", "scmax",
                                "length", "intensity", "scale")
    }
    if (addNewAdductROIs) {
        add_ROIs <- do_define_adducts(features. = feats_1, polarity = polarity)
    } else {
        add_ROIs <- matrix(nrow = 0, ncol = 8)
        colnames(iso_ROIs) <- c("mz", "mzmin", "mzmax", "scmin", "scmax",
                                "length", "intensity", "scale")
    }
    newROIs <- rbind(iso_ROIs, add_ROIs)
    ## OK, that matches the original code up to here.

    if (nrow(newROIs) == 0)
        return(feats_1)
    ## Remove ROIs that are out of mz range:
    mz_range <- range(mz)
    newROIs <- newROIs[newROIs[, "mzmin"] >= mz_range[1] &
                   newROIs[, "mzmax"] <= mz_range[2], , drop = FALSE]
    ## Remove ROIs with too low signal:
    keep_me <- logical(nrow(newROIs))
    scanindex <- as.integer(xcms:::valueCount2ScanIndex(valsPerSpect))
    for (i in 1:nrow(newROIs)) {
        vals <- .Call("getEIC", mz, int, scanindex,
                      as.double(newROIs[i, c("mzmin", "mzmax")]),
                      as.integer(newROIs[i, c("scmin", "scmax")]),
                      as.integer(length(scantime)), PACKAGE ='xcms' )
        keep_me[i] <- sum(vals$intensity, na.rm = TRUE) >= 10
    }
    newROIs <- newROIs[keep_me, , drop = FALSE]

    if (nrow(newROIs) == 0) {
        warning("No isotope or adduct ROIs for the identified features with a ",
                "valid signal found!")
        return(feats_1)
    }
    ## OK with do_predictIsotopeROIs

    ## centWave using the identified ROIs.
    roiL <- split(as.data.frame(newROIs), f = 1:nrow(newROIs))
    feats_2 <- do_detectFeatures_centWave(mz = mz, int = int,
                                          scantime = scantime,
                                          valsPerSpect = valsPerSpect, ppm = ppm,
                                          ## peakwidth = peakwidth,
                                          ## snthresh = snthresh,
                                          ## prefilter = prefilter,
                                          ## mzCenterFun = mzCenterFun,
                                          ## integrate = integrate,
                                          ## mzdiff = mzdiff, fitgauss = fitgauss,
                                          ## noise = noise,
                                          ## verboseColumns = verboseColumns,
                                          roiList = roiL,
                                          firstBaselineCheck = FALSE,
                                          roiScales = newROIs[, "scale"])
    ## Clean up of the results:
    if (nrow(feats_2) > 0) {
        ## remove NaNs
        any_na <- is.na(rowSums(feats_2[, c("mz", "mzmin", "mzmax", "rt", "rtmin",
                                            "rtmax")]))
        if (any(any_na))
            feats_2 <- feats_2[!any_na, , drop = FALSE]
        ## remove empty area
        no_area <- (feats_2[, "mzmax"] - feats_2[, "mzmin"]) == 0 ||
            (feats_2[, "rtmax"] - feats_2[, "rtmin"]) == 0
        if (any(no_area))
            feats_2 <- feats_2[!no_area, , drop = FALSE]
    }

    ## TO COMPARE.

    ## ## Compare that to the do_prodictIsotopeROIs!
    ## other_ <- xcms:::do_predictIsotopeROIs(xr, new("xcmsPeaks", res1))
    ## other_ <- lapply(other_, data.frame)
    ## other_ <- do.call(rbind, other_)
    ## rownames(other_) <- NULL
    ## checkEquals(as.data.frame(newROIs), other_)
}


