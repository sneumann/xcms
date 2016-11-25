## Functions related to the Param class and sub-classes.
#' @include DataClasses.R c.R

## Extract all slot values and put them into a list, names being the slot
## names.
.param2list <- function(x) {
    ## Get all slot names, skip those matching the provided pattern.
    sNames <- slotNames(x)
    skipSome <- grep(sNames, pattern = "^\\.")
    if(length(skipSome) > 0)
        sNames <- sNames[-skipSome]
    if(length(sNames) > 0){
        resL <- vector("list", length(sNames))
        for(i in 1:length(sNames))
            resL[[i]] <- slot(x, name = sNames[i])
        names(resL) <- sNames
        return(resL)
    }else{
          return(list())
    }
}

############################################################
## CentWaveParam

##' @inheritParams do_detectFeatures_centWave
##'
##' @return The \code{CentWaveParam} function returns a \code{CentWaveParam}
##' class instance with all of the settings specified for feature detection by
##' the centWave method.
##'
##' @rdname featureDetection-centWave
CentWaveParam <- function(ppm = 25, peakwidth = c(20, 50), snthresh = 10,
                          prefilter = c(3, 100), mzCenterFun = "wMean",
                          integrate = 1L, mzdiff = -0.001, fitgauss = FALSE,
                          noise = 0, verboseColumns = FALSE, roiList = list(),
                          firstBaselineCheck = TRUE, roiScales = numeric()) {
    return(new("CentWaveParam", ppm = ppm, peakwidth = peakwidth,
               snthresh = snthresh, prefilter = prefilter,
               mzCenterFun = mzCenterFun, integrate = integrate,
               mzdiff = mzdiff, fitgauss = fitgauss, noise = noise,
               verboseColumns = verboseColumns, roiList = roiList,
               firstBaselineCheck = firstBaselineCheck, roiScales = roiScales))
}


############################################################
## MatchedFilterParam

##' @inheritParams do_detectFeatures_matchedFilter
##'
##' @return The \code{MatchedFilterParam} function returns a
##' \code{MatchedFilterParam} class instance with all of the settings specified
##' for feature detection by the centWave method.
##'
##' @rdname featureDetection-matchedFilter
MatchedFilterParam <- function(binSize = 0.1, impute = "none",
                               baseValue = numeric(), distance = numeric(),
                               fwhm = 30, sigma = fwhm / 2.3548,
                               max = 5, snthresh = 10, steps = 2,
                               mzdiff = 0.8 - binSize * steps, index = FALSE) {
    return(new("MatchedFilterParam", binSize = binSize, impute = impute,
               baseValue = baseValue, distance = distance, fwhm = fwhm,
               sigma = sigma, max = max, snthresh = snthresh, steps = steps,
               mzdiff = mzdiff, index = index))
}
