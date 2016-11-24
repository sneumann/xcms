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

##' @description \code{CentWaveParam}: constructor function for
##' \code{CentWaveParam} classes.
##'
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
