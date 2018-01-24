#' @title Identify peaks close to calibrants' mz values
#'
#' @description Define which peaks are closest to the provided calibrants' mz
#'     values and return their index, mz value and difference between their and
#'     the calibrant's mz value.
#'
#' @note The mz values of the calibrants are first sorted. *Cave* the function
#'     requires the `peaklist` to be sorted by mz value (issue #200).
#'
#' @param peaklist `matrix` with the identified peaks of a single sample.
#'     Required columns are: `"mz"` and `intensity`
#' 
#' @param masslist `numeric` with the masses of the calibrants.
#' 
#' @param mzabs `numeric(1)` defining the absolute mz deviation.
#' 
#' @param mzppm `numeric(1)` defining the mz deviation in ppm.
#' 
#' @param neighbours `numeric(1)` with the number of neighbors.
#' 
#' @param intensity `character(1)` specifying the column in `peaklist` that
#'     contains the peaks' intensity values.
#'
#' @return A `matrix` with 3 columns:
#'     + `"pos"`: the position of the peak closest to the corresponding
#'        calibrant's mz value.
#'     + `"mass"`: the m/z of the peak.
#'     + `"dif"`: the difference between the calibrant's and peak's mz.
#'     The number of rows corresponds to the number of calibrants, unless no
#'     peak was found to be close enough to a calibrant. In the latter case
#'     the calibrant is removed from the result table.
#' 
#' @md
#' 
#' @noRd
matchpeaks <- function(peaklist, masslist, mzabs = 0.0001, mzppm = 5,
                       neighbours = 3, intensity = "into") {
    ## Dangerous, we're nowhere checking for masslist being correct.
    masslist <- masslist[order(masslist)]
    dif <- matrix(nrow = length(masslist), ncol = neighbours,
                  data = rep(1000, length(masslist) * neighbours))
    pos <- matrix(nrow = length(masslist), ncol = neighbours,
                  data = rep(0, length(masslist) * neighbours))
    mzu <- peaklist[, "mz"]
    itu <- peaklist[, intensity]
    mdif <- NA  ## smallest difference between a peak and a calibrant.
    mpos <- NA  ## position of the entry in the peak list.
    cdifs <- NA
    cposs <- NA
    lastval <- 1
    ## Loop over masslist.
    for (b in 1:length(masslist)){
        ## finding for each entry in masslist the ##neighbours closest masses
        ## in the peaklist
        a <- lastval
        sf <- FALSE ## is set to true if the neighbour-box is filled the first
        ## time for a entry in masslist
        while (a < length(mzu)){
            ## loop over peaks
            ## the currently biggest distance:
            mxd <- dif[b, max(which(abs(dif[b, ]) == max(abs(dif[b,]))))]
            ## the position of the distance
            mxp <- max(which(abs(dif[b, ]) == max(abs(dif[b,]))))
            ## If the difference between the peak and the current mass is
            ## smaller than the largest difference in dif
            if (abs(mzu[a] - masslist[b]) <= abs(mxd)) {
                dif[b, mxp] <- (mzu[a] - masslist[b])
                pos[b, mxp] <- a
                lastval <- min(pos[b, ])
                sf <- TRUE
            }else{ ## no more masses are smaller, switching to end of mzu
                if (sf) a <- length(mzu)
            }
            a <- a + 1
        }
        ## cat(min(abs(dif[b,])),"\n")
        ## checking the treshold and finding the candidate with the biggest
        ## intensity
        smalldiffs <- which(abs(dif[b, ]) <=
                            (mzabs + (mzu[pos[b, ]] / 1000000 * mzppm)))
        cdifs <- dif[b, smalldiffs]
        cposs <- pos[b, smalldiffs]
        if (length(cdifs) > 0){
            mcdi <- smalldiffs[which(itu[cposs] == max(itu[cposs]))]
            mdif[b] <- dif[b, mcdi]  ## mz - mz_masslist
            mpos[b] <- pos[b, mcdi]  ## The index of the mz in the peaklist.
        }else{
            ## No peak close.
            mdif[b] <- NA
            mpos[b] <- NA
        }
    }
    mdiffs <- mdif[!is.na(mdif)]
    mposs <- mpos[!is.na(mdif)]
    retdata <- matrix(ncol = 3, nrow = length(mdiffs),
                      data = c(mposs, mzu[mposs], mdiffs))
    colnames(retdata) = c("pos", "mass", "dif")
    retdata
}

#' @title Identify peaks close to calibrants' mz values
#'
#' @description Define which peaks are closest to the provided calibrants' mz
#'     values and return their index, mz value and difference between their and
#'     the calibrant's mz value.
#'
#' @details This function is faster than the original `matchpeaks` function and
#'     in addition does by default return the results for all calibrants' mz
#'     values, i.e. `nrow` of the result `matrix` is always equal to
#'     `length(masslist)`. To mimic the old function use `na.rm = TRUE`.
#' 
#' @note The mz values of the calibrants are first sorted.
#'
#' @param peaklist `matrix` with the identified peaks of a single sample.
#'     Required columns are: `"mz"` and `intensity`
#' 
#' @param masslist `numeric` with the masses of the calibrants.
#' 
#' @param mzabs `numeric(1)` defining the absolute mz deviation.
#' 
#' @param mzppm `numeric(1)` defining the mz deviation in ppm.
#' 
#' @param neighbours `numeric(1)` with the number of neighbors.
#' 
#' @param intensity `character(1)` specifying the column in `peaklist` that
#'     contains the peaks' intensity values.
#'
#' @param na.rm `logical(1)` whether results for calibrants for which no peak
#'     close enough were found should be excluded from the result `matrix`.
#' @return A `matrix` with 3 columns:
#'     + `"pos"`: the position of the peak closest to the corresponding
#'        calibrant's mz value.
#'     + `"mass"`: the m/z of the peak.
#'     + `"dif"`: the difference between the calibrant's and peak's mz.
#'     The number of rows matches by default (i.e. with `na.rm = FALSE`) the
#'     `length(masslist)`. The `matrix` contains `NA` values for calibrants for
#'     which no peak close to the calibrant's mz was found.
#' 
#' @md
#'
#' @author Johannes Rainer
#' 
#' @noRd
.matchpeaks2 <- function(peaklist, masslist, mzabs = 0.0001, mzppm = 5,
                         neighbours = 3, intensity = "into", na.rm = FALSE) {
    masslist <- sort(masslist)
    res_pos <- res_dif <- rep(NA_real_, length(masslist))
    peak_mz <- peaklist[, "mz"]
    peak_int <- peaklist[, "into"]
    mzppm <- mzppm / 1e6
    peak_mz_maxdiff <- mzabs + (peak_mz * mzppm) 
    for (i in seq_along(masslist)) {
        mzdiffs <- peak_mz - masslist[i]
        idx_mzdiffs <- which(abs(mzdiffs) <= peak_mz_maxdiff)
        if (length(idx_mzdiffs)) {
            the_idx <- idx_mzdiffs
            if (length(idx_mzdiffs) > 1) {
                ## The neighbours closest.
                idx_mzdiffs <- idx_mzdiffs[order(abs(mzdiffs)[idx_mzdiffs])]
                idx_mzdiffs <- idx_mzdiffs[1:min(length(idx_mzdiffs),
                                                 neighbours)]
                ## Select the one with the largest intensity.
                the_idx <- idx_mzdiffs[order(peak_int[idx_mzdiffs],
                                             decreasing = TRUE)][1]
            }
            res_pos[i] <- the_idx
            res_dif[i] <- mzdiffs[the_idx]
        }
    }
    not_na <- !is.na(res_pos)
    if (!any(not_na))
        stop("Not a single peak close to any of the specified calibrants")
    if (na.rm) {
        res <- cbind(pos = res_pos[not_na],
                     mass = peak_mz[res_pos[not_na]],
                     dif = res_dif[not_na])
    } else {
        res <- cbind(pos = res_pos,
                     mass = peak_mz[res_pos],
                     dif = res_dif)
    }
    res
}

#' @description Estimates the correction *model*.
#'
#' @param dtable `matrix` as returned by `matchpeaks` or `.matchpeaks2`.
#'
#' @param method `character(1)` defining the method that will be used to
#'     perform the calibration.
#'
#' @return A `numeric(2)` with slope and intercept of the linear model. If
#'     `method = "shift"` the slope will be `0`. Note that the first element
#'     is the slope and the second the intercept!
#'
#' @md
#'
#' @noRd
estimate <- function(dtable, method = c("linear")) {
    mdiffs <- dtable[, "dif"]
    mposs <- dtable[, "mass"]

    if (method != "shift") {
        ## complete regression
        regg <- lm(mdiffs ~ mposs)
        a <- regg$coefficients[2]
        b <- regg$coefficients[1]
    } else {
        ## only global shift
        a <- 0
        b <- mean(mdiffs)
    }
    if (method != "shift")
        c(a, b)
    else
        c(0, b)
}

#' @description Simple function to calibrate mz values using the provided
#'     parameters.
#'
#' @note `mz` values have to be increasingly sorted. At some point we might
#'     want to have a more generic implementation that supports, e.g. also
#'     loess fits.
#'
#' @param mz `numeric` with the mz values to be calibrated.
#'
#' @param method `character` defining the adjustment method. Can be either
#'     `"linear"`, `"shift"` or `"edgeshift"`.
#'
#' @param minMz `numeric(1)` with the minimal mz value. Required for
#'     `method = "edgeshift"`.
#'
#' @param maxMz `numeric(1)` with the maximal mz value. Required for
#'     `method = "edgeshift"`.
#'
#' @param slope `numeric(1)` with the slope.
#'
#' @param intercept `numeric(1)` with the intercept.
#'
#' @md
#'
#' @author Johannes Rainer
#' 
#' @noRd
.calibrate_mz <- function(mz, method = "linear", minMz, maxMz, slope = 0,
                          intercept) {
    if (is.unsorted(mz))
        stop("'mz' values have to be sorted")
    if (missing(intercept))
        stop("'intercept' is a required parameter")
    method <- match.arg(method, c("shift", "linear", "edgeshift"))
    if (method == "edgeshift") {
        if (missing(minMz) | missing(maxMz))
            stop("Parameters 'minMz' and 'maxMz' are required for method ",
                 "'edgeshift'")
        shift_lower <- mz < minMz
        shift_upper <- mz > maxMz
        mz_ <- mz
        if (any(shift_lower))
            mz_[shift_lower] <- mz[(max(which(shift_lower)) + 1)]
        if (any(shift_upper))
            mz_[shift_upper] <- mz[(min(which(shift_upper)) - 1)]
        mz <- mz - (slope * mz_ + intercept)
    } else {
        mz <- mz - (slope * mz + intercept)
    }
    mz
}
