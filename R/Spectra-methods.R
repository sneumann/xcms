#' @include Spectra.R

setMethod("show", "Spectra", function(object) {
    .show_Spectra(object, margin = "  ", print.classinfo = TRUE)
})

#' @rdname Spectra
#'
#' @section Accessing spectrum attributes:
#'
#' These methods allow to access the attributes and values of the individual
#' `Spectrum` ([Spectrum1-class] or [Spectrum2-class]) objects within the list.
#'
#' - `mz` return the m/z values of each spectrum as a `list` of `numeric`
#'   vectors.
#'
#' - `intensity` return the intensity values of each spectrum as a `list` of
#'   `numeric` vectors.
#'
#' - `rtime` return the retention time of each spectrum as a `numeric` vector
#'   with length equal to the length of `object`.
#'
#' - `precursorMz`, `precursorCharge`, `precursorIntensity`, `precScanNum`
#'   return precursor m/z values, charge, intensity and scan number for each
#'   spectrum as a `numeric` (or `integer`) vector with length equal to the
#'   length of `object`. Note that for [Spectrum1-class] objects `NA` will be
#'   returned.
#'
#' - `acquisitionNum` and `scanIndex` return the acquisition number of each
#'   spectrum and its scan index as an `integer` vector with the same length
#'   than `object`.
#'
#' - `ionCount` and `tic` return the ion count and total ion current of each
#'   spectrum.
#'
#' - `peaksCount` returns the number of peaks for each spectrum as a `integer`
#'   vector.
#'
#' - `msLevel` returns the MS level of each spectrum.
#'
#' - `collisionEnergy` returns the collision energy for each spectrum or `NA`
#'   for [Spectrum1-class] objects.
#'
#' - `polarity` returns the spectra's polarity.
#'
#' - `fromFile` returns the index from the (e.g. mzML) file the spectra where
#'   from. This applies only for spectra read using the [readMSData()] function.
#'
#' - `smoothed` whether spectra have been smoothed (i.e. processed with the
#'   [smooth()] method. Returns a `logical` of length equal to the
#'   number of spectra.
#'
#' - `isEmpty` returns `TRUE` for spectra without peak data.
#'
#' - `centroided`, `isCentroided` returns for each spectrum whether it contains
#'   *centroided* data. While `centroided` returns the internal attribute of
#'   each spectrum, `isCentroided` tries to guess whether spectra are
#'   centroided from the actual peak data.
#'
#' @md
#' 
#' @examples
#'
#' ## Extract the mz values for the individual spectra
#' mz(spl)
setMethod("mz", "Spectra", function(object) {
    lapply(object, function(z) z@mz)
})

#' @rdname Spectra
#'
#' @examples
#'
#' ## Extract the intensity values for the individual spectra
#' intensity(spl)
setMethod("intensity", "Spectra", function(object) {
    lapply(object, function(z) z@intensity)
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the retention time values for the individual spectra
#' rtime(spl)
setMethod("rtime", "Spectra", function(object) {
    vapply(object, function(z) if(length(z@rt)) z@rt else NA_real_, numeric(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the precursor m/z of each spectrum.
#' precursorMz(spl)
setMethod("precursorMz", "Spectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precursorMz)) z@precursorMz
        else NA_real_
    }, numeric(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the precursor charge of each spectrum.
#' precursorCharge(spl)
setMethod("precursorCharge", "Spectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precursorCharge)) z@precursorCharge
        else NA_integer_
    }, integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the precursor scan number for each spectrum.
#' precScanNum(spl)
setMethod("precScanNum", "Spectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precScanNum)) z@precScanNum
        else NA_integer_
    }, integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the precursor intensity of each spectrum.
#' precursorIntensity(spl)
setMethod("precursorIntensity", "Spectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precursorIntensity))
            z@precursorIntensity
        else NA_real_
    }, numeric(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the acquisition number of each spectrum.
#' acquisitionNum(spl)
setMethod("acquisitionNum", "Spectra", function(object) {
    vapply(object, function(z)
        if (length(z@acquisitionNum)) z@acquisitionNum else NA_integer_,
        integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the scan index of each spectrum.
#' scanIndex(spl)
setMethod("scanIndex", "Spectra", function(object) {
    vapply(object, function(z)
        if (length(z@scanIndex)) z@scanIndex else NA_integer_,
        integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Get the number of peaks per spectrum.
#' peaksCount(spl)
setMethod("peaksCount", "Spectra", function(object) {
    vapply(object, peaksCount, integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Get the MS level of each spectrum.
#' msLevel(spl)
setMethod("msLevel", "Spectra", function(object) {
    vapply(object, msLevel, integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Get the total ion current for each spectrum.
#' tic(spl)
setMethod("tic", "Spectra", function(object) {
    vapply(object, tic, numeric(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Get the total ion current for each spectrum.
#' ionCount(spl)
setMethod("ionCount", "Spectra", function(object) {
    vapply(object, ionCount, numeric(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the collision energy for spectrum.
#' collisionEnergy(spl)
setMethod("collisionEnergy", "Spectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@collisionEnergy)) z@collisionEnergy
        else NA_real_
    }, numeric(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Extract the file index for the spectrum.
#' fromFile(spl)
setMethod("fromFile", "Spectra", function(object) {
    vapply(object, function(z)
        if (length(z@fromFile)) z@fromFile else NA_integer_,
        integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Get the polarity for each spectrum.
#' polarity(spl)
setMethod("polarity", "Spectra", function(object) {
    vapply(object, function(z)
        if (length(z@polarity)) z@polarity else NA_integer_,
        integer(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Whether spectra are smoothed (i.e. processed with the `MSnbase::smooth`
#' ## function.
#' smoothed(spl)
setMethod("smoothed", "Spectra", function(object) {
    vapply(object, smoothed, logical(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Are spectra empty (i.e. contain no peak data)?
#' isEmpty(spl)
setMethod("isEmpty", "Spectra", function(x) {
    vapply(x, isEmpty, logical(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Do the spectra contain centroided data?
#' centroided(spl)
setMethod("centroided", "Spectra", function(object) {
    vapply(object, centroided, logical(1))
})

#' @rdname Spectra
#' 
#' @examples
#'
#' ## Do the spectra contain centroided data? Whether spectra are centroided
#' ## is estimated from the peak data.
#' isCentroided(spl)
setMethod("isCentroided", "Spectra", function(object) {
    vapply(object, isCentroided, logical(1))
})

## Still to implement:
## clean, all = FALSE
## removePeaks, t
## trimMz, mzlim
## filterMz, mz
## quantify, method = c("trapezoidation", "max", "sum", reporters, strict = FALSE)
## normalize, method = c("max", "sum", "precursor", precursorIntensity)
## bin, binSize = 1L, breaks....
## pickPeaks...
## smooth...


## Fun things:
## compareSpectra
