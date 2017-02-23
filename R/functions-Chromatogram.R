#' @include DataClasses.R
.SUPPORTED_AGG_FUN_CHROM <- c("sum", "max", "min", "mean")
names(.SUPPORTED_AGG_FUN_CHROM) <-
    c("Total ion chromatogram (TIC).", "Base peak chromatogram (BPC).",
      "Intensity representing the minimum intensity across the mz range.",
      "Intensity representing the mean intensity across the mz range.")

##' @title Validation function for Chromatogram objects
##'
##' @description This function can be used instead of the \code{validObject} to
##' check if the chromatogram is valid, without having to call the validity
##' method on all super classes.
##'
##' @param object A \code{Chromatogram} object.
##' 
##' @return \code{TRUE} if the \code{object} is valid and the error messages
##' otherwise (i.e. a \code{character}).
##' @author Johannes Rainer
##' @noRd
validChromatogram <- function(object) {
    msg <- character()
    if (length(object@rtime) != length(object@intensity))
        msg <- c(msg, "Length of 'rt' and 'intensity' have to match!")
    if (is.unsorted(object@mz))
        msg <- c(msg, "'mz' has to be increasingly ordered!")
    if (is.unsorted(object@rtime))
        msg <- c(msg, paste0("'rtime' has to be increasingly ordered!"))
    if (length(object@mz) > 0 & length(object@mz) != 2)
        msg <- c(msg, paste0("'mz' is supposed to contain the ",
                             "minimum and maximum mz values for the ",
                             "chromatogram."))
    if (length(object@filterMz) > 0 & length(object@filterMz) != 2)
        msg <- c(msg, paste0("'filterMz' is supposed to contain the ",
                             "minimum and maximum mz values of the filter",
                             " used to create the chromatogram."))
    if (length(object@parentMz) > 0 & length(object@parentMz) != 2)
        msg <- c(msg, paste0("'parentMz' is supposed to be a numeric of",
                             " length 2."))
    if (length(object@productMz) > 0 & length(object@productMz) != 2)
        msg <- c(msg, paste0("'productMz' is supposed to be a numeric of",
                             " length 2."))
    if (length(object@fromFile) > 1 | any(object@fromFile < 0))
        msg <- c(msg, paste0("'fromFile' is supposed to be a single ",
                             "positive integer!"))
    if (length(object@aggregationFun) > 1)
        msg <- c(msg, "Length of 'aggregationFun' has to be 1!")
    if (length(object@aggregationFun)) {
        if (!object@aggregationFun %in% .SUPPORTED_AGG_FUN_CHROM)
            msg <- c(msg, paste0("Invalid value for 'aggregationFun'! ",
                                 "only ",
                                 paste0("'", .SUPPORTED_AGG_FUN_CHROM,"'",
                                        collapse = ","), " are allowed!"))
    }
    if (length(msg) == 0) TRUE
    else msg
}

##' @description \code{Chromatogram}: create an instance of the
##' \code{Chromatogram} class.
##'
##' @param rtime \code{numeric} with the retention times (length has to be equal
##' to the length of \code{intensity}).
##'
##' @param intensity \code{numeric} with the intensity values (length has to be
##' equal to the length of \code{rtime}).
##'
##' @param mz \code{numeric(2)} representing the mz value range (min, max)
##' on which the chromatogram was created. This is supposed to contain the
##' \emph{real} range of mz values in contrast to the \code{filterMz} below.
##' If not applicable use \code{mzrange = c(0, 0)}.
##'
##' @param filterMz \code{numeric(2)} representing the mz value range (min,
##' max) that was used to filter the original object on mz dimension. If not
##' applicable use \code{filterMz = c(0, 0)}.
##'
##' @param parentMz \code{numeric(2)} for MRM transition tracking. Represents the
##' mz of the parent ion. See details for more information.
##' 
##' @param productMz \code{numeric(2)} for MRM transition tracking. Represents
##' the mz of the product. See details for more information.
##' 
##' @param fromFile \code{integer(1)} the index of the file within the
##' \code{\link{OnDiskMSnExp}} or \code{\link{XCMSnExp}} from which the
##' chromatogram was extracted.
##'
##' @param aggregationFun \code{character} string specifying the function that
##' was used to aggregate intensity values for the same retention time across the
##' mz range. Supported are \code{"sum"} (total ion chromatogram), \code{"max"}
##' (base peak chromatogram), \code{"min"} and \code{"mean"}.
##' 
##' @slot rtime,intensity,mz,filterMz,fromFile,aggregationFun See corresponding parameter above.
##' 
##' @rdname Chromatogram-class
Chromatogram <- function(rtime = numeric(), intensity = numeric(),
                         mz = c(0, 0), filterMz = c(0, 0),
                         parentMz = c(NA_real_, NA_real_),
                         productMz = c(NA_real_, NA_real_),
                         fromFile = integer(),
                         aggregationFun = character()) {
    return(new("Chromatogram", rtime = rtime, intensity = intensity,
               mz = range(mz), filterMz = range(filterMz),
               parentMz = range(parentMz), productMz = range(productMz),
               fromFile = as.integer(fromFile), aggregationFun = aggregationFun))
}
