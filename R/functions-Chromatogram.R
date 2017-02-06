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
    msg <- validMsg(NULL, NULL)
    if (length(object@rtime) != length(object@intensity))
        msg <- validMsg(msg, "Length of 'rt' and 'intensity' have to match!")
    if (is.unsorted(object@mzrange))
        msg <- validMsg(msg, "'mzrange' has to be increasingly ordered!")
    if (is.unsorted(object@rtime))
        msg <- validMsg(msg, paste0("'rtime' has to be increasingly ordered!"))
    if (length(object@mzrange) > 0 & length(object@mzrange) != 2)
        msg <- validMsg(msg, paste0("'mzrange' is supposed to contain the ",
                                    "minimum and maximum mz values for the ",
                                    "chromatogram."))
    if (length(object@fromFile) > 1 | any(object@fromFile < 0))
        msg <- validMsg(msg, paste0("'fromFile' is supposed to be a single ",
                                    "positive integer!"))
    if (length(object@aggregationFun) > 1)
        msg <- validMsg(msg, "Length of 'aggregationFun' has to be 1!")
    if (length(object@aggregationFun)) {
        if (!object@aggregationFun %in% .SUPPORTED_AGG_FUN_CHROM)
            msg <- validMsg(msg, paste0("Invalid value for 'aggregationFun'! ",
                                        "only ",
                                        paste0("'", .SUPPORTED_AGG_FUN_CHROM,"'",
                                              collapse = ","), " are allowed!"))
    }
    if (is.null(msg)) TRUE
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
##' @param mzrange \code{numeric(2)} representing the mz value range (min, max)
##' on which the chromatogram was created. If not applicable use
##' \code{mzrange = c(0, 0)}.
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
##' @slot rtime,intensity,mzrange,fromFile,aggregationFun See corresponding parameter above.
##' 
##' @rdname Chromatogram-class
Chromatogram <- function(rtime = numeric(), intensity = numeric(),
                         mzrange = c(0, 0), fromFile = integer(),
                         aggregationFun = character()) {
    return(new("Chromatogram", rtime = rtime, intensity = intensity,
               mzrange = range(mzrange), fromFile = as.integer(fromFile),
               aggregationFun = aggregationFun))
}
