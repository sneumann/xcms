## All methods for xcmsSet should go here.

#' @description This method updates an \emph{old} \code{xcmsSet} object to the latest
#' definition.
#' @title Update an \code{xcmsSet} object
#' @param object The \code{xcmsSet} object to update.
#' @param ... Optional additional arguments. Currently ignored.
#' @param verbose Currently ignored.
#' @return An updated \code{xcmsSet} containing all data from the input object.
#' @author Johannes Rainer
setMethod("updateObject", "xcmsSet", function(object, ..., verbose = FALSE) {
    ## Create a new empty xcmsSet and start filling it with the slot
    ## values (if present.)
    ## Define the default values.
    pks <- matrix(nrow = 0, ncol = 0)
    grps <- matrix(nrow = 0, ncol = 0)
    grpidx <- list()
    flld <- integer(0)
    phnD <- data.frame()
    theRt <- list()
    flpth <- character(0)
    prfnf <- vector("list")
    datCorr <- integer(0)
    pol <- character(0)
    prgInfo <- list()
    msl <- numeric(0)
    scnr <- numeric(0)
    prgCb <- function(progress) NULL

    ## Now replace the values with the slots... if present.
    if (.hasSlot(object, "peaks"))
        pks <- object@peaks
    if (.hasSlot(object, "groups"))
        grps <- object@groups
    if (.hasSlot(object, "groupidx"))
        grpidx <- object@groupidx
    if (.hasSlot(object, "filles"))
        flld <- object@filles
    if (.hasSlot(object, "phenoData"))
        phnD <- object@phenoData
    if (.hasSlot(object, "rt"))
        theRt <- object@rt
    if (.hasSlot(object, "filepaths"))
        flpth <- object@filepaths
    if (.hasSlot(object, "profinfo"))
        prfnf <- object@profinfo
    if (.hasSlot(object, "dataCorrection"))
        datCorr <- object@dataCorrection
    if (.hasSlot(object, "polarity"))
        pol <- object@polarity
    if (.hasSlot(object, "polarity"))
        prgInfo <- object@progressInfo
    if (.hasSlot(object, "mslevel"))
        msl <- object@mslevel
    if (.hasSlot(object, "scanrange"))
        scnr <- object@scanrange
    if (.hasSlot(object, "progressCallback"))
        prgCb <- object@progressCallback

    ## Generate the new object.
    newXs <- new("xcmsSet",
                 peaks = pks,
                 groups = grps,
                 groupidx = grpidx,
                 filled = flld,
                 phenoData = phnD,
                 rt = theRt,
                 filepaths = flpth,
                 profinfo = prfnf,
                 dataCorrection = datCorr,
                 polarity = pol,
                 progressInfo = prgInfo,
                 mslevel = msl,
                 scanrange = scnr,
                 progressCallback = prgCb
                 )
    return(newXs)
})
