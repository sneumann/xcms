netCDFStrError <- function(ncerr) {

    buflen <- 255

    .C("NetCDFStrError",
       as.integer(ncerr),
       as.integer(buflen),
       out = paste(rep(" ", buflen), collapse = ""),
       PACKAGE = "xcms")$out
}

netCDFIsFile <- function(filename) {

    ncid <- netCDFOpen(filename)
    if (!is.null(attr(ncid, "errortext")))
        return(FALSE)
    netCDFClose(ncid)

    return(TRUE)
}

netCDFOpen <- function(filename) {

    result <- .C("NetCDFOpen",
                 as.character(filename),
                 ncid = integer(1),
                 status = integer(1),
                 PACKAGE = "xcms")

    if (result$status)
        return(structure(result$status,
                         errortext = netCDFStrError(result$status)))

    return(result$ncid)
}

netCDFClose <- function(ncid) {

    result <- .C("NetCDFClose",
                 as.integer(ncid),
                 status = integer(1),
                 PACKAGE = "xcms")

    if (result$status)
        return(structure(result$status,
                         errortext = netCDFStrError(result$status)))

    result$status
}

netCDFVarID <- function(ncid, var) {

    result <- .C("NetCDFVarID",
                 as.integer(ncid),
                 as.character(var),
                 id = integer(1),
                 status = integer(1),
                 PACKAGE = "xcms")

    if (result$status)
        return(structure(result$status,
                         errortext = netCDFStrError(result$status)))

    return(result$id)
}

netCDFVarLen <- function(ncid, var) {

    if (is.character(var))
        var <- netCDFVarID(ncid, var)

    result <- .C("NetCDFVarLen",
                 as.integer(ncid),
                 as.integer(var),
                 len = integer(1),
                 status = integer(1),
                 PACKAGE = "xcms")

    if (result$status)
        return(structure(result$status,
                         errortext = netCDFStrError(result$status)))

    return(result$len)
}

netCDFVarDouble <- function(ncid, var) {

    if (is.character(var))
        var <- netCDFVarID(ncid, var)

    if (!is.null(attr(var, "errortext")))
        return(var)

    len <- netCDFVarLen(ncid, var)
    if (!is.null(attr(len, "errortext")))
        return(len)

    .C("NetCDFVarDouble",
       as.integer(ncid),
       as.integer(var),
       data = double(len),
       status = integer(1),
       PACKAGE = "xcms")$data
}

netCDFVarInt <- function(ncid, var) {

    if (is.character(var))
        var <- netCDFVarID(ncid, var)

    if (!is.null(attr(var, "errortext")))
        return(var)

    len <- netCDFVarLen(ncid, var)
    if (!is.null(attr(len, "errortext")))
        return(len)

    .C("NetCDFVarInt",
       as.integer(ncid),
       as.integer(var),
       data = integer(len),
       status = integer(1),
       PACKAGE = "xcms")$data
}

netCDFMSPoints <- function(ncid, scanIndex) {

    if (!is.integer(scanIndex)) scanIndex <- as.integer(scanIndex)

    var <- netCDFVarID(ncid, "mass_values")
    if (!is.null(attr(var, "errortext")))
        return(var)

    len <- netCDFVarLen(ncid, var)
    if (!is.null(attr(len, "errortext")))
        return(len)

    .C("NetCDFMSPoints",
       as.integer(ncid),
       as.integer(length(scanIndex)),
       scanIndex,
       as.integer(len),
       massValues = double(len),
       intensityValues = double(len),
       status = integer(1),
       PACKAGE = "xcms")[c("massValues", "intensityValues")]
}

netCDFRawData <- function(ncid) {

    rt <- netCDFVarDouble(ncid, "scan_acquisition_time")
    if (!is.null(attr(rt, "errortext")))
        stop("Could not read scan times")

    tic <- netCDFVarDouble(ncid, "total_intensity")
    if (!is.null(attr(tic, "errortext")))
        stop("Could not read total ion current")

    scanindex <- netCDFVarInt(ncid, "scan_index")
    if (!is.null(attr(scanindex, "errortext")))
        stop("Could not read scan indecies")

    pointValues <- netCDFMSPoints(ncid, scanindex)
    if (!is.null(attr(pointValues, "errortext")))
        stop("Could not read mass/intensity values")

    return(list(rt = rt, tic = tic, scanindex = scanindex,
                mz = pointValues$massValues,
                intensity = pointValues$intensityValues))
}
