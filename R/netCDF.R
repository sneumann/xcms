netCDFStrError <- function(ncerr) {

    buflen <- 255
    
    .C("NetCDFStrError",
       as.integer(ncerr),
       as.integer(buflen),
       out = paste(rep(" ", buflen), collapse = ""),
       PACKAGE = "xcms")$out
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
       DUP = FALSE, PACKAGE = "xcms")$data
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
       DUP = FALSE, PACKAGE = "xcms")$data
}
