
## .initFindPeaks <- function(all.xcms) {
##
##   assign("findPeaksMethods",
##         substr(all.xcms[grep("findPeaks.\.*", all.xcms)], start+1, 100),
##         envir=as.environment(where))
##
## }


.setXCMSOptions <- function(pkgname,xcms.opt=NA) {

    if (! any(is.na(xcms.opt))) {
        if (class(xcms.opt) != "BioCPkg")
            stop("obviously invalid package options !")

        BioC <- getOption("BioC")
        BioC$xcms <- xcms.opt
        options("BioC"=BioC)
        return()
    }

    ## add xcms specific options
    ## (not unlike what is done in 'affy')
    if (is.null(getOption("BioC"))) {
        BioC <- list()
        class(BioC) <- "BioCOptions"
        options("BioC"=BioC)
    }

    ## all findPeaks methods
    start <- nchar("findPeaks.")
    all.xcms <- ls(asNamespace(pkgname))
    findPeaks.methods <-  substr(all.xcms[grep("findPeaks\\..*", all.xcms)], start+1, 100)

    ## default for the methods
    findPeaks.method <- "matchedFilter"

    ## all groupPeaks methods
    start <- nchar("group.")
    all.xcms <- ls(asNamespace(pkgname))
    group.methods <-  substr(all.xcms[grep("group\\..*", all.xcms)], start+1, 100)

    ## default for the methods
    group.method <- "density"

    ## all groupPeaks methods
    start <- nchar("retcor.")
    all.xcms <- ls(asNamespace(pkgname))
    retcor.methods <-  substr(all.xcms[grep("retcor\\..*", all.xcms)], start+1, 100)

    ## default for the methods
    retcor.method <- "peakgroups"

    ## all fillPeaks methods
    start <- nchar("fillPeaks.")
    all.xcms <- ls(asNamespace(pkgname))
    fillPeaks.methods <-  substr(all.xcms[grep("fillPeaks\\..*", all.xcms)],
                                 start+1, 100)
    ## default method
    fillPeaks.method <- "chrom"

    ## all specDist methods
    start <- nchar("specDist.")
    all.xcms <- ls(asNamespace(pkgname))
    specDist.methods <-  substr(all.xcms[grep("specDist\\..*", all.xcms)],
                                start+1, 100)
    ## default method
    specDist.method <- "meanMZmatch"

    ## getEIC method
    getEIC.method="getEICOld"
    
    ## Sort method; see issue #180 for MSnbase
    ## sortMeth <- "auto"
    ## if (as.numeric(R.Version()$major) >= 3 & as.numeric(R.Version()$minor) >= 3)
    ##     sortMeth <- "radix"
    xcms.opt <- list(findPeaks.method = findPeaks.method,
                     findPeaks.methods = findPeaks.methods,
                     group.method = group.method,
                     group.methods = group.methods,
                     retcor.method = retcor.method,
                     retcor.methods = retcor.methods,
                     fillPeaks.method = fillPeaks.method,
                     fillPeaks.methods = fillPeaks.methods,
                     specDist.methods = specDist.methods,
                     getEIC.method = getEIC.method)
    ## No longer setting the useOriginalCode parameter as it might overwrite
    ## system wide setgings (e.g. if specified in .Rprofile

    class(xcms.opt) <- "BioCPkg"

    BioC <- getOption("BioC")
    BioC$xcms <- xcms.opt
    options("BioC"=BioC)
}
