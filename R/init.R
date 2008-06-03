
# .initFindPeaks <- function(all.xcms) {
#   
# #   assign("findPeaksMethods",
# #         substr(all.xcms[grep("findPeaks.\.*", all.xcms)], start+1, 100),
# #         envir=as.environment(where))
# #         
# }


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
  
  xcms.opt <- list(findPeaks.method=findPeaks.method, findPeaks.methods=findPeaks.methods) 

  class(xcms.opt) <- "BioCPkg"

  BioC <- getOption("BioC")
  BioC$xcms <- xcms.opt
  options("BioC"=BioC)
}
