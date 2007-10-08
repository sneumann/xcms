
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
  
  all.xcms <- ls(asNamespace(pkgname))
  find_methods <- function(type) {
    start <- nchar(type)
    substr(all.xcms[grep(paste("^", type, "\\..*", sep=""), all.xcms)], start+2, 100)
  }
  
  xcms.opt <- list(
    ## find methods for each type
    findPeaks.methods = find_methods("findPeaks"),
    removeBaseline.methods = find_methods("removeBaseline"),
    ## default for the methods
    findPeaks.method = "matchedFilter",
    removeBaseline.method = "medFilt"
  )
  
  class(xcms.opt) <- "BioCPkg"

  BioC <- getOption("BioC")
  BioC$xcms <- xcms.opt
  options("BioC"=BioC)
}
