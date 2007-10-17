
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
  
  protocols <- list(findPeaks = "matchedFilter", 
                    filterProfile = "median",
                    group = "density", 
                    retcor = "smooth", 
                    fillPeaks = "extract")
  
  xcms.opt <- list()
  for(type in names(protocols)) {
    xcms.opt[[paste(type, "methods",sep=".")]] <- find_methods(type)
    xcms.opt[[paste(type, "method", sep=".")]] <- protocols[[type]]
  }
  class(xcms.opt) <- "BioCPkg"

  BioC <- getOption("BioC")
  BioC$xcms <- xcms.opt
  options("BioC"=BioC)
}
