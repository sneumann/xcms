.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("\nThis is xcms version", packageVersion("xcms"), "\n"))
}

.onLoad <- function(libname, pkgname) {
    # require(methods)
    .setXCMSOptions(pkgname)
}

.onUnload <- function(libpath) {
    ## mzR:::rampCloseAll()  # Makes only sense if we used ramp.
    library.dynam.unload("xcms", libpath)
}
