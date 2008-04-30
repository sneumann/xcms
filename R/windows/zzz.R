
.onLoad <- function(libname, pkgname) {
    require(methods)
    .setXCMSOptions(pkgname)
}

.onUnload <- function(libpath) {
    rampCloseAll()
    if (is.loaded("NetCDFOpen"))
        library.dynam.unload("xcms", libpath)
}
