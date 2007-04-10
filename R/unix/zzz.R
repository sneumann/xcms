
.onLoad <- function(libname, pkgname) {
    require(methods)
    .setXCMSOptions(pkgname)
}

.onUnload <- function(libpath) {
    rampCloseAll()
    library.dynam.unload("xcms", libpath)
}
