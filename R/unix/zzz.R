.First.lib <- function(libname, pkgname) {
    library.dynam("xcms", pkgname, libname)
}

.Last.lib <- function(libpath) {
    rampCloseAll()
    library.dynam.unload("xcms", libpath)
}
