.First.lib <- function(libname, pkgname) {
    library.dynam("xcms", pkgname, libname)    
}

.Last.lib <- function(libpath) {
    library.dynam.unload("xcms", libpath)    
}
