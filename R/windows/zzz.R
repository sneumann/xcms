.First.lib <- function(libname, pkgname) {
    suppressWarnings(library.dynam("netcdf", pkgname, libname))
    library.dynam("xcms", pkgname, libname)    
}

.Last.lib <- function(libpath) {
    library.dynam.unload("xcms", libpath)
    library.dynam.unload("netcdf", libpath)
}
