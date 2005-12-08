.First.lib <- function(libname, pkgname) {
    dllpath <- file.path(.find.package("xcms"), "netcdfdll", "netcdf.dll")
    suppressWarnings(dyn.load(dllpath))
    library.dynam("xcms", pkgname, libname)    
}

.Last.lib <- function(libpath) {
    rampCloseAll()
    if (is.loaded("NetCDFOpen"))
        library.dynam.unload("xcms", libpath)
    dllpath <- file.path(.find.package("xcms"), "netcdfdll", "netcdf.dll")
    if (is.loaded("nc_strerror"))
        dyn.unload(dllpath)
}
