
.onLoad <- function(libname, pkgname) {
    require(methods)
    dllpath <- file.path(.find.package("xcms"), "netcdfdll", "netcdf.dll")
    suppressWarnings(dyn.load(dllpath))
    .setXCMSOptions(pkgname)
}

.onUnload <- function(libpath) {
    rampCloseAll()
    if (is.loaded("NetCDFOpen"))
        library.dynam.unload("xcms", libpath)
    dllpath <- file.path(.find.package("xcms"), "netcdfdll", "netcdf.dll")
    if (is.loaded("nc_strerror"))
        dyn.unload(dllpath)
}
