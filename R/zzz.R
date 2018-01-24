.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        paste("\nThis is xcms version", packageVersion("xcms"), "\n"))
}

.onLoad <- function(libname, pkgname) {
    # require(methods)
    .setXCMSOptions(pkgname)

    ## That below should not be really required anymore.
    ## eval(expr= ".Last" <<- function() {
    ##     if (is.loaded("mpi_initialize")){
    ##         if (mpi.comm.size(1) > 0){
    ##             mpi.close.Rslaves()
    ##         }
    ##         mpi.finalize()
    ##     }
    ## }, envir = NULL)
}

.onUnload <- function(libpath) {
    ## mzR:::rampCloseAll()  # Makes only sense if we used ramp.
    library.dynam.unload("xcms", libpath)
}
