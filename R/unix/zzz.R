.onLoad <- function(libname, pkgname) {
    # require(methods)
    .setXCMSOptions(pkgname)

    eval(expr= ".Last" <<- function() {
        if (is.loaded("mpi_initialize")){
            if (mpi.comm.size(1) > 0){
                mpi.close.Rslaves()
            }
            mpi.finalize()
        }
    }, envir = NULL)
}

.onUnload <- function(libpath) {
    mzR:::rampCloseAll()
    library.dynam.unload("xcms", libpath)
}
