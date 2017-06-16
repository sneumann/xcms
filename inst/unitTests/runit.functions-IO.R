test_mzRBackendFromFilename <- function() {
    checkEquals(xcms:::.mzRBackendFromFilename("test.mzML"), "pwiz")
    checkEquals(xcms:::.mzRBackendFromFilename("test.mzXML"), "pwiz")
    checkEquals(xcms:::.mzRBackendFromFilename("test.mzdata"), "Ramp")
    checkEquals(xcms:::.mzRBackendFromFilename("test.cdf"), "netCDF")
    checkEquals(xcms:::.mzRBackendFromFilename("test.nc"), "netCDF")
    checkEquals(xcms:::.mzRBackendFromFilename("test.bla"), NA)
    checkException(xcms:::.mzRBackendFromFilename(c(1, 2)))
}

test_mzRBackendFromFiletype <- function() {
    checkEquals(xcms:::.mzRBackendFromFiletype("mzml"), "pwiz")
    checkEquals(xcms:::.mzRBackendFromFiletype("MZXML"), "pwiz")
    checkEquals(xcms:::.mzRBackendFromFiletype("mzdata"), "Ramp")
    checkEquals(xcms:::.mzRBackendFromFiletype("netcdf"), "netCDF")
    checkEquals(xcms:::.mzRBackendFromFiletype("cdf"), "netCDF")
    checkException(xcms:::.mzRBackendFromFiletype(c(1, 2)))
    checkException(xcms:::.mzRBackendFromFiletype("bla"))
}

test_mzRBackendFromFilecontent <- function() {
    ## mzXML
    fl <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                      package = "msdata")
    checkEquals(xcms:::.mzRBackendFromFilecontent(fl), "pwiz")
    ## CDF
    fl <- system.file("cdf", "ko15.CDF", package = "msdata")
    checkEquals(xcms:::.mzRBackendFromFilecontent(fl), "netCDF")
    tmpf <- tempfile()
    file.copy(fl, tmpf)
    checkEquals(xcms:::.mzRBackendFromFilecontent(tmpf), "netCDF")
    ## mzML
    fl <- system.file("microtofq", "MM14.mzML", package = "msdata")
    checkEquals(xcms:::.mzRBackendFromFilecontent(fl), "pwiz")
    ## mzData
    fl <- system.file("iontrap", "extracted.mzData", package = "msdata")
    checkEquals(xcms:::.mzRBackendFromFilecontent(fl), "Ramp")
}
