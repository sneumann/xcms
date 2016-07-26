## All methods related to data import/export should go here.
#' @include AllGenerics.R DataClasses.R



############################################################
##
##      xcmsRaw
##
############################################################

############################################################
## write.cdf
setMethod("write.cdf", "xcmsRaw", function(object, filename) {
    require(ncdf4) || stop("Couldn't load package ncdf4 for NetCDF writing")

    scan_no <- length(object@scanindex)
    point_no <- length(object@env$mz)

    dim32bytes <- ncdim_def("_32_byte_string", "", 1:32, create_dimvar=FALSE)
    dim64bytes <- ncdim_def("_64_byte_string", "", 1:64, create_dimvar=FALSE)
    dimError   <- ncdim_def("error_number",       "", 1:1, create_dimvar=FALSE)
    dimScans   <- ncdim_def("scan_number",     "", 1:scan_no, create_dimvar=FALSE)
    dimPoints  <- ncdim_def("point_number",    "", 1:point_no, create_dimvar=FALSE)
    dimInstruments<-ncdim_def("instrument_number","",1:1, create_dimvar=FALSE)

    ## Define netCDF vars
    error_log             <- ncvar_def("error_log","", list(dim64bytes, dimError), missval=NULL, prec="char")
    scan_acquisition_time <- ncvar_def("scan_acquisition_time", "", dimScans, -1, prec="double")
    total_intensity       <- ncvar_def("total_intensity", "Arbitrary Intensity Units", dimScans, -1, prec="double")

    mass_range_min       <- ncvar_def("mass_range_min", "", dimScans, NULL, prec="double")
    mass_range_max       <- ncvar_def("mass_range_max", "", dimScans, NULL, prec="double")
    time_range_min       <- ncvar_def("time_range_min", "", dimScans, NULL, prec="double")
    time_range_max       <- ncvar_def("time_range_max", "", dimScans, NULL, prec="double")

    scan_index            <- ncvar_def("scan_index", "", dimScans, missval=-1, prec="integer")
    actual_scan_number    <- ncvar_def("actual_scan_number", "", dimScans, missval=-1, prec="integer")
    mass_values           <- ncvar_def("mass_values", "M/Z", dimPoints, -1)
    intensity_values      <- ncvar_def("intensity_values", "Arbitrary Intensity Units", dimPoints, -1)

    point_count           <- ncvar_def("point_count", "", dimScans, missval=-1, prec="integer")
    flag_count            <- ncvar_def("flag_count", "", dimScans, missval=-1, prec="integer")

    instrument_name       <- ncvar_def("instrument_name","", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_id         <- ncvar_def("instrument_id", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_mfr        <- ncvar_def("instrument_mfr", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_model      <- ncvar_def("instrument_model", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_serial_no  <- ncvar_def("instrument_serial_no", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_sw_version <- ncvar_def("instrument_sw_version", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_fw_version <- ncvar_def("instrument_fw_version", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_os_version <- ncvar_def("instrument_os_version", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_app_version<- ncvar_def("instrument_app_version", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")
    instrument_comments   <- ncvar_def("instrument_comments", "", list(dim32bytes, dimInstruments), missval=NULL, prec="char")

    ## Define netCDF definitions
    ms <- nc_create(filename, list(error_log,
                                     scan_acquisition_time,
                                     actual_scan_number,
                                     total_intensity,
                                     mass_range_min,mass_range_max,time_range_min,time_range_max,
                                     scan_index, point_count, flag_count,
                                     mass_values, intensity_values,
                                     instrument_name, instrument_id,instrument_mfr,instrument_model,
                                     instrument_serial_no,instrument_sw_version,instrument_fw_version,
                                     instrument_os_version,instrument_app_version,instrument_comments
                                     ))

    ## Add values to netCDF vars
    ncvar_put(ms, "scan_acquisition_time", object@scantime)
    ncvar_put(ms, "total_intensity", object@tic)
    ncvar_put(ms, "scan_index", object@scanindex)
    ncvar_put(ms, "actual_scan_number", object@scanindex)

    ncvar_put(ms, "time_range_min", rep(-9999, times=scan_no))
    ncvar_put(ms, "time_range_max", rep(-9999, times=scan_no))

    ncvar_put(ms, "point_count", c(object@scanindex[2:scan_no], point_no)
                 - object@scanindex[1:scan_no])
    ncvar_put(ms, "flag_count", rep(0, times=scan_no))

    mzranges <- t(sapply(1:scan_no,
                         function(scan) range(getScan(object, scan)[,"mz"])))

    ncvar_put(ms, "mass_range_min", mzranges[,1])
    ncvar_put(ms, "mass_range_max", mzranges[,2])

    ncvar_put(ms, "mass_values", object@env$mz)
    ncatt_put(ms, "mass_values", "scale_factor", 1, prec="float")

    ncvar_put(ms, "intensity_values", object@env$intensity)
    ncatt_put(ms, "intensity_values", "add_offset", 0, prec="float")
    ncatt_put(ms, "intensity_values", "scale_factor", 1, prec="float")

    ## Add ANDIMS global attributes to netCDF object
    ncatt_put(ms, 0, "dataset_completeness", "C1+C2")
    ncatt_put(ms, 0, "ms_template_revision", "1.0.1")
    ncatt_put(ms, 0, "netcdf_revision", "2.3.2")
    ncatt_put(ms, 0, "languages", "English")
    ncatt_put(ms, 0, "raw_data_mass_format", "Float")
    ncatt_put(ms, 0, "raw_data_time_format", "Short")
    ncatt_put(ms, 0, "raw_data_intensity_format", "Float")
    ncatt_put(ms, 0, "units", "Seconds")
    ncatt_put(ms, 0, "starting_scan_number", "0")
    ncatt_put(ms, 0, "global_mass_min", as.character((min(object@env$mz))))
    ncatt_put(ms, 0, "global_mass_max", as.character((max(object@env$mz))))
    ncatt_put(ms, 0, "global_intensity_min", as.character((min(object@env$intensity))))
    ncatt_put(ms, 0, "global_intensity_max", as.character((max(object@env$intensity))))
    ncatt_put(ms, 0, "calibrated_mass_min", as.character((min(object@env$intensity))))
    ncatt_put(ms, 0, "calibrated_mass_max", as.character((max(object@env$intensity))))
    ncatt_put(ms, 0, "actual_run_time_length", object@scantime[scan_no]-object@scantime[1])
    ncatt_put(ms, 0, "actual_delay_time", "1.")
    ncatt_put(ms, 0, "raw_data_uniform_sampling_flag", "0s")

    nc_close(ms)
})


