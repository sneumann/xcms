
#' Convert a `XcmsExperiment` to an `XcmsExperimentHdf5` object: export all
#' data to a HDF5 file `h5_file` and return the `XcmsExperimentHdf5`.
#'
#' @noRd
.xcms_experiment_to_hdf5 <- function(x, h5_file = character()) {
    if (!length(h5_file))
        stop("Parameter 'h5_file' is mandatory.")
    .h5_initialize_file(h5_file)
    has_chrom_peaks <- hasChromPeaks(x)
    has_features <- hasFeatures(x)
    x <- as(x, "XcmsExperimentHdf5")
    x@sample_id <- seq_along(x)
    x@hdf5_file <- h5_file
    mod_count <- 0L
    if (has_chrom_peaks) {
        ## Check if we need to change the rownames to the expected format
        ## CP<MS level><index>
        first_cpid <- rownames(x@chromPeaks)[1L]
        if (nchar(first_cpid) == 3 || substring(first_cpid, 3, 3) == "0")
            rownames(x@chromPeaks) <-
                paste0("CP", x@chromPeakData$ms_level,
                       substring(rownames(x@chromPeaks), 3))
        ## Memory-efficient export: save the data for one sample at a time. If
        ## that is too slow we could split the data and export all in one go.
        is_sample <- colnames(x@chromPeaks) == "sample"
        msl <- unique(x@chromPeakData$ms_level)

        for (i in x@sample_id) {
            idx <- unname(which(x@chromPeaks[, is_sample] == i))
            pks <- x@chromPeaks[idx, !is_sample, drop = FALSE]
            pkd <- x@chromPeakData[idx, , drop = FALSE]
            f <- factor(pkd$ms_level, levels = msl)
            pks <- split.data.frame(pks, f)
            pkd <- split.data.frame(
                pkd[, colnames(pkd) != "ms_level", drop = FALSE], f)
            names(pks) <- i
            names(pkd) <- i
            mod_count <- .h5_write_data(
                h5_file, pks, name = "chrom_peaks", ms_level = msl)
            mod_count <- .h5_write_data(
                h5_file, pkd, name = "chrom_peak_data", ms_level = msl)
        }
        slot(x, "chromPeaks", check = FALSE) <-
            x@chromPeaks[integer(), , drop = FALSE]
        slot(x, "chromPeakData", check = FALSE) <-
            x@chromPeakData[integer(), , drop = FALSE]
        slot(x, "has_chrom_peaks", check = FALSE) <- TRUE
    }
    if (has_features) {
        stop("Can not yet save feature definitions to HDF5")
        slot(x, "has_features", check = FALSE) <- TRUE
    }
    x@hdf5_mod_count <- mod_count
    x
}

#' Subset an `XcmsExperimentHdf5` object. Similar to `.subset_xcms_experiment()`
#' for `XcmsExperiment`, but optimized for `XcmsExperimentHdf5`.
#'
#' @noRd
.h5_subset_xcms_experiment <- function(x, i = integer(),
                                       keepChromPeaks = TRUE,
                                       keepAdjustedRtime = FALSE,
                                       keepFeatures = FALSE,
                                       ignoreHistory = FALSE,
                                       ...) {
    i <- i2index(i, length(x))
    if (any(i < 0)) {
        if (all(i < 0))
            i <- seq_along(x)[i]
        else stop("Mixing positive and negative indices is not supported.")
    }
    drop <- character()
    if (!keepAdjustedRtime && hasAdjustedRtime(x)) {
        svs <- unique(c(spectraVariables(x@spectra), "mz", "intensity"))
        x@spectra <- selectSpectraVariables(
            x@spectra, svs[svs != "rtime_adjusted"])
        drop <- c(drop, .PROCSTEP.RTIME.CORRECTION)
    }
    if (!keepFeatures && hasFeatures(x)) {
        stop("Subsetting with features present needs to be implemented")
        drop <- c(drop, .PROCSTEP.PEAK.GROUPING)
    }
    if (!keepChromPeaks && hasChromPeaks(x)) {
        x@has_chrom_peaks <- FALSE
        drop <- c(drop, .PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.FILLING,
                  .PROCSTEP.CALIBRATION, .PROCSTEP.PEAK.REFINEMENT)
    }
    if (!ignoreHistory && length(drop))
        x@processHistory <- dropProcessHistoriesList(
            x@processHistory, type = drop)
    x@sample_id <- x@sample_id[i]
    getMethod("[", "MsExperiment")(x, i = i)
}

#' Similar to `.xmse_merge_neighboring_peaks()` in XcmsExperiment-functions.R,
#' but this does not return the chromatograpic peaks but stores them into
#' the HDF5 file instead.
#'
#' @param x `XcmsExperimentHdf5` object with potentially multiple samples.
#'
#' @noRd
.h5_xmse_merge_neighboring_peaks <- function(x, msLevel = 1L, expandRt = 2,
                                             expandMz = 0, ppm = 10,
                                             minProp = 0.75,
                                             BPPARAM = bpparam()) {
    keep <- msLevel(spectra(x)) == msLevel
    f <- as.factor(fromFile(x)[keep])
    if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
    else rt <- rtime(spectra(x))[keep]
    ## Get the list of chromPeaks for x.
    pksl <- .h5_read_data(
        x@hdf5_file, index = x@sample_id, name = "chrom_peaks",
        ms_level = rep(msLevel, length(x@sample_id)),
        read_colnames = TRUE, read_rownames = TRUE)
    prefix <- paste0("CP", msLevel)
    cp_id <- max(c(0L, vapply(pksl, function(z)
        max(as.integer(sub(prefix, "", rownames(z)))), 1L)))
    ## Get the list of chromPeakData for x.
    pkdl <- .h5_read_data(
        x@hdf5_file, index = x@sample_id, name = "chrom_peak_data",
        ms_level = rep(msLevel, length(x@sample_id)), read_rownames = TRUE)
    ## Do refinement (in parallel)
    res <- bpmapply(
        xcms:::.merge_neighboring_peaks2,
        split(peaksData(filterMsLevel(spectra(x), msLevel = msLevel),
                        f = factor()), f), pksl, pkdl, split(rt, f),
        MoreArgs = list(expandRt = expandRt, expandMz = expandMz,
                        ppm = ppm, minProp = minProp),
        SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    ## Replace data in hdf5 for samples with changed data.
    has_merged <- which(
        vapply(res, function(z) any(is.na(rownames(z[[1L]]))), NA))
    for (i in has_merged) {
        l <- list(res[[i]]$chromPeaks)
        names(l) <- x@sample_id[i]
        .h5_write_data(h5_file = x@hdf5_file, data_list = l,
                       name = "chrom_peaks", ms_level = msLevel,
                       replace = TRUE, write_colnames = FALSE,
                       write_rownames = FALSE)
        pkd <- res[[i]]$chromPeakData
        if (!any(colnames(pkd) == "merged"))
            pkd$merged <- FALSE
        pkd$merged[grep("^CP", rownames(pkd), invert = TRUE)] <- TRUE
        l <- list(pkd)
        names(l) <- x@sample_id[i]
        .h5_write_data(h5_file = x@hdf5_file, data_list = l,
                       name = "chrom_peak_data", ms_level = msLevel,
                       replace = TRUE)
    }
    ## Report the highest CP number back.
    cp_id
}

#' Extract the `chromPeaks` `matrix` of selected samples.
#'
#'
#' @param by_sample `logical(1)` whether a `list` of `chromPeak` matrices split
#'     per sample should be returned or the *conventional* matrix with an
#'     additional column `"sample"`.
#'
#' @noRd
## .h5_chrom_peaks <- function(x, columns = character(), by_sample = TRUE) {
##     h5 <- rhdf5::H5Fopen(x@hdf5_file)
##     .h5_check_mod_count(h5, x@hdf5_mod_count)
##     grps <- .h5_dataset_names("/", h5)
##     rhdf5::H5Fclose(h5)
##     msl <- sort(as.integer(sub("ms_", "", grep("^ms_", grps, value = TRUE))))
##     ids <- rep(x@sample_id, length(msl))
##     msl <- rep(msl, each = length(x@sample_id))
##     res <- .h5_read_data(x@hdf5_file, index = ids, name = "chrom_peaks",
##                          ms_level = msl, read_colnames = TRUE,
##                          read_rownames = TRUE)
## }

################################################################################
##
##        HDF5 FUNCTIONALITY
##
################################################################################

#' Properties of the HDF5 file used for on-disk storage of xcms results:
#' - all preprocessing results are stored within the same file.
#' - storage of chrom peak detection results are organized by MS level and
#'   sample:
#'   /ms_<ms_level>/<sample id>/chrom_peaks (float array)
#'   /ms_<ms_level>/<sample id>/chrom_peaks_rownames (character array)
#'   /ms_<ms_level>/<sample id>/chrom_peaks_colnames (character array)
#'   /ms_<ms_level>/<sample id>/chrom_peak_data (list of arrays).

.h5_have_rhdf5 <- function() {
    return(requireNamespace("rhdf5", quietly = TRUE))
}

.h5_require_rhdf5 <- function() {
    if (!.h5_have_rhdf5())
        stop("Package 'rhdf5' is required for this functionality. Please ",
             "install using 'BiocManager::install(\"rhdf5\")' and try again.")
}

##  --------  READING  --------

#' Reads a single chrom peaks `matrix` from the HDF5 file. Depending on
#' `read_colnames` and `read_rownames` also the rownames and colnames are read.
#' Not reading them has performance advantages.
#'
#' @param name `character(1)` with the name of the data set to read.
#'
#' @param h5 HDF5 file handle
#'
#' @param index `integer` with the index (or indices) of the columns to read.
#'     By default, with `index = NULL` all columns are read.
#'
#' @param read_colnames `logical(1)` whether column names should also be read
#'     and set.
#'
#' @param read_rownames `logical(1)` whether row names should be read and set.
#'
#' @return numeric `matrix`
#'
#' @noRd
.h5_read_chrom_peaks <- function(name, h5, index = NULL,
                                 read_colnames = FALSE,
                                 read_rownames = FALSE) {
    d <- rhdf5::h5read(h5, name = name, index = list(NULL, index))
    if (read_rownames)
        rownames(d) <- as.vector(
            rhdf5::h5read(h5, name = paste0(name, "_rownames")))
    if (read_colnames) {
        cn <- as.vector(
            rhdf5::h5read(h5, name = paste0(name, "_colnames")))
        if (length(index))
            colnames(d) <- cn[index]
        else colnames(d) <- cn
    }
    d
}

#' Read a single chromPeakData `data.frame` from the HDF5 file. With
#' `read_rownames = TRUE` also the row names are read and set, which requires
#' an additional reading step. Note that for a `data.frame` each column
#' is stored as a separate array/DATASET type for the `data.frame` GROUP. Thus,
#' a single column can be simply read by specifying the `name` accordingly,
#' e.g. <data set name>/<column name>.
#'
#' @param name `character(1)` the name of the data set.
#'
#' @param h5 HDF5 file handle
#'
#' @param read_rownames `logical(1)` whether rownames should be read and set.
#'
#' @noRd
.h5_read_chrom_peak_data <- function(name, h5, read_rownames = FALSE, ...) {
    d <- as.data.frame(rhdf5::h5read(h5, name = name))
    if (read_rownames)
        rownames(d) <- as.vector(
            rhdf5::h5read(h5, sub("_data", "s_rownames", name)))
    d
}

#' Reads the names of the data sets of a group. This can for example be used
#' to get the column names of a `data.frame` that was saved as a GROUP of
#' DATASETs.
#'
#' @noRd
.h5_dataset_names <- function(name, h5, recursive = FALSE) {
    g <- rhdf5::H5Gopen(h5, name)
    on.exit(rhdf5::H5Gclose(g))
    rhdf5::h5ls(g, recursive = recursive, datasetinfo = FALSE)$name
}

.h5_ms_levels <- function(h5) {
    nms <- .h5_dataset_names("/", h5)
    as.integer(unique(sub("ms_", "", grep("^ms", nms, value = TRUE))))
}

#' Read selected datasets from a HDF5 file.
#'
#' @note
#'
#' Setting row and column names requires additional read steps and is hence
#' considerably slower than *just* importing the data. Also, if possible,
#' consider importing only the required column(s) using parameter `column`.
#'
#' @param h5_file `character(1)` with the HDF5 file name
#'
#' @param index `integer` with the indices/IDs of the data sets to read.
#'
#' @param name `character(1)` specifying which data should be read.
#'
#' @param ms_level `integer` with the MS level of each sample/data set that
#'     should be read. Has to have the same length than `index`.
#'
#' @param read_colnames `logical(1)` whether column names should be read and
#'     set for each `matrix`.
#'
#' @param read_rownames `logical(1)` whether row names should be read and
#'     set for each `matrix` or `data.frame`
#'
#' @param column For `name = "chrom_peak_data"`: `character(1)` allowing to
#'     select a **single** column to read. For `name = "chrom_peaks"`: `integer`
#'     with the indices of the column(s) that should be imported.
#'
#' @return `list()` with the read datasets. Will be a `list` or `numeric`
#'     matrices for `name = "chrom_peaks"` or a `list` with `data.frame`s for
#'     `name = "chrom_peak_data"`.
#'
#' @noRd
.h5_read_data <- function(h5_file = character(),
                          index = integer(),
                          name = c("chrom_peaks", "chrom_peak_data"),
                          ms_level = integer(),
                          read_colnames = FALSE,
                          read_rownames = FALSE,
                          column = NULL) {
    if (!length(index)) return(list())
    stopifnot(length(ms_level) == length(index))
    name <- match.arg(name)
    FUN <- .h5_read_chrom_peaks
    if (name == "chrom_peak_data")
        FUN <- .h5_read_chrom_peak_data
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    d <- paste0("/ms_", ms_level, "/", index, "/", name)
    if (is.character(column) && length(column) == 1L)
        d <- paste0(d, "/", column)
    lapply(d, FUN = FUN, read_colnames = read_colnames,
           read_rownames = read_rownames, index = column, h5 = h5)
}


##  --------  VALIDITY  --------

#' Compares the "mod_count" attribute from an h5 file with the expected
#' one. Throws an error if it is different.
#'
#' The `mod_count` gets incremented by any operation that writes data to the
#' HDF5 file. If the `mod_count` of the xcms result object is different than
#' the one in the HDF5 file an error is thrown.
#'
#' @param h5 file handle for the HDF5 file
#'
#' @param mod_count `integer(1)` the expected value for the modification
#'     counter.
#'
#' @return `TRUE` or throws an error.
#'
#' @noRd
.h5_check_mod_count <- function(h5, mod_count = 0L) {
    mc <- rhdf5::h5read(h5, "/header/modcount")[1L]
    if (mc != mod_count)
        stop("The HDF5 file was changed by a different process. This xcms ",
             "result object/variable is no longer valid.")
    TRUE
}

#' Checks the HDF5 file for validity.
#'
#' @param h5_file `character(1)` with the file name
#'
#' @param mod_count `integer(1)` with the expected value for the modification
#'     counter
#'
#' @return `TRUE` if the file is valid or throws an error if the file is
#'     not valid.
#'
#' @noRd
.h5_valid_file <- function(h5_file, mod_count = 0L) {
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    if (!rhdf5::H5Lexists(h5, "/header/package"))
        stop("File \"", h5_file, "\" is not in correct format.")
    res <- rhdf5::h5read(h5, "/header/package")
    if (res != "package:xcms")
        stop("File \"", h5_file, "\" is not in correct format.")
    .h5_check_mod_count(h5, mod_count = mod_count)
    TRUE
}

##  --------  WRITING  --------

.h5_compression_level <- function() 0L

#' Initializes the HDF5 file
#'
#' @noRd
.h5_initialize_file <- function(x, mod_count = 0L) {
    if (file.exists(x))
        stop("File \"", x, "\" already exists. Please choose a different name ",
             "or remove that file first.")
    h5 <- rhdf5::H5Fcreate(x)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    comp_level <- .h5_compression_level()
    rhdf5::h5createGroup(h5, "header")
    rhdf5::h5write("package:xcms", h5, "/header/package", level = comp_level)
    rhdf5::h5write(mod_count, h5, "/header/modcount", level = comp_level)
}

#' Every writing operation should increate the "mod count", i.e. the
#' count of data modifications. This function increases the mod count by +1
#' and returns this value.
#'
#' @noRd
.h5_increment_mod_count <- function(h5) {
    mc <- rhdf5::h5read(h5, "/header/modcount")[1L] + 1L
    rhdf5::h5write(mc, h5, "/header/modcount",
                   level = .h5_compression_level())
    mc
}

#' Bare writing function of the chromPeaks `matrix`
#'
#' @param x `matrix` with the chrom peak results **of a single sample**
#'
#' @param h5 HDF5 file handle.
#'
#' @param name `character(1)` with the name for the data set (e.g.
#'     `"/ms_1/1/chrom_peaks"`.
#'
#' @param level `integer(1)` with the compression level.
#'
#' @param write_colnames `logical(1)` whether to write the column names of `x`
#'     as an additional data set `paste0(name, "_colnames")`.
#'
#' @param write_rownames `logical(1)` whether to write the rownames of `x` as
#'     an additional data set `paste0(name, "_rownames")`.
#'
#' @noRd
.h5_write_chrom_peaks <- function(x, h5, name, level, write_colnames = TRUE,
                                  write_rownames = TRUE) {
    rhdf5::h5write(x, h5, name = name, level = level,
                   write.attributes = FALSE,
                   createnewfile = FALSE)
    if (write_rownames)
        rhdf5::h5write(rownames(x), h5, name = paste0(name, "_rownames"),
                       level = level, createnewfile = FALSE)
    if (write_colnames)
        rhdf5::h5write(colnames(x), h5, name = paste0(name, "_colnames"),
                       level = level, createnewfile = FALSE)
}

#' Bare writing function of the chromPeakData DATA.FRAME
#'
#' @param x `data.frame`
#'
#' @param h5 HDF5 file handle.
#'
#' @param name `character(1)` defining the name for the data set.
#'
#' @param level `integer(1)` with the compression level.
#'
#' @noRd
.h5_write_chrom_peak_data <- function(x, h5, name, level, ...) {
    rhdf5::h5writeDataset(x, h5, name, level = level,
                          DataFrameAsCompound = FALSE)
}

#' Writes data (matrix, data.frame) to a HDF5 file organized by MS level and
#' sample:
#'
#' /ms_<MS level>/<sample id>/<data set>
#'
#' Example:
#'
#' /ms_1/1/chrom_peaks
#' /ms_1/1/chrom_peak_data
#' /ms_1/2/chrom_peaks
#' /ms_1/2/chrom_peak_data
#' /ms_2/1/chromPeaks
#' ...
#'
#' @param h5_file the HDF5 file.
#'
#' @param data_list `list` of data sets (`matrix` or `data.frame`) that should
#'     be written. The names of the `list` are used as sample identifier.
#'
#' @param name `character(1)` with the name of the data set. This might be
#'     changed in future.
#'
#' @param ms_level `integer` with **same length** than `data_list` defining
#'     the MS level of the data sets.
#'
#' @param replace `logical(1)` whether eventual existing data sets should be
#'     replaced or updated. Data must be replaced if the dimensions are
#'     different.
#'
#' @noRd
.h5_write_data <- function(h5_file = character(),
                           data_list = list(),
                           name = c("chrom_peaks", "chrom_peak_data"),
                           ms_level = integer(),
                           replace = TRUE,
                           write_colnames = replace,
                           write_rownames = replace) {
    if (!length(data_list)) return(TRUE)
    stopifnot(length(ms_level) == length(data_list))
    stopifnot(!is.null(names(data_list)))
    name <- match.arg(name)
    FUN <- .h5_write_chrom_peaks
    if (name == "chrom_peak_data")
        FUN <- .h5_write_chrom_peak_data
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    comp_level <- .h5_compression_level()
    for (i in seq_along(data_list)) {
        group_ms <- paste0("/ms_", ms_level[i])
        if (!rhdf5::H5Lexists(h5, group_ms))
            rhdf5::h5createGroup(h5, group_ms)
        group_sample <- paste0(group_ms, "/", names(data_list)[i])
        if (!rhdf5::H5Lexists(h5, group_sample))
            rhdf5::h5createGroup(h5, group_sample)
        group_data <- paste0(group_sample, "/", name)
        if (replace && rhdf5::H5Lexists(h5, group_data))
            rhdf5::h5delete(h5, group_data)
        FUN(data_list[[i]], h5, group_data, comp_level,
            write_colnames = write_colnames, write_rownames = write_rownames)
    }
    .h5_increment_mod_count(h5)
}
