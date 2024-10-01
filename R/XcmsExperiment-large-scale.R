#' Very large scale data support for `XcmsExperiment`: `chromPeaks` and
#' `chromPeakData` should not be stored in memory but on disk.
#'
#' The data should be organized in the HDF5 file:
#' /header/
#' /ms_1/1/chrom_peaks
#'        /chrom_peak_data
#'       2/
#'       3/
#' /ms_2/
#'
#' So, essentially, the chromPeaks are now organized by sample index and by MS
#' level.
#'
#' FEATURE DEFINITIONS
#'
#' also a Hdf5 file, eventually the same!
#' have the mapping between features and chrom peaks as n:m two column table.
#'
#' or save the feature_index for each sample? that could then be loaded and
#' processed to extract the values/chrom peak part. that would allow again the
#' (same) chunk processing in which we load data for sets of samples and only
#' load their chrom peak matrix. the n:m matrix should however also have NA
#' integers if for one feature no chrom peak was available. Also, other
#' operations might be easier, such as gap filling as it can be done (in
#' isolation) separately per chunk.
#'
#'
#' What will the new object need:
#'
#' - extend XcmsExperiment (XcmsExperimentHdf5
#' - slot h5file with the name of the HDF5file for the chrom Peaks
#' - slot mod_count with the counter for HDF5file modification.
#' - slot with number of peaks per sample/MS level?
#' - sample_id: an integer ID for each sample that can be used for subsetting
#'   etc.
#'
#'
#' TODO: check the existing code to see how we can manage certain preprocessing
#' steps:
#'
#' * findChromPeaks
#'
#'   .mse_find_chrom_peaks_chunks (MsExperiment-functions.R); combines the
#'   peak matrices to a single matrix.
#'   Calls .mse_spectrapply_chunks (MsExperiment-functions.R) (with
#'   .mse_find_chrom_peaks_chunk (MsExperiment-functions.R): this one runs
#'   peak detection and returns them as a list.
#'
#'   Could implement a .collect_results function.
#'
#' * refineChromPeaks, MergeNeighboringPeaksParam. there are also others that
#'   might different treatment.
#'
#'   .xmse_apply_chunks (XcmsExperiment-functions.R) iterates over chunks and
#'   calls. calls .xmse_merge_neightboring_peaks on each XcmsExperiment subset
#'   (.subset_xcms_experiment). Will need a proper `[` method.
#'
#' * adjustRtime: PeakGroupsParam, gets the full raw retention times of all
#'   spectra. Question is if we can't do also this chunk-wise per file? We might
#'   need to rewrite this to enable chunk wise processing!
#'
#' * groupChromPeaks
#'
#' * fillChromPeaks

library(rhdf5)

.h5_compression_level <- function() 0L

#' Compares the "mod_count" attribute from an h5 file with the expected
#' one. Throws an error if it is different.
#'
#' The `mod_count` gets incremented by any operation that writes data to the
#' HDF5 file. If the `mod_count` of the xcms result object is different than
#' the one in the HDF5 file an error is thrown.
#'
#' @noRd
.h5_check_mod_count <- function(h5, mod_count = 0L) {
    mc <- rhdf5::h5read(h5, "/header/modcount")[1L]
    if (mc != mod_count)
        stop("The HDF5 file was changed by a different process. This xcms ",
             "result object/variable is no longer valid.")
}

.h5_increment_mod_count <- function(h5) {
    mc <- rhdf5::h5read(h5, "/header/modcount")[1L]
    rhdf5::h5write(mc + 1L, h5, "/header/modcount",
                   level = .h5_compression_level())
}

.h5_valid_file <- function(x, mod_count = 0L) {
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    if (!rhdf5::H5Lexists(h5, "/header/class"))
        stop("File ", basename(x), " is not in correct format.")
    res <- rhdf5::h5read(h5, "/header/class")
    if (res != "xcms::chromPeaks")
        stop("File ", basename(x), " is not in correct format.")
    .h5_check_mod_count(h5, mod_count = mod_count)
    character()
}

#' Initializes a h5 file
.initialize_chrom_peaks_h5_file <- function(x, mod_count = 0L) {
    h5 <- rhdf5::H5Fcreate(x)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    comp_level <- .h5_compression_level()
    rhdf5::h5createGroup(h5, "header")
    rhdf5::h5write("xcms::chromPeaks", h5, "/header/class", level = comp_level)
    rhdf5::h5write(mod_count, h5, "/header/modcount", level = comp_level)
}

#' Bare writing function of the chromPeaks MATRIX
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
    ## Note: saving rownames/colnames as attributes is not a good idea; they
    ## can not be larger than a certain size and can also only be read
    ## "as a whole"
}

#' Bare writing function of the chromPeakData DATA.FRAME
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
#'     replaced or updated. Data needs to be replaced if the dimensions are
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
    ##rhdf5::h5closeAll()
    TRUE
}

.h5_read_chrom_peaks <- function(name, h5, dimnames = TRUE) {
    d <- rhdf5::h5read(h5, name = name)
    if (dimnames) {
        attrs <- rhdf5::h5readAttributes(h5, name)
        rownames(d) <- attrs$rownames
        colnames(d) <- attrs$colnames
    }
    d
}

.h5_read_chrom_peak_data <- function(name, h5, dimnames = TRUE) {
    d <- as.data.frame(rhdf5::h5read(h5, name = name))
    if (dimnames) {
        attrs <- rhdf5::h5readAttributes(h5, sub("_data", "s", name))
        rownames(d) <- attrs$rownames
    }
    d
}

.h5_read_data <- function(h5_file = character(),
                          name = c("chrom_peaks", "chrom_peak_data"),
                          id = character(),
                          ms_level = integer(),
                          dimnames = TRUE) {
    stopifnot(length(id) == length(ms_level))
    name <- match.arg(name)
    FUN <- .h5_read_chrom_peaks
    if (name == "chrom_peak_data")
        FUN <- .h5_read_chrom_peak_data
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    nms <- paste0("/ms_", ms_level, "/", id, "/", name)
    l <- lapply(nms, FUN, h5 = h5, dimnames = dimnames)
    names(l) <- id
    l
}

################################################################################
#' Tests
#'
h5f <- "/home/jo/tmp/test.h5"
.initialize_chrom_peaks_h5_file(h5f)
h5ls(h5f)

a <- cbind(mz = c(1.12, 1.42, 13.2, 31.31),
           mzmin = c(1, 1, 13.1, 31.1),
           mzmax = c(1.13, 1.45, 13.21, 31.32),
           rt = c(12.2, 20.4, 34.1, 140.3),
           rtmin = c(10, 15.3, 30.2, 135.23),
           rtmax= c(15, 23.2, 39.2, 145.23),
           into = c(123.4, 4332, 543211, 123))
a_pd <- data.frame(ms_level = rep(1, 4), is_filled = rep(FALSE, 4))
rownames(a) <- rownames(a_pd) <- c("CP001", "CP002", "CP003", "CP004")


b <- a * 2
b_pd <- a_pd
rownames(b) <- rownames(b_pd) <- c("CP005", "CP006", "CP007", "CP008")

l <- list(a, b, a, b)
names(l) <- 1:4
l_pd <- list(a_pd, b_pd)
names(l_pd) <- 1:2

.h5_write_data(h5f, l, name = "chrom_peaks", ms_level = c(1L, 1L, 1L, 1L),
               replace = TRUE)
h5readAttributes(h5f, name = "/ms_1/1/chrom_peaks/")

.h5_write_data(h5f, l_pd, name = "chrom_peak_data",
               ms_level = c(1L, 1L), replace = TRUE)


.h5_read_data(h5f, "2", name = "chrom_peaks", ms_level = 1L)
.h5_read_data(h5f, c(1, 2), name = "chrom_peaks", ms_level = c(1L, 1L),
              dimnames = FALSE)
.h5_read_data(h5f, "2", name = "chrom_peak_data", ms_level = 1L,
              dimnames = FALSE)
.h5_read_data(h5f, c(2, 1), name = "chrom_peak_data", ms_level = c(1L, 1L),
              dimnames = TRUE)


h5 <- H5Fopen(h5f)
h5ls(h5)
h5read(h5, name = "/ms_1/2/chrom_peaks", read.attributes = TRUE)
h5readAttributes(h5, "/ms_1/2/chrom_peaks/")$colnames
h5read(h5, name = "/ms_1/2/chrom_peak_data")

h5writeDataset(a_pd, h5, name = "/test3", DataFrameAsCompound = TRUE)

.h5_check_mod_count(h5)
.h5_increment_mod_count(h5)
.h5_check_mod_count(h5)
.h5_check_mod_count(h5, 1L)

H5Fclose(h5)


library(microbenchmark)
cps <- do.call(rbind, l)
cps <- cbind(cps, sample = c(1, 1, 1, 1, 2, 2, 2, 2))
microbenchmark(
    .h5_read_chrom_peaks(h5f, "2", name = "chrom_peaks", ms_level = 1L),
    split.data.frame(cps, cps[, "sample"])[2L]
)



#' What next:
#'
#' 1) create a h5 file with the full (40GB?) data.
#'
#' Size of chromPeak matrix: 27.7GB
print(object.size(chromPeaks(chris)), unit = "GB")
#' Size of chromPeakData data.frame: 13.8GB
print(object.size(chromPeakData(chris)), unit = "GB")

f <- factor(chromPeaks(chris)[, "sample"], levels = seq_along(chris))
pkl <- split.data.frame(chromPeaks(chris), f)
pkl <- lapply(pkl, function(z) z[, colnames(z) != "sample", drop = FALSE])
h5f <- "chris.h5"
.initialize_chrom_peaks_h5_file(h5f)
h5ls(h5f)
system.time(
    .h5_write_data(h5f, pkl, name = "chrom_peaks",
                   ms_level = rep(1L, length(pkl)))
) # 199 sec
pkd <- split.data.frame(chromPeakData(chris, return.type = "data.frame"), f)
pkd <- lapply(pkd, function(z) z[, colnames(z) != "ms_level", drop = FALSE])
system.time(
    .h5_write_data(h5f, pkd, name = "chrom_peak_data",
                   ms_level = rep(1L, length(pkd)), replace = TRUE)
) # 30 sec


#' 2) subset the object to 1000QC samples and save h5 of that.
chris_qc <- chris[which(sampleData(chris)$sample_type == "Pool")]

f <- factor(chromPeaks(chris_qc)[, "sample"], levels = seq_along(chris_qc))
pkl <- split.data.frame(chromPeaks(chris_qc), f)
pkl <- lapply(pkl, function(z) z[, colnames(z) != "sample", drop = FALSE])
h5f <- "chris_qc.h5"
.initialize_chrom_peaks_h5_file(h5f)
h5ls(h5f)
system.time(
    .h5_write_data(h5f, pkl, name = "chrom_peaks",
                   ms_level = rep(1L, length(pkl)))
) # 12 sec
pkd <- split.data.frame(chromPeakData(chris_qc, return.type = "data.frame"), f)
pkd <- lapply(pkd, function(z) z[, colnames(z) != "ms_level", drop = FALSE])
system.time(
    .h5_write_data(h5f, pkd, name = "chrom_peak_data",
                   ms_level = rep(1L, length(pkd)), replace = TRUE)
) # 3 sec
h5ls(h5f)

#' 3) test timings to access chrom peaks.
#' 4) implement a chromPeaks method.
#' 5) implement a setAs method that converts XcmsExperiment to
#' XcmsExperimentHdf5.
#' 6) implement a refineChromPeaks method. See how/if that is better/faster.