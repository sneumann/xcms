#' @title Store/load contents of `MsExperiment` and `XcmsExperiment` objects
#' as/from plain text files
#'
#' @name PlainTextParam
#'
#' @export
#'
#' @family xcms result export and import formats.
#'
#' @description
#' The `PlainTextParam` class and method enable users to save/load an
#' `MsExperiment` or `XcmsExperiment` object as a collections of plain text
#' files in/from a specified folder.
#' Note that, while for all xcms results within the `XcmsExperiment` can
#' and will be exported, the full raw MS data (of the object's `Spectra` object)
#' will currently not be exported in plain text format. For `Spectra` using the
#' [MsBackendMzR()] backend, the names of the raw data files will however be
#' exported (which enables to *load* the full `Spectra` respectively
#' `MsExperiment` objects).
#'
#' For an `MsExperiment` object, the exported files include:
#'
#' - The [sampleData()] stored as a text file named *sample_data.txt*.
#'
#' - The [fileNames()] of the *Spectra* object stored in a tabular format in a
#' text file named *spectra_files.txt*.The file names will only be exported if
#' the `Spectra` object uses a [MsBackendMzR()] backend. For other backends no
#' information on raw spectra data is currently exported with `PlainTextParam`.
#'
#' - Processing queue of the `Spectra` object, ensuring that any spectra data
#' modifications are retained. It is stored in a `json` file named
#' *spectra_processing_queue.json*.
#'
#' For an `XcmsExperiment` object, the exported files are the same as those
#' for an `MsExperiment` object, with the addition of the following:
#'
#' - The [processHistory()] information of the object, stored in a `json` file
#' named *process_history.json*.
#'
#' - The chromatographic peak information obtained with [chromPeaks()] and
#' [chromPeaksData()], stored in tabular format in the text files
#' *chrom_peaks.txt* and *chrom_peak_data.txt* respectively.
#'
#' - The retention time information obtained with [adjustedRtime()] stored
#' in the text file named *rtime_adjusted.txt*.
#'
#' - The [featureDefinitions()] stored in a text file named
#' *feature_definitions.txt*. Additionally, a second file named
#' *feature_peak_index.txt* is generated to connect the features' definitions
#' with their names.
#'
#' This `param` class and method are part of the possible dispatch of the
#' generic functions `storeResults()` and `loadResults()`.
#' The folder defined in the `path` parameter will be created by calling
#' `storeResults`. If the folder already exists, previous exports in that
#' folder might get overwritten.
#'
#' @param path for `PlainTextParam` `character(1)`, defining where the files
#' are going to be stored/ should be loaded from. The default will be
#' `tempdir()`.
#'
#' @inheritParams storeResults
#'
#' @return for `PlainTextParam`: a `PlainTextParam` class. `storeResults` does
#' not return anything but saves the object to collections of different plain
#' text files to a folder.
#'
#' @author Philippine Louail, Johannes Rainer.
#'
#' @importFrom jsonlite serializeJSON write_json unserializeJSON read_json
#'
#' @examples
#' ## Load test data set of class `MsExperiment`
#' library(MsExperiment)
#' fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' pd <- data.frame(file = basename(fls),
#'                  sample = c("POOL_1", "POOL_2"),
#'                  injection_index = c(1, 19),
#'                 group = "POOL")
#' rownames(pd) <- c("1", "2")
#' mse <- readMsExperiment(fls, sampleData = pd)
#'
#' ## Define param
#' pth = file.path(tempdir(), "test")
#' param <- PlainTextParam(path = pth)
#'
#' ## Save as a collection of plain text files
#' storeResults(object = mse, param = param)
#'
#' ## Load a test data set with detected peaks, of class `XcmsExperiment`
#' faahko_sub <- loadXcmsData("faahko_sub2")
#'
#' ## Define param
#' pth = file.path(tempdir(), "test")
#' param <- PlainTextParam(path = pth)
#'
#' ## Save as a collection of plain text files
#' storeResults(object = faahko_sub, param = param)
#'
#' ## Load this saved dataset
#' faahko_load <- loadResults(param = param)
#'
NULL

#' @noRd
setClass("PlainTextParam",
         slots = c(path = "character"),
         contains = "Param",
         prototype = prototype(
             path = character()),
         validity = function(object) {
             msg <- NULL
             if (length(object@path) != 1)
                 msg <- c("'path' has to be a character string of length 1")
             msg
         })

#' @rdname PlainTextParam
#'
#' @export
PlainTextParam <- function(path = tempdir()) {
    new("PlainTextParam", path = path)
}

### methods
#' @rdname PlainTextParam
setMethod("storeResults",
          signature(object = "MsExperiment",
                    param = "PlainTextParam"),
          function(object, param){
              dir.create(path = param@path,
                         recursive = TRUE,
                         showWarnings = TRUE)
              .store_msexperiment(x = object, path = param@path)
          }
)

#' @rdname PlainTextParam
setMethod("storeResults",
          signature(object = "XcmsExperiment",
                    param = "PlainTextParam"),
          function(object, param){
              callNextMethod()
              .store_xcmsexperiment(x = object, path = param@path)
          }
)

#' @rdname PlainTextParam
setMethod("loadResults",
          signature(param = "PlainTextParam"),
          function(param){
              res <- .load_msexperiment(path = param@path)
              fl <- file.path(path, "chrom_peaks.txt")
              if (file.exists(fl))
                  res <- .load_xcmsexperiment(res, path = param@path)
              validObject(res)
              res
              }
)

#' @noRd
.store_msexperiment <- function(x, path = tempdir()) {
    .export_sample_data(as.data.frame(sampleData(x)),
                        file.path(path, "sample_data.txt"))
    .export_spectra_files(x, path = path)
    .export_spectra_processing_queue(spectra(x), path = path)
}

#' @noRd
.load_msexperiment <- function(path = character()) {
    fl <- file.path(path, "sample_data.txt")
    if (file.exists(fl))
        sd <- .import_sample_data(fl)
    else stop("No \"sample_data.txt\" file found in ", path)
    fl <- file.path(path, "spectra_files.txt")
    if (file.exists(fl)){
        sf <- .import_spectra_files(fl)
        res <- readMsExperiment(spectraFiles = sf, sampleData = sd)
        fl <- file.path(path, "spectra_processing_queue.json")
        if (file.exists(fl))
            res <- .import_processing_queue(res, fl)
    } else stop("No \"spectra_files.txt\" file found in ", path, "Spectra ",
                "data will not be restored")
}

#' Sample data
#' @noRd
.export_sample_data <- function(x, file = tempfile()) {
    write.table(x, file = file)
}

#' @noRd
.import_sample_data <- function(file = character()) {
    read.table(file)
}

#' Spectra file
#' @noRd
.export_spectra_files <- function(x, path = character()) {
    s <- spectra(x)
    if (!inherits(s@backend, "MsBackendMzR"))
        warning("Spectra data will not be exported, backend need to be of ",
                "class 'MsBackendMzR'")
    else {
        fls <- fileNames(x)
        write.table(fls, file = file.path(path, "spectra_files.txt"),
                    row.names = FALSE, col.names = FALSE)
    }
}

#' @noRd
.import_spectra_files <- function(file = character()) {
    as.character(read.table(file)[, 1L])
}

#' Processing queue
#' @param x  `Spectra` (for export) `MsExperiment` (for import)
#'
#' @noRd
.export_spectra_processing_queue <- function(x, path = character()) {
    pq <- x@processingQueue
    if (length(pq))
        write_json(serializeJSON(pq),
                   file.path(path, "spectra_processing_queue.json"))
}

#' @noRd
.import_spectra_processing_queue <- function(x, file = character()) {
    x@spectra@processingQueue <- unserializeJSON(read_json(file)[[1L]])
    x
}

#' @noRd
.store_xcmsexperiment <- function(x, path = tempdir()) {
    .export_process_history(x, path = path)
    if (hasChromPeaks(x))
        .export_chrom_peaks(x, path)
    if (hasAdjustedRtime(x))
        .export_adjusted_rtime(x, path)
    if (hasFeatures(x))
        .export_features(x, path)
}

#' @noRd
.load_xcmsexperiment <- function(x, path = character()){
    res <- as(res, "XcmsExperiment")
    res <- .import_chrom_peaks(res, path)
    fl <- file.path(path, "process_history.json")
    if (file.exists(fl))
        res <- .import_process_history(res, fl)
    else stop("No \"process_history.json\" file found in ", path)
    fl <- file.path(path, "rtime_adjusted.txt")
    if (file.exists(fl))
        res <- .import_adjusted_rtime(res, fl)
    fl <- file.path(path, "feature_definitions.txt")
    if (file.exists(fl))
        res <- .import_features(res, path)
}

#' Processing history
#' @noRd
.export_process_history <- function(x, path = character()) {
    ph <- processHistory(x)
    write_json(serializeJSON(ph), file.path(path, "process_history.json"))
}

#' @noRd
.import_process_history <- function(x, file = character()) {
    ph <- unserializeJSON(read_json(file)[[1L]])
    x@processHistory <- ph
    x
}

#' Chromatographic peaks
#' @noRd
.export_chrom_peaks <- function(x, path = character()) {
    write.table(chromPeaks(x), file = file.path(path, "chrom_peaks.txt"))
    write.table(as.data.frame(chromPeakData(x)),
                file = file.path(path, "chrom_peak_data.txt"))
}

#' @noRd
.import_chrom_peaks <- function(x, path = character()) {
    f <- file.path(path, "chrom_peaks.txt")
    pk <- as.matrix(read.table(f))
    f <- file.path(path, "chrom_peak_data.txt")
    if (!file.exists(f))
        stop("No \"chrom_peak_data.txt\" file found in ", path)
    pkd <- read.table(f)
    x@chromPeaks <- pk
    x@chromPeakData <- pkd
    x
}

#' Retention times
#' @noRd
.export_adjusted_rtime <- function(x, path = character()) {
    write.table(adjustedRtime(x), file = file.path(path, "rtime_adjusted.txt"),
                row.names = FALSE, col.names = FALSE)
}

#' @noRd
.import_adjusted_rtime <- function(x, file = character()) {
    rts <- read.table(file)[, 1L]
    x@spectra$rtime_adjusted <- as.numeric(rts)
    x
}

#' Features
#' @noRd
.export_features <- function(x, path = character()) {
    fts <- featureDefinitions(x)
    pkidx <- data.frame(
        feature_index = rep(seq_len(nrow(fts)), lengths(fts$peakidx)),
        peak_index = unlist(fts$peakidx, use.names = FALSE))
    fts$peakidx <- NA
    write.table(fts, file = file.path(path, "feature_definitions.txt"))
    write.table(pkidx, file = file.path(path, "feature_peak_index.txt"))
}

#' @noRd
.import_features <- function(x, path = character()) {
    f <- file.path(path, "feature_definitions.txt")
    fts <- read.table(f)
    f <- file.path(path, "feature_peak_index.txt")
    if (!file.exists(f))
        stop("No \"feature_peak_index.txt\" file found in ", path)
    pkidx <- read.table(f)
    fts$peakidx <- unname(split(pkidx$peak_index, pkidx$feature_index))
    x@featureDefinitions <- fts
    x
}
