#' @include XcmsExperiment.R
#' @title Store contents of `MsExperiment` and `XcmsExperiment` objects as
#' plain text files
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
#' For a `Spectra` object, the exported files include:
#'
#' - The [spectraData()] stored in a tabular format in a text file named
#' *backend_data.txt*.
#'
#' - The `processingQueueVariables`, `processing`, [processingChunkSize()] and
#' `backend` class information of the object stored in a text file named
#' *spectra_slots.txt*.
#'
#' - Processing queue of the `Spectra` object, ensuring that any spectra data
#' modifications are retained. It is stored in a `json` file named
#' *spectra_processing_queue.json*.
#'
#' Note : The Spectra object will only be exported if it uses a
#' [MsBackendMzR()] backend. Other backends are no supported as of now.
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
#' @param spectraPath for `loadResults` `character(1)`, defining the
#' absolute path where the spectra files should be imported from when loading
#' the object. The default will be set using the common file path of all the
#' spectra files when exporting.
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
#' @importFrom utils read.table write.table
#'
#' @importFrom MsExperiment MsExperiment readMsExperiment
#'
#' @importFrom MsCoreUtils common_path
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
#' faahko_load <- loadResults(object = XcmsExperiment(), param = param)
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
              write.table(as.data.frame(sampleData(object)),
                          file = file.path(param@path,
                                           "sample_data.txt"))
              ## call export of individual other objects (not MsExperiment data)
              storeResults(spectra(object), param)
              ## at some point also chromatograms, etc.
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
          signature(object = "MsExperiment",
                    param = "PlainTextParam"),
          function(object, param, spectraPath = character()){
              fl <- file.path(param@path, "sample_data.txt")
              if (file.exists(fl)){ # should i have a error if does not exist ?
                  sd <- read.table(fl)
                  rownames(sd) <- NULL #read.table force numbering of rownames
              }
              s <- loadResults(Spectra(), param)
              res <- MsExperiment(sampleData = sd, spectra = s)
              validObject(res)
              res
              }
)

#' @rdname PlainTextParam
setMethod("loadResults",
          signature(object = "XcmsExperiment",
                    param = "PlainTextParam"),
          function(object, param, spectraPath){
              res <- callNextMethod() #check if need to add spectraPath = spectraPath
              res <- .load_xcmsexperiment(res, path = param@path)
              validObject(res)
              res
          }
)

#' @rdname PlainTextParam
setMethod("storeResults", signature(object = "Spectra",
                                    param = "PlainTextParam"),
          function(object, param) {
              dir.create(path = param@path,
                         recursive = TRUE,
                         showWarnings = TRUE)
              if (!existsMethod("storeResults", c(class(object@backend)[1L],
                                                  "PlainTextParam")))
                  stop("Can not store a 'Spectra' object with backend '",
                       class(object@backend)[1L], "'")
              storeResults(object@backend, param = param)
              .export_spectra_processing_queue(object, path = param@path)
              .export_spectra_slots(object, path = param@path)
          })

#' @rdname PlainTextParam
setMethod("loadResults", signature(object = "Spectra",
                                    param = "PlainTextParam"),
          function(object, param, spectraPath = character()) {
              ## here i am NOT making a separate function for the slots
              fl  <- file.path(param@path, "spectra_slots.txt")
              if (!file.exists(fl))
                  stop("No 'spectra_slots.txt' file found in ", param@path)
              fls  <- readLines(fl)
              var_names <- sub(" =.*", "", fls)
              var_values <- sub(".* = ", "", fls)
              variables <- setNames(var_values, var_names)
              if (!existsMethod("loadResults", c(variables[["backend"]],
                                                  "PlainTextParam")))
                  stop("Can not store a 'Spectra' object with backend '",
                       variables["backend"], "'")
              b <- loadResults(object= do.call(what = variables[["backend"]],
                                               args = list()),
                          param = param, spectraPath = spectraPath) ##better way to do this ?
              s <- Spectra(b)
              s@processingQueueVariables <- unlist(strsplit(variables[["processingQueueVariables"]],
                                                     "|", fixed = TRUE))
              s@processing <- unlist(strsplit(variables[["processing"]], "|" ,
                                              fixed = TRUE))
              s@processingChunkSize <- as.numeric(variables[["processingChunkSize"]])
              fl <- file.path(param@path, "spectra_processing_queue.json")
              if (file.exists(fl))
                s <- .import_spectra_processing_queue(s, file = fl)
              s
          })



# Notes: This and the Spectra method will be moved to it's respective package
#' @rdname PlainTextParam
setMethod("storeResults", signature(object = "MsBackendMzR",
                                    param = "PlainTextParam"),
          function(object, param) {
              dir.create(path = param@path,
                         recursive = TRUE,
                         showWarnings = TRUE)
              object <-  dropNaSpectraVariables(object)
              write.table(object@spectraData,
                          file = file.path(param@path, "backend_data.txt"),
                          sep = "\t", quote = FALSE)
})

#' @rdname PlainTextParam
setMethod("loadResults", signature(object = "MsBackendMzR",
                                    param = "PlainTextParam"),
          function(object, param, spectraPath = character()) {
              b <- MsBackendMzR()
              data <- read.table(file = file.path(param@path,
                                                  "backend_data.txt"),
                                 sep = "\t", header = TRUE)
              rownames(data) <- NULL
              data <- DataFrame(data)
              b@spectraData <- data
              if (length(spectraPath) > 0){
                  old <- MsCoreUtils::common_path(dataStorage(b))
                  if (nchar(old) > 0)
                      old <- paste0(old, "/")
                  dataStorage(b) <- sub(old, spectraPath, dataStorage(b))
              }
              b
})

#' Spectra slots
#' @param x  `Spectra`
#'
#' @noRd
.export_spectra_slots <-function(x, path = character()){
    con <- file(file.path(path, "spectra_slots.txt"), open = "wt")
    on.exit(close(con))
    pq <- x@processingQueueVariables
    writeLines(paste0("processingQueueVariables = ", paste(pq, collapse = "|")),
               con = con)
    p <- x@processing
    writeLines(paste0("processing = ", paste(p,collapse = "|")), con = con)
    writeLines(paste0("processingChunkSize = ", processingChunkSize(x)),
               con = con)
    writeLines(paste0("backend = ", class(x@backend)[1L]), con = con)
}


#' Processing queue
#' @param x  `Spectra`
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
    x@processingQueue <- unserializeJSON(read_json(file)[[1L]])
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
.load_xcmsexperiment <- function(x, path = character(),
                                 spectraExport = logical()){
    x <- as(x, "XcmsExperiment")
    fl <- file.path(path, "chrom_peaks.txt")
    x <- .import_chrom_peaks(x, path)
    fl <- file.path(path, "process_history.json")
    if (file.exists(fl))
        x <- .import_process_history(x, fl)
    else stop("No \"process_history.json\" file found in ", path)
    fl <- file.path(path, "rtime_adjusted.txt")
    if (file.exists(fl))
        x <- .import_adjusted_rtime(x, fl)
    fl <- file.path(path, "feature_definitions.txt")
    if (file.exists(fl))
        x <- .import_features(x, path)
    x
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
