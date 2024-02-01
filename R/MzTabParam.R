#' @title Store contents of an `XcmsExperiment` objects as a .mztabm file
#'
#' @name MzTabParam
#'
#' @export
#'
#' @family xcms result export formats.
#'
#' @description
#' The `MzTabParam` class and method enable users to save an `XcmsExperiment`
#' object as a mzTabm format. Mainly the metadata (MTD) and Small molecule
#' feature (SMF) tables will represent the `XcmsExperiment`. The small molecule
#' summary section (SML) will be filled with `null` values as no annotation and
#' identification of compound is performed in `xcms`.
#'
#' This `param` class and method are part of the possible dispatch of the
#' generic function `storeResults`. The file will be created by calling
#' `storeResults`. If the file already exists, it will get overwritten.
#'
#' @param studyId `character` Will be both the `filename` of the object saved in
#' .mztabm format and the `mzTab-ID` in the file.
#'
#' @param polarity `character` Describes the polarity of the experiment. Two
#' inputs are possible, "positive" (default) or "negative".
#'
#' @param sampleDataColumn `character` strings corresponding to the column name(s)
#' of the `sampleData()` of the `XcmsExperiment` object with the different
#' *variables* of the experiment, for example it could be *"phenotype"*,
#' *"sample_type"*, etc...
#'
#' @param path `character(1)` Define where the file is going to be stored. The
#' default will be `tempdir()`.
#'
#' @param optionalFeatureColumns Optional columns from `featureDefinitions`
#' that should be exported too. For example it could be *"peakidx*,
#' *"npeaks"*, etc...
#'
#' @param dots `list` Correspond to any optional parameters to be passed
#' to the `featureValues` function. (e.g. parameters `method` or `value`).
#'
#' @inheritParams storeResults
#'
#' @return for `MzTabParam`: a `MzTabParam` class. `storeResults` does
#' not return anything but saves the object as a .mztabm file.
#'
#' @note
#' This function was build so that the output fit the recommendation of mztab-m
#' file format. These can be found here:
#' http://hupo-psi.github.io/mzTab/2_0-metabolomics-release/mzTab_format_specification_2_0-M_release.html
#'
#' @references
#'
#' Hoffmann N, Rein J, Sachsenberg T, Hartler J, Haug K, Mayer G, Alka O,
#' Dayalan S, Pearce JTM, Rocca-Serra P, Qi D, Eisenacher M, Perez-Riverol Y,
#' Vizcaíno JA, Salek RM, Neumann S, Jones AR. mzTab-M: A Data Standard for
#' Sharing Quantitative Results in Mass Spectrometry Metabolomics. Anal Chem.
#' 2019 Mar 5;91(5):3302-3310. doi: 10.1021/acs.analchem.8b04310. Epub 2019 Feb
#' 13. PMID: 30688441; PMCID: PMC6660005.
#'
#' @author Philippine Louail, Johannes Rainer
#'
#' @examples
#' ## Load a test data set with detected peaks, of class `XcmsExperiment`
#' test_xcms <- loadXcmsData()
#'
#' ## Define param
#' param <- MzTabParam(studyId = "test",
#'                     polarity = "positive",
#'                     sampleDataColumn = "sample_type")
#'
#' ## Save as a .mzTabm file
#' storeResults(test_xcms, param)
#'
NULL

#' @noRd
setClass("MzTabParam",
         slots = c(studyId = "character",
                   polarity = "character",
                   sampleDataColumn = "character",
                   path = "character",
                   optionalFeatureColumns = "character",
                   dots = "list"
                   ),
         contains = "Param",
         prototype = prototype(
             studyId = character(),
             polarity =  character(),
             sampleDataColumn = character(),
             path = character(),
             optionalFeatureColumns = character(),
             dots = list()
         ),
         validity = function(object) {
             msg <- NULL
             if(length(object@studyId) != 1)
                msg <- c("'studyId' has to be a character string of length 1")
             if(length(object@polarity) != 1)
                 msg <- c(msg, "'polarity' has to be a character string of length 1")
             if(length(object@sampleDataColumn) == 0)
                 msg <- c(msg, "'sampleDataColumn' cannot be empty")
             if (length(object@path) != 1)
                 msg <- c(msg, "'path' has to be a character string of length 1")
             msg
         })

#' @rdname MzTabParam
#'
#' @export
MzTabParam <- function(studyId = character(),
                       polarity = c("positive", "negative"),
                       sampleDataColumn = character(),
                       path = tempdir(),
                       optionalFeatureColumns = character(), ...) {
    polarity <- match.arg(polarity)
    new("MzTabParam", studyId = studyId, polarity = polarity,
        sampleDataColumn = sampleDataColumn, path = path,
        optionalFeatureColumns = optionalFeatureColumns, dots = list(...))
}

#' @rdname MzTabParam
setMethod("storeResults",
          signature(object = "XcmsExperiment",
                    param = "MzTabParam"),
          function(object, param){
              if(!param@sampleDataColumn %in% colnames(sampleData(object)))
                  stop("'sampleDataColumn' has to correspond to column names",
                  "of the sampleData() table")
              if(length(param@optionalFeatureColumns) != 0)
                  if(!param@optionalFeatureColumns %in% colnames(featureDefinitions(object)))
                      stop("'optionalFeatureColumns' have to correspond to",
                      "column names of the featureDefinitions() table")

              var_list <- unique(.mztab_study_variables(sampleData(object),
                                                        param@sampleDataColumn))
              fl <- file.path(param@path, paste0(param@studyId, ".mztab"))
              if (file.exists(fl)) {
                  warning("File ", basename(fl), " already exists. ",
                          "Will overwrite.")
                  file.remove(fl)
              }
              con <- file(fl, open = "at")
              on.exit(close(con))

              mtd <- .mztab_metadata(object, study_id = param@studyId,
                                     polarity = param@polarity,
                                     col_phenotype = param@sampleDataColumn)
              .mztab_export(mtd, con)
              writeLines("", con)

              sml <- .mztab_small_molecule_summary(n_sample = length(object),
                                                   var_list = var_list)
              .mztab_export(sml, con)
              writeLines("", con)

              smf <- do.call(.mztab_small_molecule_feature,
                             c(list(object = object,
                                    opt_columns = param@optionalFeatureColumns),
                                    param@dots))
              .mztab_export(smf, con)
          })


### Helper functions

#' @description
#'
#' Create the metadata `matrix` (MTD). Use the .MTD static object as a basis.
#'
#' @noRd
.mztab_metadata <- function(object, study_id, polarity, col_phenotype) {
    n_sample <- length(object)
    seq_sample <- seq_len(n_sample)
    base <- .MTD
    base[base == "replace_id"]  <- study_id
    msrun <- cbind(
        name = paste0("ms_run[", seq_sample, "]-location"),
        value = unlist(fileNames(object)),
        order = .prefix_zero(seq_sample)
    )
    if (polarity == "positive") {
        pol <- cbind(
            name = paste0("ms_run[", seq_sample, "]-scan_polarity[1]"),
            value = "[MS, MS:1000130, positive scan, ]",
            order = .prefix_zero(seq_sample))
    } else {
        pol <- cbind(
            name = paste0("ms_run[", seq_sample, "]-scan_polarity[1]"),
            value = "[MS, MS:1000129, negative scan, ]",
            order = .prefix_zero(seq_sample))
    }
    msrun <- rbind(msrun, pol)
    msrun <- msrun[order(msrun[, "order"]), c("name", "value")]

    assay <- cbind(
        name = paste("assay", seq_sample),
        value = paste0("assay[", seq_sample, "]"),
        order = .prefix_zero(seq_sample)
    )
    assay_ref <- cbind(
        name = paste0("assay[", seq_sample, "]-ms_run_ref"),
        value = paste0("ms_run[", seq_sample, "]"),
        order = .prefix_zero(seq_sample)
    )
    assay <- rbind(assay, assay_ref)
    assay <- assay[order(assay[, "order"]), c("name", "value")]

    var <- .mztab_study_variable_entries(sampleData(object), col_phenotype)

    mtd <- rbind(
        base[1:4, ],
        msrun,
        assay,
        var,
        base[5:nrow(base), ]
    )
    cbind(id = "MTD", mtd)
}

#' @description
#'
#' Create the *empty* small molecule summary (SML) `matrix`.
#' Use the .SML static object as a basis
#'
#' @noRd
#'
#' @examples
#'
#' .mztab_small_molecule_summary(5, c(1, 2, 3))
.mztab_small_molecule_summary <- function(n_sample, var_list) {
    sml <- c(.SML, paste0("abundance_assay[", seq_len(n_sample) ,"]"),
             paste0("abundance_study_variable[",
                    seq_len(length(var_list)) ,"]"),
             paste0("abundance_variation_study_variable[",
                    seq_len(length(var_list)) ,"]"))
    rbind(sml, c("SML", "1", rep("null", length(sml) - 2L)))
}

#' @description
#'
#' Create the small molecule feature (SMF) `matrix` One row is one feature
#' defined in xcms. Use the .SMF static object.
#' object as a basis.
#'
#' @noRd
#'
#' @param object `XcmsExperiment`.
#'
#' @param opt_columns `character` defining optional columns in
#'     `featureDefinitions` that should be exported too.
#'
#' @param ... optional parameters for the `featureValues` call.
.mztab_small_molecule_feature <- function(object, opt_columns = character(),
                                          ...) {
    fts <- featureDefinitions(object)
    smf <- matrix(data = "null", ncol = length(.SMF), nrow = nrow(fts),
                  dimnames = list(character(), .SMF))
    smf[, "SFH"] <- "SMF"
    smf[, "SMF_ID"] <- seq_len(nrow(fts))
    smf[, "exp_mass_to_charge"] <- as.character(fts$mzmed, digits = 18)
    smf[, "retention_time_in_seconds"] <- as.character(fts$rtmed, digits = 15)
    smf[, "retention_time_in_seconds_start"] <- as.character(fts$rtmin,
                                                             digits = 15)
    smf[, "retention_time_in_seconds_end"] <- as.character(fts$rtmax,
                                                           digits = 15)

    fvals <- featureValues(object, ...)
    colnames(fvals) <- paste0("abundance_assay[", seq_len(ncol(fvals)), "]")
    fvals <- apply(fvals, 2L, as.character, digits = 15)
    fvals[is.na(fvals)] <- "null"

    smf <- cbind(smf, fvals)

    if (length(opt_columns)) {
        tmp <- do.call(
            cbind, lapply(opt_columns,
                          function(z) {
                              if (z == "peakidx") {
                                  vapply(fts[, z], paste0, character(1),
                                         collapse = "| ")
                              } else as.character(fts[, z], digits = 15)
                          }))
        colnames(tmp) <- paste0("opt_", opt_columns)
        smf <- cbind(smf, tmp)
    }
    smf <- rbind(colnames(smf), smf)
    unname(smf)
}

#' @description
#'
#' Define the `character` vector with all study variables.
#'
#' @param x `data.frame` with sample annotations
#'
#' @param variable `character` with the column name(s) containing the study
#'     variables.
#'
#' @return `character` with the study variables, being a concatenation of
#'     the column name `"|"` and the value of the variable.
#'
#' @noRd
#'
#' @examples
#'
#' x <- data.frame(sex = c("male", "female", "female", "male", "male"),
#'                 group = c("case", "case", "control", "case", "control"))
#'
#' .mztab_study_variables(x, variable = c("sex", "group"))
#'
#' .mztab_study_variables(x, "sex")
.mztab_study_variables <- function(x = data.frame(), variable = character(),
                                   sep = ":") {
    do.call(cbind, lapply(variable, function(z) paste0(z, sep, x[, z])))
}

#' @description
#'
#' Create a `matrix` with the *study_variable* metadata content based.
#'
#' @param x `data.frame` representing the `sampleData` of the `MsExperiment`.
#'
#' @param variable `character` defining the columns in `x` containing the
#'     sample variables.
#'
#' @return `character` `matrix` (two columns) with the result.
#'
#' @noRd
#'
#' @examples
#'
#' .mztab_study_variable_entries(x, "sex")
#'
#' .mztab_study_variable_entries(x, c("group", "sex"))
.mztab_study_variable_entries <- function(x, variable, sep = "| ") {
    svar <- .mztab_study_variables(x, variable)
    unique_svar <- unique(as.vector(svar))
    res <- matrix(character(), ncol = 2, nrow = 0,
                  dimnames = list(character(), c("name", "value")))
    for (i in seq_along(unique_svar)) {
        idx <- which(svar == unique_svar[i], arr.ind = TRUE)
        res <- rbind(
            res,
            matrix(ncol = 2,
            c(paste0("study_variable[", i, "]"),
              paste0("study_variable[", i, "]-assay_refs"),
              paste0(unique_svar[i]),
              paste0("assay[", idx[, "row"], "]", collapse = sep)
              )))
    }
    res
}

#' @noRd
.prefix_zero <- function(x) {
    sprintf(paste0("%0", ceiling(log10(max(x) + 1)), "d"), x)
}

#' @noRd
.mztab_export <- function(x, con) {
    x <- apply(x, 1L, paste0, collapse = "\t")
    writeLines(x, con)
}

### Static objects
.MTD <- cbind(
    name = c("mzTab-version",
             "mzTab-ID",
             "software[1]",
             "quantification_method",
             "cv[1]-label",
             "cv[1]-full_name",
             "cv[1]-version",
             "cv[1]-uri",
             "database[1]",
             "database[1]-prefix",
             "database[1]-version",
             "database[1]-uri",
             "small_molecule-quantification_unit",
             "small_molecule_feature-quantification_unit",
             "id_confidence_measure[1]"),
    value = c("2.0.0-M",
              "replace_id",
              "[MS, MS:4711, xcms, 3.1.1]",
              "[MS, MS:1001834, LC-MS label-free quantitation analysis]",
              "MS",
              "PSI-MS controlled vocabulary",
              "4.1.138",
              "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo",
              "[,, \"no database\", null ]",
              "null",
              "Unknown",
              "null",
              "null",
              "[MS, MS:1001841, LC-MS feature volume, ]",
              "null"))

.SML <- c("SMH", "SML_ID","SMF_ID_REFS", "database_identifier",
          "chemical_formula", "smiles", "inchi", "chemical_name",
          "uri", "theoretical_neutral_mass", "adduct_ions", "reliability",
          "best_id_confidence_measure", "best_id_confidence_value")

.SMF <- c("SFH", "SMF_ID","SME_ID_REFS", "SME_ID_REF_ambiguity_code",
          "adduct_ion", "isotopomer", "exp_mass_to_charge", "charge",
          "retention_time_in_seconds", "retention_time_in_seconds_start",
          "retention_time_in_seconds_end")
