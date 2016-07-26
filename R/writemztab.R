## TODO:
## study-variable , separated -> DONE
## assay[1]-assay[2]-assay[3]-assay[4]-assay[5]-assay[6]

## link sample - assay:
##
## MTD	sample[1]-description	KO15
## MTD	assay[1]-sample_ref	sample[1]

## utility function, combining different length objects into a dataframe
## padding short columns with NA
rbind.ragged <- function(x, y) {
    x <- as.data.frame(x) 
    y <- as.data.frame(y) 
    colnames(x) <- seq(1:ncol(x))
    colnames(y) <- seq(1:ncol(y))
    rbind.fill(x,y)
}

cvTerm <- function(CV, accession, name, value) {
    paste("[", paste(CV, accession, name, value, sep=", "), "]", sep="")    
}
#cvTerm("MS", "MS:1000443", "Mass Analyzer Type", "Orbitrap")

mzFileType <- function(paths) {
    result <- character(length(paths))

    idx <- grepl("([Cc][Dd][Ff]$)|([Nn][Cc]$)", paths)
    result[idx] <- "netCDF"

    idx <- grepl("([Mm][Zz])?[Xx][Mm][Ll]$", paths)
    result[idx] <- "mzXML"

    idx <- grepl("[Mm][Zz][Dd][Aa][Tt][Aa]$", paths)
    result[idx] <- "mzData"

    idx <- grepl("[Mm][Zz][Mm][Ll]$", paths)
    result[idx] <- cvTerm("MS", "MS:1000584", "Proteomics Standards Inititative mzML file format", "mzML file")

    result
}

## paths <- c("bla.cdf", "foo.nc", "blub.mzdata", "blub.mzXML", "blub.XmL",
##            "foo/bar.bar/baz.mzML", "mzData/careful.mzML")
## cbind(paths, mzFileType(paths))


mzTabHeader <- function(mztab, version, mode, type, description, xset) {
    runs <- filepaths(xset)
    names(runs) <- paste("ms_run[", 1:length(runs), "]-location", sep="")


    samples <- paste("sample[", 1:length(runs), "]", sep="")
    names(samples) <- paste("assay[", 1:length(runs), "]-sample_ref", sep="")
    
    sampleDesc <- sampnames(xset)
    names(sampleDesc) <- paste("sample[", 1:length(runs), "]-description", sep="")    
    filetypes <- mzFileType(runs)
    names(runs) <- paste("ms_run[", 1:length(filetypes), "]-format", sep="")

    assays <- paste("ms_run[", seq(along=runs), "]", sep="")
    names(assays) <- paste("assay[", seq(along=runs), "]-ms_run_ref", sep="")

    variableAssays <- unlist(tapply(seq(along=sampclass(xset)), sampclass(xset), function(x)
                                    paste(paste("assay[",x,"]", sep=""), collapse=",")))
    names(variableAssays) <- paste("study_variable[", seq(along=variableAssays), "]-assay_refs", sep="")
    
    variableDescriptions <- unique(as.character(sampclass(xset)))
    names(variableDescriptions) <- paste("study_variable[", seq(along=variableDescriptions), "]-description", sep="")
    
    mztab <- rbind.ragged(mztab, mzTabAddComment("Meta data section"))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD",
                                                  c("mzTab-version"=version,
                                                    "mzTab-mode"=mode,
                                                    "mzTab-type"=type,
                                                    "description"=description)))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", runs))

    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", samples))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", sampleDesc))
    
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", assays))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", variableAssays))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", variableDescriptions))
}

mzTabAddComment <- function(comments) {
    cbind.data.frame("COM", comments, stringsAsFactors=FALSE)
}

mzTabAddTagValue <- function(section, values) {
    cbind.data.frame(section, names(values), values, stringsAsFactors=FALSE)
}

mzTabAddValues <- function(mztab, headers, section, values) {
    h <- cbind.data.frame(headers, t(names(values)), stringsAsFactors=FALSE)
    v <- cbind.data.frame(section, values, stringsAsFactors=FALSE)
    
    mztab <- rbind.ragged(mztab, h)
    mztab <- rbind.ragged(mztab, v)
}

mzTabAddSME <- function(mztab, xset) {
    runs <- seq(along=sampnames(xset))
    variables <- seq(along=levels(sampclass(xset)))
    
    idHeaders <- c("identifier", "description", "chemical_formula",
                   "smiles", "inchi_key", "database", "database_version")

    searchHeaders1 <- c("search_engine", "best_search_engine_score")
   
    searchHeaders2 <- paste("search_engine_score_ms_run[", runs, "]", sep="")
    
    searchHeaders3 <- c("reliability", "modifications")

    featureHeaders <- c("charge", "adduct_ion", "exp_mass_to_charge",
                        "calc_mass_to_charge", "calc_neutral_mass", "retention_time",
                        "retention_time_window", "uri", "spectra_ref")

    abundanceAssayHeaders <- paste("smallmolecule_abundance_assay[", runs, "]", sep="")
    
    
    abundanceVariableHeaders <- unlist(lapply(variables, FUN=function(v) c(paste("smallmolecule_abundance_study_variable[", v,"]", sep=""),
                                  paste("smallmolecule_abundance_stddev_study_variable[", v,"]", sep=""),
                                  paste("smallmolecule_abundance_std_error_study_variable[", v,"]", sep=""))))

    optHeaders <- "opt_global_cv_isotopic_mass_trace"

    headers <- c(idHeaders,
                 searchHeaders1, searchHeaders2, searchHeaders3,
                 featureHeaders, abundanceAssayHeaders, abundanceVariableHeaders, optHeaders)

    g <- groups(xset)
    v <- groupval(xset, value="into")
    
    result <-  as.data.frame(matrix(character(0), ncol=length(headers), nrow=nrow(g)))
    colnames(result) <- headers

    # Calculate median/mean per study variable
    #variableAssays <- unlist(tapply(seq(along=sampclass(xset)), sampclass(xset), function(x)
    #                         paste(paste("assay[",x,"]", sep=""), collapse=",")))
    #names(variableAssays) <- paste("study_variable[", seq(along=variableAssays), "]-assay_refs", sep="")
    
    
    result[,"retention_time"] <- g[,"rtmed"]
    result[,"exp_mass_to_charge"] <- g[,"mzmed"]
    result[, grepl("smallmolecule_abundance_assay", colnames(result))] <- v
    
    mztab <- mzTabAddValues(mztab, "SEH", "SME", result)
    
}

writeMzTab <- function(object, filename) {
    write.table(object, file=filename,
                row.names=FALSE, col.names=FALSE,
                quote=TRUE, sep="\t", na="\"\"")

}

########################
## Example for faahKO
##

if (FALSE) {
    library(Rcpp) ## for rbind.fill
    library(plyr) ## for rbind.fill
    library(faahKO)

    if(! exists("xs")) {
        xs <- group(faahko)
    }

    mzt <- data.frame(character(0))
    mzt <- xcms:::mzTabHeader(mzt,
                       version="1.1.0", mode="Complete", type="Quantification",
                       description="faahKO",
                       xset=xs)
    mzt <- xcms:::mzTabAddSME(mzt, xs)
    ##mzt
    
    xcms:::writeMzTab(mzt, "faahKO.mzTab")
}

#############################
if (FALSE) {
    library(MSnbase)
    m <- readMzTabData("faahKO.mzTab")
}



       
