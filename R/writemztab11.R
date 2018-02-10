## TODO:
## study-variable , separated -> DONE
## assay[1]-assay[2]-assay[3]-assay[4]-assay[5]-assay[6]
##
## Validation through Validator by Nils Hoffmann
## https://github.com/nilshoffmann/jmzTab-m/
##
## java -jar  /home/sneumann/src/jmzTab-m/cli/target/jmztabm-cli-1.0.0-SNAPSHOT.jar -check inFile=faahKO.mzTab -level Warn
##
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

cvTerm <- function(CV, accession, name, value="") {
    paste("[", paste(CV, accession, name, value, sep=", "), "]", sep="")    
}
#cvTerm("MS", "MS:1000443", "Mass Analyzer Type", "Orbitrap")

## Similar to rowSum
rowSd <- function(x, na.rm) apply(x, 1, sd, na.rm=na.rm)    

mzFileType <- function(paths) {
    result <- character(length(paths))
    result <- "Unknown"
    
    idx <- grepl("([Cc][Dd][Ff]$)|([Nn][Cc]$)", paths)
    result[idx] <- cvTerm("EDAM", "format:3650", "netCDF")

    idx <- grepl("([Mm][Zz])?[Xx][Mm][Ll]$", paths)
    result[idx] <- cvTerm("EDAM", "format:3654", "mzXML")

    idx <- grepl("[Mm][Zz][Dd][Aa][Tt][Aa]$", paths)
    result[idx] <- cvTerm("MS", "MS:1000564","PSI mzData format")

    idx <- grepl("[Mm][Zz][Mm][Ll]$", paths)
    result[idx] <- cvTerm("MS", "MS:1000584", "Proteomics Standards Inititative mzML file format", "mzML file")

    result
}

## paths <- c("bla.cdf", "foo.nc", "blub.mzdata", "blub.mzXML", "blub.XmL",
##            "foo/bar.bar/baz.mzML", "mzData/careful.mzML")
## cbind(paths, mzFileType(paths))


mzTabHeader <- function(mztab, version,
                        id="The ID of the mzTab file.",
                        title="The file’s human readable title.",
                        description="The file’s human readable description.",
                        sample_processing=c("[,,extraction,]", "[,,derivatisation,]"),
                        software=c("[MS, MS:1002205, ProteoWizard msconvert, 4711 ]",
                                   "[MS, MS:1001582, XCMS, 2.99.6 ]"),
                        xset,
                        value) {
    runs <- paste("file:", filepaths(xset), sep="/")
    names(runs) <- paste("ms_run[", 1:length(runs), "]-location", sep="")


    samples <- paste("ms_run[", 1:length(runs), "]", sep="")
    names(samples) <- paste("assay[", 1:length(runs), "]-ms_run_ref", sep="")
    
    sampleDesc <- sampnames(xset)
    names(sampleDesc) <- paste("sample[", 1:length(runs), "]-description", sep="")    
    
    filetypes <- mzFileType(runs)
    names(filetypes) <- paste("ms_run[", 1:length(filetypes), "]-format", sep="")

    msruns <- paste("ms_run[", seq(along=runs), "]", sep="")
    names(msruns) <- paste("assay[", seq(along=runs), "]-ms_run_ref", sep="")

    ms_runs <- sampnames(xset)
    names(ms_runs) <- paste("ms_run[", seq(along=runs), "]", sep="")

    variableAssays <- unlist(tapply(seq(along=sampclass(xset)), sampclass(xset), function(x)
                                    paste(paste("assay[",x,"]", sep=""), collapse=", ")))
    names(variableAssays) <- paste("study_variable[", seq(along=variableAssays), "]-assay_refs", sep="")
    
    variableDescriptions <- unique(as.character(sampclass(xset)))
    names(variableDescriptions) <- paste("study_variable[", seq(along=variableDescriptions), "]-description", sep="")
    
    mztab <- rbind.ragged(mztab, mzTabAddComment("Meta data section"))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD",
                                                  c("mzTab-version"=version,
                                                    "mzTab-ID"=id,
                                                    "title"=title,
                                                    "description"=description)))

    names(sample_processing) <- paste("sample_processing[", seq(along=sample_processing), "]", sep="")
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", sample_processing))

    names(software) <- paste("software[", seq(along=software), "]", sep="")
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", software))

    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", value)) # XXXvalue

    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", ms_runs))
 
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", samples))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", sampleDesc))
    
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", runs))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", msruns))

    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", variableDescriptions))
    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", variableAssays))

    mztab <- rbind.ragged(mztab, mzTabAddTagValue("MTD", filetypes))    

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

mzTabAddSML <- function(mztab, xset, value) {
    runs <- seq(along=sampnames(xset))
    variables <- seq(along=levels(sampclass(xset)))

    idHeaders <- c("SML_ID", "SMF_ID_REFS",
                   "database_identifier", "chemical_formula", "smiles", "inchi", "chemical_name", "uri",
                   "theoretical_neutral_mass", "exp_mass_to_charge",
                   "retention_time", "adduct_ions", "reliability", "best_id_confidence_measure", "best_id_confidence_value")

#, "abundance_assay[1]", "abundance_assay[2]", "abundance_assay[3]", "abundance_assay[4]", "abundance_assay[5]", "abundance_assay[6]#", "abundance_study_variable[1]", "coeffvar_study_variable[1]", "abundance_study_variable[2]", "coeffvar_study_variable[2]"
    
    ## idHeaders <- c("SML_ID", "SMF_ID_REFS", "database_identifier", "chemical_formula",
    ##                "theoretical_neutral_mass", "exp_mass_to_charge",
    ##                "retention_time", "adduct_ions", "reliability", "uri",
    ##                "best_search_engine", "best_smallmolecule_id_confidence_measure[1_n]")
                   
    abundanceAssayHeaders <- paste("abundance_assay[", runs, "]", sep="")    
    
    abundanceVariableHeaders <- unlist(lapply(variables, FUN=function(v) {
        c(paste("abundance_study_variable[", v,"]", sep=""),
          paste("abundance_coeffvar_study_variable[", v,"]", sep=""))
    }))

    headers <- c(idHeaders,
                 abundanceAssayHeaders, abundanceVariableHeaders)

    g <- groups(xset)
    SMF_ID_REFS <- sapply(xs@groupidx, function(x) paste(x, collapse="|"))
    v <- groupval(xset, value="into") ## STN: Hardcoded, needs fixing
    
    result <-  as.data.frame(matrix(character(0), ncol=length(headers), nrow=nrow(g)))
    colnames(result) <- headers

    # Calculate median/mean per study variable
    #variableAssays <- unlist(tapply(seq(along=sampclass(xset)), sampclass(xset), function(x)
    #                         paste(paste("assay[",x,"]", sep=""), collapse=",")))
    #names(variableAssays) <- paste("study_variable[", seq(along=variableAssays), "]-assay_refs", sep="")
    
    
    result[,"SML_ID"] <- seq(1,nrow(v))
    result[,"SMF_ID_REFS"] <- SMF_ID_REFS 
    result[,"retention_time"] <- g[,"rtmed"]
    result[,"exp_mass_to_charge"] <- g[,"mzmed"] 
    result[,"reliability"] <- 4 ## "unknown compound"
    result[, grepl("abundance_assay", colnames(result))] <- v

    
    cl <- sampclass(xset)    
    means <- sapply(unique(cl), function(g) rowMeans(v[,cl==g,drop=FALSE], na.rm=TRUE))
    sds <- sapply(unique(cl), function(g) rowSd(v[,cl==g,drop=FALSE], na.rm=TRUE))
    coeffs <- means / sds

    result[, grepl("abundance_study_variable", colnames(result))] <- means
    result[, grepl("abundance_coeffvar_study_variable", colnames(result))] <- coeffs
        
    mztab <- mzTabAddValues(mztab, "SMH", "SML", result)
    
}

mzTabAddSMF <- function(mztab, xset, value) {
    runs <- seq(along=sampnames(xset))
    variables <- seq(along=levels(sampclass(xset)))
    
    idHeaders <- c("SMF_ID", "SME_ID_REF_Ambiguity_code")

    featureHeaders <- c("charge", "adduct_ion", "exp_mass_to_charge",
                        "retention_time", "retention_time_start", "retention_time_end")

    abundanceAssayHeaders <- paste("quant_assay[", runs, "]", sep="")
    
    headers <- c(idHeaders,
                 featureHeaders, abundanceAssayHeaders)

    g <- groups(xset)
    
    v <- groupval(xset, value=value)
    
    result <-  as.data.frame(matrix(character(0), ncol=length(headers), nrow=nrow(g)))
    colnames(result) <- headers

    # Calculate median/mean per study variable
    #variableAssays <- unlist(tapply(seq(along=sampclass(xset)), sampclass(xset), function(x)
    #                         paste(paste("assay[",x,"]", sep=""), collapse=",")))
    #names(variableAssays) <- paste("study_variable[", seq(along=variableAssays), "]-assay_refs", sep="")
    
    
    result[,"SMF_ID"] <- seq(1,nrow(v))
    result[,"retention_time"] <- g[,"rtmed"]
    result[,"retention_time_start"] <- g[,"rtmin"]
    result[,"retention_time_end"] <- g[,"rtmax"]
    result[,"exp_mass_to_charge"] <- g[,"mzmed"]
    result[, grepl("quant_assay", colnames(result))] <- v
    
    mztab <- mzTabAddValues(mztab, "SEH", "SMF", result)    
}


writeMzTab <- function(mtd, sml, filename) {
    write.table(mtd, file=filename,
                row.names=FALSE, col.names=FALSE,
                quote=FALSE, sep="\t", na="")
    write.table(sml, file=filename, append=TRUE,
                row.names=FALSE, col.names=FALSE,
                quote=FALSE, sep="\t", na="null")

}

########################
## Example for faahKO
##

if (TRUE) {
    library(Rcpp) ## for rbind.fill
    library(plyr) ## for rbind.fill
    library(faahKO)

    if(! exists("xs")) {
        xs <- group(faahko)
    }

    values <- c("quantification_method"="[MS,MS:1002019,label-free raw feature quantitation,]",
                "small_molecule-quantification_unit"="[PRIDE, PRIDE:0000395, Ratio, ]",
                "small_molecule-identification_reliability"="[PRIDE, PRIDE:0000395, Ratio, ]" ## DUMMY!
                )

    mztmtd <- data.frame(character(0))
    mztmtd <- mzTabHeader(mztmtd,
                       version="1.0.99", description="faahKO",
                       xset=xs,
                       value=values)

    mztsml <- data.frame(character(0))
    mztsml <- mzTabAddSML(mztsml, xs, value=value) # needs to be done
#    mzt <- mzTabAddSMF(mzt, xs, value=value) # needs to be done
    ##mzt
    
    writeMzTab(mztmtd, mztsml, "faahKO.mzTab")
}

#############################
if (FALSE) {
    library(MSnbase)
    m <- readMzTabData("faahKO.mzTab")
}



       
