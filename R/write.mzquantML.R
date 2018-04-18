
setMethod("write.mzQuantML", "xcmsSet", function(object, filename) {

    require(XML) || stop("We need library(XML) to write mzQuantML")

    mzq <- buildMzq(object)

    ## the sink() workaround seems to be needed for proper indenting.
    sink(file=filename)
    cat(saveXML(mzq)) ##, prefix = '<?xml version="1.0"?>\n'))
    sink(NULL)
})

verify.mzQuantML <- function(filename,
                             xsdfilename=system.file('xsd/mzQuantML_1_0_0.xsd', package = "xcms")) {
    xsd = xmlTreeParse(xsdfilename, isSchema =TRUE, useInternal = TRUE)
    doc = xmlInternalTreeParse(filename)
    xmlSchemaValidate(xsd, doc)
}

buildCVlist <- function() {

    CvList <- newXMLNode("CvList")

    CVs <- list(
        newXMLNode("Cv", attrs=c(
                             id="PSI-MS",
                             fullName="Proteomics Standards Initiative Mass Spectrometry Ontology",
                             version="2.29.0",
                             uri="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo")),


        newXMLNode("Cv", attrs=c(
                             id="UO",
                             fullName="Unit Ontology",
                             version="1.20",
                             uri="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo")))

    addChildren(CvList, kids=CVs)
}

buildAuditCollection <- function() {
    newXMLNode("AuditCollection",
               .children=list(
                   newXMLNode("Person",
                              attrs=c(
                                  firstName="Steffen",
                                  lastName="Neumann",
                                  id="PERS_STN"),
                              .children=list(newXMLNode("Affiliation",
                                  attrs=c(organization_ref="ORG_IPB")))),
                   newXMLNode("Organization",
                              attrs=c(
                                  name="Leibniz Institute of Plant Biochemistry",
                                  id="ORG_IPB"))
                   ))
}

buildCvParams <- function(cvparams) {
    lapply(cvparams, function(cv) {
        newXMLNode("cvParam", attrs=cv)
    })
}

buildAnalysisSummary <- function() {
    newXMLNode("AnalysisSummary",
               .children=buildCvParams(
                   list(c(accession="MS:1001834", cvRef="PSI-MS", name="LC-MS label-free quantitation analysis"),
                        c(accession="MS:1002019", cvRef="PSI-MS", name="label-free raw feature quantitation",  value="false"))
                   ))
}

buildInputFiles <- function(xs) {
    f <- filepaths(xs)
    rfgs <- lapply(1:length(f), function(i) {
        x <- f[i]
        name <- basename(x)
        newXMLNode("RawFilesGroup",
                   attrs=c(id=paste("rfg_", i, sep="")),
                   .children=newXMLNode("RawFile",
                       attrs=c(location=x, name=name, id=paste("rf_", i, sep=""))))
    })
    newXMLNode("InputFiles", .children=rfgs)
}

buildSoftwareList <- function() {
    newXMLNode("SoftwareList",
               .children=list(
                   newXMLNode("Software",
                              attrs=c(
                                  id="xcms",
                                  version=packageVersion("xcms")),
                              .children=buildCvParams(
                                  list(c(accession="MS:1001830", cvRef="PSI-MS", name="XCMS")))
                              )))
}

buildDataProcessingList <- function() {
    newXMLNode("DataProcessingList",
               .children=newXMLNode("DataProcessing",
                   attrs=c(order="1",
                       software_ref="xcms", id="DP1"),
                   .children=newXMLNode("ProcessingMethod",
                   attrs=c(order="1"))))
}

buildAssayList <- function(xs) {
    f <- filepaths(xs)
    assays <- lapply(1:length(f), function(i) {
        x <- f[i]
        name <- basename(x)
        newXMLNode("Assay",
                   attrs=c(
                       rawFilesGroup_ref=paste("rfg_", i, sep=""),
                       name=x,
                       id=paste("assay", i, sep="_")),
                   .children=list(newXMLNode("Label",
                       .children=newXMLNode("Modification",
                           buildCvParams(list(c(accession="MS:1002038", cvRef="PSI-MS", name="unlabeled sample")))
                           ))))
    })
    newXMLNode("AssayList",
               attrs=c(id="AssayList_1"),
               .children=assays)
}

buildStudyVariableList <- function(xs) {
    ## These variables are given by sampclass(xs)
    sc <- sampclass(xs)
    svlevels <- levels(sc)
    StudyVariables <- lapply (seq(1,length(svlevels)), function(i) {
        sv <- svlevels[i]
        ##        assays <- paste (sampnames(xs)[sampclass(xs) == sv], collapse="-")
        assay_ref <- paste("assay", which(sc == sv), sep="_", collapse=" ")

        newXMLNode("StudyVariable",
                   attrs=c(
                       name=svlevels[i],
                       id=paste("SV_COLLAPSED", i, sep="_")),
                   .children=list(newXMLNode("Assay_refs", assay_ref)))
    })

    ## iff phenoData has >1 columns, we also add the StudyVariables
    ## taken from the phenoData columns individually
    ## The variables are obtained via phenoData(xs)
    ## and written to mzQuantML in a column-wise fashion
    ##
    p <- cbind(phenoData(xs), "dummy")

    phenoVariables <- list()
    if (ncol(p)>1) {
        phenoVariables <- unlist(lapply(seq(1, ncol(p)), function(j) {
            sc <- p[,j]
            svlevels <- levels(sc)
            StudyVariables <- lapply (seq(1,length(svlevels)), function(i) {
                sv <- svlevels[i]
                assay_ref <- paste("assay", which(sc == sv), sep="_", collapse=" ")

                newXMLNode("StudyVariable",
                           attrs=c(
                               name=svlevels[i],
                               id=paste("SV_PHENODATA", j, i, sep="_")),
                           .children=list(newXMLNode("Assay_refs", assay_ref)))
            })
        }))
    }

    newXMLNode("StudyVariableList",
               .children=c(StudyVariables, phenoVariables))

}


buildFeatureList <- function(xs, parent) {
    sn <- sampnames(xs)
    pks <- cbind(id=paste("Feature", 1:nrow(peaks(xs)), sep="_"), peaks(xs)[, c("mz", "rt", "sample")])

    Features <- apply(pks, MARGIN=1, FUN=function(p) {
        newXMLNode("Feature", attrs=c(p[c("id", "rt", "mz")], charge="null"))
    })

    snum <- pks[,"sample"]
    unique_snum <- unique(snum)
    idxList <- lapply(unique_snum, function(x) {
        which(snum==x)
    })

    dummy <- lapply(1:length(idxList), function(i) {

        parent$addNode(newXMLNode("FeatureList",
                                  attrs=c(id=paste("FeatureList_", i, sep=""),
                                      rawFilesGroup_ref=paste("rfg_", i, sep="")),
                                  .children=Features[idxList[[i]]]))
    })
}

buildSmallMoleculeList <- function(xs) {


    feature_refs <- apply(groupval(xs, value="index"), MARGIN=1, FUN=function(x) {paste("Feature", x, sep="_", collapse=" ")})
    data <- cbind(groupnames(xs),
                  feature_refs,
                  groupval(xs, value="into"))
    naidx <- is.na(data)
    if (any(naidx)) {
        warning("xcmsSet still contains NA values, filling zero in. Better use fillPeaks() !")
        data[naidx] <- 0
    }


    SmallMolecules <- apply(data[,1:2], 1, FUN=function(x) {
        newXMLNode("SmallMolecule",
                   attrs=c(id=x[1]),
                   .children=list(newXMLNode("Feature_refs", x[2])))
    })

    DataType <- newXMLNode("DataType",
                           .children=buildCvParams(list(c(accession="MS:1001840",
                               cvRef="PSI-MS",
                               name="LC-MS feature intensity"))))

    assay_ref <- paste("assay", 1:length(sampnames(xs)), sep="_", collapse=" ")
    ColumnIndex <- newXMLNode("ColumnIndex", assay_ref)

    DataMatrix <- newXMLNode("DataMatrix",
                             .children=apply(data, MARGIN=1, FUN=function(row) {
                                 newXMLNode("Row",
                                            paste(row[c(-1,-2)], collapse=" "),
                                            attrs=c(object_ref=row[1]))
                             }))

    AssayQuantLayer <- newXMLNode("AssayQuantLayer",
                                  attrs=c(id="xset_1"),
                                  .children=list(DataType,
                                      ColumnIndex,
                                      DataMatrix))

    newXMLNode("SmallMoleculeList",
               attrs=c(id="SML_1"),
               .children=c(SmallMolecules,
                   AssayQuantLayer))
}

buildMzq <- function(xs) {
    mzqVersion="1.0.0"
    schemaLocation="http://psidev.info/psi/pi/mzQuantML/1.0.0 ../../../schema/mzQuantML_1_0_0.xsd"

    mzq = xmlTree(tag="MzQuantML",
        attrs=c(
            version=mzqVersion,
            id="mzq-generated-from-xcms",
            creationDate="2012-11-26T15:13:08.330Z",
            "xsi:schemaLocation"=schemaLocation),
        namespaces = c(
            "http://psidev.info/psi/pi/mzQuantML/1.0.0",
            xsi="http://www.w3.org/2001/XMLSchema-instance")
        )

    mzq$addNode(buildCVlist())
    mzq$addNode(buildAuditCollection())
    mzq$addNode(buildAnalysisSummary())
    mzq$addNode(buildInputFiles(xs))
    mzq$addNode(buildSoftwareList())
    mzq$addNode(buildDataProcessingList())
    mzq$addNode(buildAssayList(xs))
    mzq$addNode(buildStudyVariableList(xs))
    mzq$addNode(buildSmallMoleculeList(xs))

    buildFeatureList(xs, parent=mzq)

    mzq$closeTag()
}
