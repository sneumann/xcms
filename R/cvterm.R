#' S4 class to store, access and link the PSI-MS style cvTerms
#'
#' This class can hold OBO terms, and in particular is used to set and retrieve
#' strings in the form [STATO,STATO:0000169,fold change,,] as used in e.g. the mzTab format.
#' See here for an example: https://github.com/HUPO-PSI/mzTab/blob/master/examples/labelfree_SQI.mzTab
#'
#' \describe{
#'    \item{ontology}{ID of the ontology that contains the term.}
#'
#'    \item{accession}{Accession of the term including Ontology prefix.}
#'
#'    \item{name}{Human-readable label of the ontology term.}
#' 
#'    \item{value}{An optional free-text value of the term.}
#'  }
#' @name cvTern-class
#' @rdname cvTerm-class
#### ' @exportClass cvTerm # NOTYETEXPORTED

.cvTerm <- setClass("cvTerm",
                   representation=list(ontology="character",
                              accession="character",
                              name="character",
                              value="character"),
                   prototype(ontology="",accession="",name="",value="")
                   )

cvTerm <- function(ontology="", accession="", name="", value="") 
    .cvTerm(ontology=ontology, accession=accession, name=name, value=value) 


setOldClass("cvTerm", S4Class = "cvTerm")

as.character.cvTerm <- function(object, ...) {
    paste("[",
          object@ontology, ",",
          object@accession, ",",
          object@name, ",",
          object@value, ",",
          "]", sep="")}

setAs("cvTerm", "character", function(from) as.character.cvTerm(from))

setMethod("show",
          "cvTerm",
          function(object) {
              cat(as(object, "character"), "\n")})

if (FALSE) {
    
    cv <- cvTerm(ontology="ontology",
                 accession="cv ref",
                 name="name",
                 value="value")
    
    isS4(cv)
    show(cv)

    cvTerm(ontology="STATO", accession="STATO_0000169", name="fold change", value="")

    cvTerm(ontology="STATO", accession="STATO_0000169", name="fold change")

    cvTerm("STATO", "STATO_0000169", "fold change")
    

    dr <- data.frame(name=c("A", "B", "C"),
                     tvalue=c(3,6,9),
                     pvalue=c(0.1,0.2,0.3),
                     ko16=c(11,22,33))

    attr(dr$pvalue, "cvTerms") <- list(cvTerm("STATO", "STATO:0000169", "fold change"))
    attr(dr$ko16, "cvTerms") <- list(cvTerm("MS", "MS:1001841", "LC-MS feature volume"))
    

}


