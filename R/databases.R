
.ppm<-function(Mr, Mm){ ## Mr Mz real Mm Mz Measured
    ppm<-abs((10^6)*(Mr-Mm)/Mm) ##abs for positive number
    return(ppm)
}

ppmDev<-function(Mr, ppmE=5){
    error<-(ppmE/10^6)*Mr ## just take the ppm as a percentage
    deviation<-c(Mr+error, Mr-error)## 1 is high 2 is low
    return(deviation)
}

read.mascot<-function(file, type="csv"){
    ## Experimental
    if(!file.exists(file))
        stop("File doesnt exist")
    if(type=="csv"){
        lines<-readLines(file) ##check to see if file is ~ normal
        if (grep("Protein hit", lines)) {
            x<-lines[-seq(grep("Protein hit", lines)+1)]
            myDF <-read.csv(textConnection(x), header=T)
        } else {
            stop("Noncompatable csv file\n")
        }
    } else if (type=="xml"){
        stop("XML files are currently not supported, support comming soon")
    }
    return(myDF)
}

# KeggSearch <- function(object, DBsearchMS) {
#
#     libname <- 'KEGGSOAP'
#     KEGG.status <- try(require(libname, character.only = TRUE, quietly=TRUE))
#
#     if (class(KEGG.status) == "try-error")
#         stop("Couldn't load KEGGSOAP\n")
#
#     if(class(object)=="xcmsSet"){
#         groupmat<-groups(object)
#         neutralmass <- groupmat[,"mzmed"] + ifelse(DBsearchMS < 0, 1, -1)
#
#         KEGGcmpd<-array(0, dim=length(neutralmass))
#         for(i in 1:length(neutralmass)){
#             KEGGcmpd[i]<-search.compounds.by.mass(neutralmass, range(ppmDev(neutralmass[i], DBsearchMS)))[1]
#         }
#         return(KEGGcmpd)
#     }
# }
