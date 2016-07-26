## Methods for xcmsFragments.
#' @include functions-xcmsFragments.R

############################################################
## show
setMethod("show", "xcmsFragments", .xcmsFragments.show)

############################################################
## collect
setMethod("collect", "xcmsFragments", .xcmsFragments.collect)

############################################################
## plotTree
setMethod("plotTree", "xcmsFragments", .xcmsFragments.plotTree)

############################################################
## hasMSn
setMethod("hasMSn", "xcmsFragments", .xcmsFragments.hasMSn)

############################################################
## findneutral
setMethod("findneutral", "xcmsFragments", function(object, find, ppmE=25, print=TRUE) {
    find<-range(ppmDev(Mr=find, ppmE))
    spectra<-unique(object@peaks[,"MSnParentPeakID"])
    found<-matrix(ncol=10)

    for (i in 1:length(spectra)){
        if(spectra[i] > 0){
            losses<-diff(sort(object@peaks[object@peaks[,"MSnParentPeakID"] ==spectra[i],"mz"] ))
            if(length(losses) > 0){
                if(length(which(losses > find[1] & losses < find[2])) > 0){
                    idx<-which(object@peaks[,"MSnParentPeakID"] == spectra[i])
                    PrecursorMZ<-object@peaks[object@peaks[idx,"MSnParentPeakID"],"mz"]
                    CE<-object@peaks[object@peaks[idx,"MSnParentPeakID"],"CollisionEnergy"]
                    dat<-object@peaks[idx,c("MSnParentPeakID", "msLevel", "rt", "mz",
                                            "intensity", "Sample", "GroupPeakMSn")]
                    found<-rbind(found, cbind(NeutralLoss=c(losses,0),  PrecursorMz=PrecursorMZ,
                                              dat[order(dat[,"mz"]),], CollisionEnergy=CE))
                }
            }
        }
    }
    if(nrow(found) >1){
        found<-found[2:nrow(found),]
    } else{
        cat("Nothing found\n")
        return(0)
    }
    if (print == TRUE){
        cat("We looked for", find[2], " to", find[1], "and found:\n")
        print(found)
    }
    return(found)
})

############################################################
## findMZ
setMethod("findMZ", "xcmsFragments", function(object, find, ppmE=25, print=TRUE) {
    find<-range(ppmDev(Mr=find, ppmE))
    fragMZ<-0

    found<- which(object@peaks[,"mz"] > find[1] & object@peaks[,"mz"] < find[2] & object@peaks[,"msLevel"] >1)
    if(length(found) <1){
        cat("nothing was found\n")
        return(0)
    }

    PrecursorMZ<-object@peaks[object@peaks[found,"MSnParentPeakID"],"mz"]
    CE<-object@peaks[object@peaks[found,"MSnParentPeakID"],"CollisionEnergy"]
    foundFrag<-cbind(PrecursorMz=PrecursorMZ, object@peaks[found,c("MSnParentPeakID", "msLevel",
                     "rt", "mz", "intensity", "Sample", "GroupPeakMSn")], CollisionEnergy=CE)

    if (print == TRUE){
        cat("We looked for", find[2], " to", find[1], "and found:\n")
        print(foundFrag)
    }
    return(foundFrag)
})

