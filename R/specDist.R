######################################### wrapper

setGeneric("specDist", function(object, ...) standardGeneric("specDist"))
setMethod("specDist", signature(object="xcmsSet"), function(object,peakIDs1,peakIDs2,  method=getOption("BioC")$xcms$specDist.method, mzabs=0.001, mzppm=10,symmetric=FALSE,...) {
   # stop("!!! xsSet")
    a <- peakIDs1
    a <- peakIDs2 ## causes a stop if one of the variables is not defined (wie schwachsinnig ist das denn???? )
    peaks <- object@peaks
    peakTable1 <- peaks[peakIDs1,c("mz","into")]
    peakTable2 <- peaks[peakIDs2,c("mz","into")]

    method <- match.arg(method, getOption("BioC")$xcms$specDist.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("specDist", method, sep=".")
    distance <- do.call(method, alist<-list(peakTable1, peakTable2, ...))
    distance
    })

setMethod("specDist", signature(object="xsAnnotate"), function(object,object2=NA, peakIDs1=NA,peakIDs2=NA,PSpec1=NA,PSpec2=NA, method=getOption("BioC")$xcms$specDist.method, 
                                                               mzabs=0.001, mzppm=10,symmetric=FALSE,...) {
    if (is.na(peakIDs1[1])&is.na(peakIDs2[1])&is.na(PSpec1)) stop("No index given. You must define at least either two vectors of peak-IDs (peakIDs1, peakIDs2) or the indicies for two Pseudospectra (PSpec1, PSpec2 )!")
    if (class(object2)!="xsAnnotate") {
        object2<-object
        }
    peaks1 <- object@peaks
    peaks2 <- object2@peaks
    method <- match.arg(method, getOption("BioC")$xcms$specDist.methods)
        if (is.na(method))
            stop("unknown method : ", method)
        method <- paste("specDist", method, sep=".")

    if (!is.na(peakIDs1[1]) & !is.na(peakIDs2[1])){
        peakTable1 <- peaks1[peakIDs1,c("mz","into")]
        peakTable2 <- peaks2[peakIDs2,c("mz","into")]
        distance <- do.call(method, alist<-list(peakTable1, peakTable2, ...))
        return(distance)
    }else{
        if (!is.na(PSpec1)&!is.na(PSpec2)) {
            if (length(object@pspectra[[PSpec1]])>1) { peakTable1 <- cbind(object@pspectra[[PSpec1]], object@peaks[object@pspectra[[PSpec1]],])
                }else{ peakTable1 <- cbind(object@pspectra[[PSpec1]], t(as.matrix(object@peaks[object@pspectra[[PSpec1]],])))}
            colnames(peakTable1) <- c("oPeak",colnames(object@peaks))

            if (length(object2@pspectra[[PSpec2]])>1) { peakTable2 <- cbind(object2@pspectra[[PSpec2]], object2@peaks[object2@pspectra[[PSpec2]],])
                }else{ peakTable2 <- cbind(object2@pspectra[[PSpec2]], t(as.matrix(object2@peaks[object2@pspectra[[PSpec2]],])))}
            colnames(peakTable2) <- c("oPeak",colnames(object2@peaks))    
            distance <- do.call(method, alist<-list(peakTable1, peakTable2, mzabs=mzabs, mzppm=mzppm,...))
            return(distance)
            }
        }
    })

### specDist - functions
setGeneric("specDist.meanMZmatch", function(peakTable1, peakTable2, ...) standardGeneric("specDist.meanMZmatch"))
setMethod("specDist.meanMZmatch", signature(peakTable1="matrix", peakTable2="matrix"), function(peakTable1, peakTable2, matchdist=1, matchrate=1, mzabs=0.001, mzppm=10,symmetric=TRUE,...){
       # stop("test")
        mz1 <- peakTable1[,"mz"]
        mz2 <- peakTable2[,"mz"]
        len1<-length(mz1)
        len2<-length(mz2)
        dn<-NULL
        idx<-minifm(peakTable1[,"mz"],peakTable2[,"mz"],mzabs=mzabs,mzppm=mzppm, symmetric=symmetric)

        for (a in 1:length(idx[[1]])) 
            dn[a] <- abs(mz1[idx[[1]][a]]-mz2[idx[[2]][a]])
       
        if (length(dn)>0){
           #cat(length(dn))
            return( (sum(dn) / length(dn)*matchdist) + ((len1+len2)/(length(dn)*2))*matchrate )
            }else{
            return(NA)
            }
        }
    )

minifm <- function(v1, v2, mzabs,mzppm,symmetric)
    {
    tol=mzabs + (max(c(v1,v2))/1000000)*mzppm
    m<-fastMatch(v1,v2,tol=tol,symmetric=symmetric);
    idx1<-NA
    idx2<-NA
    idpos<-1
    for (a in 1:length(m)){
        if (is.numeric(m[[a]])) {
            mns <- abs(v2[m[[a]]]-v1[a])
            if (min(mns)<(mzabs + (v1[a]/1000000)*mzppm)) {
                idx1[idpos]<-a
                idx2[idpos]<-m[[a]][which(mns == min(mns))]
                idpos <- idpos + 1
                }
            }
        }
    list(idx1=idx1, idx2=idx2)
    }

setGeneric("specDist.cosinus", function(peakTable1, peakTable2, ...) standardGeneric("specDist.cosinus"))
setMethod("specDist.cosinus", signature(peakTable1="matrix", peakTable2="matrix"), function(peakTable1, peakTable2, mzExp=0.6, intExp=3, nPdiff=2, 
                                                                                            nPmin=8, mzabs=0.001, mzppm=10,symmetric=FALSE,...){
    len1 <- nrow(peakTable1)
    len2 <- nrow(peakTable2)
    if ((len1<len2 & len2/len1 > nPdiff)|(len1>len2 & len1/len2 > nPdiff)) {
        return(NA)
        }else{
        idx<-minifm(peakTable1[,"mz"],peakTable2[,"mz"],mzabs=mzabs,mzppm=mzppm, symmetric=symmetric)
        int1<-NA
        int2<-NA
        if (!is.na(idx$idx1[1]) ) {
            for (a in 1:length(idx$idx1)){
                int1[a] <- peakTable1[idx$idx1[a],"into"]^intExp * peakTable1[idx$idx1[a],"mz"]^mzExp
                int2[a] <- peakTable2[idx$idx2[a],"into"]^intExp * peakTable2[idx$idx2[a],"mz"]^mzExp
                }
    
        if (length(int1)>nPmin){
                N <- length(int1)
                Nu <- len2-N
                Nl <- len1-N
                Fd <- ( sum(int1 * int2) ^ 2 ) / ( sum(int1^2) * sum(int2^2))
                summe <- 1
                if (N>1){
                    summe <- 0
                    for (a in 2:N) {
                        kl <- ( (int1[a] / int1[a-1]) * (int2[a-1] / int2[a]) )
                        if (kl>=1) {
                            h <- 1 
                            }else {
                            h <- -1
                            }
                        summe <- summe + kl^h
                        }
                    }
                Fr <- 1/N * summe
                dist <- ((Nu + Nl)* Fd + N * Fr) / (Nu + N + Nl)
                return(dist)
                }else{
                return(NA)
                }
            }else{
            return(NA)
            }
        }
     })


setGeneric("specDist.peakCount", function(peakTable1, peakTable2, ...) standardGeneric("specDist.peakCount"))
setMethod("specDist.peakCount", signature(peakTable1="matrix", peakTable2="matrix"), function(peakTable1, peakTable2, mzabs=0.001, mzppm=10,symmetric=FALSE,...){
    idx<-minifm(peakTable1[,"mz"],peakTable2[,"mz"],mzabs=mzabs,mzppm=mzppm,symmetric)
    return(length(idx$idx1))
        })
