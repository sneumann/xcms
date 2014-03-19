setGeneric("AutoLockMass", function(object) standardGeneric("AutoLockMass"))

setMethod("AutoLockMass", "xcmsRaw", function(object) {
    if(length(grep("xml|mzData|mzXML|mzML", object@filepath, ignore.case=TRUE)) >= 1){
        tempFreq<-diff(which(diff(object@scantime) == 0))-1
        idx <- which(tempFreq != floor(mean(tempFreq))) ## only needed for newer lockmass signal
        if(is.nan(mean(tempFreq)) ){
            dn<-density(diff(object@scantime))
            lockMassScans <- quantile(dn$x, .75) ## hopefully always correct (?)
            inx<-which(diff(object@scantime) >= lockMassScans) ## these seems to be some of the new files
            return(inx)
        }else if(all(tempFreq == mean(tempFreq)) ){
            freqLock<-mean(tempFreq)
        } else if(all(idx == which(tempFreq != floor(mean(tempFreq) )) )){
            ## for the newer mzML and mzXML not sure why the change?
            ## This means that there is only one gap :( ??
            stop("This file is different from the normally seen files and requires special programming\n
                        This functionality has not been implemented yet\n ")
            ## these files seem to come either from newer MS units or/and msconvert ....
        } else {
            freqLock<-mean(tempFreq)
            warning("\nLock mass frequency wasn't detected correctly", immediate.=TRUE)
        }

        if(diff(object@scantime[1:5])[1] == 0 ){
            start<-1
        } else{
            start<-freqLock
        }
        return(makeacqNum(object, freqLock, start))

    } else if(length(grep("cdf", object@filepath, ignore.case=TRUE)) >= 1){
        ## check to see if we have the X02.CDF files around
        ## These files should be the lock mass channel
        file02<-list.files(gsub("01.CDF", "02.CDF", object@filepath), recursive=T)
        if(length(file02)> 0){
            xr<-xcmsRaw(file02)
            lockMass<-sapply(xr@scantime, function(x, object){
                which.min(abs(object@scantime - x))
            }, object)
            return(lockMass)
        } else {
            ## we couldn't find the files so lets try to find them automatically
            hr <- hist(diff(object@scantime), breaks=4, plot=FALSE)
            if(length(hr$counts) > 2){
                idx<-which(hr$counts == 0)
                ## could have something here about which way the plot is ie cor R is - or +
                                        # if(cor(hr$mids, hr$counts) < 0){
                inx<-which(diff(object@scantime) >= hr$mids[(max(idx))])
                                        # } else {
                                        #       inx<-which()
                                        # }
            }else if(length(hr$counts) == 2){
                inx<-which(diff(object@scantime) >= hr$mids[2])
            } else {
                stop("File appears to have been run without lock mass\n ")
            }
            if(length(inx) <= 1){
                warning("\nLock mass frequency wasn't detected", immediate.=TRUE)
                return(0)
            }
            ## above we're looking for scantimes that are much longer than the normal scan times
            tempFreq<-diff(inx)-1
            if(all(tempFreq == median(tempFreq)) ){
                freqLock<-median(tempFreq)
            }else{
                freqLock<-median(tempFreq)
                warning("Lock mass frequency wasn't detected correctly\n", immediate.=TRUE)
            }

            if(inx[1] == 0 || inx[1] == 1){
                start<-1
            }else{
                start<-freqLock
            }
                                        #return(inx)
            return(makeacqNum(object, freqLock, start))
        }
    } else{
        stop("Couldn't detect file type\n")
    }
})

setGeneric("makeacqNum", function(object, freq, start=1) standardGeneric("makeacqNum"))

setMethod("makeacqNum", "xcmsRaw", function(object, freq, start=1) {

    freq<-freq+1 ##nessary for the start at +1 and others since 1st scan is +1

    acqNum<-numeric()
    fo<-seq(from=start, to=length(object@scanindex), by=freq)
    for(i in fo){
        acqNum<-c(acqNum, i,i+1)
    }
    return(acqNum)
})


setGeneric("stitch", function(object, lockMass, ...) standardGeneric("stitch"))

setMethod("stitch", "xcmsRaw", function(object, lockMass) {
    if(length(grep("xml|mzData", object@filepath, ignore.case=TRUE)) >= 1){
        type<-stitch.xml
    } else if(length(grep("cdf", object@filepath, ignore.case=TRUE)) >= 1){
        ## lets check to see if lockMass is one scan or two
        if(any(diff(lockMass) == 0)){
            type<-stitch.netCDF.new
        }else {
            type<-stitch.netCDF
        }
    } else{
        stop("Unknown stitch method \n")
    }

    invisible(do.call(type, list(object, lockMass)))
})

setGeneric("stitch.xml", function(object, lockMass) standardGeneric("stitch.xml"))

setMethod("stitch.xml", "xcmsRaw", function(object, lockMass) {

    ob<-new("xcmsRaw")
    ob@env$mz<-object@env$mz
    ob@env$intensity<-object@env$intensity
    ob@scanindex<-object@scanindex
    ob@scantime<-object@scantime

    ob@acquisitionNum<-1:length(ob@scanindex)
    ob@filepath<-object@filepath
    ob@mzrange<-range(ob@env$mz)
    ob@profmethod<-object@profmethod
    ob@tic<-object@tic
    ob@profparam<-list()

    arr<-array(dim=c(2,max(diff(ob@scanindex)), length(ob@scanindex)))
    if(lockMass[1] == 1){
        lockMass<-lockMass[3:length(lockMass)]
    }
    lockMass<-matrix(lockMass, ncol=2, byrow=TRUE)
    if((lockMass[nrow(lockMass),2]+2) > length(ob@scanindex)){
        lockMass<-lockMass[1:(nrow(lockMass)-1),]
    }

    for(i in 1:(length(ob@scanindex)-1)){
        if(any(i == lockMass[,1])){
            arr[1,,i] <-c(object@env$mz[(object@scanindex[(i-1)]+1):object@scanindex[i]],
                          rep(NA, (max(diff(object@scanindex))-
                                   length((object@scanindex[(i-1)]+1):object@scanindex[i])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[(i-1)]+1):object@scanindex[i]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[(i-1)]+1):object@scanindex[i])) ))

        } else if(any(i == lockMass[,2])){
            arr[1,,i] <-c(object@env$mz[(object@scanindex[i+1]+1):object@scanindex[(i+2)]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i+1]+1):object@scanindex[(i+2)])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[i+1]+1):object@scanindex[(i+2)]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i+1]+1):object@scanindex[(i+2)])) ))

        } else{
            arr[1,,i] <-c(object@env$mz[(object@scanindex[i]+1):object@scanindex[i+1]],
                          rep(NA, (max(diff(object@scanindex))-
                                   length((object@scanindex[i]+1):object@scanindex[i+1])) ))

            arr[2,,i] <-c(object@env$intensity[(object@scanindex[i]+1):object@scanindex[i+1]],
                          rep(NA, (max(diff(object@scanindex)) -
                                   length((object@scanindex[i]+1):object@scanindex[i+1])) ))
        }
        ## mz is in 1; Intensity is in 2
        ##remake scanindex
        if(i == 1){
            ob@scanindex[i]<-as.integer(0)

        }else if(i == length(ob@scanindex)-1){
            ob@scanindex[i]<-as.integer(length(na.omit(arr[1,,(i-1)]))+ob@scanindex[(i-1)])
            ob@scanindex[i+1]<-as.integer(length(na.omit(arr[1,,i]))+ob@scanindex[i])
                                        #			ob@scanindex[i+1]<-as.integer(length(ob@env$mz))
        }else{
            ob@scanindex[i]<-as.integer(length(na.omit(arr[1,,(i-1)]))+ob@scanindex[(i-1)])
        }
    }

    NAidx<-is.na(arr[1,,])
    ob@env$mz<-as.numeric(arr[1,,][!NAidx])
    ob@env$intensity<-as.numeric(arr[2,,][!NAidx])
    ob<-remakeTIC(ob) ## remake TIC

    return(ob)
})

setGeneric("stitch.netCDF", function(object, lockMass) standardGeneric("stitch.netCDF"))

setMethod("stitch.netCDF", "xcmsRaw", function(object, lockMass) {
    if(length(lockMass) == 0 | all(lockMass == 0)){
        return(object)
    }

    ob<-new("xcmsRaw")

    ob@filepath<-object@filepath
    ob@mzrange<-range(ob@env$mz)
    ob@profmethod<-object@profmethod
    ob@profparam<-list()

    arr<-array(dim=c(2,max(diff(object@scanindex)), (length(object@scanindex)+length(lockMass)) ))
    ob@scanindex <- integer(length=length(arr[1,1,]))
    ob@acquisitionNum<-1:length(ob@scanindex)

    if(lockMass[1] == 1){
        lockMass<-lockMass[3:length(lockMass)]
    }
    lockMass<-matrix(lockMass, ncol=2, byrow=TRUE)
    if((lockMass[nrow(lockMass),2]+2) > length(ob@scanindex)){
        lockMass<-lockMass[1:(nrow(lockMass)-1),]
    } ## remove the last lock mass scan if it's at the end of the run

    add<-0
    arrMax<-length(arr[1,,1])
    scanIx<-integer(length(arr[1,1,]))
    for(i in 1:length(object@scanindex)){
                                        #		if((i+add) > length(object@scanindex)){
                                        #			break
                                        #		}
        scan<-getScan(object, i)
        arr[1,,i+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
        arr[2,,i+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
        scanIx[(i+add)+1]<- (scanIx[(i+add)])+nrow(scan)
        ##proably going to need a cut at the end of scanIx +1 problem

        if(any(i == lockMass[,1])){
            arr[1,,i+1+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
            arr[2,,i+1+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
            scanIx[(i+add)+2]<- (scanIx[1+i+add])+nrow(scan)

            scan<-getScan(object, i+1)
            arr[1,,i+2+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
            arr[2,,i+2+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
            scanIx[(i+add)+3]<- (scanIx[(i+2)+add])+nrow(scan)

            add<-add+2
        }
    }

    NAidx<-is.na(arr[1,,])
    ob@env$mz<-as.numeric(arr[1,,][!NAidx])
    ob@env$intensity<-as.numeric(arr[2,,][!NAidx])
    ##remake scanindex
                                        #	scanInx<- as.integer(apply(arr[1,,], 2, function(x){
                                        #		inx<-is.na(x)
                                        #		length(x[!inx]) ## need to add these length together
                                        #	}))
    ob@scanindex<-as.integer(scanIx)
    ob@scantime <- sapply(1:length(ob@scanindex), function(x, time){
        time*x
    }, mean(diff(object@scantime)))
    ob<-remakeTIC(ob) ## remake TIC

    return(ob)
})

setGeneric("stitch.netCDF.new", function(object, lockMass) standardGeneric("stitch.netCDF.new"))

setMethod("stitch.netCDF.new", "xcmsRaw", function(object, lockMass) {
    if(length(lockMass) == 0 | all(lockMass == 0)){
        return(object)
    }

    ob<-new("xcmsRaw")

    ob@filepath<-object@filepath
    ob@mzrange<-range(ob@env$mz)
    ob@profmethod<-object@profmethod
    ob@profparam<-list()

    arr<-array(dim=c(2,max(diff(object@scanindex)), (length(object@scanindex)+length(lockMass)) ))
    ob@scanindex <- integer(length=length(arr[1,1,]))
    ob@acquisitionNum<-1:length(ob@scanindex)

    if(lockMass[1] == 1){
        lockMass<-lockMass[2:length(lockMass)]
    }
                                        # lockMass<-matrix(lockMass, ncol=2, byrow=TRUE)

    add<-0
    arrMax<-length(arr[1,,1])
    scanIx<-integer(length(arr[1,1,]))
    for(i in 1:length(object@scanindex)){
        scan<-getScan(object, i)
        arr[1,,i+add]<- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
        arr[2,,i+add]<- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
        scanIx[(i+add)+1]<- (scanIx[(i+add)])+nrow(scan)

        if(any(i == lockMass)){
            arr[1,,i+1+add]  <- c(scan[,"mz"], rep(NA, (arrMax-length(scan[,"mz"])) ))
            arr[2,,i+1+add]  <- c(scan[,"intensity"], rep(NA, (arrMax-length(scan[,"intensity"])) ))
            scanIx[(i+add)+2]<- (scanIx[1+i+add])+nrow(scan)

            add<-add+1
            ## for the moment lets be dirty and add the scan before
            ## upgrade later to 1/2 and 1/2 from each scan
        }
    }

    NAidx<-is.na(arr[1,,])
    ob@env$mz<-as.numeric(arr[1,,][!NAidx])
    ob@env$intensity<-as.numeric(arr[2,,][!NAidx])
    ## above is to remove any NA buffers from the array

    ob@scanindex<-as.integer(scanIx)
    ob@scantime <- sapply(1:length(ob@scanindex), function(x, time){
        time*x
    }, mean(diff(object@scantime))) ## remake the scantime vector
    ob<-remakeTIC(ob) ## remake TIC

    return(ob)
})

remakeTIC<-function(object){
    for(i in 1:length(object@scanindex)){
        object@tic[i]<-sum(getScan(object, i)[,"intensity"])
    }
    return(object)
}
