setGeneric("AutoLockMass", function(object) standardGeneric("AutoLockMass"))

setMethod("AutoLockMass", "xcmsRaw", function(object) {
	tempFreq<-diff(which(diff(object@scantime) == 0))-1
	if(all(tempFreq == mean(tempFreq)) ){
		freqLock<-mean(tempFreq)
	}else{
		freqLock<-mean(tempFreq)
		warning("Lock mass frequency wasn't detected correctly\n", immediate.=TRUE)
	}

	if(diff(object@scantime[1:5])[1] == 0 ){
		start<-1
	}else{
		start<-freqLock
	}
	return(makeacqNum(object, freqLock, start))
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


setGeneric("stitch", function(object, lockMass) standardGeneric("stitch"))

setMethod("stitch", "xcmsRaw", function(object, lockMass) {
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
				rep(NA, (max(diff(object@scanindex))-length((object@scanindex[(i-1)]+1):object@scanindex[i])) ))

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
				rep(NA, (max(diff(object@scanindex))-length((object@scanindex[i]+1):object@scanindex[i+1])) ))

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

	return(ob)
})


