
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
		require(stats)
		if(file.exists("mz.txt")){
			unlink("mz.txt")
		}
		if(file.exists("intensity.txt")){
			unlink("intensity.txt")
		}

		hold<-0
		opps<-0
		scanindex<-list()

		for(i in 1:length(object@acquisitionNum)){
			percent<-(i/length(object@acquisitionNum))*100
			#cat(paste(round(percent), "%\r"))
			if( i == 1 || i == 2 ){ 
				scan<-getScan(object, i)
				cat(scan[,"mz"], " ", append=TRUE, file="mz.txt")
				cat(scan[,"intensity"], " ", append=TRUE, file="intensity.txt")
				next # since lock mass is 1 & 2 there is no gap to fill :)
			}else if(i == (hold+1)){ 
				next 
			}else if(all(lockMass != i ) || i >= (length(object@acquisitionNum) - 2) ){
				#Get all the other scans into the text file as well
				scan<-getScan(object, i)
				cat(scan[,"mz"], " ", append=TRUE, file="mz.txt")
				cat(scan[,"intensity"], " ", append=TRUE, file="intensity.txt")
				next
			}else if(any(lockMass == i)){
				hold<-i
				cat(paste(basename(object@filepath), " :\t", round(percent), "%\r"))

				scan<-getScan(object, i-1)
				scanB<-getScan(object, i+2)
				cat(scan[,"intensity"], " ", scanB[,"intensity"], " ", append=TRUE, file="intensity.txt")
				cat(scan[,"mz"], " ", scanB[,"mz"], " ",append=TRUE, file="mz.txt")
				##Fill the gap with the scan before and jsut after lock Mass scans 
			}
		}
		cat("\n")

		##Now all loops are finished we can add the mz.txt and Intensity.txt 
		##to the object and make the scanindex :)
		gc()
		if(file.exists("mz.txt") && file.exists("intensity.txt")){
			#cat(paste("Reading m/z & intensity Vectors...\n", sep=""))
			mz<-scan("mz.txt", what="numeric")
			intensity<-scan("intensity.txt", what="numeric")
		}else{
			stop("ERROR: Hard Drive error - m/z and Intensity Vectors not found\n")
		}
		mz<-as.numeric(mz)
		intensity<-as.numeric(intensity)

		scanIdx<-0
		for(k in 1:length(mz)){
			if( (k+2) == (length(mz))){
				break
			} else if(mz[k+1] > mz[k]){
				#cat(k, " <-k \n")
				if(mz[k+2] < mz[k+1] && mz[k] < mz[k+2]){
					stop("Error at ", k, " : value greater than element k\n")
				}else{
					next
				}
			}else if(mz[k+1] < mz[k]){
				scanIdx<-c(scanIdx, k)
	#			cat("N--> ",n," <--N \n")
			}
		}
		
		ob<-new("xcmsRaw")
		ob@scanindex<-as.integer(scanIdx)
		ob@env$mz<-mz
		ob@env$intensity<-intensity
		ob@acquisitionNum<-1:length(scanIdx)
		ob@filepath<-object@filepath
		ob@mzrange<-range(mz)
		ob@profmethod<-object@profmethod
		ob@tic<-object@tic
		ob@scantime<-object@scantime
		ob@profparam<-list()
		
		rm(object,scanIdx, mz, intensity)
		gc()
		cat("\n")
		return(ob)
})

