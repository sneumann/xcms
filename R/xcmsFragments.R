require(methods) || stop("Couldnt load package methods")

setClass("xcmsFragments", representation(peaks = "matrix",
                                   MS2spec = "list",
                                   specinfo = "matrix"
                                        #, pipeline = "xcmsRawPipeline"
                                         ),
         prototype(peaks = matrix(nrow = 0, ncol = 6),
                          MS2spec=NULL,
                          specinfo=NULL
                                        #, pipeline = new("xcmsRawPipeline")
                   ))

xcmsFragments <- function(xs = NULL, ...) {
    object <- new("xcmsFragments")
    object <- collect(object,xs,...)
    object
}

.xcmsFragments.show <- function(object) {
        if(object@peaks){
            cat("An \"xcmsFragments\" object with ",nrow(object@peaks),
                " peaks in",length(unique(object@peaks[,"rt"])),"Spectra\n")
            cat("From Level 1 to",max(object@peaks[,"msLevel"]),
                "Number of Samples: ",length(unique(object@peaks[,"Sample"])),".\n")

            for (a in unique(object@peaks[,"Sample"])){
                cat("\nSample",a,":\n")
                for (b in unique(object@peaks[which(object@peaks[,"Sample"]==a),"msLevel"])) {
                    cat("   ",
                        length(which(object@peaks[which(object@peaks[,"Sample"]==a),"msLevel"] == b)),
                            "Peaks in Level",b,"\n")
                }
            }
       } else {
            cat("An \"xcmsFragments\" object with", length(object@specinfo[,"ref"]), " collected mass spectra\n\n")
            cat("Mass range:", paste(round(range(object@specinfo[,"preMZ"]), 4), collapse = "-"), "m/z\n")
            paste(cat("Collision Energy range:", paste(range(object@specinfo[,"CollisionEnergy"])), "V"), sep="")

        }
        memsize <- object.size(object)
        cat("\nMemory usage:", signif(memsize/2^20, 3), "MB\n")
    }
setMethod("show", "xcmsFragments", .xcmsFragments.show)


.xcmsFragments.collect <- function(object,xs,xraw=NULL,compMethod="floor", snthresh=20, mzgap=.2, uniq=TRUE) {
    ## This is called after findPeaks()
    ## e.g. during xcmsFragments() constructor
    ## This creates the association between the xcmsSet peaks
    ## and the MSn scans in the raw files.

    if (class(xs)=="xcmsSet") {
        ms1peaks <- peaks(xs)
    } else if (class(xs)=="matrix") {
        ms1peaks <- xs
    } else {
        stop("xs is neither xcmsSet nor matrix")
    }

    numAloneSpecs<-0 ## msnSpecs without ms1-parentspecs
    numMs1Peaks<- length(ms1peaks[,"mz"])
    npPeakID<-1:numMs1Peaks
    npMSnParentPeakID<-rep(0,numMs1Peaks)
    npMsLevel<-rep(1,numMs1Peaks)
    npMz<-ms1peaks[,"mz"]
    npRt<-ms1peaks[,"rt"]
    npIntensity<-ms1peaks[,"into"]
    npSample<- ms1peaks[,"sample"]

    PeakNr <- numMs1Peaks ## PeakNr+1 is the beginning peakindex for msn-spectra

    if (class(xs)=="xcmsSet"){ ## a matrix referrs to only one xr, the (so i hope) given one as parameter
	paths <- length(xs@filepaths)
	}else{paths=1}

    ## looking for every Sample-xcmsRaw (only if xs is given)
    for (NumXcmsPath in 1:paths){
        if (class(xs)=="xcmsSet"){
		xcmsRawPath <- xs@filepaths[NumXcmsPath]
        	xr <- xcmsRaw(xcmsRawPath, includeMSn = TRUE)
		}else{xr <- xraw}

        for (z in 1:length(xr@msnScanindex)){ ## looking for every msn
            if (z<length(xr@msnScanindex)) {
                GoTo<- xr@msnScanindex[z+1]
            } else {
                GoTo<- length(xr@env$msnMz)
            }
            ## the last value in the msnScanidex is the start of the last msnscan,
            ## this reaches till the end of the mz-vector
            ## now find out from which ID this spec comes from:
            if (compMethod=="floor") {pMZ<- which(floor(npMz) == floor(xr@msnPrecursorMz[z]))}
            if (compMethod=="round") {pMZ<- which(round(npMz) == round(xr@msnPrecursorMz[z]))}
            if (compMethod=="none")  {pMZ<- which(      npMz  == xr@msnPrecursorMz[z])}

            ActualParentPeakID <- 0
            biggestRt <- 0

            for (i in pMZ) { ## which rt of a n-1 scan is the biggest? taking those peakID
                if ((npMsLevel[i] == (xr@msnLevel[z]-1)) & (npSample[i] == NumXcmsPath)) {
                    if (npRt[i]<xr@msnRt[z]) {
                        if (npRt[i]>=biggestRt) {
                            biggestRt<- npRt[i]
                            ActualParentPeakID <- i
                        }
                    }
                }
            }
            if (ActualParentPeakID==0) {
                numAloneSpecs<-numAloneSpecs+1
            } else {
                MzTable = new("matrix", ncol=2,nrow=length(xr@env$msnMz[(xr@msnScanindex[z]+1):GoTo]),
                              data=c(xr@env$msnMz[(xr@msnScanindex[z]+1):GoTo],
                              xr@env$msnIntensity[(xr@msnScanindex[z]+1):GoTo]))
                colnames(MzTable)<- c("mz","intensity")

                ## calling the mini-peakpick
                npeaks <- specPeaks(MzTable, sn=snthresh, mzgap=mzgap)

                if (nrow(npeaks) >0) {
                    for (numPeaks in 1:nrow(npeaks)) {
                        ## for every picked peaks in the PeakPicked-List
                        PeakNr<- PeakNr+1
                        npPeakID[PeakNr]<- PeakNr # increasing peakid
                        npMSnParentPeakID[PeakNr]<- ActualParentPeakID
                        npMsLevel[PeakNr]<- xr@msnLevel[z]
                        npRt[PeakNr]<- xr@msnRt[z]
                        npMz[PeakNr]<- npeaks[numPeaks,"mz"]
                        npIntensity[PeakNr]<- npeaks[numPeaks,"intensity"]
                        npSample[PeakNr]<- NumXcmsPath
                        
                    }
                }
            }
        }
    }

    fragmentColnames <- c("peakID", "MSnParentPeakID","msLevel","rt", "mz",
                          "intensity", "Sample","GroupPeakMSn")

    ## later this is TRUE if the MS1-Peak is Part of a Group and has MSNs behind
    npGroupPeakMSn <- rep(FALSE,length(npSample))
    object@peaks  <- new("matrix", nrow = length(npMz), ncol = length(fragmentColnames),
                         data=c(npPeakID,npMSnParentPeakID,npMsLevel,npRt,npMz,npIntensity,
                         npSample,npGroupPeakMSn))
    colnames(object@peaks) <- fragmentColnames
    cat(length(npPeakID),"Peaks picked,",numAloneSpecs,"MSn-Specs ignored.\n")
    object
}
setGeneric("collect", function(object, ...) standardGeneric("collect"))
setMethod("collect", "xcmsFragments", .xcmsFragments.collect)

.xcmsFragments.plotTree <- function(object, mzRange=range(object@peaks[,"mz"]),
                                    rtRange=range(object@peaks[,"rt"]),
                                    xcmsSetPeakID=NULL, xcmsFragmentPeakID=NULL,
                                    textOnly=FALSE) {
    if (!textOnly) {
            require(Rgraphviz) || {
        	warning("Rgraphviz was not found, fallback to text mode.")
		textOnly=TRUE
		}
        }
    gm<-NULL
    nodeNames<-NULL

    ## recursive method for building the tree
    miniPlotTree <- function(ms1mass,level,peak,gm){
        if (level>1 |
            ((object@peaks[peak,"rt"]>=rtRange[1])&
             (object@peaks[peak,"rt"]<=rtRange[2])&
             (object@peaks[peak,"mz"]>=mzRange[1])&
             (object@peaks[peak,"mz"]<=mzRange[2]))){ ## this is the user-defined mz/rt range
            spc<-"";
            if (textOnly){
                if (object@peaks[peak,"msLevel"]>1){
                    for (b in 2:object@peaks[peak,"msLevel"]) {
                        spc<-paste(spc,"    ")
                    } ## indention of the line depends on the depth of the node in the tree
                }
                cat("|",spc,
                    "#", object@peaks[peak,"peakID"],
                    "_M", ms1mass,
                    "T", object@peaks[peak,"rt"],
                    "F", object@peaks[peak,"mz"],
                    ## "L", object@peaks[peak,"msLevel"],
                    "\n",sep="")
                } else {
                if (object@peaks[peak,"msLevel"]>1) {
                    gm <-rbind(gm,c(object@peaks[peak,"MSnParentPeakID"],
                                    object@peaks[peak,"peakID"]))
                }
            }

            for (b in which((object@peaks[,"msLevel"]>level)
                            & (object@peaks[,"MSnParentPeakID"] == peak) )) {
                gm <- miniPlotTree(ms1mass,object@peaks[b,"msLevel"],b,gm)
            }
        }
        return(gm)
    }

    if (!is.null(xcmsFragmentPeakID)) {
        ## i have only one peak (with all subpeaks) to plot
        gm <- miniPlotTree(object@peaks[xcmsFragmentPeakID,"mz"],
                           object@peaks[xcmsFragmentPeakID,"msLevel"],
                           xcmsFragmentPeakID,gm)
        } else {
            if (!is.null(xcmsSetPeakID)) {
                ## i have only one peak (with all subpeaks) to plot
                gm <- miniPlotTree(object@peaks[xcmsSetPeakID,"mz"],1,xcmsSetPeakID,gm)
            } else {
                for (a in which(object@peaks[,"msLevel"]==1)) {
                    gm <- miniPlotTree(object@peaks[a,"mz"],1,a,gm)
                }
            }
        }
    nodeNames=NULL
    nodeLabels=NULL

    for (a in object@peaks[,"peakID"]) {
        nodeNames <-c(nodeNames, paste(object@peaks[a,"peakID"]))
        nodeLabels <- c(nodeLabels, paste(round(object@peaks[a,"rt"], digits=3),"\\\n",
                                          round(object@peaks[a,"mz"], digits=3),"\\\n"))
    }

    nodes = list(label = nodeNames)
    nodes$label[nodeNames]  = nodeLabels
    if (!is.null(gm)) {
        if (!textOnly){
            colnames(gm) = c("fromID","toID")
            V <- as.character(unique(as.vector(gm)))
            OG <- new("graphNEL",nodes=V)
            edgemode(OG) <- "directed"
            OG <- addEdge(from=as.character(gm[,"fromID"]),
                          to=as.character(gm[,"toID"]),
                          graph=OG,weights=1)
            globA =list(node=list(shape="circle",width=10)  )
            Rgraphviz:::plot(OG,nodeAttrs=nodes, attrs=globA)
        }
    }
}
setGeneric("plotTree", function(object, ...) standardGeneric("plotTree"))
setMethod("plotTree", "xcmsFragments", .xcmsFragments.plotTree)

.xcmsFragments.hasMSn <- function(object,  xcmsSetPeakID) {
    ## Quick check whether MSn exists for some peak in xcmsSet
    any(object@peaks[,"MSnParentPeakID"] == xcmsSetPeakID)
}
setGeneric("hasMSn", function(object, ...) standardGeneric("hasMSn"))
setMethod("hasMSn", "xcmsFragments", .xcmsFragments.hasMSn)


xcmsFragments.makeXS <- function(xs,xf,FUNS,filename){
	## gets an xcmsSet and the corresponding xcmsFragments
	## the xs is a grouped, RTcorrected, Regrouped xs, the xf is made from this
	## returns the xs with #Samples noew rows containing the hamming-distance to the mean
	
	GetMSnVector <- function(object, xcmsSetPeakID) {
	## Returns a vector which contains all peakIDs of the msnTree with the Tree-parentmass xcmsSetPeakID		
		treewalk <- function(peakIDs,parentID){
			lnPeaks <- which((object@peaks[,"MSnParentPeakID"] == parentID))
			pIDs <- c(peakIDs, lnPeaks)
			for (a in lnPeaks){
				pIDs <- treewalk(pIDs,object@peaks[a,"peakID"])
				}
			pIDs
			}
		pIDs=NULL
		pIDs <- treewalk(pIDs,xcmsSetPeakID)	
		pIDs
		}

	## first step: Getting information: which grouped peak in which sample has a ms2Spec hanging on it
	cat("starting ")
	MSNinfo=matrix(nrow=length(groupidx(xs)), ncol=length(sampnames(xs)), 
		       data=rep(FALSE,(length(groupidx(xs))*length(sampnames(xs)))) ) ## Table(NumSamples * NumGroups)
	for (a in 1:length(groupidx(xs))){
		groV <- groupval(xs,"medret","index")[a,]
		for (b in 1:length(groV)){
			hmsn<-hasMSn(xf,groV[b])
			xf@peaks[groV[b],"GroupPeakMSn"] <- hmsn
			MSNinfo[a,xf@peaks[groV[b],"Sample"]]<-
			MSNinfo[a,xf@peaks[groV[b],"Sample"]] || hmsn
			}
		}
	
	samples <- unique(xf@peaks[,"Sample"])
	nsamp<-length(samples)
	x<-c(1:length(groupidx(xs)))
	AllMSn <- sapply(x,FUN=function(x){sum(MSNinfo[x,])}) == ncol(MSNinfo)
	## AllMSn: wich group is lucky to have MSNs in ALL samples?
	## Processing those groups:
	## first: creating a xs with the msnPeaks from all samples:
	distances=matrix(ncol=nsamp*length(FUNS),data=rep(0,(length(AllMSn)*nsamp*length(FUNS))))
	for (g in which(AllMSn)){
		##peakIDs <- groupidx(xs)[[g]]
		peakIDs <- groupval(xs,"medret","index")[g,]
		Nmz=NULL
		Nrt=NULL
		Nit=NULL
		Nsa=NULL
		for (a in 1:length(peakIDs)){ ## getting the MSns for all Samples the current Group
			MsnPeaks <- GetMSnVector(xf,peakIDs[a]) 
			Nmz <- c(Nmz,xf@peaks[MsnPeaks,"mz"])
			Nrt <- c(Nrt,xf@peaks[MsnPeaks,"rt"])
			Nit <- c(Nit,xf@peaks[MsnPeaks,"intensity"])
			Nsa <- c(Nsa,xf@peaks[MsnPeaks,"Sample"])
			}
		xsn <- new("xcmsSet")
		coln<-c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","maxo","sn","sample")
		xsn@peaks=new("matrix", ncol=length(coln), data=c(Nmz,Nmz,Nmz,Nrt,Nrt,Nrt,Nit,Nit,Nit,Nsa))
		colnames(xsn@peaks) <- coln
		phenoData(xsn) <- samples
		sampnames(xsn) <- samples
		sampclass(xsn) <- samples
		## The temporary-xs is ready now
		## apppying the given function to the object:
		#stop("does the xsSets have RTs?")
		dists=NULL
		for (fu in 1:length(FUNS)){
			dists <- c(dists,FUNS[[fu]](xsn))
			}
		distances[g,1:(nsamp*length(FUNS))] <- dists
		}#cat("loopend\n")
		
	## adding the information to groupval(grouped_object)
	cat("Groups:", length(AllMSn),"complete:",length(which(AllMSn)))
	#stop("hope it gets until here")
	write.csv(distances,file=filename)
	distances
	}
#setGeneric("makeXS", function(object, ...) standardGeneric("makeXS"))
#setMethod("makeXS", "xcmsFragments", .xcmsFragments.makeXS)

getXS<-function(xs,xf,g)
	{
	GetMSnVector <- function(object, xcmsSetPeakID) {
	## Returns a vector which contains all peakIDs of the msnTree with the Tree-parentmass xcmsSetPeakID		
		treewalk <- function(peakIDs,parentID){
			lnPeaks <- which((object@peaks[,"MSnParentPeakID"] == parentID))
			pIDs <- c(peakIDs, lnPeaks)
			for (a in lnPeaks){
				pIDs <- treewalk(pIDs,object@peaks[a,"peakID"])
				}
			pIDs
			}
		pIDs=NULL
		pIDs <- treewalk(pIDs,xcmsSetPeakID)	
		pIDs
		}

	## first step: Getting information: which grouped peak in which sample has a ms2Spec hanging on it
	samples <- unique(xf@peaks[,"Sample"])
	nsamp<-length(samples)
	x<-c(1:length(groupidx(xs)))
	peakIDs <- groupval(xs,"medret","index")[g,]
		Nmz=NULL
		Nrt=NULL
		Nit=NULL
		Nsa=NULL
		for (a in 1:length(peakIDs)){ ## getting the MSns for all Samples the current Group
			MsnPeaks <- GetMSnVector(xf,peakIDs[a]) 
			Nmz <- c(Nmz,xf@peaks[MsnPeaks,"mz"])
			Nrt <- c(Nrt,xf@peaks[MsnPeaks,"rt"])
			Nit <- c(Nit,xf@peaks[MsnPeaks,"intensity"])
			Nsa <- c(Nsa,xf@peaks[MsnPeaks,"Sample"])
			}
		#stop("RTs filled?")
		xsn <- new("xcmsSet")
		coln<-c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","maxo","sn","sample")
		xsn@peaks=new("matrix", ncol=length(coln), data=c(Nmz,Nmz,Nmz,Nrt,Nrt,Nrt,Nit,Nit,Nit,Nsa))
		colnames(xsn@peaks) <- coln
		phenoData(xsn) <- samples
		sampnames(xsn) <- samples
	xsn
	}


if (!isGeneric("simMatrix") )
    setGeneric("simMatrix", function(object, ...) standardGeneric("simMatrix"))

setMethod( "simMatrix", "xcmsFragments", function(object, ppmval=500, ...) {
##experimental
    if(length(object@MS2spec) >=1){
	stop("Undefined object slot\nPlease collect spectra first!")
    }
    ##takes very long time and needs some sorting but good info returned!!
    ##make into C++ code ??
    spectra<-object@MS2spec
    Allscores<-matrix(0, nrow=length(spectra))
    for (i in 1:dim(object@specinfo)[1] ){
	cat(object@specinfo[i,"preMZ"], " ")
	if(dim(spectra[[i]])[1] <= 1) {
	    scores<-rep(1, length(spectra))
	} else {
	    for (j in 1:dim(object@env$specinfo)[1]){
#		cat(i, ":i ", j, ":j || ")
		if(dim(spectra[[j]])[1] <=1 ){
		    score<-1
		} else {
		    neutralLoss1<-sort(diff(sort(spectra[[i]][,"mz"], decreasing=T) ))
		    neutralLoss2<-sort(diff(sort(spectra[[j]][,"mz"], decreasing=T) ))
		    score<-score_fun(neutralLoss1, neutralLoss2, ppmval)
		}
		if (j == 1){
		    scores<-score
		} else {
		    scores<-c(scores, score)
		}
	    }
	}
	if(i == 1){
	    Allscores<-scores
	}else{
	    Allscores<-rbind(Allscores,scores)
	}
    }
    colnames(Allscores)<-paste(object@specinfo[,"rtmin"], "/", object@specinfo[,"preMZ"], sep="")
    rownames(Allscores)<-paste(object@specinfo[,"rtmin"], "/", object@specinfo[,"preMZ"], sep="")
    return(Allscores)
 })

if (!isGeneric("spectralClust") )
    setGeneric("spectralClust", function(object, ...) standardGeneric("spectralClust"))

setMethod( "spectralClust", "xcmsFragments", function(object, ...){
##experimental
    clust<-data.frame()
    info<-object@specinfo
    spectra<-object@MS2spec
    for(i in 1:length(info[,"ref"]) ){
	if(dim(spectra[[i]])[1] ==1 ) {
	    i<-i+1
	}
	neutralLoss<-sort(diff(spectra[[i]][,"mz"] ))
	ind<-neutralLoss > 1.5
	neutralLOSS_I<-neutralLoss[ind][1]
	neutralLOSS_II<-neutralLoss[ind][2]
	small_I<-min(spectra[[i]][,"mz"])
	small_II<-sort(spectra[[i]][,"mz"])[2]
	
	if (dim(clust)[1] == 0){
	    clust <- matrix(c(small_I, small_II, neutralLOSS_I, neutralLOSS_II), ncol=4)
	} else {
	    clust <-rbind(clust, c(small_I, small_II, neutralLOSS_I, neutralLOSS_II))
	}
    }
    
    return(clust)

})

if (!isGeneric("findneutral") )
    setGeneric("findneutral", function(object, ...) standardGeneric("findneutral"))

setMethod("findneutral", "xcmsFragments", function(object, find, ppmE=10, print=TRUE, ...) {
    find<-ppmDev(Mr=find, ppmE) #gets the deviation window for a mass [1] is top [2] is min
    neutral<-0
    found<-0
    spectra<-object@env$MS2spec
    for (i in 1:dim(object@env$specinfo)[1] ){
	#cat(object@env$specinfo[i,"preMZ"], " ")
	if(dim(spectra[[i]])[1] >= 2) {
	    neutral<-sort(abs(diff(sort(spectra[[i]][,"mz"])) ))
	    if(any (neutral < find[1] & neutral > find[2])== TRUE){
		if(is.interger(found)){
		    found<-object@env$specinfo[i, c("ref", "preMZ", "rtmin", "rtmax")]
		}else{
		    found<-rbind(found, object@env$specinfo[i, c("ref", "preMZ", "rtmin", "rtmax")] )
		}
	    }
	}
    }
    if (print == TRUE){
	cat("We looked for", find[2], " to", find[1], "and found:\n")
	print(found)
    }
    return(found)
})


if (!isGeneric("findMZ") )
    setGeneric("findMZ", function(object, ...) standardGeneric("findMZ"))

setMethod("findMZ", "xcmsFragments", function(object, find, ppmE=10, print=TRUE, ...) {
    find<-ppmDev(Mr=find, ppmE) #gets the deviation window for a mass [1] is top [2] is min
    fragMZ<-0
    found<-0
    spectra<-object@env$MS2spec
    for (i in 1:dim(object@env$specinfo)[1] ){
	#cat(object@env$specinfo[i,"preMZ"], " ")
	if(dim(spectra[[i]])[1] >= 2) {
	    fragMZ<-sort(spectra[[i]][,"mz"])
	    if(any (fragMZ < find[1] & fragMZl > find[2])== TRUE){
		if(is.interger(found)){
		    found<-object@env$specinfo[i, c("ref", "preMZ", "rtmin", "rtmax")]
		}else{
		    found<-rbind(found, object@env$specinfo[i, c("ref", "preMZ", "rtmin", "rtmax")] )
		}
	    }
	}
    }
    if (print == TRUE){
	cat("We looked for", find[2], " to", find[1], "and found:\n")
	print(found)
    }
    return(found)
})

if (!isGeneric("searchMetlin") )
    setGeneric("searchMetlin", function(object, ...) standardGeneric("searchMetlin"))

setMethod( "searchMetlin", "xcmsFragments", function(object, ppmfrag=10, ppmMZ= 5, file, MS1data=FALSE, metXML="metlin", ...) {
    if(metXML=="metlin"){
        metlinfile<-"http://metlin.scripps.edu/download/MSMS.XML"
        metlinMS<-"http://metlin.scripps.edu/download/MS.XML"
    }else{
        metlinfile<-metXML
    }
    if(MS1data==TRUE){
        met.xml<-read.metlin(metlinfile, MS1=TRUE, MSxml=metlinMS)
    }else{
        met.xml<-read.metlin(metlinfile, MS1=FALSE)
    }
    #prob<-confidence(metlinfile, ppmMZ, ppmfrag, mode=method)
    
    Na<-c(23, 22.98976)
    H<-c(1, 1.00794)
    H.n<-c(-1, -1.00794)
    Cl<-c(-35, -35.453) # hard coded values of molecular weight of elements
    e<-c(0, 0.0005)
    elements<-rbind(H, Na, H.n, Cl, e)
    CEref<-c(0, 10, 20, 40)

    check=FALSE
    for(i in 1:dim(object@specinfo)[1] ){
	mz.diff<-object@specinfo[i, "preMZ"] - round(as.numeric(met.xml[,"preMZ"]), 0)
	Index<-0
	if(any(mz.diff) == elements["H",1] ){ # check to see if what type of ionisation it is and if it's an adduct
	    deviate<-ppmDev(object@specinfo[i, "AccMZ"]-elements["H",2], ppmMZ)
	    exp.mode<-c("+", "")
	    Index<-which(as.numeric(met.xml[,"preMZ"]) < deviate[1] & as.numeric(met.xml[,"preMZ"]) > deviate[2])
	} else if (any(mz.diff) == elements["Na",1]){
	    deviate<-ppmDev(object@specinfo[i, "AccMZ"]-(elements["Na",2]+elements["e",2]), ppmMZ)
	    exp.mode<-c("+","Na")
	    Index<-which(as.numeric(met.xml[,"preMZ"]) < deviate[1] & as.numeric(met.xml[,"preMZ"]) > deviate[2])
	} else if (any(mz.diff) == elements["H.n",1]){
	    deviate<-ppmDev(object@specinfo[i, "AccMZ"]+elements["H.n",2], ppmMZ)
	    exp.mode<-c("-", "")
	    Index<-which(as.numeric(met.xml[,"preMZ"]) < deviate[1] & as.numeric(met.xml[,"preMZ"]) > deviate[2])
	} else if (any(mz.diff) == elements["Cl",1]){
	    deviate<-ppmDev(object@specinfo[i, "AccMZ"]+elements["Cl",2], ppmMZ)
	    exp.mode<-c("-", "Cl")
	    Index<-which(as.numeric(met.xml[,"preMZ"]) < deviate[1] & as.numeric(met.xml[,"preMZ"]) > deviate[2])
	} ## index the DB by accurate mass
	met<-met.xml[Index,]
        CE<-object@specinfo[i,"CollisionEnergy"]
	if(dim(met)[1] ==0){ ## If it doesn't exist go on.
	    next 
	}

	uni.met<-unique(met[,"name"]) ## This comes out as a factor, clean it?
	for (j in 1:length(uni.met)){
	    cat(" ", object@specinfo[i,"preMZ"])
	    nameIndex<- which(met[,"name"] == uni.met[j])
	    IonIndex<-which(met[, "mode"] == exp.mode[1] & met[, "adduct"] == exp.mode[2])
            if(length(IonIndex) == 0) {
		IonIndex<-which(met[, "mode"] == "*" &  met[, "adduct"] == "M")
	    }
            if (any(CEref == CE)){ ## Do we have the same collision energy ?
                CeIndex<-which(met[,"collisionEnergy"] == CE)
            }else{
                greaterThan<-which(CEref > CE)[1]
                LessThan<-which(CEref < CE)
                LessThan<-LessThan[length(LessThan)] 
                if((CE-CEref[LessThan]) > (CEref[greaterThan] -CE)){
                    CeIndex<-which(met[,"collisionEnergy"] == CEref[greaterThan])
                }else if ((CE-CEref[LessThan]) < (CEref[greaterThan] -CE)) {
                    CeIndex<-which(met[,"collisionEnergy"] == CEref[LessThan])
                }else{
                    CeIndex<-which(met[,"collisionEnergy"] == CEref[greaterThan])
                }
            }
            SpecIndex<-overlap(nameIndex, IonIndex, CeIndex)
            object@MS2spec[[i]][, "mz"]<-as.matrix(object@MS2spec[[i]][, "mz"])
	    if(dim(met[SpecIndex,])[1] == 0){
		##This should only happen when we pick up a molecule that doesn't
		## corrospond to the picked up ionisation or adduct that we think it is
		next ## so jump ahead
	    }else if(met[SpecIndex,"frag.MZ"][1] == 0 ){ ##only used if MS1 matching as well
		if (j ==1 ){
		    dist<-c(i, j, object@specinfo[i,"AccMZ"], object@specinfo[i,"rtmin"], object@specinfo[i, "rtmax"], object@specinfo[i, "CollisionEnergy"], met[SpecIndex, "collisionEnergy"][1], 0.0,  met[SpecIndex,"preMZ"][1])
		    name<-met[SpecIndex, "name"][1]
		    mode<-met[SpecIndex, "mode"][1]
		    adduct<-met[SpecIndex, "adduct"][1]
		}else{
		    dist<-rbind(dist, c(i, j, object@specinfo[i,"AccMZ"], object@specinfo[i,"rtmin"], object@specinfo[i, "rtmax"], object@specinfo[i, "CollisionEnergy"], met[SpecIndex, "collisionEnergy"][1], 0.0,  met[SpecIndex,"preMZ"][1]))
		    name<-c(name, met[SpecIndex, "name"][1])
		    mode<-c(mode, met[SpecIndex, "mode"][1])
		    adduct<-c(adduct, met[SpecIndex, "adduct"][1])
		}	
	    } else { ## Just MS/MS matching no MS^n yet!!
                #metSpec<-as.numeric(met[SpecIndex,"frag.MZ"])
                #metSpec<-as.numeric(cbind(metSpec, met[SpecIndex, "int"]))
                #metSpec<-data.frame(as.numeric(met[SpecIndex,"frag.MZ"]),  as.numeric(met[SpecIndex, "int"]))
                #colnames(metSpec)<-c("mz", "intensity")
                #metSpec<-deisotopeNclean(metSpec)
		if(j ==1){ ## if it doesn't exist make it
		    score<-score_fun(as.numeric(met[SpecIndex,"frag.MZ"]), object@MS2spec[[i]][,"mz"], ppmval=ppmfrag)
		    if(length(file)){
			eicdir<-paste(file, "_spectra", sep="")
			if(!file.exists(eicdir)){
			    dir.create(eicdir)
			}
			if(capabilities("png")){
			    png(file.path(eicdir, paste(i,"-",j, ".png", sep="")), width=1024, height=768)
			} else{
			    pdf(file.path(eicdir, "%03d.pdf"), width=1024/72, height=768/72, onefile=FALSE)
			}
		    }
		    plot.metlin(met[SpecIndex, c("frag.MZ", "int")], object@MS2spec[[i]], i, j, MZlabel=object@specinfo[i,"preMZ"], ... )
		    if(length(file)){
			dev.off()
		    }

		    #score<-dnorm(score, prob[[2]], prob[[1]])
		    dist<-c(i, j, object@specinfo[i,"AccMZ"], object@specinfo[i,"rtmin"], object@specinfo[i, "rtmax"], object@specinfo[i, "CollisionEnergy"], met[SpecIndex, "collisionEnergy"][1], score, met[SpecIndex,"preMZ"][1], similar(sort(met[SpecIndex,"frag.MZ"]), sort(object@MS2spec[[i]][,"mz"]), ppmfrag), distance(sort(met[SpecIndex,"frag.MZ"]), sort(object@MS2spec[[i]][,"mz"]), ppmfrag), length(met[SpecIndex,"preMZ"]))

		    name<-met[SpecIndex, "name"][1]
		    mode<-met[SpecIndex, "mode"][1]
		    adduct<-met[SpecIndex, "adduct"][1]
		}else{ ## when it does exist add to it (has to be a better way?!?!)
		    score<-score_fun(as.numeric(met[SpecIndex,"frag.MZ"]), object@MS2spec[[i]][,"mz"], ppmval=ppmfrag)
		    
		    if(length(file)){
			eicdir<-paste(file, "_spectra", sep="")
			if(!file.exists(eicdir)){
			    dir.create(eicdir)
			}
			if(capabilities("png")){
			    png(file.path(eicdir, paste(i,"-",j, ".png", sep="")), width=1024, height=768)
			} else{
			    pdf(file.path(eicdir, "%03d.pdf"), width=1024/72, height=768/72, onefile=FALSE)
			}
		    }
		    plot.metlin(met[SpecIndex, c("frag.MZ", "int")], object@MS2spec[[i]], i, j, MZlabel=object@specinfo[i,"preMZ"], ... )
		    if(length(file)){
			dev.off()
		    }
		    #score<-dnorm(score, prob[[2]], prob[[1]])
		    dist<-rbind(dist,c(i, j, object@specinfo[i, "AccMZ"], object@specinfo[i, "rtmin"], object@specinfo[i, "rtmax"], object@specinfo[i, "CollisionEnergy"], met[SpecIndex, "collisionEnergy"][1], score, met[SpecIndex,"preMZ"][1], similar(sort(met[SpecIndex,"frag.MZ"]), sort(object@MS2spec[[i]][,"mz"]), ppmfrag), distance(sort(met[SpecIndex,"frag.MZ"]), sort(object@MS2spec[[i]][,"mz"]), ppmfrag), length(met[SpecIndex,"preMZ"]) ) )
		    name<-c(name, met[SpecIndex, "name"][1])
		    mode<-c(mode, met[SpecIndex, "mode"][1])
		    adduct<-c(adduct, met[SpecIndex, "adduct"][1])
		}
	    }
	}
	if(check == FALSE) {
	    all.dist<-dist
	    all.name<-name
	    all.mode<-mode
	    all.adduct<-adduct
	    check <-TRUE
	} else{
	    all.dist<-rbind(all.dist, dist)
	    all.name<-c(all.name, name)
	    all.mode<-c(all.mode, mode)
	    all.adduct<-c(all.adduct, adduct)
	}
    }

    if (check == FALSE) {
	cat("\nNothing match the masses in metlin, try a larger ppm range. \n")
    }else{
	dist.df<-data.frame(all.dist, all.name, all.mode, all.adduct, stringsAsFactors=FALSE, row.names=1:length(all.name))
	colnames(dist.df)<-c("A", "B", "Precursor Ion", "rtmin", "rtmax", "CollisionEnergy experiment", "CollisionEnergy Reference", "Percentage Match", "Metlin Mass", "# matching", "# non-matching", "Total # Ref ion", "Metlin ID Name",  "Ionization", "Adduct")
	#rownames(dist.df)<-rep("", dim(dist.df)[1])
	cat("\nDone", "\n")
	
	file.tsv<-paste(file, ".tsv", sep="")
	write.table(dist.df, file=file.tsv, sep="\t", row.names=FALSE)
	invisible(dist.df)
    }
})


plot.metlin<-function(MetSpec, ExpSpec, placeA, placeB, MZlabel,col=c("red", "blue"), neg=TRUE){
    ExpSpec[,"intensity"]<-ExpSpec[,"intensity"]/max(ExpSpec[,"intensity"])*100
    maxMZ<-max(ExpSpec[,1], MetSpec[,1]) ##I know could be better with "mz" but not the same sort later
    ##ind<-which(MetSpec[,"int"] <= 90)
    ##MetSpec[ind, "int"]<-MetSpec[ind, "int"]/max(MetSpec[ind,"int"])*75
    if(neg == FALSE){
        op <- par(mfrow = c(2, 1),pty = "m", adj=0.5) # 1 x 2 pictures on one plot
    
        plot(ExpSpec, type="h", col=col[1], ylab="Intensity", xlab="", main=paste("MS/MS spectra for " , MZlabel , sep=""), sub=paste(placeA, " - ", placeB, sep=""), xlim=c(0, (maxMZ+(maxMZ/2))))
        abline(0,0)
        plot(MetSpec, type="h", col=col[2], xlab="M/Z", ylab="", main="Metlin Reference Spectrum", xlim=c(0, (maxMZ+(maxMZ/2) )) )
        abline(0,0)
        par(op)
    }else{
        plot(ExpSpec[,1:2],type="h", col=col[1], ylim=c(-max(c(ExpSpec[,2], MetSpec[,2])), max(c(ExpSpec[,2], MetSpec[,2]))), main=paste("MS/MS spectra for " , MZlabel , sep=""), sub=paste(placeA, " - ", placeB, sep=""), ylab="Intensity",xlab="m/z")
        points(specA[,1], -(MetSpec[,2]), type="h", col=col[2])
        abline(0,0)
        legend("bottomright", paste("Reference", sep=""), col="blue", pch=16)
        legend("topright", paste("Experimental", sep=""), col="red",pch=16)

        
        for(i in 1:nrow(ExpSpec)){
            text(ExpSpec[i,1], ExpSpec[i,2], round(ExpSpec[i,1], 1))
        }
        for(j in 1:nrow(MetSpec)){
            text(MetSpec[j,1], -(MetSpec[j,2]), round(MetSpec[j,1], 1))
        }
    }
}

if (!isGeneric("simSearch") )
  setGeneric("simSearch", function(object,...) standardGeneric("simSearch"))

setMethod( "simSearch", "xcmsFragments", function(object, ppmfrag=20, percent=50, file, fullReport=FALSE, ...) {
    metlinfile<-"http://metlin.scripps.edu/download/MSMS.XML"
    spectra<-metlinToList(metlinfile)
    cat("Data converted\nProcessing data...\n")
    check=FALSE
    for(i in 19:length(object@MS2spec)){
	if(!is.matrix(object@MS2spec[[i]])) {
	    next ## go on if not neutral loss capible 
	}
        if(dim(object@MS2spec[[i]])[1] <= 1 ){
            next ## go on cus no neutral losses!!
        }
#cat(" i ->", i, "\t")
	neutralExp<-sort(abs(diff(sort(as.matrix(object@MS2spec[[i]])[,"mz"] )) ))
        cat(paste(object@specinfo[i,"preMZ"], " ", sep=""))
	for(j in 1:length(spectra)){
	    if(dim(spectra[[j]])[1] < 1 ){
	        next ## if we have something which doesn't have a neutral loss go on
	    }
#cat(" j->", j)
            if(length(spectra[[j]][,"frag.MZ"]) <= 1){ ##some things are getting by don't know why!!
                next
            }
	    neutralMet<-sort(abs(diff(sort(spectra[[j]][,"frag.MZ"])) ))##make neutral losses
            
            NeutScore<-score_fun(neutralMet, neutralExp, ppmfrag)
            FragScore<-score_fun(spectra[[j]][,"frag.MZ"], object@MS2spec[[i]][,"mz"], ppmfrag)
            
	    if(NeutScore >= percent | FragScore >= percent){
                #cat(paste(" .", sep=""))
                neutralFreq<-ms2Freq(c(neutralMet, neutralExp), ppmError=ppmfrag)
                fragmentFreq<-ms2Freq(c(spectra[[j]][,"frag.MZ"], object@MS2spec[[i]][,"mz"]), ppmError=ppmfrag)
                indNeutral<-which.max(neutralFreq[[2]])
                indFrag<-which.max(fragmentFreq[[2]])
                commonNeutral<-neutralFreq[[1]][indNeutral]
                commonFrag<-fragmentFreq[[1]][indFrag]
                

	        if(check==FALSE){
                    result<-c(round(object@specinfo[i, "AccMZ"], 4), object@specinfo[i, "rtmin"], object@specinfo[i, "rtmax"], object@specinfo[i, "CollisionEnergy"], FragScore, NeutScore, commonNeutral, commonFrag, spectra[[j]][1,"name"], spectra[[j]][1,"preMZ"], spectra[[j]][1,"collisionEnergy"])
                    check<-TRUE
                 }else {
                    result<-rbind(result, c(round(object@specinfo[i, "AccMZ"],4), object@specinfo[i, "rtmin"], object@specinfo[i, "rtmax"], object@specinfo[i, "CollisionEnergy"],FragScore, NeutScore, commonNeutral, commonFrag, spectra[[j]][1,"name"], spectra[[j]][1,"preMZ"], spectra[[j]][1,"collisionEnergy"] ))
	        }
	    }
	}
        
    }
    
    colnames(result)<-c("m/z", "rtmin", "rtmax", "Experiment Collision Energy", "Fragment Score", "Neutral Score", "Common Neutral loss", "Common Fragment", "Compound Name", "Metlin Mass", "Collision Energy")
    rownames(result)<-rep("", dim(result)[1])
    cat(paste("Searching done\n", "Grouping...", sep=""))
    RmzUnique<-unique(round(as.numeric(result[,"m/z"]), 1) ) ##Make the groups
    for(k in 1:length(RmzUnique)){
        cat(" ", k)
        index<-which(round(as.numeric(result[,"m/z"]),1) == RmzUnique[k])
        
        comNeutFreq<-ms2Freq(as.numeric(result[index, "Common Neutral loss"]), ppmfrag)
        FreqIndexNeut<-order(comNeutFreq[[2]])
        TotalCommonLoss<-comNeutFreq[[1]][FreqIndexNeut][1:5]

        comFragFreq<-ms2Freq(as.numeric(result[index, "Common Fragment"]), ppmfrag)
        FreqIndexFrag<-order(comFragFreq[[2]])
        TotalCommonFrag<-comFragFreq[[1]][FreqIndexFrag][1:5]

        if(!exists("CommonIonResult")){
            CommonIonResult<-cbind(round(as.numeric(result[index, "m/z"][1]),1), paste(TotalCommonFrag, collapse=" :"), paste(TotalCommonLoss, collapse=" :"))
        } else {
            CommonIonResult<-rbind(CommonIonResult, cbind(round(as.numeric(result[index, "m/z"][1]),1), paste(TotalCommonFrag, collapse=" :"), paste(TotalCommonLoss, collapse=" :")))
        }
    }
    colnames(CommonIonResult)<-c("m/z", "Top 5 Common Fragment", "Top 5 common Neutral Losses")
    write.table(CommonIonResult, file=paste(file, ".tsv", sep=""), sep="\t", row.names=FALSE)

    cat("\n")
    if (fullReport == TRUE){
        write.table(result, file=paste(file, "_fullReport.tsv", sep=""), sep="\t", row.names=FALSE)
    }
    invisible(result)
})

score_fun<-function(ref, exp, ppmval){
    ref<-sort(ref)
    exp<-sort(exp)
#    Sscore<-similar(ref, exp, ppm=ppmval)# / max(c(ref, exp))
##error when the ref or exp is 2X the length of the ref or exp :(
    Sscore<-0
    for(i in 1:length(ref)){
      for(j in 1:length(exp)){
        if(ppm(ref[i], exp[j]) < ppmval){
          Sscore<-Sscore+1
        }
      }
    }
    Dscore<-distance(ref, exp, ppm=ppmval)# / max(c(ref, exp)) ## get the scores for the spectra
    lenR<-length(ref)
    lenE<-length(exp) ##find the lengths of the spectra

    score<-Sscore-Dscore ##Calculate the unequiolavated score
    topScore<-min(lenE, lenR) - (max(lenE, lenR) - min(lenR, lenE)) ## This is 100% match (similarity - distance )
    ZeroScore<-0 - max(lenE, lenR)## This is the 0% (similarity - distance)

    score<- 100*abs((score-ZeroScore)/(ZeroScore -topScore) )## And finally equiovalated to a percentage
    invisible(round(score,2))
}

ms2Freq<-function(fragments, ppmError=10){
    freq<-0
    observedFragments<-0
    ObservedFreq<-0
    for(i in 1:length(fragments)){
        for(j in 1:length(fragments)){
            if(ppm(fragments[j], fragments[i]) <= ppmError){
                freq<-freq+1
            }
            if(j == length(fragments)){
                observedFragments<-c(observedFragments, fragments[i])
                ObservedFreq<-c(ObservedFreq, freq)
                freq<-0
            }
        }
    }
    x<-list()
    x[[1]]<-observedFragments
    x[[2]]<-ObservedFreq
    invisible(x)
}


