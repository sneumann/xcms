require(methods) || stop("Couldn't load package methods")

setClass("xcmsFragments", representation(peaks = "matrix"
                                        #, pipeline = "xcmsRawPipeline"
                                         ),
         prototype(peaks = matrix(nrow = 0, ncol = 6)
                                        #, pipeline = new("xcmsRawPipeline")
                   ))

xcmsFragments <- function(xs = NULL, ...) {
    object <- new("xcmsFragments")
    object <- collect(object,xs,...)
    object
}

.xcmsFragments.collect <- function(object, xs,compMethod="floor", snthresh=20, mzgap=.2, uniq=TRUE) {
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

    ## looking for every Sample-xcmsRaw
    for (NumXcmsPath in 1:length(xs@filepaths)){
        xcmsRawPath <- xs@filepaths[NumXcmsPath]
        xr <- xcmsRaw(xcmsRawPath, includeMSn = TRUE)

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

                if (nrow(npeaks) > 0) {
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
            plot(OG,nodeAttrs=nodes, attrs=globA)
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

.xcmsFragments.show <- function(object) {
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

        memsize <- object.size(object)
        cat("\nMemory usage:", signif(memsize/2^20, 3), "MB\n")
    }
setMethod("show", "xcmsFragments", .xcmsFragments.show)
