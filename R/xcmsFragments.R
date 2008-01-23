require(methods) || stop("Couldn't load package methods")

setClass("xcmsFragments", representation(peaks = "matrix"
                                         #, pipeline = "xcmsRawPipeline"
                                         ),
         prototype(peaks = matrix(nrow = 0, ncol = 6)
                   #, pipeline = new("xcmsRawPipeline")
                   ))

xcmsFragments <- function(xs = NULL,
                          ..., pipeline = NULL) {
    object <- new("xcmsFragments")
#    pipeline(object) <- pipeline
    object
}

.xcmsFragments.collect <- function(object, xs, xr, compMethod="floor", snthresh=20, mzgap=.2, uniq=TRUE) {
    ## This is called after findPeaks()
    ## e.g. during xcmsSet() constructor
    ## This creates the association between the xcmsSet peaks
    ## and the MSn scans in the raw files.

    if (class(xs)=="xcmsSet") {
        ms1peaks <- peaks(xs)
    } else if (class(xs)=="matrix") {
        ms1peaks <- xs
    } else {
        stop("xs is neither xcmsSet nor matrix")
    }

    numMs1Peaks<- length(ms1peaks[,"mz"])
    npPeakID<-1:numMs1Peaks
    npMSnParentPeakID<-rep(0,numMs1Peaks)
    npMsLevel<-rep(1,numMs1Peaks)
    npMz<-ms1peaks[,"mz"]
    npRt<-ms1peaks[,"rt"]
    npIntensity<-ms1peaks[,"into"]

    PeakNr <- numMs1Peaks ## PeakNr+1 is the beginning peakindex for msn-spectra
     ## now to do: insert the msndata from the xcmsraw-object into the table
    for (z in 1:length(xr@msnScanindex)){ ## looking for every msn-spectrum
        if (z<length(xr@msnScanindex)){
            GoTo<- xr@msnScanindex[z+1]
            }else{
            GoTo<- length(xr@env$msnMz)
            }            ## the last value in the msnScanidex is the start of the last msnscan, this reaches til the end of the mz-vector

        ## now i have to find out from which ID this spec comes from:
        if (compMethod=="floor") {pMZ<- which(floor(npMz) == xr@msnPrecursorMz[z])}
        if (compMethod=="round") {pMZ<- which(round(npMz) == xr@msnPrecursorMz[z])}
        if (compMethod=="none")  {pMZ<- which(      npMz  == xr@msnPrecursorMz[z])}
        ActualParentPeakID<- 0
        biggestRt<- 0
        for (i in pMZ) ## which rt of a n-1 scan is the biggest? taking those peakID
            if (npMsLevel[i] == (xr@msnLevel[z]-1))
                if (npRt[i]<xr@msnRt[z])
                    if (npRt[i]>=biggestRt)
                        {
                        biggestRt<- npRt[i]
                        ActualParentPeakID <- i
                        }
        if (ActualParentPeakID==0){
        cat("Warning: the Scan",z,"L",xr@msnLevel[z],"mz=",xr@msnPrecursorMz[z], "Has no ParentPeak!\n")  ## was ist wenn der peak in der tabelle nicht existiert?
            }else{
                    #cat("msnScan",z,"PpeakiD",ActualParentPeakID,"\n")
                    MzTable = new("matrix", ncol=2,nrow=length(xr@env$msnMz[(xr@msnScanindex[z]+1):GoTo]),
                    data=c(xr@env$msnMz[(xr@msnScanindex[z]+1):GoTo],xr@env$msnIntensity[(xr@msnScanindex[z]+1):GoTo])  )
                    colnames(MzTable)<- c("mz","intensity")
                    npeaks <- specPeaks(MzTable, sn=snthresh, mzgap=mzgap) ## calling the mini-peakpick
                    if (nrow(npeaks)>0){
                        for (numPeaks in 1:nrow(npeaks)){ ## for every picked peaks in the PeakPicked-List
                            PeakNr<- PeakNr+1
                            npPeakID[PeakNr]<- PeakNr # increasing peakid
                            npMSnParentPeakID[PeakNr]<- ActualParentPeakID
                            npMsLevel[PeakNr]<- xr@msnLevel[z]
                            npRt[PeakNr]<- xr@msnRt[z]
                            npMz[PeakNr]<- npeaks[numPeaks,"mz"]
                            npIntensity[PeakNr]<- npeaks[numPeaks,"intensity"]
                            }
                        }
                    }
         }

#         ## Window and centroid from initial MS1 peak
#         mz <- c(peak$mzmin, peak$mz, peak$mzmax)
#         rt <- c(peak$rtmin, peak$rt, peak$rtmax)
#         msLevel <- 2
#
#
#         ## the getMSn* should go to xcmsRaw
#         ##
#
#         while (msLevel <= maxMsLevel) {
#             ## Select the "best fitting/most beatiful" child
#             fragmentScanId <- getMSnSpecId (xr, msLevel, peaks$mz, peak$rt)
#
#             if (!is.null(fragmentScanID)) {
#                 warning("No MSn found for peak", peak)
#                 break
#             }
#
#             ## Let's say Scan Nr. 17 in the xcmsRaw
#             fragmentScanId <- c(17)
#
#             ## Call "mini-Peakpicker on MSn spectra
#             ## remove small, noisy raw data
#             fragmentPeaks <- specPeaks(, sn=snthresh, mzgap=mzgap)
#             ## assemble xcmsFragments row
#             ## newpeaks[i, ] <- c(newpeakid, peak$peakID,  xr$msnPrecursorScan[17],
#             ##                    msLevel, fragmentPeaks$mz, xr$msnRt[17], fragmentPeaks$intensity)
#
#             ## Collect children in next level
#             msLevel <- msLevel +1
#         }
#
#         rbind(object@peaks, newpeaks)
#     }
#  stop("jetze müsste ferddich sein")
  fragmentColnames <- c("peakID", "MSnParentPeakID","msLevel","rt", "mz", "intensity")
  object@peaks  <- new("matrix", nrow = length(npMz), ncol = length(fragmentColnames),data=c(npPeakID,npMSnParentPeakID,npMsLevel,npRt,npMz,npIntensity))
  colnames(object@peaks) <- fragmentColnames
  object
  }
setGeneric("collect", function(object, ...) standardGeneric("collect"))
setMethod("collect", "xcmsFragments", .xcmsFragments.collect)

.xcmsFragments.plotTree <- function(object, mzRange=range(object@peaks[,"mz"]), rtRange=range(object@peaks[,"rt"]),
                                    xcmsSetPeakID=NULL, xcmsFragmentPeakID=NULL) {
    range(object@peaks[,"mz"])
    miniPlotTree <- function(ms1mass,level,peak){
        if (level>1 |
            ((object@peaks[peak,"rt"]>=rtRange[1])&
            (object@peaks[peak,"rt"]<=rtRange[2])&
            (object@peaks[peak,"mz"]>=mzRange[1])&
            (object@peaks[peak,"mz"]<=mzRange[2]))){ ## this is the mz/rt range the user can define
            spc<-"";
            if (object@peaks[peak,"msLevel"]>1){
                for (b in 2:object@peaks[peak,"msLevel"]) {spc<-paste(spc,"    ")} ## indention of the line depends on the depth of the node in the tree
                    }
            cat("|",spc,"#",object@peaks[peak,"peakID"],"_M",ms1mass,"T",object@peaks[peak,"rt"],"F",object@peaks[peak,"mz"],"L",object@peaks[peak,"msLevel"],"\n",sep="")
            for (b in which((object@peaks[,"msLevel"]>level) & (object@peaks[,"MSnParentPeakID"] == peak) )){
                miniPlotTree(ms1mass,object@peaks[b,"msLevel"],b)
                }
            }
        }
    if (!is.null(xcmsFragmentPeakID)) {
        miniPlotTree(object@peaks[xcmsFragmentPeakID,"mz"],object@peaks[xcmsFragmentPeakID,"msLevel"],xcmsFragmentPeakID) ## i have only one peak (with all subpeaks) to plot
        }else{
        if (!is.null(xcmsSetPeakID)) { ## i have only one peak (with all subpeaks) to plot
            miniPlotTree(object@peaks[xcmsSetPeakID,"mz"],1,xcmsSetPeakID)
            }else{
                for (a in which(object@peaks[,"msLevel"]==1)){
                miniPlotTree(object@peaks[a,"mz"],1,a)
                }
            }
        }
    }

setGeneric("plotTree", function(object, ...) standardGeneric("plotTree"))
setMethod("plotTree", "xcmsFragments", .xcmsFragments.plotTree)

xcmsFragments.hasMSn <- function(object,  xcmsSetPeakID) {
 ## Quick check whether MSn exists for some peak in xcmsSet
 ## return msnSpecID what? there are more than 1 msnspecids for one peak!
 return(any(object@peaks[,"msnParentPeakID"] == xcmsSetPeakID))
}

# .xcmsFragments.distance <- function(object, msnSpecID1, msnSpecID2) {
#     ## A spectrum distance measure
#     warning("NYI")
#     return (17)
# }
#
 .xcmsFragments.show <- function(object) {

     cat("An \"xcmsFragments\" object with ",nrow(object@peaks)," peaks in",length(unique(object@peaks[,"rt"])),"Spectra\n")
     cat("From Level 1 to",max(object@peaks[,"msLevel"]),"\n")
#    show(object@pipeline)
     cat("\n")
     memsize <- object.size(object)
     cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
}
#
setMethod("show", "xcmsFragments", .xcmsFragments.show)
