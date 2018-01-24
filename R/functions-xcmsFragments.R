## Functions for xcmsFragments
#' @include DataClasses.R

############################################################
## xcmsFragments constructor.
xcmsFragments <- function(xs = NULL, ...) {
    object <- new("xcmsFragments")
    object <- collect(object,xs,...)
    object
}

.xcmsFragments.show <- function(object) {
    if(length(object@peaks)){
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
    npCollisionEnergy<-rep(0, numMs1Peaks)

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

        npCollisionEnergy<-xr@msnCollisionEnergy

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
                        npPeakID[PeakNr]<- PeakNr ## increasing peakid
                        npMSnParentPeakID[PeakNr]<- ActualParentPeakID
                        npMsLevel[PeakNr]<- xr@msnLevel[z]
                        npRt[PeakNr]<- xr@msnRt[z]
                        npMz[PeakNr]<- npeaks[numPeaks,"mz"]
                        npIntensity[PeakNr]<- npeaks[numPeaks,"intensity"]
                        npSample[PeakNr]<- NumXcmsPath
                        npCollisionEnergy[PeakNr]<-xr@msnCollisionEnergy[z]
                    }
                }
            }
        }
    }

    fragmentColnames <- c("peakID", "MSnParentPeakID","msLevel","rt", "mz",
                          "intensity", "Sample","GroupPeakMSn", "CollisionEnergy")

    ## later this is TRUE if the MS1-Peak is Part of a Group and has MSNs behind
    npGroupPeakMSn <- rep(FALSE,length(npSample))
    object@peaks  <- new("matrix", nrow = length(npMz), ncol = length(fragmentColnames),
                         data=c(npPeakID,npMSnParentPeakID,npMsLevel,npRt,npMz,npIntensity,
                         npSample,npGroupPeakMSn, npCollisionEnergy))
    colnames(object@peaks) <- fragmentColnames
    cat(length(npPeakID),"Peaks picked,",numAloneSpecs,"MSn-Specs ignored.\n")
    object
}

.xcmsFragments.plotTree <- function(object, mzRange=range(object@peaks[,"mz"]),
                                    rtRange=range(object@peaks[,"rt"]),
                                    xcmsSetPeakID=NULL, xcmsFragmentPeakID=NULL,
                                    textOnly=FALSE) {
    libname <- "Rgraphviz"

    if (!textOnly) {
        (require(libname,character.only=TRUE,quietly=TRUE)) || {
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
                    "##", object@peaks[peak,"peakID"],
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
            Rgraphviz::plot(OG,nodeAttrs=nodes, attrs=globA)
        }
    }
}

.xcmsFragments.hasMSn <- function(object,  xcmsSetPeakID) {
    ## Quick check whether MSn exists for some peak in xcmsSet
    any(object@peaks[,"MSnParentPeakID"] == xcmsSetPeakID)
}

xcmsFragments.makeXS <- function(xs,xf,FUNS,filename){
    ## gets an xcmsSet and the corresponding xcmsFragments
    ## the xs is a grouped, RTcorrected, Regrouped xs, the xf is made from this
    ## returns the xs with ##Samples now rows containing the hamming-distance to the mean

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
    ## AllMSn: which group is lucky to have MSNs in ALL samples?
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
        ##stop("does the xsSets have RTs?")
        dists=NULL
        for (fu in 1:length(FUNS)){
            dists <- c(dists,FUNS[[fu]](xsn))
        }
        distances[g,1:(nsamp*length(FUNS))] <- dists
    }##cat("loopend\n")

    ## adding the information to groupval(grouped_object)
    cat("Groups:", length(AllMSn),"complete:",length(which(AllMSn)))
    ##stop("hope it gets until here")
    write.csv(distances,file=filename)
    distances
}

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
    ##stop("RTs filled?")
    xsn <- new("xcmsSet")
    coln<-c("mz","mzmin","mzmax","rt","rtmin","rtmax","into","maxo","sn","sample")
    xsn@peaks=new("matrix", ncol=length(coln), data=c(Nmz,Nmz,Nmz,Nrt,Nrt,Nrt,Nit,Nit,Nit,Nsa))
    colnames(xsn@peaks) <- coln
    phenoData(xsn) <- samples
    sampnames(xsn) <- samples
    xsn
}
