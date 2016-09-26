## On the long run it would be nice to have all generics in here.
## Alphabetically ordered.

## A
setGeneric("absent", function(object, class, minfrac) standardGeneric("absent"))
setGeneric("AutoLockMass", function(object) standardGeneric("AutoLockMass"))

## C
setGeneric("calibrate", function(object, ...) standardGeneric("calibrate"))
setGeneric("collect", function(object, ...) standardGeneric("collect"))

## D
setGeneric("deepCopy", function(object) standardGeneric("deepCopy"))
setGeneric("diffreport", function(object, ...) standardGeneric("diffreport"))

## F
setGeneric("filepaths", function(object) standardGeneric("filepaths"))
setGeneric("filepaths<-", function(object, value) standardGeneric("filepaths<-"))
setGeneric("fillPeaks.chrom", function(object, ...)
    standardGeneric("fillPeaks.chrom"))
setGeneric("fillPeaks.MSW", function(object, ...)
    standardGeneric("fillPeaks.MSW"))
setGeneric("fillPeaks", function(object, ...) standardGeneric("fillPeaks"))
setGeneric("findMZ", function(object, find, ppmE=25, print=TRUE)
    standardGeneric("findMZ"))
setGeneric("findmzROI", function(object, ...) standardGeneric("findmzROI"))
setGeneric("findneutral", function(object, find, ppmE=25, print=TRUE)
    standardGeneric("findneutral"))
setGeneric("findKalmanROI", function(object, ...)
    standardGeneric("findKalmanROI"))
setGeneric("findPeaks", function(object, ...) standardGeneric("findPeaks"))
setGeneric("findPeaks.centWave", function(object, ...)
    standardGeneric("findPeaks.centWave"))
setGeneric("findPeaks.addPredictedIsotopeFeatures", function(object, ...)
    standardGeneric("findPeaks.addPredictedIsotopeFeatures"))
setGeneric("findPeaks.centWaveWithPredictedIsotopeROIs", function(object, ...)
    standardGeneric("findPeaks.centWaveWithPredictedIsotopeROIs"))
setGeneric("findPeaks.massifquant", function(object, ...)
    standardGeneric("findPeaks.massifquant"))
setGeneric("findPeaks.matchedFilter", function(object, ...)
    standardGeneric("findPeaks.matchedFilter"))
setGeneric("findPeaks.MSW", function(object, ...)
    standardGeneric("findPeaks.MSW"))
setGeneric("findPeaks.MS1", function(object, ...)
    standardGeneric("findPeaks.MS1"))


## G
setGeneric("getEIC", function(object, ...) standardGeneric("getEIC"))
setGeneric("getMsnScan", function(object, ...) standardGeneric("getMsnScan"))
setGeneric("getPeaks", function(object, ...) standardGeneric("getPeaks"))
setGeneric("getScan", function(object, ...) standardGeneric("getScan"))
setGeneric("getSpec", function(object, ...) standardGeneric("getSpec"))
setGeneric("getXcmsRaw", function(object, ...) standardGeneric("getXcmsRaw"))
## There's too many "group" methods here...
setGeneric("group.density", function(object, ...) standardGeneric("group.density"))
setGeneric("group.mzClust", function(object, ...) standardGeneric("group.mzClust"))
setGeneric("group.nearest", function(object, ...) standardGeneric("group.nearest"))
setGeneric("group", function(object, ...) standardGeneric("group"))
setGeneric("groupidx", function(object) standardGeneric("groupidx"))
setGeneric("groupidx<-", function(object, value) standardGeneric("groupidx<-"))
setGeneric("groupnames", function(object, ...) standardGeneric("groupnames"))
setGeneric("groups", function(object) standardGeneric("groups"))
setGeneric("groups<-", function(object, value) standardGeneric("groups<-"))
setGeneric("groupval", function(object, ...) standardGeneric("groupval"))

## H
setGeneric("hasMSn", function(object, ...) standardGeneric("hasMSn"))

## I
setGeneric("image", function(x, ...) standardGeneric("image"))
setGeneric("isCentroided", function(object, ...) standardGeneric("isCentroided"))

## L
setGeneric("levelplot", function(x, data, ...) standardGeneric("levelplot"))
setGeneric("loadRaw", function(object, ...) standardGeneric("loadRaw"))

## M
setGeneric("makeacqNum", function(object, freq, start=1) standardGeneric("makeacqNum"))
setGeneric("mslevel", function(object, ...) standardGeneric("mslevel"))
setGeneric("mslevel<-", function(object, value) standardGeneric("mslevel<-"))
setGeneric("msnparent2ms", function(object, ...) standardGeneric("msnparent2ms"))
setGeneric("msn2ms", function(object, ...) standardGeneric("msn2ms"))
setGeneric("mzrange", function(object) standardGeneric("mzrange"))

## P
setGeneric("peaks<-", function(object, value) standardGeneric("peaks<-"))
setGeneric("peakTable", function(object, ...) standardGeneric("peakTable"))
setGeneric("plotChrom", function(object, ...) standardGeneric("plotChrom"))
setGeneric("plotEIC", function(object, ...) standardGeneric("plotEIC"))
setGeneric("plotPeaks", function(object, ...) standardGeneric("plotPeaks"))
setGeneric("plotRaw", function(object, ...) standardGeneric("plotRaw"))
setGeneric("plotrt", function(object, ...) standardGeneric("plotrt"))
setGeneric("plotScan", function(object, ...) standardGeneric("plotScan"))
setGeneric("plotSpec", function(object, ...) standardGeneric("plotSpec"))
setGeneric("plotSurf", function(object, ...) standardGeneric("plotSurf"))
setGeneric("plotTIC", function(object, ...) standardGeneric("plotTIC"))
setGeneric("plotTree", function(object, ...) standardGeneric("plotTree"))
setGeneric("present", function(object, class, minfrac) standardGeneric("present"))
setGeneric("profinfo", function(object) standardGeneric("profinfo"))
setGeneric("profinfo<-", function(object, value) standardGeneric("profinfo<-"))
setGeneric("profMedFilt", function(object, ...) standardGeneric("profMedFilt"))
setGeneric("profMethod", function(object) standardGeneric("profMethod"))
setGeneric("profMethod<-", function(object, value) standardGeneric("profMethod<-"))
setGeneric("profRange", function(object, ...) standardGeneric("profRange"))
setGeneric("profStep", function(object) standardGeneric("profStep"))
setGeneric("profStep<-", function(object, value) standardGeneric("profStep<-"))
setGeneric("profStepPad<-", function(object, value) standardGeneric("profStepPad<-"))
setGeneric("profMz", function(object) standardGeneric("profMz"))
setGeneric("progressCallback", function(object) standardGeneric("progressCallback"))
setGeneric("progressCallback<-", function(object, value) standardGeneric("progressCallback<-"))
setGeneric("progressInfoUpdate", function(object) standardGeneric("progressInfoUpdate"))

## R
setGeneric("rawEIC", function(object, ...) standardGeneric("rawEIC"))
setGeneric("rawMat", function(object, ...) standardGeneric("rawMat"))
setGeneric("rawMZ", function(object, ...) standardGeneric("rawMZ"))
setGeneric("retcor", function(object, ...) standardGeneric("retcor"))
setGeneric("retcor.peakgroups", function(object, ...) standardGeneric("retcor.peakgroups"))
setGeneric("retcor.obiwarp", function(object, ...) standardGeneric("retcor.obiwarp"))
setGeneric("revMz", function(object, ...) standardGeneric("revMz"))
setGeneric("rtrange", function(object) standardGeneric("rtrange"))

## S
setGeneric("sampclass", function(object) standardGeneric("sampclass"))
setGeneric("sampclass<-", function(object, value) standardGeneric("sampclass<-"))
setGeneric("sampnames", function(object) standardGeneric("sampnames"))
setGeneric("sampnames<-", function(object, value) standardGeneric("sampnames<-"))
setGeneric("scanrange", function(object, ...) standardGeneric("scanrange"))
setGeneric("scanrange<-", function(object, value) standardGeneric("scanrange<-"))
setGeneric("sortMz", function(object, ...) standardGeneric("sortMz"))
setGeneric("specDist", function(object, ...) standardGeneric("specDist"))
setGeneric("specDist.meanMZmatch",
           function(peakTable1, peakTable2, matchdist=1, matchrate=1,
                    mzabs=0.001, mzppm=10, symmetric=TRUE)
               standardGeneric("specDist.meanMZmatch"))
setGeneric("specDist.cosine",
           function(peakTable1, peakTable2, mzabs = 0.001, mzppm = 10,
                    mzExp = 0.6, intExp = 3, nPdiff = 2, nPmin = 8,
                    symmetric = FALSE)
               standardGeneric("specDist.cosine"))
setGeneric("specDist.peakCount",
           function(peakTable1, peakTable2, mzabs=0.001, mzppm=10,symmetric=FALSE)
               standardGeneric("specDist.peakCount"))
setGeneric("stitch", function(object, lockMass, ...) standardGeneric("stitch"))
setGeneric("stitch.xml", function(object, lockMass) standardGeneric("stitch.xml"))
setGeneric("stitch.netCDF", function(object, lockMass) standardGeneric("stitch.netCDF"))
setGeneric("stitch.netCDF.new", function(object, lockMass) standardGeneric("stitch.netCDF.new"))


## W
setGeneric("write.cdf", function(object, ...) standardGeneric("write.cdf"))
setGeneric("write.mzdata", function(object, ...) standardGeneric("write.mzdata"))
setGeneric("write.mzQuantML", function(object, ...) standardGeneric("write.mzQuantML"))

## X
setGeneric("xcmsSource", function(object, ...) standardGeneric("xcmsSource"))

