# Internal page for hidden aliases

For S4 methods that require a documentation entry but only clutter the
index.

## Usage

``` r
# S4 method for class 'XcmsExperiment'
show(object)

# S4 method for class 'XcmsExperimentHdf5'
show(object)

# S4 method for class 'XcmsExperimentHdf5,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'XcmsExperimentHdf5'
filterMsLevel(object, msLevel. = uniqueMsLevels(object))

# S4 method for class 'XcmsExperimentHdf5'
filterIsolationWindow(object, mz = numeric())

# S4 method for class 'XcmsExperimentHdf5'
filterRt(object, rt, msLevel. = uniqueMsLevels(object))

# S4 method for class 'XcmsExperimentHdf5'
filterMzRange(object, mz = numeric(), msLevel. = uniqueMsLevels(object))

# S4 method for class 'XcmsExperimentHdf5,Param'
findChromPeaks(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  add = FALSE,
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperimentHdf5'
dropChromPeaks(object, keepAdjustedRtime = FALSE)

# S4 method for class 'XcmsExperimentHdf5'
hasChromPeaks(object, msLevel = integer())

# S4 method for class 'XcmsExperimentHdf5'
chromPeaks(object) <- value

# S4 method for class 'XcmsExperimentHdf5'
chromPeaks(
  object,
  rt = numeric(),
  mz = numeric(),
  ppm = 0,
  msLevel = integer(),
  type = c("any", "within", "apex_within"),
  isFilledColumn = FALSE,
  columns = character(),
  bySample = FALSE
)

# S4 method for class 'XcmsExperimentHdf5'
chromPeakData(object) <- value

# S4 method for class 'XcmsExperimentHdf5,MergeNeighboringPeaksParam'
refineChromPeaks(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperimentHdf5,CleanPeaksParam'
refineChromPeaks(object, param = CleanPeaksParam(), msLevel = 1L)

# S4 method for class 'XcmsExperimentHdf5'
hasFilledChromPeaks(object)

# S4 method for class 'XcmsExperimentHdf5,ChromPeakAreaParam'
fillChromPeaks(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperimentHdf5'
dropFilledChromPeaks(object)

# S4 method for class 'XcmsExperimentHdf5'
manualChromPeaks(
  object,
  chromPeaks = matrix(numeric()),
  samples = seq_along(object),
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperimentHdf5,BetaDistributionParam'
chromPeakSummary(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperimentHdf5'
dropAdjustedRtime(object)

# S4 method for class 'XcmsExperimentHdf5'
hasFeatures(object, msLevel = integer())

# S4 method for class 'XcmsExperimentHdf5'
dropFeatureDefinitions(object, keepAdjustedRtime = FALSE)

# S4 method for class 'XcmsExperimentHdf5,Param'
groupChromPeaks(object, param, msLevel = 1L, add = FALSE)

# S4 method for class 'XcmsExperimentHdf5'
featureArea(
  object,
  mzmin = min,
  mzmax = max,
  rtmin = min,
  rtmax = max,
  features = character(),
  msLevel = 1L,
  minMzWidthPpm = 0
)

# S4 method for class 'XcmsExperimentHdf5'
featureDefinitions(object) <- value

# S4 method for class 'XcmsExperimentHdf5'
featureDefinitions(
  object,
  mz = numeric(),
  rt = numeric(),
  ppm = 0,
  type = c("any", "within", "apex_within"),
  msLevel = integer()
)

# S4 method for class 'XcmsExperimentHdf5'
featureValues(
  object,
  method = c("medret", "maxint", "sum"),
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = NA_real_,
  msLevel = integer()
)

# S4 method for class 'XcmsExperimentHdf5'
manualFeatures(object, peakIdx = list(), msLevel = 1L)

# S4 method for class 'XcmsExperimentHdf5'
chromatogram(
  object,
  rt = matrix(nrow = 0, ncol = 2),
  mz = matrix(nrow = 0, ncol = 2),
  aggregationFun = "sum",
  msLevel = 1L,
  chunkSize = 2L,
  isolationWindowTargetMz = NULL,
  return.type = c("XChromatograms", "MChromatograms"),
  include = character(),
  chromPeaks = c("apex_within", "any", "none"),
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperimentHdf5'
chromPeakSpectra(
  object,
  method = c("all", "closest_rt", "closest_mz", "largest_tic", "largest_bpi"),
  msLevel = 2L,
  expandRt = 0,
  expandMz = 0,
  ppm = 0,
  skipFilled = FALSE,
  peaks = character(),
  chromPeakColumns = c("rt", "mz"),
  return.type = c("Spectra", "List"),
  ...
)

# S4 method for class 'XcmsExperimentHdf5'
featureSpectra(
  object,
  msLevel = 2L,
  expandRt = 0,
  expandMz = 0,
  ppm = 0,
  skipFilled = FALSE,
  return.type = c("Spectra", "List"),
  features = character(),
  method = c("all", "closest_rt", "closest_mz", "largest_tic", "largest_bpi"),
  chromPeakColumns = c("rt", "mz"),
  featureColumns = c("rtmed", "mzmed"),
  ...
)

# S4 method for class 'XcmsExperimentHdf5'
featureChromatograms(
  object,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "max",
  features = character(),
  return.type = "XChromatograms",
  chunkSize = 2L,
  mzmin = min,
  mzmax = max,
  rtmin = min,
  rtmax = max,
  ...,
  progressbar = TRUE,
  BPPARAM = bpparam()
)

fixedRt(object)

fixedMz(object)

# S4 method for class 'MsFeatureData'
show(object)

# S4 method for class 'MsFeatureData'
hasAdjustedRtime(object)

# S4 method for class 'MsFeatureData'
hasFeatures(object, msLevel = integer())

# S4 method for class 'MsFeatureData'
hasChromPeaks(object, msLevel = integer())

# S4 method for class 'MsFeatureData'
adjustedRtime(object)

# S4 method for class 'MsFeatureData'
adjustedRtime(object) <- value

# S4 method for class 'MsFeatureData'
dropAdjustedRtime(object, rtraw)

# S4 method for class 'MsFeatureData'
featureDefinitions(object, msLevel = integer())

# S4 method for class 'MsFeatureData'
featureDefinitions(object) <- value

# S4 method for class 'MsFeatureData'
dropFeatureDefinitions(object, dropAdjustedRtime = FALSE)

# S4 method for class 'MsFeatureData'
chromPeaks(object)

# S4 method for class 'MsFeatureData'
chromPeaks(object) <- value

# S4 method for class 'MsFeatureData'
dropChromPeaks(object)

# S4 method for class 'MsFeatureData'
chromPeakData(object, columns = character())

# S4 method for class 'MsFeatureData'
chromPeakData(object) <- value

# S4 method for class 'OnDiskMSnExp'
hasAdjustedRtime(object)

# S4 method for class 'CentWaveParam'
ppm(object)

# S4 method for class 'CentWaveParam'
ppm(object) <- value

# S4 method for class 'CentWaveParam'
peakwidth(object)

# S4 method for class 'CentWaveParam'
peakwidth(object) <- value

# S4 method for class 'CentWaveParam'
snthresh(object)

# S4 method for class 'CentWaveParam'
snthresh(object) <- value

# S4 method for class 'CentWaveParam'
prefilter(object)

# S4 method for class 'CentWaveParam'
prefilter(object) <- value

# S4 method for class 'CentWaveParam'
mzCenterFun(object)

# S4 method for class 'CentWaveParam'
mzCenterFun(object) <- value

# S4 method for class 'CentWaveParam'
integrate(f)

# S4 method for class 'CentWaveParam'
integrate(object) <- value

# S4 method for class 'CentWaveParam'
mzdiff(object)

# S4 method for class 'CentWaveParam'
mzdiff(object) <- value

# S4 method for class 'CentWaveParam'
fitgauss(object)

# S4 method for class 'CentWaveParam'
fitgauss(object) <- value

# S4 method for class 'CentWaveParam'
noise(object)

# S4 method for class 'CentWaveParam'
noise(object) <- value

# S4 method for class 'CentWaveParam'
verboseColumns(object)

# S4 method for class 'CentWaveParam'
verboseColumns(object) <- value

# S4 method for class 'CentWaveParam'
roiList(object)

# S4 method for class 'CentWaveParam'
roiList(object) <- value

# S4 method for class 'CentWaveParam'
firstBaselineCheck(object)

# S4 method for class 'CentWaveParam'
firstBaselineCheck(object) <- value

# S4 method for class 'CentWaveParam'
roiScales(object)

# S4 method for class 'CentWaveParam'
roiScales(object) <- value

# S4 method for class 'MatchedFilterParam'
binSize(object)

# S4 method for class 'MatchedFilterParam'
binSize(object) <- value

# S4 method for class 'MatchedFilterParam'
impute(object)

# S4 method for class 'MatchedFilterParam'
impute(object) <- value

# S4 method for class 'MatchedFilterParam'
baseValue(object)

# S4 method for class 'MatchedFilterParam'
baseValue(object) <- value

# S4 method for class 'MatchedFilterParam'
distance(object)

# S4 method for class 'MatchedFilterParam'
distance(object) <- value

# S4 method for class 'MatchedFilterParam'
fwhm(object)

# S4 method for class 'MatchedFilterParam'
fwhm(object) <- value

# S4 method for class 'MatchedFilterParam'
sigma(object)

# S4 method for class 'MatchedFilterParam'
sigma(object) <- value

# S4 method for class 'MatchedFilterParam'
max(x)

# S4 method for class 'MatchedFilterParam'
max(object) <- value

# S4 method for class 'MatchedFilterParam'
snthresh(object)

# S4 method for class 'MatchedFilterParam'
snthresh(object) <- value

# S4 method for class 'MatchedFilterParam'
steps(object)

# S4 method for class 'MatchedFilterParam'
steps(object) <- value

# S4 method for class 'MatchedFilterParam'
mzdiff(object)

# S4 method for class 'MatchedFilterParam'
mzdiff(object) <- value

# S4 method for class 'MatchedFilterParam'
index(object)

# S4 method for class 'MatchedFilterParam'
index(object) <- value

# S4 method for class 'MassifquantParam'
ppm(object)

# S4 method for class 'MassifquantParam'
ppm(object) <- value

# S4 method for class 'MassifquantParam'
peakwidth(object)

# S4 method for class 'MassifquantParam'
peakwidth(object) <- value

# S4 method for class 'MassifquantParam'
snthresh(object)

# S4 method for class 'MassifquantParam'
snthresh(object) <- value

# S4 method for class 'MassifquantParam'
prefilter(object)

# S4 method for class 'MassifquantParam'
prefilter(object) <- value

# S4 method for class 'MassifquantParam'
mzCenterFun(object)

# S4 method for class 'MassifquantParam'
mzCenterFun(object) <- value

# S4 method for class 'MassifquantParam'
integrate(f)

# S4 method for class 'MassifquantParam'
integrate(object) <- value

# S4 method for class 'MassifquantParam'
mzdiff(object)

# S4 method for class 'MassifquantParam'
mzdiff(object) <- value

# S4 method for class 'MassifquantParam'
fitgauss(object)

# S4 method for class 'MassifquantParam'
fitgauss(object) <- value

# S4 method for class 'MassifquantParam'
noise(object)

# S4 method for class 'MassifquantParam'
noise(object) <- value

# S4 method for class 'MassifquantParam'
verboseColumns(object)

# S4 method for class 'MassifquantParam'
verboseColumns(object) <- value

# S4 method for class 'MassifquantParam'
criticalValue(object)

# S4 method for class 'MassifquantParam'
criticalValue(object) <- value

# S4 method for class 'MassifquantParam'
consecMissedLimit(object)

# S4 method for class 'MassifquantParam'
consecMissedLimit(object) <- value

# S4 method for class 'MassifquantParam'
unions(object)

# S4 method for class 'MassifquantParam'
unions(object) <- value

# S4 method for class 'MassifquantParam'
checkBack(object)

# S4 method for class 'MassifquantParam'
checkBack(object) <- value

# S4 method for class 'MassifquantParam'
withWave(object)

# S4 method for class 'MassifquantParam'
withWave(object) <- value

# S4 method for class 'MSWParam'
snthresh(object)

# S4 method for class 'MSWParam'
snthresh(object) <- value

# S4 method for class 'MSWParam'
verboseColumns(object)

# S4 method for class 'MSWParam'
verboseColumns(object) <- value

# S4 method for class 'MSWParam'
scales(object)

# S4 method for class 'MSWParam'
scales(object) <- value

# S4 method for class 'MSWParam'
nearbyPeak(object)

# S4 method for class 'MSWParam'
nearbyPeak(object) <- value

# S4 method for class 'MSWParam'
peakScaleRange(object)

# S4 method for class 'MSWParam'
peakScaleRange(object) <- value

# S4 method for class 'MSWParam'
ampTh(object)

# S4 method for class 'MSWParam'
ampTh(object) <- value

# S4 method for class 'MSWParam'
minNoiseLevel(object)

# S4 method for class 'MSWParam'
minNoiseLevel(object) <- value

# S4 method for class 'MSWParam'
ridgeLength(object)

# S4 method for class 'MSWParam'
ridgeLength(object) <- value

# S4 method for class 'MSWParam'
peakThr(object)

# S4 method for class 'MSWParam'
peakThr(object) <- value

# S4 method for class 'MSWParam'
tuneIn(object)

# S4 method for class 'MSWParam'
tuneIn(object) <- value

# S4 method for class 'MSWParam'
addParams(object)

# S4 method for class 'MSWParam'
addParams(object) <- value

# S4 method for class 'CentWavePredIsoParam'
snthreshIsoROIs(object)

# S4 method for class 'CentWavePredIsoParam'
snthreshIsoROIs(object) <- value

# S4 method for class 'CentWavePredIsoParam'
maxCharge(object)

# S4 method for class 'CentWavePredIsoParam'
maxCharge(object) <- value

# S4 method for class 'CentWavePredIsoParam'
maxIso(object)

# S4 method for class 'CentWavePredIsoParam'
maxIso(object) <- value

# S4 method for class 'CentWavePredIsoParam'
mzIntervalExtension(object)

# S4 method for class 'CentWavePredIsoParam'
mzIntervalExtension(object) <- value

# S4 method for class 'CentWavePredIsoParam'
polarity(object)

# S4 method for class 'CentWavePredIsoParam'
polarity(object) <- value

# S4 method for class 'PeakDensityParam'
sampleGroups(object)

# S4 method for class 'PeakDensityParam'
sampleGroups(object) <- value

# S4 method for class 'PeakDensityParam'
bw(object)

# S4 method for class 'PeakDensityParam'
bw(object) <- value

# S4 method for class 'PeakDensityParam'
minFraction(object)

# S4 method for class 'PeakDensityParam'
minFraction(object) <- value

# S4 method for class 'PeakDensityParam'
minSamples(object)

# S4 method for class 'PeakDensityParam'
minSamples(object) <- value

# S4 method for class 'PeakDensityParam'
binSize(object)

# S4 method for class 'PeakDensityParam'
binSize(object) <- value

# S4 method for class 'PeakDensityParam'
maxFeatures(object)

# S4 method for class 'PeakDensityParam'
maxFeatures(object) <- value

# S4 method for class 'PeakDensityParam'
ppm(object)

# S4 method for class 'MzClustParam'
sampleGroups(object)

# S4 method for class 'MzClustParam'
sampleGroups(object) <- value

# S4 method for class 'MzClustParam'
ppm(object)

# S4 method for class 'MzClustParam'
ppm(object) <- value

# S4 method for class 'MzClustParam'
absMz(object)

# S4 method for class 'MzClustParam'
absMz(object) <- value

# S4 method for class 'MzClustParam'
minFraction(object)

# S4 method for class 'MzClustParam'
minFraction(object) <- value

# S4 method for class 'MzClustParam'
minSamples(object)

# S4 method for class 'MzClustParam'
minSamples(object) <- value

# S4 method for class 'NearestPeaksParam'
sampleGroups(object)

# S4 method for class 'NearestPeaksParam'
sampleGroups(object) <- value

# S4 method for class 'NearestPeaksParam'
mzVsRtBalance(object)

# S4 method for class 'NearestPeaksParam'
mzVsRtBalance(object) <- value

# S4 method for class 'NearestPeaksParam'
absMz(object)

# S4 method for class 'NearestPeaksParam'
absMz(object) <- value

# S4 method for class 'NearestPeaksParam'
absRt(object)

# S4 method for class 'NearestPeaksParam'
absRt(object) <- value

# S4 method for class 'NearestPeaksParam'
kNN(object)

# S4 method for class 'NearestPeaksParam'
kNN(object) <- value

# S4 method for class 'PeakGroupsParam'
minFraction(object)

# S4 method for class 'PeakGroupsParam'
minFraction(object) <- value

# S4 method for class 'PeakGroupsParam'
extraPeaks(object)

# S4 method for class 'PeakGroupsParam'
extraPeaks(object) <- value

# S4 method for class 'PeakGroupsParam'
smooth(x)

# S4 method for class 'PeakGroupsParam'
smooth(object) <- value

# S4 method for class 'PeakGroupsParam'
span(object)

# S4 method for class 'PeakGroupsParam'
span(object) <- value

# S4 method for class 'PeakGroupsParam'
family(object)

# S4 method for class 'PeakGroupsParam'
family(object) <- value

# S4 method for class 'PeakGroupsParam'
peakGroupsMatrix(object)

# S4 method for class 'PeakGroupsParam'
peakGroupsMatrix(object) <- value

# S4 method for class 'PeakGroupsParam'
subset(x)

# S4 method for class 'PeakGroupsParam'
subset(object) <- value

# S4 method for class 'PeakGroupsParam'
subsetAdjust(object)

# S4 method for class 'PeakGroupsParam'
subsetAdjust(object) <- value

# S4 method for class 'ObiwarpParam'
binSize(object)

# S4 method for class 'ObiwarpParam'
centerSample(object)

# S4 method for class 'ObiwarpParam'
centerSample(object) <- value

# S4 method for class 'ObiwarpParam'
response(object)

# S4 method for class 'ObiwarpParam'
response(object) <- value

# S4 method for class 'ObiwarpParam'
distFun(object)

# S4 method for class 'ObiwarpParam'
distFun(object) <- value

# S4 method for class 'ObiwarpParam'
gapInit(object)

# S4 method for class 'ObiwarpParam'
gapInit(object) <- value

# S4 method for class 'ObiwarpParam'
gapExtend(object)

# S4 method for class 'ObiwarpParam'
gapExtend(object) <- value

# S4 method for class 'ObiwarpParam'
factorDiag(object)

# S4 method for class 'ObiwarpParam'
factorDiag(object) <- value

# S4 method for class 'ObiwarpParam'
factorGap(object)

# S4 method for class 'ObiwarpParam'
factorGap(object) <- value

# S4 method for class 'ObiwarpParam'
localAlignment(object)

# S4 method for class 'ObiwarpParam'
localAlignment(object) <- value

# S4 method for class 'ObiwarpParam'
initPenalty(object)

# S4 method for class 'ObiwarpParam'
initPenalty(object) <- value

# S4 method for class 'ObiwarpParam'
subset(x)

# S4 method for class 'ObiwarpParam'
subset(object) <- value

# S4 method for class 'ObiwarpParam'
subsetAdjust(object)

# S4 method for class 'ObiwarpParam'
subsetAdjust(object) <- value

# S4 method for class 'FillChromPeaksParam'
expandMz(object)

# S4 method for class 'FillChromPeaksParam'
expandMz(object) <- value

# S4 method for class 'FillChromPeaksParam'
expandRt(object)

# S4 method for class 'FillChromPeaksParam'
expandRt(object) <- value

# S4 method for class 'FillChromPeaksParam'
ppm(object)

# S4 method for class 'FillChromPeaksParam'
ppm(object) <- value

# S4 method for class 'ProcessHistory'
show(object)

# S4 method for class 'XProcessHistory'
show(object)

# S4 method for class 'XCMSnExp'
show(object)

# S4 method for class 'XcmsResult,PeakGroupsParam'
adjustRtimePeakGroups(object, param = PeakGroupsParam(), msLevel = 1L)

# S4 method for class 'XChromatogram'
show(object)

# S4 method for class 'XChromatograms'
show(object)

# S4 method for class 'xcmsEIC'
show(object)

# S4 method for class 'xcmsFragments'
show(object)

# S4 method for class 'xcmsPeaks'
show(object)

# S4 method for class 'xcmsRaw'
show(object)

# S4 method for class 'xcmsSet'
show(object)
```

## Value

Not applicable
