# Grouping FTICR-MS data with xcms

## Introduction

This document describes how to use
*[xcms](https://bioconductor.org/packages/3.24/xcms)* for the analysis
of direct injection mass spec data, including peak detection,
calibration and correspondence (grouping of peaks across samples).

## Peak detection

Prior to any other analysis step, peaks have to be identified in the
mass spec data. In contrast to the typical metabolomics workflow, in
which peaks are identified in the chromatographic (time) dimension, in
direct injection mass spec data sets peaks are identified in the m/z
dimension. *[xcms](https://bioconductor.org/packages/3.24/xcms)* uses
functionality from the *MassSpecWavelet* package to identify such peaks.

Below we load the required packages. For information on the parallel
processing setup please see the *BiocParallel* vignette.

[`library`](https://rdrr.io/r/base/library.html)`(`[`MSnbase`](https://lgatto.github.io/MSnbase)`)`` `[`library`](https://rdrr.io/r/base/library.html)`(`[`xcms`](https://github.com/sneumann/xcms)`)`` `[`library`](https://rdrr.io/r/base/library.html)`(`[`MassSpecWavelet`](https://github.com/zeehio/MassSpecWavelet)`)`` `[`library`](https://rdrr.io/r/base/library.html)`(`[`MsDataHub`](https://rformassspectrometry.github.io/MsDataHub)`)`` `[`register`](https://rdrr.io/pkg/BiocParallel/man/register.html)`(`[`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)`(``)``)`

In this documentation we use an example data set from the
`r Biocpkg("MsDataHub")` package. Assuming that
*[MsDataHub](https://bioconductor.org/packages/3.24/MsDataHub)* is
installed, it will obtain the files from
<https://doi.org/10.5281/zenodo.18494293> and load the data set. We
create also a `data.frame` describing the experimental setup based on
the file names.

`## We're using 2 samples per condition`` ``## from https://doi.org/10.5281/zenodo.18494293`` ``mzML_files`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`` `` `[`HAM004_641fE_14.11.07..Exp1.extracted.mzML`](https://rformassspectrometry.github.io/MsDataHub/reference/FTICR.html)`(``)``,`` `` `[`HAM004_641fE_14.11.07..Exp2.extracted.mzML`](https://rformassspectrometry.github.io/MsDataHub/reference/FTICR.html)`(``)``,`` `` `[`HAM005_641fE_14.11.07..Exp1.extracted.mzML`](https://rformassspectrometry.github.io/MsDataHub/reference/FTICR.html)`(``)``,`` `` `[`HAM005_641fE_14.11.07..Exp2.extracted.mzML`](https://rformassspectrometry.github.io/MsDataHub/reference/FTICR.html)`(``)`` ``)`` `` ``## Create a data.frame assigning samples to sample groups, i.e. ham4 and ham5.`` ``grp`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(``"ham4"``, ``"ham4"``, ``"ham5"``, ``"ham5"``)`` ``pd`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``filename ``=`` `[`basename`](https://rdrr.io/r/base/basename.html)`(``mzML_files``)``, sample_group ``=`` ``grp``)`` `` ``## Load the data.`` ``ham_raw`` ``<-`` `[`readMSData`](https://lgatto.github.io/MSnbase/reference/readMSData.html)`(``files ``=`` ``mzML_files``,`` `` pdata ``=`` `[`AnnotatedDataFrame`](https://rdrr.io/pkg/Biobase/man/class.AnnotatedDataFrame.html)`(``pd``)``,`` `` mode ``=`` ``"onDisk"``)`

    ## Warning in .testReadMSDataInput(environment()): Reading different file formats in.
    ## This is untested and you are welcome to try it out.
    ## Please report back!

The data files are from *direct injection* mass spectrometry
experiments, i.e. we have only a single spectrum available for each
sample and no retention times.

`## Only a single spectrum with an *artificial* retention time is available`` ``## for each sample`` `[`rtime`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)`(``ham_raw``)`

    ## F1.S1 F2.S1 F3.S1 F4.S1 
    ##    -1    -1    -1    -1

Peaks are identified within each spectrum using the *mass spec wavelet*
method.

`## Define the parameters for the peak detection`` ``msw`` ``<-`` `[`MSWParam`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md)`(``scales ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``4``, ``9``)``, nearbyPeak ``=`` ``TRUE``, winSize.noise ``=`` ``500``,`` `` SNR.method ``=`` ``"data.mean"``, snthresh ``=`` ``10``)`` `` ``ham_prep`` ``<-`` `[`findChromPeaks`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)`(``ham_raw``, param ``=`` ``msw``)`` `` `[`head`](https://rdrr.io/r/utils/head.html)`(`[`chromPeaks`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``ham_prep``)``)`

    ##            mz    mzmin    mzmax rt rtmin rtmax    into     maxo       sn intf
    ## CP01 403.2367 403.2279 403.2447 -1    -1    -1 4735258 372259.4 22.97534   NA
    ## CP02 409.1845 409.1747 409.1936 -1    -1    -1 4158404 310572.1 20.61382   NA
    ## CP03 413.2677 413.2585 413.2769 -1    -1    -1 6099006 435462.6 27.21723   NA
    ## CP04 423.2363 423.2266 423.2459 -1    -1    -1 2708391 174252.7 14.74527   NA
    ## CP05 427.2681 427.2574 427.2779 -1    -1    -1 6302089 461385.6 32.50050   NA
    ## CP06 437.2375 437.2254 437.2488 -1    -1    -1 7523070 517917.6 34.37645   NA
    ##           maxf sample
    ## CP01  814693.1      1
    ## CP02  732119.9      1
    ## CP03 1018994.8      1
    ## CP04  435858.5      1
    ## CP05 1125644.3      1
    ## CP06 1282906.5      1

## Calibration

The `calibrate` method can be used to correct the m/z values of
identified peaks. The currently implemented method requires identified
peaks and a list of m/z values for known calibrants. The identified
peaks m/z values are then adjusted based on the differences between the
calibrants’ m/z values and the m/z values of the closest peaks (within a
user defined permitted maximal distance). Note that this method does
presently only calibrate identified peaks, but not the original m/z
values in the spectra.

Below we demonstrate the `calibrate` method on one of the data files
with artificially defined calibration m/z values. We first subset the
data set to the first data file, extract the m/z values of 3 peaks and
modify the values slightly.

`## Subset to the first file.`` ``first_file`` ``<-`` `[`filterFile`](https://lgatto.github.io/MSnbase/reference/MSnExp-class.html)`(``ham_prep``, file ``=`` ``1``)`` `` ``## Extract 3 m/z values`` ``calib_mz`` ``<-`` `[`chromPeaks`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``first_file``)``[`[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``4``, ``7``)``, ``"mz"``]`` ``calib_mz`` ``<-`` ``calib_mz`` ``+`` ``0.00001`` ``*`` `[`runif`](https://rdrr.io/r/stats/Uniform.html)`(``1``, ``0``, ``0.4``)`` ``*`` ``calib_mz`` ``+`` ``0.0001`

Next we calibrate the data set using the previously defined *artificial*
calibrants. We are using the `"edgeshift"` method for calibration that
adjusts all peaks within the range of the m/z values of the calibrants
using a linear interpolation and shifts all chromatographic peaks
outside of that range by a constant factor (the difference between the
lowest respectively largest calibrant m/z with the closest peak’s m/z).
Note that in a *real* use case, the m/z values would obviously represent
known m/z of calibrants and would not be defined on the actual data.

`## Set-up the parameter class for the calibration`` ``prm`` ``<-`` `[`CalibrantMassParam`](https://sneumann.github.io/xcms/reference/calibrate-calibrant-mass.md)`(``mz ``=`` ``calib_mz``, method ``=`` ``"edgeshift"``,`` `` mzabs ``=`` ``0.0001``, mzppm ``=`` ``5``)`` ``first_file_calibrated`` ``<-`` `[`calibrate`](https://sneumann.github.io/xcms/reference/calibrate.md)`(``first_file``, param ``=`` ``prm``)`

To evaluate the calibration we plot below the difference between the
adjusted and raw m/z values (y-axis) against the raw m/z values.

`diffs`` ``<-`` `[`chromPeaks`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``first_file_calibrated``)``[``, ``"mz"``]`` ``-`` `` `[`chromPeaks`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``first_file``)``[``, ``"mz"``]`` `` `[`plot`](https://rdrr.io/r/base/plot.html)`(``x ``=`` `[`chromPeaks`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``first_file``)``[``, ``"mz"``]``, xlab ``=`` `[`expression`](https://rdrr.io/r/base/expression.html)`(``m``/``z``[``raw``]``)``,`` `` y ``=`` ``diffs``, ylab ``=`` `[`expression`](https://rdrr.io/r/base/expression.html)`(``m``/``z``[``calibrated``]`` ``-`` ``m``/``z``[``raw``]``)``)`

![](xcms-direct-injection_files/figure-html/calibrationresult-1.png)

## Correspondence

Correspondence aims to group peaks across samples to define the
*features* (ions with the same m/z values). Peaks from single spectrum,
direct injection MS experiments can be grouped with the *MZclust*
method. Below we perform the correspondence analysis with the
`groupChromPeaks` method using default settings.

`## Using default settings but define sample group assignment`` ``mzc_prm`` ``<-`` `[`MzClustParam`](https://sneumann.github.io/xcms/reference/groupChromPeaks.md)`(``sampleGroups ``=`` ``ham_prep``$``sample_group``)`` ``ham_prep`` ``<-`` `[`groupChromPeaks`](https://sneumann.github.io/xcms/reference/groupChromPeaks.md)`(``ham_prep``, param ``=`` ``mzc_prm``)`

Getting an overview of the performed processings:

`ham_prep`

    ## MSn experiment data ("XCMSnExp")
    ## Object size in memory: 0.04 Mb
    ## - - - Spectra data - - -
    ##  MS level(s): 1 
    ##  Number of spectra: 4 
    ##  MSn retention times: -1:59 - -1:59 minutes
    ## - - - Processing information - - -
    ## Data loaded [Tue Sep  1 08:32:55 2026] 
    ##  MSnbase version: 2.39.5 
    ## - - - Meta data  - - -
    ## phenoData
    ##   rowNames: 1 2 3 4
    ##   varLabels: filename sample_group
    ##   varMetadata: labelDescription
    ## Loaded from:
    ##   [1] 187c3767766c_10386...  [4] 187c44b3d682_10392
    ##   Use 'fileNames(.)' to see all files.
    ## protocolData: none
    ## featureData
    ##   featureNames: F1.S1 F2.S1 F3.S1 F4.S1
    ##   fvarLabels: fileIdx spIdx ... spectrum (36 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## - - - xcms preprocessing - - -
    ## Chromatographic peak detection:
    ##  method: MSW 
    ##  38 peaks identified in 4 samples.
    ##  On average 9.5 chromatographic peaks per sample.
    ## Correspondence:
    ##  method: mzClust 
    ##  20 features identified.
    ##  Median mz range of features: 9.1553e-05
    ##  Median rt range of features: 0

The peak group information, i.e. the *feature* definitions can be
accessed with the `featureDefinitions` method.

[`featureDefinitions`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``ham_prep``)`

    ## DataFrame with 20 rows and 10 columns
    ##          mzmed     mzmin     mzmax     rtmed     rtmin     rtmax    npeaks
    ##      <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    ## FT01   402.285   402.285   402.286        -1        -1        -1         2
    ## FT02   403.237   403.237   403.237        -1        -1        -1         4
    ## FT03   405.109   405.109   405.109        -1        -1        -1         2
    ## FT04   409.184   409.184   409.185        -1        -1        -1         2
    ## FT05   410.144   410.144   410.145        -1        -1        -1         2
    ## ...        ...       ...       ...       ...       ...       ...       ...
    ## FT16   437.238   437.238   437.238        -1        -1        -1         2
    ## FT17   438.240   438.240   438.240        -1        -1        -1         2
    ## FT18   439.151   439.151   439.151        -1        -1        -1         2
    ## FT19   441.130   441.130   441.131        -1        -1        -1         2
    ## FT20   445.293   445.292   445.293        -1        -1        -1         2
    ##           ham4      ham5     peakidx
    ##      <numeric> <numeric>      <list>
    ## FT01         0         2       16,28
    ## FT02         2         2 17,29,1,...
    ## FT03         0         2       18,30
    ## FT04         2         0        10,2
    ## FT05         0         2       19,31
    ## ...        ...       ...         ...
    ## FT16         2         0        6,13
    ## FT17         2         0        7,14
    ## FT18         0         2       26,37
    ## FT19         0         2       38,27
    ## FT20         2         0        15,8

Plotting the raw data for direct injection samples involves a little
more processing than for LC/GC-MS data in which we can simply use the
`chromatogram` method to extract the data. Below we extract the
m/z-intensity pairs for the peaks associated with the first feature. We
thus first identify the peaks for that feature and define their m/z
values range. Using this range we can subsequently use the `filterMz`
function to sub-set the full data set to the signal associated with the
feature’s peaks. On that object we can then call the `mz` and
`intensity` functions to extract the data.

`## Get the peaks belonging to the first feature`` ``pks`` ``<-`` `[`chromPeaks`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``ham_prep``)``[`[`featureDefinitions`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)`(``ham_prep``)``$``peakidx``[[``1``]``]``, ``]`` `` ``## Define the m/z range`` ``mzr`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`min`](https://rdrr.io/r/base/Extremes.html)`(``pks``[``, ``"mzmin"``]``)`` ``-`` ``0.001``, `[`max`](https://rdrr.io/r/base/Extremes.html)`(``pks``[``, ``"mzmax"``]``)`` ``+`` ``0.001``)`` `` ``## Subset the object to the m/z range`` ``ham_prep_sub`` ``<-`` `[`filterMz`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)`(``ham_prep``, mz ``=`` ``mzr``)`` `` ``## Extract the mz and intensity values`` ``mzs`` ``<-`` `[`mz`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)`(``ham_prep_sub``, bySample ``=`` ``TRUE``)`` ``ints`` ``<-`` `[`intensity`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)`(``ham_prep_sub``, bySample ``=`` ``TRUE``)`` `` ``## Plot the data`` `[`plot`](https://rdrr.io/r/base/plot.html)`(``3``, ``3``, pch ``=`` ``NA``, xlim ``=`` `[`range`](https://rdrr.io/r/base/range.html)`(``mzs``)``, ylim ``=`` `[`range`](https://rdrr.io/r/base/range.html)`(``ints``)``, main ``=`` ``"FT01"``,`` `` xlab ``=`` ``"m/z"``, ylab ``=`` ``"intensity"``)`` ``## Define colors`` ``cols`` ``<-`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``"#ff000080"``, `[`length`](https://rdrr.io/r/base/length.html)`(``mzs``)``)`` ``cols``[``ham_prep_sub``$``sample_group`` ``==`` ``"ham5"``]`` ``<-`` ``"#0000ff80"`` ``tmp`` ``<-`` `[`mapply`](https://rdrr.io/r/base/mapply.html)`(``mzs``, ``ints``, ``cols``, FUN ``=`` ``function``(``x``, ``y``, ``col``)`` ``{`` `` `[`points`](https://rdrr.io/r/graphics/points.html)`(``x``, ``y``, col ``=`` ``col``, type ``=`` ``"l"``)`` ``}``)`

![](xcms-direct-injection_files/figure-html/feature1-1.png)

To access the actual intensity values of each feature in each sample the
`featureValue` method can be used. The setting `value = "into"` tells
the function to return the integrated signal for each peak (one
representative peak) per sample.

`feat_vals`` ``<-`` `[`featureValues`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)`(``ham_prep``, value ``=`` ``"into"``)`` `[`head`](https://rdrr.io/r/utils/head.html)`(``feat_vals``)`

    ##      187c3767766c_10386 187c6a5a17f9_10387 187c1e5daedb_10391
    ## FT01                 NA                 NA            4095293
    ## FT02            4735258            6202418            4811391
    ## FT03                 NA                 NA            2982453
    ## FT04            4158404            5004546                 NA
    ## FT05                 NA                 NA            2872023
    ## FT06            6099006            4950642                 NA
    ##      187c44b3d682_10392
    ## FT01            4804763
    ## FT02            2581183
    ## FT03            2268984
    ## FT04                 NA
    ## FT05            2133219
    ## FT06                 NA

`NA` is reported for features in samples for which no peak was
identified at the feature’s m/z value. In some instances there might
still be a signal at the feature’s position in the raw data files, but
the peak detection failed to identify a peak. For these cases signal can
be recovered using the `fillChromPeaks` method that integrates all raw
signal at the feature’s location. If there is no signal at that location
an `NA` is reported.

`ham_prep`` ``<-`` `[`fillChromPeaks`](https://sneumann.github.io/xcms/reference/fillChromPeaks.md)`(``ham_prep``, param ``=`` `[`FillChromPeaksParam`](https://sneumann.github.io/xcms/reference/fillChromPeaks.md)`(``)``)`` `` `[`head`](https://rdrr.io/r/utils/head.html)`(`[`featureValues`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)`(``ham_prep``, value ``=`` ``"into"``)``)`

    ##      187c3767766c_10386 187c6a5a17f9_10387 187c1e5daedb_10391
    ## FT01           768754.0          1230140.4            4095293
    ## FT02          4735257.5          6202417.6            4811391
    ## FT03           652566.6           374109.9            2982453
    ## FT04          4158404.5          5004546.3            1221031
    ## FT05           652201.1           403448.4            2872023
    ## FT06          6099006.3          4950641.7            1573988
    ##      187c44b3d682_10392
    ## FT01          4804762.5
    ## FT02          2581183.1
    ## FT03          2268984.5
    ## FT04          1241294.4
    ## FT05          2133219.4
    ## FT06           977694.5

## Further analysis

Further analysis, i.e. detection of features/metabolites with
significantly different abundances, or PCA analyses can be performed on
the feature matrix using functionality from other R packages, such as
*[limma](https://bioconductor.org/packages/3.24/limma)*.

## Session information

[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] MsDataHub_1.13.1       MassSpecWavelet_1.79.2 xcms_4.11.3           
    ##  [4] BiocParallel_1.47.0    MSnbase_2.39.5         S4Vectors_0.51.9      
    ##  [7] Biobase_2.73.2         BiocGenerics_0.59.12   generics_0.1.4        
    ## [10] mzR_2.47.0             Rcpp_1.1.2             BiocStyle_2.41.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
    ##   [3] MultiAssayExperiment_1.39.0 magrittr_2.0.5             
    ##   [5] farver_2.1.2                MALDIquant_1.22.3          
    ##   [7] rmarkdown_2.31              fs_2.1.0                   
    ##   [9] ragg_1.5.2                  vctrs_0.7.3                
    ##  [11] memoise_2.0.1               BiocBaseUtils_1.15.1       
    ##  [13] htmltools_0.5.9             S4Arrays_1.13.0            
    ##  [15] progress_1.2.3              AnnotationHub_4.3.2        
    ##  [17] curl_8.0.0                  signal_1.8-1               
    ##  [19] SparseArray_1.13.2          mzID_1.51.0                
    ##  [21] sass_0.4.10                 bslib_0.12.0               
    ##  [23] htmlwidgets_1.6.4           desc_1.4.3                 
    ##  [25] plyr_1.8.9                  httr2_1.3.0                
    ##  [27] impute_1.87.0               cachem_1.1.0               
    ##  [29] igraph_2.3.3                lifecycle_1.0.5            
    ##  [31] iterators_1.0.14            pkgconfig_2.0.3            
    ##  [33] Matrix_1.7-6                R6_2.6.1                   
    ##  [35] fastmap_1.2.0               MatrixGenerics_1.25.0      
    ##  [37] clue_0.3-68                 digest_0.6.39              
    ##  [39] pcaMethods_2.5.0            AnnotationDbi_1.75.2       
    ##  [41] ExperimentHub_3.3.2         textshaping_1.0.5          
    ##  [43] GenomicRanges_1.65.1        RSQLite_3.53.3             
    ##  [45] filelock_1.0.3              Spectra_1.23.3             
    ##  [47] httr_1.4.8                  abind_1.4-8                
    ##  [49] compiler_4.6.1              withr_3.0.3                
    ##  [51] bit64_4.8.4                 doParallel_1.0.17          
    ##  [53] S7_0.2.2                    PTMods_1.1.0               
    ##  [55] DBI_1.3.0                   Chromatograms_1.3.3        
    ##  [57] MASS_7.3-66                 MsExperiment_1.15.0        
    ##  [59] rappdirs_0.3.4              DelayedArray_0.39.6        
    ##  [61] tools_4.6.1                 PSMatch_1.17.0             
    ##  [63] otel_0.2.0                  glue_1.8.1                 
    ##  [65] QFeatures_1.23.1            grid_4.6.1                 
    ##  [67] cluster_2.1.8.3             reshape2_1.4.5             
    ##  [69] gtable_0.3.6                preprocessCore_1.75.0      
    ##  [71] tidyr_1.3.2                 data.table_1.18.6.1        
    ##  [73] hms_1.1.4                   MetaboCoreUtils_1.21.1     
    ##  [75] XVector_0.53.0              BiocVersion_3.24.0         
    ##  [77] foreach_1.5.2               pillar_1.11.1              
    ##  [79] stringr_1.6.0               limma_3.69.4               
    ##  [81] dplyr_1.2.1                 BiocFileCache_3.3.0        
    ##  [83] lattice_0.23-1              bit_4.6.0                  
    ##  [85] tidyselect_1.2.1            Biostrings_2.81.6          
    ##  [87] knitr_1.51                  bookdown_0.48              
    ##  [89] IRanges_2.47.5              Seqinfo_1.3.2              
    ##  [91] ProtGenerics_1.45.0         SummarizedExperiment_1.43.0
    ##  [93] xfun_0.60                   statmod_1.5.2              
    ##  [95] matrixStats_1.5.0           stringi_1.8.9              
    ##  [97] lazyeval_0.2.3              yaml_2.3.12                
    ##  [99] evaluate_1.0.5              codetools_0.2-20           
    ## [101] MsCoreUtils_1.25.4          tibble_3.3.1               
    ## [103] BiocManager_1.30.27         cli_3.6.6                  
    ## [105] affyio_1.83.0               systemfonts_1.3.2          
    ## [107] jquerylib_0.1.4             dbplyr_2.6.0               
    ## [109] png_0.1-9                   XML_3.99-0.24              
    ## [111] parallel_4.6.1              pkgdown_2.2.1.9000         
    ## [113] ggplot2_4.0.3               blob_1.3.0                 
    ## [115] prettyunits_1.2.0           AnnotationFilter_1.37.0    
    ## [117] MsFeatures_1.21.0           scales_1.4.0               
    ## [119] affy_1.91.0                 ncdf4_1.24                 
    ## [121] purrr_1.2.2                 crayon_1.5.3               
    ## [123] rlang_1.3.0                 KEGGREST_1.53.6            
    ## [125] vsn_3.81.0
