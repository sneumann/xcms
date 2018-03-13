## unit tests will not be done if RUnit is not available
if(require("RUnit", quietly=TRUE)) {

    ## --- Setup ---
    
    pkg <- "xcms" # <-- Change to package name!
    if(Sys.getenv("RCMDCHECK") == "FALSE") {
        ## Path to unit tests for standalone running under Makefile (not R CMD check)
        ## PKG/tests/../inst/unitTests
        path <- file.path(getwd(), "..", "inst", "unitTests")
    } else {
        ## Path to unit tests for R CMD check
        ## PKG.Rcheck/tests/../PKG/unitTests
        path <- system.file(package=pkg, "unitTests")
    }
    cat("\nRunning unit tests\n")
    print(list(pkg=pkg, getwd=getwd(), pathToUnitTests=path))

    library(package=pkg, character.only=TRUE)
    library(package="faahKO", character.only=TRUE)

    attr(faahko, "filepaths") <- sapply(
        as.list(basename(attr(faahko, "filepaths"))),
        function(x) system.file("cdf", if (length(grep("ko",x)) > 0) "KO" else  "WT" ,x, package = "faahKO"))

    library(BiocParallel)
    if (.Platform$OS.type == "unix") {
        prm <- MulticoreParam()
    } else {
        prm <- SnowParam()
    }
    register(bpstart(prm))
    
    ## Create some objects we can re-use in different tests:
    ## Needed in runit.XCMSnExp.R
    faahko_3_files <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                        system.file('cdf/KO/ko16.CDF', package = "faahKO"),
                        system.file('cdf/KO/ko18.CDF', package = "faahKO"))
    
    ## An xcmsRaw for the first file:
    faahko_xr_1 <- xcmsRaw(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
                           profstep = 0)
    faahko_od <- readMSData(faahko_3_files, mode = "onDisk")

    faahko_xod <- findChromPeaks(faahko_od, param = CentWaveParam(noise = 10000,
                                                                  snthresh = 40))
    faahko_xs <- xcmsSet(faahko_3_files, profparam = list(step = 0),
                         method = "centWave", noise = 10000, snthresh = 40)

    ## faahko_xod <- findChromPeaks(faahko_od, param = CentWaveParam(noise = 5000))
    ## faahko_xs <- xcmsSet(faahko_3_files, profparam = list(step = 0),
    ##                      method = "centWave", noise = 5000)
    ## Doing also the retention time correction etc
    od_x <- faahko_od
    xod_x <- faahko_xod
    xod_xg <- groupChromPeaks(
        xod_x, param = PeakDensityParam(
                   sampleGroups = rep(1, length(fileNames(xod_x)))))
    xod_xgr <- adjustRtime(xod_xg, param = PeakGroupsParam(span = 0.4))
    xod_xgrg <- groupChromPeaks(
        xod_xgr, param = PeakDensityParam(
                     sampleGroups = rep(1, length(fileNames(xod_x)))))

    xod_r <- adjustRtime(as(od_x, "XCMSnExp"), param = ObiwarpParam())
    
    faahko_grouped_filled <- fillPeaks(group(faahko))
    faahko_grouped_retcor_filled <- fillPeaks(group(retcor(group(
        updateObject(faahko)))))
    
    ## microtofq
    library(msdata)
    microtofq_fs <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
                      system.file("microtofq/MM8.mzML", package = "msdata"))
    microtofq_xr <- xcmsRaw(microtofq_fs[1], profstep = 0)
    microtofq_od <- readMSData(microtofq_fs, mode = "onDisk")

    ## Direct injection data:
    fticrf <- list.files(system.file("fticr", package = "msdata"),
                         recursive = TRUE, full.names = TRUE)
    fticr <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
    fticr_xod <- findChromPeaks(fticr, MSWParam(scales = c(1, 7),
                                                peakThr = 80000, ampTh = 0.005,
                                                SNR.method = "data.mean",
                                                winSize.noise = 500))
    fticr_xs <- xcmsSet(method="MSW", files=fticrf[1:2], scales=c(1,7),
                        SNR.method='data.mean' , winSize.noise=500,
                        peakThr=80000,  amp.Th=0.005)
    
    ## microtofq_xod <- findChromPeaks(microtofq_od, param = MSWParam())
    ## If desired, load the name space to allow testing of private functions
    ## if (is.element(pkg, loadedNamespaces()))
    ##     attach(loadNamespace(pkg), name=paste("namespace", pkg, sep=":"), pos=3)
    ##
    ## or simply call PKG:::myPrivateFunction() in tests

    ## --- Testing ---

    ## Define tests
    testSuite <- defineTestSuite(name=paste(pkg, "unit testing"),
                                 dirs=path)
    ## Run
    tests <- runTestSuite(testSuite)

    ## Default report name
    pathReport <- file.path(path, "report")

    ## Report to stdout; disable text files
    cat("------------------- UNIT TEST SUMMARY ---------------------\n\n")
    printTextProtocol(tests, showDetails=FALSE)
    ## printTextProtocol(tests, showDetails=FALSE,
    ##                   fileName=paste(pathReport, "Summary.txt", sep=""))
    ## printTextProtocol(tests, showDetails=TRUE,
    ##                   fileName=paste(pathReport, ".txt", sep=""))

    ## ## Report to HTML file
    ## printHTMLProtocol(tests, fileName=paste(pathReport, ".html", sep=""))

    ## Return stop() to cause R CMD check stop in case of
    ##  - failures i.e. FALSE to unit tests or
    ##  - errors i.e. R errors
    tmp <- getErrors(tests)
    if(tmp$nFail > 0 | tmp$nErr > 0) {
        stop(paste("\n\nunit testing failed (#test failures: ", tmp$nFail,
                   ", #R errors: ",  tmp$nErr, ")\n\n", sep=""))
    }
} else {
    warning("cannot run unit tests -- package RUnit is not available")
}

