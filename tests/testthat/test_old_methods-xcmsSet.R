test_that("diffreport etc works", {
    g <- group(faahko)
    f <- fillPeaks(g)
    d <- diffreport(f)
    ## Fake xcmsSet with 1 sample in class 1
    sampclass(f) <- c(1,2,2,3,3,3,3,3,3,3,3,3)
    d <- diffreport(f, class1=1, class2=2)

    ## Test anova
    d <- diffreport(f)
    sampclass(f) <- rep(1:4,3)
    d <- diffreport(f)
    d$anova
})

test_that("findPeaks.MSW works", {
    ## We do expect an error if we have multiple spectra (issue #237)
    expect_error(findPeaks.MSW(faahko_xr_1))
})

test_that("fillPeaks,xcmsSet and filled flag works", {
    xsg <- group(faahko)
    xsgf <- fillPeaks(xsg, method = "chrom")

    expect_equal(nrow(peaks(xsg)) + length(xsgf@filled), nrow(peaks(xsgf)))
})

test_that("fillPeaks,xcmsSet columns works", {
    xsg <- group(faahko)
    xsg <- group(faahko_xs)
    peaks(xsg) <- cbind(peaks(xsg), anotherColumn=4711)

    oldCnames <- colnames(peaks(xsg))
    xsgf <- fillPeaks(xsg) # parallel disabled: , nSlaves=2)

    newCnames <- colnames(peaks(xsgf))
    expect_equal(oldCnames, newCnames)

    ## Check dims if nothing to do
    oldDims <- dim(peaks(xsgf))
    xsgf2 <- fillPeaks(xsgf) # parallel disabled: , nSlaves=2)
    newDims <- dim(peaks(xsgf2))
    expect_equal(oldDims, newDims)

    ## Case where only some samples have NA values
    xsg <- group(faahko_xs, minfrac=1)
    xsgf <- fillPeaks(xsg) # parallel disabled: , nSlaves=2)
    sampclass(xsgf) <- c(rep("KO", 1), rep("WT", 2))
    xsgf <- group(xsgf, minfrac=1)
    xsgf <- fillPeaks(xsgf) # parallel disabled: , nSlaves=2)
})

test_that("getEIC,xcmsSet works", {
    ## xset <- fillPeaks(group(faahko))
    xset <- faahko_grouped_filled
    ## xset <- faahko_grouped_filled
    e <- getEIC(xset, sampleidx = c(1,2), groupidx = c(1,2), rtrange=200)
    expect_equal(sampnames(e), c("ko15", "ko16"))
    ## plot(e)
    ## Reproduce issue #92
    e <- getEIC(xset, sampleidx = c(5, 9), groupidx = c(1, 2), rtrange = 200)
    expect_equal(sampnames(e), sampnames(xset)[c(5, 9)])
    ## Compare with raw data.
    rtr <- matrix(c(2876, 2932), nrow = 1)
    mzr <- matrix(c(200.1, 200.1), nrow = 1)
    e <- getEIC(xset, sampleidx = c(1), mzrange = mzr,
                rtrange = rtr, rt = "raw")
    ## Read the raw data of file 1:
    xr <- xcmsRaw(filepaths(xset)[1], profstep = profStep(xset))
    e_2 <- getEIC(xr, mzrange = mzr, rtrange = rtr, step = 0.1)
    expect_equal(e_2@eic[[1]][[1]], e@eic[[1]][[1]])
    ## Check what happens if we select another -> issue #92
    e <- getEIC(xset, sampleidx = c(5, 9), mzrange = mzr,
                rtrange = rtr, rt = "raw")
    ## sample 5
    xr <- xcmsRaw(filepaths(xset)[5], profstep = profStep(xset))
    e_2 <- getEIC(xr, mzrange = mzr, rtrange = rtr, step = 0.1)
    expect_equal(e_2@eic[[1]][[1]], e@eic[[1]][[1]])
    ## sample 9
    xr <- xcmsRaw(filepaths(xset)[9], profstep = profStep(xset))
    e_2 <- getEIC(xr, mzrange = mzr, rtrange = rtr, step = 0.1)
    expect_equal(e_2@eic[[1]][[1]], e@eic[[2]][[1]])
})

test_that("getEIC,xcmsSet works after retcor", {
    ## xset <- fillPeaks(group(retcor(group(faahko))))
    xset <- faahko_grouped_retcor_filled
    ## xset <- faahko_processed
    opt.warn <- options("warn")$warn
    options("warn" = 2) ## turns warning into errors
    e <- getEIC(xset, sampleidx=c(1,2), groupidx=c(1,2),
                rt="corrected", rtrange=200)
    options("warn" = opt.warn)
    plot(e)
})

test_that("getXcmsRaw,xcmsSet works", {

    xset <- faahko_grouped_retcor_filled
    
    ## get the first as raw data file.
    xr <- getXcmsRaw(xset, sampleidx = 1)
    ## apply the rt correction
    xr <- getXcmsRaw(xset, sampleidx = 1, rt="raw")
    ## get the first 4
    xr <- getXcmsRaw(xset, sampleidx = 1:4)
    ## check if the settings are translated correctly
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(file, step = 2)
    xr <- getXcmsRaw(xs)
    xr_orig <- xcmsRaw(file, profstep = 2)
    expect_equal(xr, xr_orig)
    ## check the prof step:
    expect_equal(profStep(xs), 2)
    ## check all *new* methods for xcmsSet
    expect_equal(mslevel(xs), mslevel(xr))
    expect_equal(profMethod(xs), profMethod(xr))
    expect_equal(profStep(xs), profStep(xr))
    ## scanrange for the xcmsSet is NULL which means we're reading all data from
    ## the raw data files, while the one of the xcmsRaw is always
    ## (1, lenght(scans)).
    ## expect_equal(scanrange(xs), scanrange(xr))
    expect_equal(profinfo(xs), profinfo(xr))
    ## testing alternative scan range.
    xr2 <- getXcmsRaw(xs, scanrange = c(5, 100))
    expect_equal(scanrange(xr2), c(1, 96))
    xr2_orig <- xcmsRaw(file, scanrange = c(5, 100), profstep = 2)
    expect_equal(xr2, xr2_orig)
    ## This scanrange is expected to be from 1 to length(xr@scantime)
    expect_equal(scanrange(xr2), c(1, length(xr2@scantime)))
    ## Test xcmsSet with scanrange:
    xs <- xcmsSet(file, step = 2, scanrange = c(5, 100))
    expect_equal(scanrange(xs), c(5, 100))
    ## BUT: if we extract the xcmsRaw from this xcmsSet object we will get
    ## (1, 96) instead of (5, 100).
    xr <- getXcmsRaw(xs)
    expect_equal(scanrange(xr), c(1, length(xr@scantime)))
    xr_2 <- xcmsRaw(file, scanrange = c(5, 100), profstep = 2)
    expect_equal(xr, xr_2)
})

test_that("getXcmsRaw,xcmsSet works, issue #44", {
    ## Subset to two files.
    xsetRaw <- updateObject(faahko)
    xsetRaw <- xsetRaw[, 1:2]

    ## First sample is reference, i.e. no rt adjustment performed
    xs <- retcor(group(xsetRaw), method = "obiwarp", center = 1)
    ## Second is corrected, first is center:
    expect_identical(xs@rt$raw[[1]], xs@rt$corrected[[1]])
    expect_true(!all(xs@rt$raw[[2]] == xs@rt$corrected[[2]]))

    ## Now, if we get the second raw file we expect to get the raw times.
    expect_true(all(xs@rt$corrected[[1]] == xs@rt$raw[[1]])) # TRUE, wouldn't do correction.
    expect_false(all(xs@rt$corrected[[2]] == xs@rt$raw[[2]])) ## FALSE, would do correction.

    ## We get the raw data for the second file; this one was corrected and
    ## thus it's retention time is expected to be different from raw.
    xr2 <- getXcmsRaw(xs, sampleidx = 2, rt = "corrected")
    expect_identical(xr2@scantime, xs@rt$corrected[[2]])

    expect_false(all(xr2@scantime == xs@rt$raw[[2]]))  ## That should be FALSE!
})

test_that("group,GroupDensity doesn't fail with OnePeak", {
    xs <- faahko
    p <- peaks(xs)
    peaks(xs) <- p[1,,drop=FALSE]
    g <- group(xs, minsamp=1, minfrac=0.001, method="density")
})

test_that("group,xcmsSet Nearest works", {
    xs <- faahko
    p <- peaks(xs)
    g <- group(xs, method="nearest")
    expect_equal(range(unlist(g@groupidx))[1],  1)
})

test_that("retcor,xcmsSet obiwarp doesn't fail", {
    faahko_sub <- faahko[,1:2]
    xr <- retcor(faahko_sub, method="obiwarp", profStep = 10)
    xr <- retcor(faahko_sub, method="obiwarp", localAlignment=1, profStep = 10)
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="cor", profStep = 10)
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="cor_opt", profStep = 10)
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="cov", profStep = 10)
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="euc", profStep = 10)
    xr <- retcor(faahko_sub, method="obiwarp", distFunc="prd", profStep = 10)
})

test_that("sampclass,xcmsSet works", {
    library(faahKO)
    xset <- faahko
    ## grouping the peaks
    xset <- group(xset, method="density")
    ## reversing the order of the classes.
    xset.revorder <- xset
    sampclass(xset.revorder) <- c(rep("WT", 6), rep("KO", 6))
    xset.revorder <- group(xset.revorder, method="density")
    ## check if we get what we want:
    expect_equal(groups(xset)[, "KO"], groups(xset.revorder)[, "WT"])

    ## repeat that but submitting already a factor
    xset.revorder.f <- xset
    sampclass(xset.revorder.f) <- factor(c(rep("WT", 6), rep("KO", 6)))
    xset.revorder.f <- group(xset.revorder.f, method="density")
    expect_equal(groups(xset)[, "KO"], groups(xset.revorder.f)[, "WT"])

    ## next: pheno data contains a column class with a factor.
    pd <- data.frame(class=factor(c(rep("WT", 6), rep("KO", 6))))
    xset.pheno <- xset
    phenoData(xset.pheno) <- pd
    xset.pheno <- group(xset.pheno, method="density")
    expect_equal(groups(xset)[, "KO"], groups(xset.pheno)[, "WT"])

    ## next checking what happens if we submit a multi-column data.frame
    ## to sampclass<-
    pd <- data.frame(dummy=rep(c("a", "b"), 6), group=c(rep("KO", 6), rep("WT", 6)),
                     pair=c(1:6, 1:6))
    phenoData(xset) <- pd
    ## No class column, so we're returning the interaction.
    expect_equal(sampclass(xset), interaction(pd, drop=TRUE))
    ## now we're going to submit 2 columns of the pd
    sampclass(xset) <- pd[, c("group", "pair")]
    expect_equal(as.character(sampclass(xset)),
                 as.character(interaction(pd[, c("group", "pair"), drop=TRUE])))
})

test_that("sampnames,xcmsSet and sampclass,xcmsSet work", {
    file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
    xs <- xcmsSet(files=file, snames=c("A"), method = "centWave", noise = 10000)
    expect_equal(sampnames(xs), "A")
    expect_equal(rownames(phenoData(xs)), "A")
    xs <- xcmsSet(files=file, sclass=c("A"), method = "centWave", noise = 10000)
    expect_equal(as.character(sampclass(xs)), "A")
})

test_that("[,xcmsSet works", {
    .compare <- function(xs, xsub, idx){
        ## sampclass?
        expect_equal(as.character(sampclass(xsub)),
                     as.character(sampclass(xs)[idx]))
        ## sampnames?
        expect_equal(sampnames(xsub),
                     sampnames(xs)[idx])
        ## filepaths?
        expect_equal(filepaths(xsub), filepaths(xs)[idx])
        ## peaks
        peaks.all <- peaks(xs)
        peaks.sub <- peaks(xsub)
        cat("Comparing peaks...")
        for(i in 1:length(idx)){
            ## peaks
            expect_equal(peaks.sub[peaks.sub[, "sample"]==i, "into"],
                         peaks.all[peaks.all[, "sample"]==idx[i], "into"])
            ## rt, raw
            expect_equal(xs@rt$raw[[idx[i]]], xsub@rt$raw[[i]])
            ## rt, corrected
            expect_equal(xs@rt$corrected[[idx[i]]], xsub@rt$corrected[[i]])
        }
        cat("OK\n")
        if(length(xs@groupidx) > 0){
            cat("Comparing peak groups...")
            ## check the groups and groupidx using the groupval method.
            into.all <- groupval(xs, value="into", method="maxint")
            into.sub <- groupval(xsub, value="into", method="maxint")
            for(i in 1:length(idx)){
                expect_equal(into.all[, idx[i]], into.sub[, i])
            }
            cat("OK\n")
        }
    }

    xset <- faahko
    idx <- 8
    xsub <- xset[, idx]
    .compare(xset, xsub, idx)
    ## repeating the same using different ordering.
    idx <- c(8, 1, 5)
    xsub <- xset[, idx]
    .compare(xset, xsub, idx)

    ## performing a grouping, retention time correction and second grouping
    ##data(faahko, package="faahKO")
    xset <- group(faahko, method="density")
    xset <- retcor(xset, method="loess", family="symmetric")
    xset <- group(xset)
    ##xset <- fillPeaks(xset)
    xsub <- xset[, idx]
    .compare(xset, xsub, idx)

    ## subset by names
    idx <- c("ko15", "wt19", "ko21", "wt22")
    xsub <- xset[, idx]
    intidx <- match(idx, sampnames(xset))
    .compare(xset, xsub, intidx)

    ## and at last using logical.
    idx <- rep(FALSE, 12)
    idx[c(2, 8)] <- TRUE
    xsub <- xset[, idx]
    intidx <- which(idx)
    .compare(xset, xsub, intidx)
    ## testing some errors...
    expect_error(xset[1, ])
    expect_error(xset[, 20])
    expect_error(xset[, "not there"])

    fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
            system.file('cdf/KO/ko16.CDF', package = "faahKO"),
            system.file('cdf/KO/ko18.CDF', package = "faahKO"),
            system.file('cdf/KO/ko19.CDF', package = "faahKO"))
    xs_1 <- xcmsSet(fs, profparam = list(step = 0), method = "centWave",
                    noise = 10000, snthresh = 50)
    ## Testing subsetting with .processHistory:
    xsub <- xs_1[, c(2, 3)]
    ph <- xs_1@.processHistory[c(2, 3)]
    ph[[1]]@fileIndex <- 1L
    ph[[2]]@fileIndex <- 2L
    expect_equal(xsub@.processHistory, ph)
    ## Reverse ordering:
    xsub <- xs_1[, c(3, 2)]
    ph <- xs_1@.processHistory[[3]]
    ph@fileIndex <- 1L
    expect_equal(.getProcessHistory(xsub, fileIndex = 1), list(ph))
    ph <- xs_1@.processHistory[[2]]
    ph@fileIndex <- 2L
    expect_equal(.getProcessHistory(xsub, fileIndex = 2), list(ph))

    ## Add fake ProcessHistory before and after the real ones.
    ph <- xs_1@.processHistory
    ph <- c(list(ProcessHistory(fileIndex = 1:4)), ph,
            list(ProcessHistory(fileIndex = 1:4)))
    xs_1@.processHistory <- ph
    xsub <- xs_1[, c(2, 3)]
    ph <- xs_1@.processHistory[c(1, 3, 4, 6)]
    ph[[1]]@fileIndex <- 1L:2L
    ph[[2]]@fileIndex <- 1L
    ph[[3]]@fileIndex <- 2L
    ph[[4]]@fileIndex <- 1L:2L
    expect_equal(xsub@.processHistory, ph)
})

test_that("sampclass,xcmsSet works with unused groups", {
    classes <- sampclass(faahko)
    levels(classes) <- c(levels(classes), "Leftover")
    xs <- faahko
    sampclass(xs) <- classes
    xsg <- group(xs)
    expect_equal(sampclass(xs), sampclass(xsg))
})

test_that("updateObject,xcmsSet works", {
    data(xs)

    newXs <- updateObject(xs)
    expect_identical(xs@peaks, newXs@peaks)
    ## xs might be an old object.
    if (!.hasSlot(xs, "scanrange"))
        expect_identical(newXs@scanrange, numeric(0))
    expect_identical(xs@groups, newXs@groups)
})
