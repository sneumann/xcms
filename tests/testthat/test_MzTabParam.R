faahko <- loadXcmsData("faahko_sub2")
faahko <- groupChromPeaks(
    faahko, PeakDensityParam(sampleGroups = rep(1, length(faahko))))

xmse_full <- loadXcmsData()

test_that(".mztab_metadata works", {
    mtd <- .mztab_metadata(xmse_full, study_id = "test", polarity = "negative",
                         col_phenotype = "sample_type")

    expect_false(all(mtd != "test"))
    expect_true(is.matrix(mtd))
    expect_true(is.character(mtd))
    expect_equal(mtd[, "id"], rep("MTD", nrow(mtd)))
    expect_equal(ncol(mtd), 3)
    expect_equal(length(grep("assay ", mtd[, 2])), length(xmse_full))
    expect_equal(length(grep("negative ", mtd[, 3])), length(xmse_full))
})

test_that(".mztab_small_molecule_feature works", {
    a <- .mztab_small_molecule_feature(faahko)
    expect_true(is.matrix(a))
    expect_true(is.character(a))
    expect_equal(a[1, 1], "SFH")
    expect_equal(nrow(a), nrow(featureDefinitions(faahko)) + 1)
    expect_equal(length(grep("abundance_assay", a[1, ])), length(faahko))

    b <- .mztab_small_molecule_feature(
        faahko, opt_columns = c("peakidx", "ms_level"),
        method = "sum", value = "maxo")
    expect_equal(nrow(a), nrow(b))
    expect_true(ncol(b) > ncol(a))
    expect_true(any(b[1, ] == "opt_peakidx"))
    expect_true(any(b[1, ] == "opt_ms_level"))
    ## Check that optional parameters are correctly passed
    a2 <- a[2:nrow(a), a[1, ] == "abundance_assay[2]"]
    a2[a2 == "null"] <- NA
    a2 <- as.numeric(a2)
    b2 <- b[2:nrow(b), b[1, ] == "abundance_assay[2]"]
    b2[b2 == "null"] <- NA
    b2 <- as.numeric(b2)
    expect_equal(is.na(a2), is.na(b2))
    expect_true(all(a2[!is.na(a2)] != b2[!is.na(b2)]))
    ## Check that we don't loose information through import/export
    expect_equal(a2, unname(featureValues(faahko)[, 2]))
    expect_equal(b2, unname(featureValues(
        faahko, method = "sum", value = "maxo")[, 2]))
})

test_that(".mztab_study_variables works", {
    x <- data.frame(a = 1:3, b = c(1, 1, 2))
    res <- .mztab_study_variables(x, variable = "a")
    expect_true(is.character(res))
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 1L)
    expect_equal(res[, 1L], c("a:1", "a:2", "a:3"))

    res <- .mztab_study_variables(x, variable = c("a", "b"))
    expect_true(is.character(res))
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 2L)
    expect_equal(res[, 1L], c("a:1", "a:2", "a:3"))
    expect_equal(res[, 2L], c("b:1", "b:1", "b:2"))
})

test_that(".mztab_study_variable_entries works", {
    x <- data.frame(a = 1:3, b = c(1, 1, 2))
    res <- .mztab_study_variable_entries(x, variable = "a")
    expect_true(is.character(res))
    expect_true(is.matrix(res))
    expect_true(ncol(res) == 2L)
    expect_equal(unname(res[1, 1]), "study_variable[1]")
    expect_equal(unname(res[2, 1]), "study_variable[1]-assay_refs")
    expect_equal(unname(res[1, 2]), "a:1")
    expect_equal(unname(res[2, 2]), "assay[1]")

    res2 <- .mztab_study_variable_entries(x, variable = c("a", "b"))
    expect_true(length(res2) > length(res))
})

test_that("storeResults,MzTabParam works", {
    d <- tempdir()
    p <- MzTabParam(studyId = "test_study", path = d,
                    sampleDataColumn = "sample_index",
                    optionalFeatureColumns = "peakidx")
    storeResults(faahko, p)
    expect_true(file.exists(file.path(d, "test_study.mztab")))
    res <- readLines(file.path(d, "test_study.mztab"))
    expect_true(length(res) > 0L)
    expect_true(length(grep("^MTD", res)) > 0)
    expect_true(length(grep("^SML", res)) > 0)
    expect_true(length(grep("^SMF", res)) > 0)
    ## Check for empty lines
    expect_true(length(grep(c("^MTD|SML|SMF"), res, invert = TRUE)) == 2)
})
