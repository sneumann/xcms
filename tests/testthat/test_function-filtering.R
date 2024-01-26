xmse_full <- loadXcmsData("xmse")

test_that("RsdFilter", {
    # Test default parameters
    filter <- RsdFilter()
    expect_is(filter, "RsdFilter")
    expect_equal(filter@threshold, 0.3)
    expect_equal(filter@qcIndex, integer())
    expect_equal(filter@na.rm, TRUE)
    expect_equal(filter@mad, FALSE)

    #test with XcmsExperiment object
    filter <- RsdFilter(qcIndex = sampleData(xmse_full)$sample_type == "QC")
    filtered_xmse <- filterFeatures(xmse_full, filter)
    expect_lte(nrow(featureDefinitions(filtered_xmse)),
               nrow(featureDefinitions(xmse_full)))

    #test error message when index too small of too large
    filter_1 <- RsdFilter(qcIndex = integer())
    filter_2 <- RsdFilter(qcIndex = c(sampleData(xmse_full)$sample_type == "QC",
                                      TRUE))
    filter_3 <- RsdFilter(qcIndex = c(2, 15))
    expect_error(filterFeatures(xmse_full, filter_1))
    expect_error(filterFeatures(xmse_full, filter_2))
    expect_error(filterFeatures(xmse_full, filter_3))

    #test in SummarizedExperiment object
    res <- quantify(xmse_full)
    filter <- RsdFilter(qcIndex = res$sample_type == "QC")
    filtered_res <- filterFeatures(res, filter)
    expect_lte(nrow(filtered_res), nrow(res))

    #test same amount of features in in both object type
    expect_equal(nrow(featureDefinitions(filtered_xmse)),
               nrow(filtered_res))
})

test_that("DratioFilter", {
    # Test default parameters
    filter <- DratioFilter()
    expect_is(filter, "DratioFilter")
    expect_equal(filter@threshold, 0.5)
    expect_equal(filter@qcIndex, integer())
    expect_equal(filter@studyIndex, integer())
    expect_equal(filter@na.rm, TRUE)
    expect_equal(filter@mad, FALSE)

    #test with XcmsExperiment object
    filter <- DratioFilter(
        qcIndex = sampleData(xmse_full)$sample_type == "QC",
        studyIndex = sampleData(xmse_full)$sample_type == "study")
    filtered_xmse <- filterFeatures(xmse_full, filter)
    expect_lte(nrow(featureDefinitions(filtered_xmse)),
               nrow(featureDefinitions(xmse_full)))

    #test in SummarizedExperiment object
    res <- quantify(xmse_full)
    filter <- DratioFilter(qcIndex = res$sample_type == "QC",
                           studyIndex = res$sample_type == "study")
    filtered_res <- filterFeatures(res, filter)
    expect_lte(nrow(filtered_res), nrow(res))

    #test same amount of features in in both object type
    expect_equal(nrow(featureDefinitions(filtered_xmse)),
                 nrow(filtered_res))
})

test_that("PercentMissingFilter", {
    # Test default parameters
    filter <- PercentMissingFilter()
    expect_is(filter, "PercentMissingFilter")
    expect_equal(filter@threshold, 30)
    expect_equal(filter@f, character())

    #test with XcmsExperiment object
    filter <- PercentMissingFilter(f = sampleData(xmse_full)$sample_type)
    filtered_xmse <- filterFeatures(xmse_full, filter)
    expect_lte(nrow(featureDefinitions(filtered_xmse)),
               nrow(featureDefinitions(xmse_full)))

    #test in SummarizedExperiment object
    res <- quantify(xmse_full)
    filter <- PercentMissingFilter(f = res$sample_type)
    filtered_res <- filterFeatures(res, filter)
    expect_lte(nrow(filtered_res), nrow(res))

    #test same amount of features in in both object type
    expect_equal(nrow(featureDefinitions(filtered_xmse)),
                 nrow(filtered_res))
})

test_that("BlankFlag", {
    # Test default parameters
    filter <- BlankFlag()
    expect_is(filter, "BlankFlag")
    expect_equal(filter@threshold, 2)
    expect_equal(filter@na.rm, TRUE)
    expect_equal(filter@qcIndex, integer())
    expect_equal(filter@blankIndex, integer())

    #test with XcmsExperiment object - using study sample as blanks
    filter <- BlankFlag(
        qcIndex = sampleData(xmse_full)$sample_type == "QC",
        blankIndex = sampleData(xmse_full)$sample_type == "study")
    filtered_xmse <- filterFeatures(xmse_full, filter)
    expect_true("possible_contaminants" %in%
                    colnames(featureDefinitions(filtered_xmse)))

    #test in SummarizedExperiment object - using study sample as blanks
    res <- quantify(xmse_full)
    filter <- BlankFlag(qcIndex = res$sample_type == "QC",
                        blankIndex = res$sample_type == "study")
    filtered_res <- filterFeatures(res, filter)
    expect_true("possible_contaminants" %in% colnames(rowData(filtered_res)))

    expect_equal(featureDefinitions(filtered_xmse)$possible_contaminants,
                 rowData(filtered_res)$possible_contaminants)
})
