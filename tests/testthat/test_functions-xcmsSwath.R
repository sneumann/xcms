test_that("chromatogram correlation works", {
  
  # create
  rtime1 <- seq(5,25, step = 1)
  intensity1 <- dnorm(rtime1, mean=14, sd=1.0)*200
  
  rtime2 <- seq(5.5, 25.5, step = 1)
  intensity2 <- dnorm(rtime2, mean=14, sd=1.0)*500
  
  ## bogus chromatograms
  ms1chrom <- new("Chromatogram",
                  rtime = rtime1,
                  intensity = intensity1)
  
  ms2chrom <- new("Chromatogram",
                  rtime = rtime2,
                  intensity = intensity2)
  
  # check that retention times are the same
  expect_equal(rtime(.align_chromatogram(ms1chrom, ms2chrom)), rtime(ms1chrom))
  
  # check that correlation with NA values fails
  expect_equal(.correlate_chromatogram(ms1chrom, ms2chrom, use = "everything"), NA)
  
  # check correlation value is correct (Pearson, standard setting)
  expect_equal(.correlate_chromatogram(ms1chrom, ms2chrom), 0.9956445)
  
  # check correlation value is correct (Kendall)
  expect_equal(.correlate_chromatogram(ms1chrom, ms2chrom, method = "kendall"), 0.997264)
  
  # check correlation value is correct (Spearman)
  expect_equal(.correlate_chromatogram(ms1chrom, ms2chrom, method = "spearman"), 0.999622)
  
})