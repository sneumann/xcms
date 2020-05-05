test_that("specDist works", {
    data(faahko, package="faahKO")
    dist1 <- specDist(faahko, which(faahko@peaks[,"sample"] == 1),
                      which(faahko@peaks[,"sample"] == 2),
                      method = "meanMZmatch")
    dist2 <- specDist(faahko, which(faahko@peaks[,"sample"] == 3),
                      which(faahko@peaks[,"sample"] == 4),
                      method = "cosine")
    dist3 <- specDist(faahko, which(faahko@peaks[,"sample"] == 5),
                      which(faahko@peaks[,"sample"] == 6),
                      method = "peakCount")
    expect_equal(round(dist1, 3), 3.159)
    expect_equal(round(dist2, 3), 3.183)
    expect_equal(dist3, 129)
})

