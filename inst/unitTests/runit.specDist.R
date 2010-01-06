test.specDist <- function() {
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 dist1 <- specDist(xs,c(1:50),c(50:100), method="meanMZmatch")
 dist2 <- specDist(an, PSpec1=1,PSpec2=1, method="cosinus")
 checkEquals(dist1, 50.5)
 checkEquals(round(dist2,3), 0.966)    
}
