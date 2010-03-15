test.xcmsRawms123 <- function() {

  filename <- system.file('iontrap/extracted.mzData', package = "msdata")

  x1 <- xcmsRaw(filename, includeMSn=TRUE)
  x2 <- xcmsRaw(filename, includeMSn=TRUE, mslevel=2)
  x3 <- xcmsRaw(filename, includeMSn=TRUE, mslevel=3)

  checkTrue(length(x1@env$msnMz) == length(x2@env$mz) + length(x3@env$mz))

  checkTrue(all(x1@msnLevel[1:6]==2))
  
  checkTrue(all(x1@msnScanindex[1:6] == x2@scanindex[1:6]))
  
}

test.xcmsSetms2mf <- function() {

  filename <- system.file('iontrap/DhexD4_2.mzData', package = "msdata")
  xs2 <- xcmsSet(filename, mslevel=2)
}

test.xcmsSetms2cw <- function() {

  filename <- system.file('iontrap/DhexD4_2.mzData', package = "msdata")
  xs2 <- xcmsSet(filename, method="centWave", mslevel=2)

}
