test.write.mzQuantML <- function() {
    xsg <- group(faahko)

    mzqFile <- paste(tempdir(), "faahKO.mzq.xml", sep="/")
    write.mzQuantML(xsg, mzqFile)

    v <- verify.mzQuantML(filename=mzqFile)

    checkTrue(v$status == "0")    
}

