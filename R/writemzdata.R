
setMethod("write.mzdata", "xcmsRaw", function(object, filename) {

    require(caTools) || stop("We need library(caTools) to encode base64 values in mzData")

    sink(file=filename)

    cat('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    cat('<mzData xsi:noNamespaceSchemaLocation="http://psidev.sourceforge.net/ms/xml/mzdata/mzdata.xsd" version="1.05" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" accessionNumber="1">\n',
        '<cvLookup cvLabel="PSI" fullName="The PSI Ontology" version="1.4.0" address="http://psidev.sourceforge.net/ontology/"/>\n')
    cat ('<spectrumList count="',length(object@scanindex),'">\n', sep="")

    mslevel <- 1; # Write MS1 first
    offset <- 0;
    polarities <- c("Negative", "Positive", "Unknown")

    if (length(object@scanindex) > 0) {
        for (id in 1:length(object@scanindex)) {

            cat ('<spectrum id="', id, '">\n', sep="")
            cat ('<spectrumDesc>\n')
            cat ('<spectrumSettings>\n')
            cat ('<spectrumInstrument msLevel="1">\n')

            cat('<cvParam cvLabel="PSI" accession="PSI:1000039" ',
                'name="TimeInSeconds" value="', object@scantime[id], '"/>\n', sep="")

            if (!is.null(object@polarity) & !is.na(object@polarity[id])) {
                cat('<cvParam cvLabel="PSI" accession="PSI:1000037" ',
                    'name="Polarity" value="', polarities[object@polarity[id]], '"/>\n', sep="")
            }

            cat ("</spectrumInstrument>\n")
            cat ("</spectrumSettings>\n")
            cat ('</spectrumDesc>\n')

            target <- new("raw")
            peaks <- getScan(object, id)

            if (is.unsorted(peaks[,"mz"])) { ## fix "bad" scans
                peaks <- peaks[order(peaks[,"mz"]),]
            }

            cat ('<mzArrayBinary>\n')
            cat ('<data precision="32" endian="big" length="', nrow(peaks), '">', sep="")
            cat (base64encode(writeBin(peaks[,"mz"], con=target, size=4, endian="big")))
            cat ('</data>\n')
            cat ('</mzArrayBinary>\n')

            cat ('<intenArrayBinary>\n')
            cat ('<data precision="32" endian="big" length="', nrow(peaks), '">', sep="")
            cat (base64encode(writeBin(peaks[,"intensity"], con=target, size=4, endian="big")))
            cat ('</data>\n')
            cat ('</intenArrayBinary>\n')

            cat ('</spectrum>\n')

            offset <- id
        }
    }

    mslevel <- object@msnLevel
    if (length(object@msnScanindex) >0 ) {
        for (id in seq(1, length(object@msnScanindex))) {

            cat ('<spectrum id="', id=id+offset, '">\n', sep="")
            cat ('<spectrumDesc>\n')
            cat ('<spectrumSettings>\n')
            cat ('<spectrumInstrument msLevel="', msLevel=mslevel[id], '">\n', sep="")

            cat('<cvParam cvLabel="PSI" accession="PSI:1000039" ',
                'name="TimeInSeconds" value="', object@msnRt[id], '"/>\n', sep="")

            if (!is.null(object@polarity) & !is.na(object@polarity[id])) {
                cat('<cvParam cvLabel="PSI" accession="PSI:1000037" ',
                    'name="Polarity" value="', polarities[object@polarity[id]], '"/>\n', sep="")
            }

            cat ('</spectrumInstrument>\n')
            cat ('</spectrumSettings>\n')

            if (mslevel[id] > 1) { ## only for mslevel >1, here always true
                cat ('<precursorList count="1">\n')
                cat ('<precursor msLevel="', mslevel[id]-1, '" ',
                     'spectrumRef="', max(0, object@msnPrecursorScan[id], na.rm=T), '">\n', sep="")
                cat('<ionSelection>\n')
                cat('<cvParam cvLabel="PSI" accession="PSI:1000040"',
                    ' name="m/z" value="', object@msnPrecursorMz[id], '"/>\n', sep="")
                cat('</ionSelection>\n')

                cat('<activation>\n')
                cat('<cvParam cvLabel="PSI" accession="PSI:1000045" ',
                    'name="Collision Energy" value="', object@msnCollisionEnergy[id], '"/>\n', sep="")
                cat('</activation>\n')

                cat ('</precursor>\n')
                cat ('</precursorList>\n')

            }
            cat ('</spectrumDesc>\n')

            target <- new("raw")
            peaks <- getMsnScan(object, id)

            if (is.unsorted(peaks[,"mz"])) { ## fix "bad" scans
                peaks <- peaks[order(peaks[,"mz"]),]
            }

            cat ('<mzArrayBinary>\n')
            cat ('<data precision="32" endian="big" length="', nrow(peaks), '">', sep="")
            cat (base64encode(writeBin(peaks[,"mz"], con=target, size=4, endian="big")))
            cat ('</data>\n')
            cat ('</mzArrayBinary>\n')


            cat ('<intenArrayBinary>\n')
            cat ('<data precision="32" endian="big" length="', nrow(peaks), '">', sep="")
            cat (base64encode(writeBin(peaks[,"intensity"], con=target, size=4, endian="big")))
            cat ('</data>\n')
            cat ('</intenArrayBinary>\n')

            cat ('</spectrum>\n')

        }
    }

    cat ('</spectrumList>\n')
    cat ('</mzData>\n')

    sink(NULL)
})

