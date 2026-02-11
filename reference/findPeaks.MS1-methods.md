# Collecting MS1 precursor peaks

Collecting Tandem MS or MS\$^n\$ Mass Spectrometry precursor peaks as
annotated in XML raw file

## Methods

- object = "xcmsRaw":

  ` findPeaks.MS1(object) `

## Details

Some mass spectrometers can acquire MS1 and MS2 (or MS\$^n\$ scans)
quasi simultaneously, e.g. in data dependent tandem MS or DDIT mode.

Since xcmsFragments attaches *all* MS\$^n\$ peaks to MS1 peaks in
xcmsSet, it is important that findPeaks and xcmsSet do not miss any MS1
precursor peak.

To be sure that all MS1 precursor peaks are in an xcmsSet, findPeaks.MS1
does not do an actual peak picking, but simply uses the annotation
stored in mzXML, mzData or mzML raw files.

This relies on the following XML tags:

mzData:
` <spectrum id="463"> <spectrumInstrument msLevel="2"> <cvParam cvLabel="psi" accession="PSI:1000039" name="TimeInSeconds" value="92.7743"/> </spectrumInstrument> <precursor msLevel="1" spectrumRef="461"> <cvParam cvLabel="psi" accession="PSI:1000040" name="MassToChargeRatio" value="462.091"/> <cvParam cvLabel="psi" accession="PSI:1000042" name="Intensity" value="366.674"/> </precursor> </spectrum> `

mzXML:
` <scan num="17" msLevel="2" retentionTime="PT1.5224S"> <precursorMz precursorIntensity="125245">220.1828003</precursorMz> </scan> `

Several mzXML and mzData converters are known to create incomplete
files, either without intensities (they will be set to 0) or without the
precursor retention time (then a reasonably close rt will be chosen.
NYI).

## Arguments

- object:

  `xcmsRaw` object

## Value

A matrix with columns:

- mz, mzmin, mzmax:

  annotated MS1 precursor selection mass

- rt, rtmin, rtmax:

  annotated MS1 precursor retention time

- into, maxo, sn:

  annotated MS1 precursor intensity

## Author

Steffen Neumann, <sneumann@ipb-halle.de>

## See also

[`findPeaks-methods`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md)
[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
