\name{diffreport-methods}
\docType{methods}
\alias{diffreport}
\alias{diffreport,xcmsSet-method}
\title{Create report of analyte differences}
\description{
  Create a report showing the most significant differences between
  two sets of samples. Optionally create extracted ion chromatograms
  for the most significant differences.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{diffreport(object, class1 = levels(sampclass(object))[1],
                   class2 = levels(sampclass(object))[2],
                   filebase = character(), eicmax = 0, eicwidth = 200,
                   sortpval = TRUE, classeic = c(class1,class2),
                   value=c("into","maxo","intb"), metlin = FALSE,
				   h=480,w=640, mzdec=2, ...)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{class1}{
    character vector with the first set of sample classes to be
    compared
  }
  \item{class2}{
    character vector with the second set of sample classes to be
    compared
  }
  \item{filebase}{
    base file name to save report, \code{.tsv} file and \code{_eic}
    will be appended to this name for the tabular report and EIC
    directory, respectively. if blank nothing will be saved
  }
  \item{eicmax}{
    number of the most significantly different analytes to create
    EICs for
  }
  \item{eicwidth}{
    width (in seconds) of EICs produced
  }
  \item{sortpval}{
    logical indicating whether the reports should be sorted by
    p-value
  }
  \item{classeic}{
    character vector with the sample classes to include in the EICs
  }
  \item{value}{
    intensity values to be used for the diffreport. \cr
    If \code{value="into"}, integrated peak intensities are used. \cr
    If \code{value="maxo"}, maximum peak intensities are used. \cr
    If \code{value="intb"}, baseline corrected integrated peak intensities are used (only available if peak detection was done by \code{\link{findPeaks.centWave}}).
  }
  \item{metlin}{
    mass uncertainty to use for generating link to Metlin metabolite
    database. the sign of the uncertainty indicates negative or
    positive mode data for M+H or M-H calculation. a value of FALSE
    or 0 removes the column
  }
  \item{h}{
    Numeric variable for the height of the eic and boxplots that are printed out.
  }
  \item{w}{
    Numeric variable for the width of the eic and boxplots print out made.
  }
  \item{mzdec}{
    Number of decimal places of title m/z values in the eic plot.
  }
  \item{...}{
    optional arguments to be passed to \code{\link{mt.teststat}}
  }
}
\details{
  This method handles creation of summary reports with statistics
  about which analytes were most significantly different between
  two sets of samples. It computes Welch's two-sample t-statistic
  for each analyte and ranks them by p-value. It returns a summary
  report that can optionally be written out to a tab-separated file.

  Additionally, it does all the heavy lifting involved in creating
  superimposed extracted ion chromatograms for a given number of
  analytes.  It does so by reading the raw data files associated
  with the samples of interest one at a time. As it does so, it
  prints the name of the sample it is currently reading. Depending
  on the number and size of the samples, this process can take a
  long time.

  If a base file name is provided, the report (see Value section)
  will be saved to a tab separated file. If EICs are generated,
  they will be saved as 640x480 PNG files in a newly created
  subdirectory. However this parameter can be changed with the
  commands arguments. The numbered file names correspond to the
  rows in the report.

  Chromatographic traces in the EICs are colored and labeled by
  their sample class. Sample classes take their color from the
  current palette. The color a sample class is assigned is dependent
  its order in the \code{xcmsSet} object, not the order given in
  the class arguments. Thus \code{levels(sampclass(object))[1]}
  would use color \code{palette()[1]} and so on. In that way, sample
  classes maintain the same color across any number of different
  generated reports.

  When there are multiple sample classes, xcms will produce boxplots of the
  different classes and will generate a single anova p-value statistic.
  Like the eic's the plot number corresponds to the row number in the
  report.
}

\value{
  A data frame with the following columns:

  \item{fold}{
    mean fold change (always greater than 1, see \code{tstat} for
    which set of sample classes was higher)
  }
  \item{tstat}{
    Welch's two sample t-statistic, positive for analytes having
    greater intensity in \code{class2}, negative for analytes having
    greater intensity in \code{class1}
  }
  \item{pvalue}{p-value of t-statistic}
  \item{anova}{p-value of the anova statistic if there are multiple classes}
  \item{mzmed}{median m/z of peaks in the group}
  \item{mzmin}{minimum m/z of peaks in the group}
  \item{mzmax}{maximum m/z of peaks in the group}
  \item{rtmed}{median retention time of peaks in the group}
  \item{rtmin}{minimum retention time of peaks in the group}
  \item{rtmax}{maximum retention time of peaks in the group}
  \item{npeaks}{number of peaks assigned to the group}
  \item{Sample Classes}{
    number samples from each sample class represented in the group
  }
  \item{metlin}{A URL to metlin for that mass}
  \item{...}{one column for every sample class}
  \item{Sample Names}{integrated intensity value for every sample}
  \item{...}{one column for every sample}
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{mt.teststat}},
  \code{\link{palette}}
}
\keyword{methods}
\keyword{file}
