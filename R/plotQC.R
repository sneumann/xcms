#' Plot m/z and RT deviations for QC purposes without external reference data
#'
#' Use "democracy" to determine the average m/z and RT deviations
#' for a grouped xcmsSet, and dependency on sample or absolute m/z
#'
#' plotQC() is a warpper to create a set of diagnostic plots.
#' For the m/z deviations, the median of all m/z withon one group are assumed.
#'
#' @usage plotQC(object, sampNames, sampColors, sampOrder, what)
#'
#' @aliases plotQC
#'
#' @param object     A grouped \code{\link{xcmsSet}}
#' @param sampNames  Override sample names (e.g. with simplified names)
#' @param sampColors Provide a set of colors (default: monochrome ?)
#' @param sampOrder  Override the order of samples, e.g. to bring them in order
#'                   of measurement to detect time drift
#' @param what       A vector of which QC plots to generate.
#'  "mzdevhist": histogram of mz deviations. Should be gaussian shaped. If it is multimodal, then some peaks seem to have a systematically higher m/z deviation
#'  "rtdevhist":  histogram of RT deviations. Should be gaussian shaped. If it is multimodal, then some peaks seem to have a systematically higher RT deviation
#'  "mzdevmass": Shows whether m/z deviations are absolute m/z dependent, could indicate miscalibration
#'  "mzdevtime": Shows whether m/z deviations are RT dependent, could indicate instrument drift
#'  "mzdevsample": median mz deviation for each sample, indicates outliers
#'  "rtdevsample": median RT deviation for each sample, indicates outliers
#'
#'
#' @return List with four matrices, each of dimension features * samples:
#'  "mz":    median mz deviation for each sample
#'  "mzdev": median mz deviation for each sample
#'  "rt":    median RT deviation for each sample
#'  "rtdev": median RT deviation for each sample
#'
#' @examples
#' library(faahKO)
#' xsg <- group(faahko)
#'
#' plotQC(xsg, what="mzdevhist")
#' plotQC(xsg, what="rtdevhist")
#' plotQC(xsg, what="mzdevmass")
#' plotQC(xsg, what="mzdevtime")
#' plotQC(xsg, what="mzdevsample")
#' plotQC(xsg, what="rtdevsample")
#'
#' @author Michael Wenk, Michael Wenk <michael.wenk@@student.uni-halle.de>

plotQC <- function(object,
                   sampNames = NULL,
                   sampColors = NULL,
                   sampOrder=NULL,
                   what=c("mzdevhist",
                          "rtdevhist",
                          "mzdevmass",
                          "mzdevtime",
                          "mzdevsample",
                          "rtdevsample")) {

  if (inherits(object, "xcmsSet")) {
    deviations <- list(mzs=groupval(object, value = "mz"),
                       rts=groupval(object, value = "rt"))
    deviations$mzdev <- deviations$mzs - groups(object)[,"mzmed"]
    deviations$rtdev <- deviations$rts - groups(object)[,"rtmed"]

  } else if (inherits(object, "XCMSnExp")) {
    deviations <- list(mzs=featureValues(object, value = "mz"),
                       rts=featureValues(object, value = "rt"))
    deviations$mzdev <- deviations$mzs - featureDefinitions(object)[,"mzmed"]
    deviations$rtdev <- deviations$rts - featureDefinitions(object)[,"rtmed"]
  } else {
    stop("object not xcmsSet nor XCMSnExp")
  }

  if (missing(sampNames) || is.null(sampNames)) {
    if (class(object)=="xcmsSet") {
      sampNames <- sampnames(object)
    } else if (class(object)=="XCMSnExp") {
      sampNames <- sampleNames(object)
    } else {
      stop("object not xcmsSet nor XCMSnExp")
    }
  }

  n <- ncol(deviations$mzs)
  if (missing(sampColors) || is.null(sampColors)) {
    sampColors <- rainbow(n)
  }

  if (missing(sampOrder) || is.null(sampOrder)) {
    sampOrder <- 1:n
  }

  ## Plot histograms of deviation
  if ("mzdevhist" %in% what) {
    hist(deviations$mzdev, breaks=100,
         ylab = "Number of Peaks",
         xlab = "m/z Deviation",
         main = "m/z Deviation")
  }

  if ("rtdevhist" %in% what) {
    hist(deviations$rtdev, breaks=100, ylab = "Number of Peaks",
         xlab = "Retention Time Deviation",
         main = "Retention Time Deviation")
  }

  if ("mzdevmass" %in% what) {
    ## Plot mass deviation depending on absolute mass
    # Add extra space to right of plot area; change clipping to figure
    # par(mar=c(6, 8, 8.1, 10), xpd=TRUE)
    # par(oma = c(0,0,0,10))
    plot(x = as.vector(deviations$mzs),
         y = deviations$mzdev,
         col=sampColors, pch = ".", type = "p",
         main = "m/z Deviation vs. m/z",
         xlab = "m/z", ylab = "m/z deviation")

    for(i in 1:n) {
      data <- na.omit(data.frame(mzs = deviations$mzs[,i], mzdev = deviations$mzdev[,i]))
      lo <- loess(formula = mzdev ~ mzs, data = data)
      currSampleValues <- deviations$mzs[!is.na(deviations$mzs[,i]),i]
      currSampleValues <- currSampleValues[order(currSampleValues)]
      lines(currSampleValues, predict(lo), col=sampColors[i], lwd = 2)
    }
    legend("topright", legend=sampNames, col=sampColors, pch = c(1), cex=0.35,  ncol=4, title = "samples")
  }


  if ("mzdevtime" %in% what) {
    ## Plot mass deviation depending on retention time
    ## Add extra space to right of plot area; change clipping to figure
    ##par(mar=c(5.1, 5.1, 8.1, 8.1), xpd=TRUE)
    plot(x = as.vector(deviations$rts), y = deviations$mzdev, pch = ".", type = "p",
         col=sampColors, main = "m/z deviation vs. retention time",
         xlab = "Retention Time", ylab = "m/z deviation")

    for(i in 1:n){
      data <- na.omit(data.frame(rts = deviations$rts[,i], mzdev = deviations$mzdev[,i]))
      lo <- loess(formula = mzdev ~ rts, data = data)
      currSampleValues <- deviations$rts[!is.na(deviations$rts[,i]),i]
      currSampleValues <- currSampleValues[order(currSampleValues)]
      lines(currSampleValues, predict(lo), col=sampColors[i], lwd = 1)
    }
    legend("topright", legend=sampNames, col=sampColors, pch = c(1), cex=0.35,  ncol=4, title = "samples")
  }

  ## still to come: median deviations per sample,
  ## to detect corrupt samples

  if ("mzdevsample" %in% what) {
    barplot(apply(deviations$mzdev[,sampOrder], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE)),
            col = sampColors, xlab = "", ylab = "m/z Deviation",
            names.arg = sampNames, las = 2)
  }
  if ("rtdevsample" %in% what) {
    barplot(apply(deviations$rtdev[,sampOrder], MARGIN=2, FUN=function(x) median(x, na.rm=TRUE)),
            col = sampColors, xlab = "", ylab = "Retention Time Deviation",
            names.arg = sampNames, las = 2)
  }
  invisible(deviations)
}
