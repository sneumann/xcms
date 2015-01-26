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
#' @return No return value
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
    
  if (missing(sampNames) || is.null(sampNames)) {
      sampNames <- sampnames(object)
  }
  
  n <- length(sampclass(object))  
  if (missing(sampColors) || is.null(sampColors)) {
    sampColors <- rainbow(n)
  }

  if (missing(sampOrder) || is.null(sampOrder)) {
    sampOrder <- 1:n
  }

  deviations <- getdeviations(object, sampNames)
  deviation_mzs <- deviations$mzdev
  deviation_rts <- deviations$rtdev
  mzs <- deviations$mzs
  rts <- deviations$rts  
  
  ## Plot histograms of deviation
  if ("mzdevhist" %in% what) {
      hist(deviation_mzs, breaks=100,
           ylab = "Number of Peaks",
           xlab = "m/z Deviation",
           main = "m/z Deviation")
  }

  if ("rtdevhist" %in% what) {
      hist(deviation_rts, breaks=100, ylab = "Number of Peaks",
           xlab = "Retention Time Deviation",
           main = "Retention Time Deviation")
  }

  if ("mzdevmass" %in% what) {
  ## Plot mass deviation depending on absolute mass
  # Add extra space to right of plot area; change clipping to figure
  # par(mar=c(6, 8, 8.1, 10), xpd=TRUE)
  # par(oma = c(0,0,0,10))
  plot(x = as.vector(mzs), 
       y = deviation_mzs, 
       col=sampColors, pch = ".", type = "p", 
       main = "m/z Deviation vs. m/z", 
       xlab = "m/z", ylab = "m/z deviation")
  
  for(i in 1:n) {
    data <- na.omit(data.frame(mzs = mzs[,i], deviation_mzs = deviation_mzs[,i]))
    lo <- loess(formula = deviation_mzs ~ mzs, data = data)
    currSampleValues <- mzs[!is.na(mzs[,i]),i]
    currSampleValues <- currSampleValues[order(currSampleValues)]
    lines(currSampleValues, predict(lo), col=sampColors[i], lwd = 2)
  }
  legend("topright", legend=sampNames, col=sampColors, pch = c(1), cex=0.35,  ncol=4, title = "samples")
}  
  

    if ("mzdevtime" %in% what) {
        ## Plot mass deviation depending on retention time
        ## Add extra space to right of plot area; change clipping to figure
        ##par(mar=c(5.1, 5.1, 8.1, 8.1), xpd=TRUE)
  plot(x = as.vector(rts), y = deviation_mzs, pch = ".", type = "p", 
       col=sampColors, main = "m/z deviation vs. retention time", 
       xlab = "Retention Time", ylab = "m/z deviation")
  
  for(i in 1:n){
    data <- na.omit(data.frame(rts = rts[,i], deviation_mzs = deviation_mzs[,i]))
    lo <- loess(formula = deviation_mzs ~ rts, data = data)
    currSampleValues <- rts[!is.na(rts[,i]),i]
    currSampleValues <- currSampleValues[order(currSampleValues)]
    lines(currSampleValues, predict(lo), col=sampColors[i], lwd = 1)
  }
  legend("topright", legend=sampNames, col=sampColors, pch = c(1), cex=0.35,  ncol=4, title = "samples")
}

  ## still to come: median deviations per sample, 
  ## to detect corrupt samples

  if ("mzdevsample" %in% what) {
      barplot(apply(deviation_mzs[,sampOrder], MAR=2, FUN=function(x) median(x, na.rm=TRUE)),
              col = sampColors, xlab = "", ylab = "m/z Deviation", 
              names.arg = sampNames, las = 2)
  }
  if ("rtdevsample" %in% what) {      
      barplot(apply(deviation_rts[,sampOrder], MAR=2, FUN=function(x) median(x, na.rm=TRUE)),
              col = sampColors, xlab = "", ylab = "Retention Time Deviation", 
              names.arg = sampNames, las = 2)
  }
}


getdeviations <- function(object, sampNames = NULL) {
  if (missing(sampNames) || is.null(sampNames)) {
      sampNames <- sampnames(object)
  }
  n <- length(sampnames(object))  
  p <- peaks(object)
  
  ## Get a matrix of all mz and another one of all RT in all samples
  pidx <- groupval(object, value="index")
  colnames(pidx) <- sampNames
  mzs <- p[pidx, "mz"]
  dim(mzs) <- c((length(mzs)/n),n)

  rts <- p[pidx, "rt"]
  dim(rts) <- c((length(rts)/n),n)
  
  ## Calculate deviation between median mz (or RT) and each observed mz (or RT)
  result <- list(mzs=mzs,
                 rts=rts,
                 mzdev=mzs - groups(object)[,"mzmed"],
                 rtdev=rts - groups(object)[,"rtmed"])
}



## plotQC(xsg, what="mzdevhist")
## plotQC(xsg, what="rtdevhist")
## plotQC(xsg, what="mzdevmass")
## plotQC(xsg, what="mzdevtime")
## plotQC(xsg, what="mzdevsample")
## plotQC(xsg, what="rtdevsample")
