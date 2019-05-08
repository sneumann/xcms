## function "align" EICs
# are the MS1 chroms, y are the MS2
.align_chromatogram <- function(x, y, na.value = NA) {

  # linear interpolation of y on x
  y_aligned <- approx(rtime(y), intensity(y), rtime(x))
  
  # deal with NA
  if(!is.na(na.value)) {
    y_aligned[is.na()] <- na.value
  }
  
  # correct rtime and int in chromatogram object
  y@rtime <- y_aligned$x
  y@intensity <- y_aligned$y
  
  # return
  return(y)
  
}

# function to correlate EICs
# are the MS1 chroms, y are the MS2
# standard for use is "pairwise.complete.obs", which uses only pairs for which both values are given
.correlate_chromatogram <- function(x, y,
                                    na.value = NA,
                                    use = "pairwise.complete.obs",
                                    method = c("pearson", "kendall", "spearman")) {
  
  # check if x and y are on the same time scale
  if(!all(rtime(x) == rtime(y))) {
    y <- .align_chromatogram(x, y,
                             na.value = na.value)
  }
  
  # correlate values
  corCoeff <- cor(intensity(x),
                  intensity(y),
                  use = use,
                  method = method)
  
  # return
  return(corCoeff)
  
}

## wrapper function for multiple chromatograms
# are the MS1 chroms, y are the MS2
.correlate_chromatograms <-function(x, y) {
  
  # check if same amount of chroms are in x and y
  if(!ncol(x) == ncol(y)) {
    stop("Number of samples different between MS1 and MS2 chroms")
  }
  
  # make matrix to store correlation values
  corMatrix <- matrix(0, nrow = nrow(y), ncol = ncol(y))
  
  # iterate and make correlation
  # TODO
  # use biocParallel backend???
  
  # return
  return(corMatrix)
  
}