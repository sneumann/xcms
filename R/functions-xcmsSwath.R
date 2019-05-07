## function "align" EICs
.alignEics <- function(x, y) {

  # linear interpolation of y on x
  y_aligned <- approx(rtime(y), intensity(y), rtime(x))
  
  # deal with NA
  # TODO NA to 0
  
  # correct rtime and int in chromatogram object
  y@rtime <- y_aligned$x
  y@intensity <- y_aligned$y
  
  # return
  return(y)
  
}

# function to correlate EICs
.correlateEics <- function(x, y) {
  
  # check if x and y are on the same time scale
  if(!rtime(x) == rtime(y)) {
    y <- .alignEics(x, y)
  }
  
  # correlate values
  corValues <- cor.test(intensity(x), intensity(y))
  
  # return
  return(list(corValues$estimate, corValues$p.value))
  
}
