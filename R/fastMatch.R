
fastMatch <- function(x,y,tol=0.001) {

    if (any(is.na(y)))
        stop("NA's are not allowed in y !\n")
    
    ok <- !(is.na(x))
    ans <- order(x)
    keep <- seq_along(ok)[ok]
    xidx <- ans[ans %in% keep]
    xs <- x[xidx]
    
    yidx <- order(y)
    ys <- y[yidx] 

    if (!is.double(xs)) 
      xs <- as.double(xs)
    if (!is.double(ys)) 
      ys<- as.double(ys)
    
    if (!is.integer(xidx)) 
      xidx <- as.integer(xidx)
    if (!is.integer(yidx)) 
      yidx <- as.integer(yidx)  

  .Call("fastMatch", xs, ys, xidx, yidx, as.integer(length(x)), as.double(tol) , PACKAGE ='CAMERA' )
  
}
