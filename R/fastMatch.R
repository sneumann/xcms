
fastMatch <- function(x,y,tol=0.001, symmetric=FALSE) {

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

    fm <- .Call("fastMatch", xs, ys, xidx, yidx, as.integer(length(x)), as.double(tol) , PACKAGE ='xcms' )
    fm2 <- vector("list", length=length(fm))
    ##stop("!")
    if (symmetric){
        for (a in 1:length(fm)) {
            if (!is.null(fm[[a]][1])){
                tmp<-NULL
                for (b in 1:length(fm[[a]])){
                    if ((abs(x[a]-y[fm[[a]]][b]) == min(abs(x[a]-y[fm[[a]][b]]),
                            abs(x[a]  -y[fm[[a]][b]-1]),
                            abs(x[a]  -y[fm[[a]][b]+1]),
                            abs(x[a-1]-y[fm[[a]][b]]),
                            abs(x[a+1]-y[fm[[a]][b]]), na.rm=T)
                         )) {
                        tmp<-c(tmp, fm[[a]][b])
                    }
                }
                fm2[[a]]<-tmp
            }
        }
    }else {
        fm2 <- fm}
    fm2
}
