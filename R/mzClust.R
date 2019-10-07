mzClust_hclust <- function(x, eppm, eabs)
{
    N <- length(x)
    d <- dist(x)
    g <- .C("R_mzClust_hclust",
            x = as.double(x),
            num = N,
            d = as.double(d),
            g = integer(length = N),
            eppm = as.double(eppm),
            eabs = as.double(eabs), NAOK = TRUE)$g
    return(g)
}

mzClustGeneric <- function(p,sampclass=NULL,
                           mzppm = 20,
                           mzabs = 0,
                           minsamp = 1,
                           minfrac=0.5)
{
    makeBin <- function(pos)
    {
        if(pos > numpeaks)
            return(list(pos=pos,bin=c(-1)))

        bin <- pord[pos]
        pos <- pos+1
        basepeak <- p[bin[1],1]
        error_range <- c(basepeak,
                         basepeak * error_window + basepeak + 2 * mzabs)
        while(pos < numpeaks && p[pord[pos], 1] <= error_range[2]) {
            bin <- c(bin, pord[pos])
            pos <- pos + 1
        }

        if(pos %% (numpeaks%/%100+1) == 0) {
            cat(format(((pos-1)/numpeaks*100),digits=1,nsmall=2)," ")
            flush.console()
        }

        lst <- list(pos=pos,bin=bin)
        lst
    }
    meanDeviationOverLimit <- function(bin)
    {
        bin_mz <- p[bin,1]
        m <- mean(bin_mz)
        error_range <- c(m-ppm_error*m-mzabs, ppm_error*m+m+mzabs)
        if(length(bin_mz[(bin_mz > error_range[2]) |
                         (bin_mz < error_range[1])]) > 0 ) {
            return(TRUE)
        } else { FALSE }
    }
    bin2output <- function(bin)
    {

        gcount <- integer(length(classnum))
        if(length(gcount) != 0) {
            for(i in seq(along = bin)) {
                class_idx <- sampclass[p[bin[i],2]]
                gcount[class_idx] <- gcount[class_idx] + 1
            }
        }
        ## make sure, that if no classes given, 'any' is false
        if(length(bin) < minsamp || (!any(gcount >= classnum*minfrac)
                                     && length(gcount)>0))
            return(list())
        groupvec <- c(rep(NA,4+length(gcount)))
        groupvec[1] <- mean(p[bin,1])
        groupvec[2:3] <- range(p[bin,1])
        groupvec[4] <- length(bin)
        sorted <- order(p[bin,1])
        grp_members <- bin[sorted]
        groupvec[4+seq(along = gcount)] <- gcount
        lst <- list(stat=groupvec,members=grp_members)
        lst
    }
    ppm_error <- mzppm / 1000000
    error_window <- 2 * ppm_error

    ## numeric version of classlabel
    if(is.null(sampclass)) {
        classnum <- integer(0)
        classnames <- seq(along=classnum)
    } else {
        classnames <- levels(sampclass)
        sampclass <- as.vector(unclass(sampclass))

        classnum <- integer(max(sampclass))
    }

    for (i in seq(along = classnum))
        classnum[i] <- sum(sampclass == i)

    numpeaks <- nrow(p)

    groupmat <- matrix(nrow = 512, ncol = 4 + length(classnum))
    groupindex <- vector("list", 512)

    pord <- order(p[,1])
    pos <- c(1)
    binNumber <- 1
    newbin <- makeBin(pos)
    binA <- newbin$bin
    pos <- newbin$pos
    while(1) {
        if (binNumber +4 > nrow(groupmat)) {
            groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat),
                                               ncol = ncol(groupmat)))
            groupindex <- c(groupindex, vector("list", length(groupindex)))
        }
        newbin <- makeBin(pos)
        binB <- newbin$bin
        pos <- newbin$pos



        if(binB[1] < 0) {
            out <- bin2output(binA)
            if(length(out) != 0) {
                groupmat[binNumber,] <- out$stat
                groupindex[[binNumber]] <- out$members
                binNumber <- binNumber + 1
            }
            break
        }
        max_binA <- max(p[binA,1])
        min_binB <- min(p[binB,1])

        binclust <- 0
        if(max_binA + max_binA*error_window+2*mzabs >= min_binB &&
           min_binB - min_binB*error_window -2*mzabs <= max_binA) {
            binC <- c(binA,binB)
            binclust <- 1
        } else {
            if(meanDeviationOverLimit(binA)) {
                binC <- binA
                binclust <- 2
            }
        }

        ## case: not in range or mean deviation over limit
        ## perform hierarchical clustering
        if(binclust != 0) {
            groups <- mzClust_hclust(p[binC,1],ppm_error,mzabs)

            last_group <- groups[which.max(p[binC,1])]
            binA <- binC[which(groups == last_group)]
            if(max(groups) >1) {
                for(c in 1:max(groups)) {
                    if(c == last_group) { next }
                    tmp_grp <- which(groups == c)
                    tmp_c <- binC[tmp_grp]
                    out <- bin2output(tmp_c)
                    if(length(out) != 0) {
                        groupmat[binNumber,] <- out$stat
                        groupindex[[binNumber]] <- out$members
                        binNumber <- binNumber + 1
                    }
                }
            }
        }

        if(binclust != 1) {
            out <- bin2output(binA)
            if(length(out) != 0) {
                groupmat[binNumber,] <- out$stat
                groupindex[[binNumber]] <- out$members
                binNumber <- binNumber + 1
            }
            binA <- binB
        }


    }
    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "npeaks", classnames)

    binNumber <- binNumber - 1
    groupmat <- groupmat[seq(length = binNumber),]
    groupindex <- groupindex[seq(length = binNumber)]
    cat("\n")
    flush.console()
    return(list(mat=groupmat,idx=groupindex))
}
