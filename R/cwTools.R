cwt <- function (ms, scales = 1, wavelet = "mexh") 
{ ## modified from package MassSpecWavelet 
    if (wavelet == "mexh") {
        psi_xval <- seq(-8, 8, length = 1024)
        psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * 
            exp(-psi_xval^2/2)
    }
    else if (is.matrix(wavelet)) {
        if (nrow(wavelet) == 2) {
            psi_xval <- wavelet[1, ]
            psi <- wavelet[2, ]
        }
        else if (ncol(wavelet) == 2) {
            psi_xval <- wavelet[, 1]
            psi <- wavelet[, 2]
        }
        else {
            stop("Unsupported wavelet format!")
        }
    }
    else {
        stop("Unsupported wavelet!")
    }
    oldLen <- length(ms)
    ms <- extendNBase(ms, nLevel = NULL, base = 2)
    len <- length(ms)
    nbscales <- length(scales)
    wCoefs <- NULL
    psi_xval <- psi_xval - psi_xval[1]
    dxval <- psi_xval[2]
    xmax <- psi_xval[length(psi_xval)]
    for (i in 1:length(scales)) {
        scale.i <- scales[i]
        f <- rep(0, len)
        j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
        if (length(j) == 1) 
            j <- c(1, 1)
        lenWave <- length(j)
        f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
        if (length(f) > len) 
            {i<-i-1;break;}   #  stop(paste("scale", scale.i, "is too large!"))
        wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
        wCoefs.i <- c(wCoefs.i[(len - floor(lenWave/2) + 1):len], 
            wCoefs.i[1:(len - floor(lenWave/2))])
        wCoefs <- cbind(wCoefs, wCoefs.i)
    }
    if (i < 1) return(NA)
    scales <- scales[1:i]
    if (length(scales) == 1) 
        wCoefs <- matrix(wCoefs, ncol = 1)
    colnames(wCoefs) <- scales
    wCoefs <- wCoefs[1:oldLen, , drop = FALSE]
    wCoefs
}

extendNBase <- function(x, nLevel=1, base=2, ...) 
{ ## from package MassSpecWavelet 
	if (!is.matrix(x)) x <- matrix(x, ncol=1)	
	
	nR <- nrow(x)
	if (is.null(nLevel)) {
		nR1 <- nextn(nR, base)		
	} else {
		nR1 <- ceiling(nR / base^nLevel) * base^nLevel		
	}
	if (nR != nR1) {
		x <- extendLength(x, addLength=nR1-nR, ...)
	}
	x
}

extendLength <-
function(x, addLength=NULL, method=c('reflection', 'open', 'circular'), direction=c('right', 'left', 'both')) 
{       ## from package MassSpecWavelet
	if (is.null(addLength)) stop('Please provide the length to be added!')
	if (!is.matrix(x)) x <- matrix(x, ncol=1)	
	method <- match.arg(method)
	direction <- match.arg(direction)
	
	nR <- nrow(x)
	nR1 <- nR + addLength
	if (direction == 'both') {
		left <- right <- addLength
	} else if (direction == 'right') {
		left <- 0
		right <- addLength
	} else if (direction == 'left') {
		left <- addLength
		right <- 0
	}
	
	if (right > 0) {
		x <- switch(method,
			reflection =rbind(x, x[nR:(2 * nR - nR1 + 1), , drop=FALSE]),
			open = rbind(x, matrix(rep(x[nR,], addLength), ncol=ncol(x), byrow=TRUE)),
			circular = rbind(x, x[1:(nR1 - nR),, drop=FALSE]))
	}

	if (left > 0) {
		x <- switch(method,
			reflection =rbind(x[addLength:1, , drop=FALSE], x),
			open = rbind(matrix(rep(x[1,], addLength), ncol=ncol(x), byrow=TRUE), x),
			circular = rbind(x[(2 * nR - nR1 + 1):nR,, drop=FALSE], x))
	}
	if (ncol(x) == 1)  x <- as.vector(x)
	
	x
}

getLocalMaximumCWT <-
function(wCoefs, minWinSize=5, amp.Th=0) 
{        ## from package MassSpecWavelet 
	localMax <- NULL
	scales <- as.numeric(colnames(wCoefs))
	
	for (i in 1:length(scales)) {
		scale.i <- scales[i]
		winSize.i <- scale.i * 2 + 1
		if (winSize.i < minWinSize) {
			winSize.i <- minWinSize
		} 
		temp <- localMaximum(wCoefs[,i], winSize.i)
		localMax <- cbind(localMax, temp)
	}
	# Set the values less than peak threshold as 0
	localMax[wCoefs < amp.Th] <- 0
	colnames(localMax) <- colnames(wCoefs)
	rownames(localMax) <- rownames(wCoefs)
	localMax
}

localMaximum <-
function (x, winSize = 5) 
{   ## from package MassSpecWavelet
	len <- length(x)
	rNum <- ceiling(len/winSize)
	
	## Transform the vector as a matrix with column length equals winSize
	##		and find the maximum position at each row.
	y <- matrix(c(x, rep(x[len], rNum * winSize - len)), nrow=winSize)
	y.maxInd <- apply(y, 2, which.max)
	## Only keep the maximum value larger than the boundary values
	selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
	
	## keep the result
	localMax <- rep(0, len)
	localMax[(selInd-1) * winSize + y.maxInd[selInd]] <- 1
	
	## Shift the vector with winSize/2 and do the same operation
	shift <- floor(winSize/2)
	rNum <- ceiling((len + shift)/winSize)	
	y <- matrix(c(rep(x[1], shift), x, rep(x[len], rNum * winSize - len - shift)), nrow=winSize)
	y.maxInd <- apply(y, 2, which.max)
	## Only keep the maximum value larger than the boundary values
	selInd <- which(apply(y, 2, function(x) max(x) > x[1] & max(x) > x[winSize]))
	localMax[(selInd-1) * winSize + y.maxInd[selInd] - shift] <- 1
	
	## Check whether there is some local maxima have in between distance less than winSize
	maxInd <- which(localMax > 0)
	selInd <- which(diff(maxInd) < winSize)
	if (length(selInd) > 0) {
		selMaxInd1 <- maxInd[selInd]
		selMaxInd2 <- maxInd[selInd + 1]
		temp <- x[selMaxInd1] - x[selMaxInd2]
		localMax[selMaxInd1[temp <= 0]] <- 0
		localMax[selMaxInd2[temp > 0]] <- 0
	}
	
	localMax
}

getRidge <-
function(localMax, iInit=ncol(localMax), step=-1, iFinal=1, minWinSize=5, gapTh=3, skip=NULL) 
{    ## modified from package MassSpecWavelet

	scales <- as.numeric(colnames(localMax))
	if (is.null(scales))  scales <- 1:ncol(localMax)
	
	maxInd_curr <- which(localMax[, iInit] > 0)
	nMz <- nrow(localMax)
	
	if (is.null(skip))	{
		skip <- iInit + 1
	}

	## Identify all the peak pathes from the coarse level to detail levels (high column to low column)
	## Only consider the shortest path
	if ( ncol(localMax) > 1 ) colInd <- seq(iInit+step, iFinal, step)
	 else colInd <- 1
	ridgeList <- sapply(maxInd_curr, as.list)
	names(ridgeList) <- maxInd_curr
	peakStatus <- sapply(rep(0, length(maxInd_curr)), as.list)
	names(peakStatus) <- maxInd_curr
	
	## orphanRidgeList keep the ridges disconnected at certain scale level
	## Changed by Pan Du 05/11/06
	orphanRidgeList <- NULL
	orphanRidgeName <- NULL
	nLevel <- length(colInd)

	for (j in 1:nLevel) {
		col.j <- colInd[j]
		scale.j <- scales[col.j]

		if (colInd[j] == skip) {
			oldname <- names(ridgeList)
			ridgeList <- sapply(ridgeList, function(x) c(x, x[length(x)]))
			#peakStatus <- sapply(peakStatus, function(x) c(x, x[length(x)]))
			names(ridgeList) <- oldname
			#names(peakStatus) <- oldname
			next
		}
		
		if (length(maxInd_curr) == 0) {
			maxInd_curr <- which(localMax[, col.j] > 0)
			next
		}

		## The slide window size is proportional to the CWT scale
		winSize.j <- scale.j * 2 + 1
		if (winSize.j < minWinSize) {
			winSize.j <- minWinSize
		}

		selPeak.j <- NULL
		remove.j <- NULL
		for (k in 1:length(maxInd_curr)) {
			ind.k <- maxInd_curr[k]
			start.k <- ifelse(ind.k-winSize.j < 1, 1, ind.k-winSize.j)
			end.k <- ifelse(ind.k+winSize.j > nMz, nMz, ind.k+winSize.j)
			ind.curr <- which(localMax[start.k:end.k, col.j] > 0) + start.k - 1
			#ind.curr <- which(localMax[, col.j] > 0)
			if (length(ind.curr) == 0) {
				status.k <- peakStatus[[as.character(ind.k)]]
				if (status.k > gapTh & scale.j >= 2) {
					temp <- ridgeList[[as.character(ind.k)]]
					orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp)-status.k)]))
					orphanRidgeName <- c(orphanRidgeName, paste(col.j + status.k + 1, ind.k, sep='_'))
					remove.j <- c(remove.j, as.character(ind.k))  
					next
				} else {
					ind.curr <- ind.k
					peakStatus[[as.character(ind.k)]] <- status.k + 1
				}
			} else {
				peakStatus[[as.character(ind.k)]] <- 0
				if (length(ind.curr) >= 2)  ind.curr <- ind.curr[which.min(abs(ind.curr - ind.k))]
			}
			ridgeList[[as.character(ind.k)]] <- c(ridgeList[[as.character(ind.k)]], ind.curr)
			selPeak.j <- c(selPeak.j, ind.curr)
		}
		## Remove the disconnected lines from the currrent list
		if (length(remove.j) > 0) {
			removeInd <- which(names(ridgeList) %in% remove.j)
			ridgeList <- ridgeList[-removeInd]
			peakStatus <- peakStatus[-removeInd]
		}

		## Check for duplicated selected peaks and only keep the one with the longest path.
		dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
		if (length(dupPeak.j) > 0) {
			removeInd <- NULL
			for (dupPeak.jk in dupPeak.j) {
				selInd <- which(selPeak.j == dupPeak.jk)
				selLen <- sapply(ridgeList[selInd], length)
				removeInd.jk <- which.max(selLen)
				removeInd <- c(removeInd, selInd[-removeInd.jk])
				orphanRidgeList <- c(orphanRidgeList, ridgeList[removeInd.jk])
				orphanRidgeName <- c(orphanRidgeName, paste(col.j, selPeak.j[removeInd.jk], sep='_'))
			}
			selPeak.j <- selPeak.j[-removeInd]
			ridgeList <- ridgeList[-removeInd]
			peakStatus <- peakStatus[-removeInd]
		}
		
		## Update the names of the ridgeList as the new selected peaks
		#if (scale.j >= 2) {
			names(ridgeList) <- selPeak.j
			names(peakStatus) <- selPeak.j			
		#}
		
		## If the level is larger than 3, expand the peak list by including other unselected peaks at that level
		if (scale.j >= 2) {
			maxInd_next <- which(localMax[, col.j] > 0)
			unSelPeak.j <- maxInd_next[!(maxInd_next %in% selPeak.j)]
			newPeak.j <- sapply(unSelPeak.j, as.list)
			names(newPeak.j) <- unSelPeak.j
			## Update ridgeList
			ridgeList <- c(ridgeList, newPeak.j)
			maxInd_curr <- c(selPeak.j, unSelPeak.j)
			## Update peakStatus
			newPeakStatus <- sapply(rep(0, length(newPeak.j)), as.list)
			names(newPeakStatus) <- newPeak.j
			peakStatus <- c(peakStatus, newPeakStatus)
		} else {
			maxInd_curr <- selPeak.j
		}
	}
	## Attach the peak level at the beginning of the ridge names
	names(ridgeList) <- paste(1, names(ridgeList), sep='_')
	names(orphanRidgeList) <- orphanRidgeName
	## Combine ridgeList and orphanRidgeList
	ridgeList <- c(ridgeList, orphanRidgeList)
	
	## Reverse the order as from the low level to high level.
	ridgeList <- lapply(ridgeList, rev)
	## order the ridgeList in increasing order
	#ord <- order(selPeak.j)
	#ridgeList <- ridgeList[ord]
	
	## Remove possible duplicated ridges
	ridgeList <- ridgeList[!duplicated(names(ridgeList))]

	attr(ridgeList, 'class') <- 'ridgeList'
	attr(ridgeList, 'scales') <- scales
	ridgeList
}

rmean <- function (X,width = min(length(X), 20))
{ ## modified running() from package gtools
    width <- min(length(X),width)
    n <- length(X)
    from <- 1:n
    to <- pmin((1:n) + width - 1, n)
    elements <- apply(cbind(from, to), 1, function(x) seq(x[1], x[2]))
    len <- sapply(elements, length)
    skip <- (len < width)
    elements <- elements[!skip]
    funct <- function(which, what, fun) mean(what[which])
    sapply(elements, funct, what = X, fun = mean)
}

odd <- function (x) x != as.integer(x/2) * 2;

gauss <- function(x, h, mu, sigma){
    h*exp(-(x-mu)^2/(2*sigma^2))
}

fitGauss <- function(td,d,pgauss=NA) {
 if (!any(is.na(pgauss))) { mu <- pgauss$mu; sigma <- pgauss$sigma;h <- pgauss$h }
 fit <- try(nls(d ~ SSgauss(td,mu,sigma,h)), silent = TRUE)
 if (class(fit) == "try-error") 
    fit <- try(nls(d ~ SSgauss(td,mu,sigma,h),algo='port'), silent = TRUE)
 if (class(fit) == "try-error")  return(rep(NA,3))

 as.data.frame(t(fit$m$getPars()))
}

wnoisedet <- function(x,lev,num) {
  cnt <- 0; cpos <- NULL; fpos <- NULL
  for (k in 1:length(x)) {
    if (x[k] > lev) {cnt <- cnt + 1; cpos <- c(cpos,k)}
    else {
      if (cnt >= num) fpos <- c(fpos,cpos) 
      cnt <- 0
      cpos <- NULL
      }   
  }
  fpos
} 

noiseQuant <- function(d) {
 quantile(d,1/3)
}

remove.NA.features <- function(peaks) {
  rmi <- which(is.na(peaks[,"mu"]))
  if (length(rmi) >0) return(peaks[-rmi,]) else return(peaks)
}

joinOverlappingFeatures <- function(td,d,scantime,scan.range,peaks,maxGaussOverlap=0.5,maxGaussErr=0.3) {  
  peaks <- as.data.frame(peaks)
  gausspeaks <- which(!is.na(peaks[,"mu"])); Ngp <- length(gausspeaks)
  newpeaks <- peaks
  if (Ngp > 1) {
    comb <- which(upper.tri(matrix(0,Ngp,Ngp)),arr=TRUE)
    remove <- overlap <- rep(FALSE,Ngp)
    for (k in 1:dim(comb)[1]) {
        p1 <- gausspeaks[comb[k,1]]; p2 <- gausspeaks[comb[k,2]]
        overlap[k] <- gaussCoverage(xlim=scan.range,h1=peaks[p1,"h"],mu1=peaks[p1,"mu"],s1=peaks[p1,"sigma"], 
            h2=peaks[p2,"h"],mu2=peaks[p2,"mu"],s2=peaks[p2,"sigma"]) >= maxGaussOverlap
    }
  
    if (any(overlap)) {
      if (all(overlap)) {
        ## simple case, use border maxima
        newmin <- min(peaks[,"lmin"])
        newmax <- max(peaks[,"lmax"])
        newpeaks <- peaks[1,]
        newpeaks["scpos"] <- -1; newpeaks["scmin"] <- -1; newpeaks["scmax"] <- -1; ## not defined after join
        newpeaks["maxo"] <- max(peaks[,"maxo"])  ; newpeaks["sn"] <- max(peaks[,"sn"])
        newpeaks["lmin"] <- newmin   ; newpeaks["lmax"] <-  newmax 
        newpeaks["rtmin"] <- scantime[td[newmin]]   ; newpeaks["rtmax"] <- scantime[td[newmax]]     
        newpeaks["rt"] <- weighted.mean(peaks[,"rt"],w=peaks[,"maxo"])
        ## re-fit gaussian
        md <- max(d[newmin:newmax]);d1 <- d[newmin:newmax]/md;
        pgauss <- fitGauss(td[newmin:newmax],d[newmin:newmax],pgauss = list(mu=td[newmin] + (td[newmax]-td[newmin])/2,sigma=td[newmax]-td[newmin],h=max(peaks[,"h"])))
        if (!any(is.na(pgauss)) && all(pgauss > 0)) {
              newpeaks["mu"] <- pgauss$mu; newpeaks["sigma"] <- pgauss$sigma; newpeaks["h"] <- pgauss$h;  
              newpeaks["egauss"] <- sqrt((1/length(td[newmin:newmax])) * sum(((d1-gauss(td[newmin:newmax],pgauss$h/md,pgauss$mu,pgauss$sigma))^2)))
        } else stop('Panic: re-fit after join (all) failed.')
      } else {
        ## check each combination, join if overlapping, mark smaller peak as deleted
        overlaps <- which(overlap)
        newpeaks <- peaks
        for (o in overlaps) {
          p1 <- gausspeaks[comb[o,1]]; p2 <- gausspeaks[comb[o,2]]
          jp <- peaks[p1,]
          newmin <- min(peaks[c(p1,p2),"lmin"])
          newmax <- max(peaks[c(p1,p2),"lmax"])
          jp["scpos"] <- -1; jp["scmin"] <- -1; jp["scmax"] <- -1; ## not defined after join
          jp["maxo"] <- max(peaks[c(p1,p2),"maxo"])  
          jp["sn"] <- max(peaks[c(p1,p2),"sn"])
          jp["lmin"] <- newmin; jp["lmax"] <-  newmax 
          jp["rtmin"] <- scantime[td[newmin]]   ; jp["rtmax"] <- scantime[td[newmax]]     
          jp["rt"] <- weighted.mean(peaks[c(p1,p2),"rt"],w=peaks[c(p1,p2),"maxo"])
          if (peaks[p1,"h"] < peaks[p2,"h"]) { remove[p1] <- TRUE; newpeaks[p2,] <- jp } else {remove[p2] <- TRUE; newpeaks[p1,] <- jp}
        }
        ## collect remaining peaks
        newpeaks <- newpeaks[-which(remove),]
        ## re-fit joined peaks
        jpi <- which(newpeaks[,"scpos"] == -1) ## joined peaks
        for (j in jpi) {
          newmin <- newpeaks[j,"lmin"]; newmax <- newpeaks[j,"lmax"]
          md <- max(d[newmin:newmax]);d1 <- d[newmin:newmax]/md;
          pgauss <- fitGauss(td[newmin:newmax],d[newmin:newmax],pgauss = list(mu=td[newmin] + (td[newmax]-td[newmin])/2,sigma=td[newmax]-td[newmin],h=newpeaks[j,"h"]))
          if (!any(is.na(pgauss)) && all(pgauss > 0)) {
                newpeaks[j,"mu"] <- pgauss$mu; newpeaks[j,"sigma"] <- pgauss$sigma; newpeaks[j,"h"] <- pgauss$h;  
                newpeaks[j,"egauss"] <- sqrt((1/length(td[newmin:newmax])) * sum(((d1-gauss(td[newmin:newmax],pgauss$h/md,pgauss$mu,pgauss$sigma))^2)))
          } else stop('Panic: re-fit after join (2) failed.')
          if (newpeaks[j,"egauss"] > maxGaussErr) stop('Panic: re-fit after join (2) failed, err > maxGaussErr')
        }
      }
    } # any overlap
  }
 as.matrix(newpeaks)
}

descendMinTol <- function(d,startpos,maxDescOutlier) {
  l <- startpos[1]; r <- startpos[2]; outl <- 0; N <- length(d)
  # left
  while ((l > 1) && (d[l] > 0) && outl <= maxDescOutlier) {
      if (outl > 0) vpos <- opos else vpos <- l
      if (d[l-1] > d[vpos]) outl <- outl + 1 else outl <- 0
      if (outl == 1) opos <- l
      l <- l -1
  }
  if (outl > 0) l <- l + outl
  # right
  outl <- 0;
  while ((r < N) && (d[r] > 0) && outl <= maxDescOutlier) {
      if (outl > 0) vpos <- opos else vpos <- r
      if (d[r+1] > d[vpos]) outl <- outl + 1 else outl <- 0
      if (outl == 1) opos <- r
      r <- r + 1
  } 
  if (outl > 0) r <- r - outl
  c(l,r)
}

cent <- function(x) {
  floor(length(x)/2)
}

gaussCoverage <- function(xlim,h1,mu1,s1,h2,mu2,s2) { 
  overlap <- NA
  by = 0.05
  a <- s2^2 - s1^2
  c <- -( 2 * s1^2 * s2^2 * (log(h1) - log(h2)) + (s1^2 * mu2^2) - (s2^2 * mu1^2) )
  b <- ((2 * s1^2 *mu2) - (2 * s2^2 * mu1))
  D <- b^2 - (a*c)
  if (a==0) {S1 <- -c/b; S2 <- NA
   } else if ((D < 0) || ((b^2 - (4*a*c)) < 0)) {S1 <- S2 <- NA
     } else {
        S1 <- (-b + sqrt(b^2 - (4*a*c))) / (2*a)
        S2 <- (-b - sqrt(b^2 - (4*a*c))) / (2*a)
        if (S2 < S1) {tmp<-S1; S1<-S2; S2<-tmp}
      }  
  if (!is.na(S1)) if (S1 < xlim[1] || S1 > xlim[2]) S1 <- NA
  if (!is.na(S2)) if (S2 < xlim[1] || S2 > xlim[2]) S2 <- NA

  x <- seq(xlim[1],xlim[2],by=by)
  vsmall <- min(sum(gauss(x,h1,mu1,s1)), sum(gauss(x,h2,mu2,s2)))
  
  if (!is.na(S1) && !is.na(S2)) {
    x0 <- seq(xlim[1],S1,by=by)
    xo <- seq(S1,S2,by=by)
    x1 <- seq(S2,xlim[2],by=by)
    if (gauss(x0[cent(x0)],h1,mu1,s1) < gauss(x0[cent(x0)],h2,mu2,s2)) ov1 <- sum(gauss(x0,h1,mu1,s1))                    else ov1 <- sum(gauss(x0,h2,mu2,s2))
    if (gauss(xo[cent(xo)],h1,mu1,s1) < gauss(xo[cent(xo)],h2,mu2,s2)) ov <- sum(gauss(xo,h1,mu1,s1))                    else ov <- sum(gauss(xo,h2,mu2,s2)) 
    if (gauss(x1[cent(x1)],h1,mu1,s1) < gauss(x1[cent(x1)],h2,mu2,s2)) ov2 <- sum(gauss(x1,h1,mu1,s1))                    else ov2 <- sum(gauss(x1,h2,mu2,s2)) 
    overlap <- ov1 + ov + ov2
  } else
  if (is.na(S1) && is.na(S2)) { # no overlap -> intergrate smaller function
     if (gauss(x[cent(x)],h1,mu1,s1) < gauss(x[cent(x)],h2,mu2,s2)) overlap <- sum(gauss(x,h1,mu1,s1))                    else overlap <- sum(gauss(x,h2,mu2,s2))
  } else
  if (!is.na(S1) || !is.na(S2)) {
    if (is.na(S1)) S0 <- S2 else S0 <- S1
    x0 <- seq(xlim[1],S0,by=by)
    x1 <- seq(S0,xlim[2],by=by)
    if (gauss(x0[cent(x0)],h1,mu1,s1) < gauss(x0[cent(x0)],h2,mu2,s2)) ov1 <- sum(gauss(x0,h1,mu1,s1))                    else ov1 <- sum(gauss(x0,h2,mu2,s2))
    if (gauss(x1[cent(x1)],h1,mu1,s1) < gauss(x1[cent(x1)],h2,mu2,s2)) ov2 <- sum(gauss(x1,h1,mu1,s1))                    else ov2 <- sum(gauss(x1,h2,mu2,s2))
    overlap <- ov1 + ov2
  }
  
  overlap / vsmall
}
