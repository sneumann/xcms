matchpeaks <- function(peaklist, masslist, mzabs=0.0001, mzppm=5, neighbours=3, intensity="into"){
    masslist <- masslist[order(masslist)]
    dif <- matrix(nrow=length(masslist), ncol=neighbours, data=rep(1000,length(masslist)*neighbours))
    pos <- matrix(nrow=length(masslist), ncol=neighbours, data=rep(0,length(masslist)*neighbours))
    mzu <- peaklist[,"mz"]
    itu <- peaklist[,intensity]
    mdif <- NA
    mpos <- NA
    cdifs <- NA
    cposs <- NA
    lastval<-1
    for (b in 1:length(masslist)){
        ## finding for each entry in masslist the ##neighbours closest masses in the spectra
        a <- lastval
        sf <- FALSE ## is set to true if the neighbour-box is filled the first time for a entry in masslist
        while (a < length(mzu)){
            mxd <- dif[b,max(which(abs(dif[b,])==max(abs(dif[b,]))))] ## the currently biggest distance
            mxp <-       max(which(abs(dif[b,])==max(abs(dif[b,])))) ## the position of the distance
            if (abs(mzu[a]-masslist[b]) <= abs(mxd)){
                dif[b,mxp] <- (mzu[a]-masslist[b])
                pos[b,mxp] <- a
                lastval <- min(pos[b,])
                sf <- TRUE
            }else{ ## no more masses are smaller, switching to end of mzu
                if (sf) a <- length(mzu)
            }
            a <- a + 1
        }
        ## cat(min(abs(dif[b,])),"\n")
        ## checking the treshold and finding the candidate with the biggest intensity
        smalldiffs <- which(abs(dif[b,]) <= (mzabs+(mzu[pos[b,]]/1000000*mzppm)))
        cdifs <- dif[b,smalldiffs]
        cposs <- pos[b,smalldiffs]
        if (length(cdifs)>0){
            mcdi <- smalldiffs[which(itu[cposs] == max(itu[cposs]))]
            mdif[b] <- dif[b,mcdi]
            mpos[b] <- pos[b,mcdi]
        }else{
            mdif[b]<-NA
            mpos[b]<-NA
        }
    }
    mdiffs <- mdif[!is.na(mdif)]
    mposs <- mpos[!is.na(mdif)]
    retdata <- matrix(ncol=3, nrow=length(mdiffs), data=c(mposs,mzu[mposs],mdiffs))
    colnames(retdata) = c("pos","mass","dif")
    retdata
}

estimate <- function(dtable,method="linear"){
    mdiffs <- dtable[,"dif"]
    mposs <- dtable[,"mass"]

    if (method!="shift"){
        ## complete regression
        regg <- lm(mdiffs ~ mposs)
        a <- regg$coefficients[2]
        b <- regg$coefficients[1]
    }else{
        ## only global shift
        a <- 0
        b <- mean(mdiffs)
    }
    if (method != "shift") return (c(a,b))
    else return(c(0,b))
}
