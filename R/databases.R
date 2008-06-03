read.metlinMS<- function(xml){
    reading<-readLines(xml)
    xml.mat<-matrix(nrow=length(reading),ncol=3)
    name<-grep("name", reading)
    Pmz<-grep("mass", reading)
    nameVAL<-reading[name]
    PmzVal<-reading[Pmz]

    pattern1<-"\t{3}<mass>(.*)</mass>"
    pattern2<-"\t{3}<name>(.*)</name>"
    PmzVal<-gsub(pattern1, "\\1", PmzVal, perl=T)
    nameVAL<-gsub(pattern2, "\\1", nameVAL, perl=T)

    nameVAL.correct<-nameVAL[2:length(nameVAL)]
    PmzVal.correct<-as.numeric(PmzVal[2:length(PmzVal)])
    met.mat<-cbind(nameVAL.correct, PmzVal.correct)
    colnames(met.mat)<-c("name", "MZ")
    metMS.df<-as.data.frame(met.mat, stringsAsFactors=FALSE)
    metMS.df[,"MZ"]<-as.numeric(metMS.df[,"MZ"])
    metMS.df<-metMS.df[order(metMS.df[,"MZ"]), ]

    return(metMS.df)
}

distance<-function(met, xcm, ppmval, matrix=FALSE){
    l.met<-length(met)
    l.xcm<-length(xcm)
    d<-array(0, dim=c(l.met+1, l.xcm+1))

    d[,1] <- 1:(l.met+1)
    d[1,] <- 1:(l.xcm+1)
    d[1,1] <- 0

    for(i in 2:(l.met+1)){
	for(j in 2:(l.xcm+1)){
		if(ppm(met[i-1], xcm[j-1]) <= ppmval ){ ## ppm cal use
			cost<- 0
			#subcost<-0
			} else {
			cost<- 1
			#subcost<-holdsub
		}
		d[i,j]<- min(d[i-1,j] +cost, ##inserting peak
			     d[i,j-1] +cost, ##deleteing peak
			     d[i-1,j-1] + cost) ##
	}
    }
    if(matrix == TRUE){
        return(d) ##check print whole matrix
    }else{
	return(d[l.met+1, l.xcm+1])
    }
}

similar<-function(met, xcm, ppmval, matrix=FALSE){
    l.met<-length(met)
    l.xcm<-length(xcm)
    d<-array(0, dim=c(l.met+1, l.xcm+1)) #we can cheat and use an AoA:)

    d[,1] <- 1:(l.met+1)
    d[1,] <- 1:(l.xcm+1)
    d[1,1]<-max(l.met,l.xcm) ##Put the max simlarity at the start

    for(i in 2:(l.met+1)){
	for(j in 2:(l.xcm+1)){
		if(ppm(met[i-1], xcm[j-1]) <= ppmval ){ ## ppm cal use
			cost<- 0
			#subcost<-0
			} else {
			cost<- 1
			#subcost<-holdsub
		}
		d[i,j]<- max(d[i-1,j] -cost, ##inserting peak
			     d[i,j-1] -cost, ##deleteing peak
			     d[i-1,j-1] - cost) ## match
	}
    }
    if(matrix == TRUE){
        return(d) ##check print whole matrix
    }else{
	return(max(l.met,l.xcm) - d[l.met+1, l.xcm+1]) ##This give number of similar masses :)
    }
}

ppm<-function(Mr, Mm){ ## Mr Mz real Mm Mz Measured
    ppm<-abs((10^6)*(Mr-Mm)/Mm) ##abs for positive number
    return(ppm)
}

ppmDev<-function(Mr, ppmE=5){
    error<-(ppmE/10^6)*Mr ## just take the ppm as a percentage
    deviation<-c(Mr+error, Mr-error)## 1 is high 2 is low
    return(deviation)
}

read.mascot<-function(file, type="csv"){
    ## Experimental
    if(!file.exists(file))
	stop("File doesnt exist")
    if(type=="csv"){
	lines<-readLines(file) ##check to see if file is ~ normal
	if (grep("Protein hit", lines)) {
		x<-lines[-seq(grep("Protein hit", lines)+1)]
		myDF <-read.csv(textConnection(x), header=T)
	} else {
		stop("Noncompatable csv file\n")
	}
    } else if (type=="xml"){
	stop("XML files are currently not supported, support comming soon")
    }
    return(myDF)
}

KeggSearch <- function(metabo, write=FALSE) {
    ##Experimental
    KEGG<-"http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=100&dbkey=kegg&keywords="
    Mascot<-read.mascot(masfile, type)
    result<-vector()
    for(i in 1:dim(metabo)[1]){
	KEGG<-paste(KEGG,metabo[i,"name"], sep="")
	KEGG<-readLines(url(KEGG),warn=FALSE) ## Don't tell me that EOF was incomplete
	for(j in 1:dim(Mascot)[1,])
		xover<-agrep(Mascot[j,"prot_name"], KEGG) ##Do the best we can grep isn't always happy
		if(xover){
			result[j]<-paste(Mascot[j,"prot_name"], "and",
                                         metabo[i,"name"], "have been found to be in the same pathway.", sep="")
			rm(xover)
		}
    }
	if(write){
		cat(result, file="PathwayMatch.txt", sep="\n")
	}else {
		cat(result, sep="\n")
	}
}

