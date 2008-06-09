read.metlin<-function(xml, MS1=TRUE, MSxml) {
#Parsing the METLIN XML File
    #cat(paste("Reading Metlin...", sep="")) #if metlin gets really big then add this so user knows what's happening :)
    reading<-readLines(xml)
    xml.mat<-matrix(nrow=length(reading),ncol=3)
    name<-grep("name", reading)# OK it's messy like hell but I didn't want to add the XML package for XCMS
    mz<-grep("mz", reading) #Anyway the XML package sucks!
    int<-grep("intensity", reading)
    Pmz<-grep("mass", reading)
    mode<-grep("mode", reading)
    adduct<-grep("adduct", reading)
    CE<-grep("collisionE", reading)

    mzVAL<-reading[mz]
    intVAL<-reading[int]
    nameVAL<-reading[name]
    PmzVal<-reading[Pmz]
    adductVAL<-reading[adduct]
    modeVAL<-reading[mode]
    CEVAL<-reading[CE]

    pattern<-"\t{3}<mz>(.*)</mz>"# "\t{3}(?:<mz>|<name>(?????)(?:</mz>|</name>))"
    pattern1<-"\t{3}<mass>(.*)</mass>"
    pattern2<-"\t{3}<name>(.*)</name>"
    pattern3<-"\t{3}<mode>(.*)</mode>"
    pattern4<-"\t{3}<adduct>(.*)</adduct>"
    pattern5<-"\t{3}<intensity>(.*)</intensity>"
    pattern6<-"\t{3}<collisionE>(.*)</collisionE>"

    mzVAL<-gsub(pattern, "\\1", mzVAL, perl=T)
    intVAL<-gsub(pattern5, "\\1", intVAL, perl=T)
    PmzVal<-gsub(pattern1, "\\1", PmzVal, perl=T)
    nameVAL<-gsub(pattern2, "\\1", nameVAL, perl=T)
    modeVAL<-gsub(pattern3, "\\1", modeVAL, perl=T)
    adductVAL<-gsub(pattern4, "\\1", adductVAL, perl=T)
    CEVAL<-gsub(pattern6, "\\1", CEVAL, perl=T)

    mzVAL.correct<-as.numeric(mzVAL[2:length(mzVAL)])
    intVAL.correct<-intVAL[2:length(intVAL)]
    nameVAL.correct<-nameVAL[2:length(nameVAL)]
    PmzVal.correct<-as.numeric(PmzVal[2:length(PmzVal)])
    mode.correct<-modeVAL[2:length(modeVAL)]
    adduct.correct<-adductVAL[2:length(adductVAL)]
    CE.correct<-CEVAL[2:length(CEVAL)]

    met.mat<-cbind(nameVAL.correct, mzVAL.correct, intVAL.correct, PmzVal.correct, mode.correct, adduct.correct, CE.correct)

    met.mat<-as.data.frame(met.mat, stringsAsFactors=FALSE)
    colnames(met.mat)<-c("name", "frag.MZ", "int", "preMZ", "mode", "adduct", "collisionEnergy")
    met.mat[,"frag.MZ"]<-as.numeric(met.mat[,"frag.MZ"])
    met.mat[,"int"]<-as.numeric(met.mat[,"int"])
    met.mat[,"preMZ"]<-as.numeric(met.mat[,"preMZ"])
    met.mat[,"collisionEnergy"]<-as.numeric(met.mat[,"collisionEnergy"])
    met.mat<-met.mat[order(met.mat[,"preMZ"]), ]
    

    if(MS1==TRUE){
	MS1.df<-suppressWarnings(read.metlinMS(MSxml) ) ##Need to check why warnings?
	MS2.name<-unique(met.mat[,"name"])
	logi<-as.logical(match(MS1.df[,"name"], MS2.name, nomatch=0))
	MS1.df<-data.frame(name=MS1.df[!logi,"name"], frag.MZ=0 , int=0, preMZ=MS1.df[!logi,"MZ"], mode="*", adduct="M", collisionEnergy=0) 
## add fake values to the MS1 data
	met.mat<-rbind(MS1.df, met.mat)
    }

    return(met.mat)
}

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

