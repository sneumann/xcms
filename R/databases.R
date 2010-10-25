read.metlin<-function(xml, MS1=TRUE, MSxml) {
	xmlOpen<-url(xml)
    load(xmlOpen)##This creates MS2_data
	close(xmlOpen)
	met.mat<-MS2_data[[1]]
	for(i in 1:length(MS2_data[[2]])){
		idx<-which(met.mat[,"adduct"] == MS2_data[[2]][i])
		met.mat[idx,"adduct"]<-names(MS2_data[[2]])[i]
	}
	MS1.df<-xcms:::read.metlinMS("http://metlin.scripps.edu/download/MS.Rda")
	UniMolid<-unique(met.mat[,"molid"])
	met.mat$name<-rep("", nrow(met.mat))
	cat("loading...\n")
	for(i in 1:length(UniMolid)){
		MS2idx<-which(met.mat[,"molid"] == UniMolid[i])
		MS1idx<-which(MS1.df[,"molid"]  == UniMolid[i])
		
		met.mat$name[MS2idx]<-MS1.df[MS1idx,"name"]
		cat(round((i/length(UniMolid)))*100, "%\r")
	}
	cat("\n")

	colnames(met.mat)<-c("molid", "PreMZ", "mode", "collisionEnergy", "adduct", "frag.MZ", "int", "name")
	met.mat<-met.mat[order(met.mat[,"PreMZ"]), ]

    if(MS1==TRUE){
        MS2.name<-unique(met.mat[,"name"])
        logi<-as.logical(match(MS1.df[,"name"], MS2.name, nomatch=0))
        MS1.df<-data.frame(molid=MS1.df[!logi,"molid"],  PreMZ=MS1.df[!logi,"MZ"], mode="*", collisionEnergy=0, 
 							adduct="M",frag.MZ=0 , int=0, name=MS1.df[!logi,"name"]) 
## add fake values to the MS1 data
        met.mat<-rbind(MS1.df, met.mat)
    }
    return(met.mat)
}

read.metlinMS<- function(xml){
    xmlOpen<-url(xml)
	load(xmlOpen) ##this create ans
	close(xmlOpen)
    return(ans)
}

metlinToList<-function(metlinfile){
    met.xml<-read.metlin(metlinfile, MS1=FALSE)
    cat("Converting Data type\n")
    ref<-1
    spectra<-list()
    met.xml<-met.xml[order(met.xml[,"name"]),]
    names<-unique(met.xml[,"name"])
    for(i in 1:length(names)){
        idx<-which(met.xml[,"name"] == names[i])
        CE<-unique(met.xml[idx,"collisionEnergy"])
        tempMetMat<-met.xml[idx,]
        for(j in 1:length(CE)){
            ind<-which(tempMetMat[,"collisionEnergy"] == CE[j])
            ion<-unique(tempMetMat[ind,"mode"])
            foo<-tempMetMat[ind,]
            for(k in 1:length(ion)){
                ix<-which(foo[,"mode"]== ion[k])
                MetSpec<-foo[ix,]
                spectra[[ref]]<-MetSpec
                ref<-ref+1
            }
        }
        percent<-(i/length(names))*100
        cat(paste("\r", round(percent), "% Done", sep=""))
    }
    cat("\n")
    invisible(spectra)
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

KeggSearch <- function(object, DBsearchMS) {
	
    libname <- 'KEGGSOAP'
    KEGG.status <- try(require(libname, character.only = TRUE, quietly=TRUE))
    
    if (class(KEGG.status) == "try-error")
        stop("Couldn't load KEGGSOAP\n")
	
	if(class(object)=="xcmsSet"){
		groupmat<-groups(object)
		neutralmass <- groupmat[,"mzmed"] + ifelse(DBsearchMS < 0, 1, -1)

		KEGGcmpd<-array(0, dim=length(neutralmass))
		for(i in 1:length(neutralmass)){
			KEGGcmpd[i]<-search.compounds.by.mass(neutralmass, range(ppmDev(neutralmass[i], DBsearchMS)))[1]
		}
		return(KEGGcmpd)
	}
}

metlinMS1search<-function(object, DBsearchMS){
	met<-read.metlinMS("http://metlin.scripps.edu/download/MS.Rda")
	if(class(object) =="xcmsSet"){
		groupmat<-groups(object)
		neutralmass <- groupmat[,"mzmed"] + ifelse(DBsearchMS < 0, 1, -1)
	} else if(class(object) =="xcmsPeaks"){
		neutralmass<-object[,"mz"] + ifelse(DBsearchMS < 0, 1, -1)
	}else{
		cat("Method requires an xcmsPeaks or xcmsSet object\n")
		return(0)
	}
	
	ans<-matrix(0,ncol=2,nrow=length(neutralmass))
	colnames(ans)<-c("name", "mass")
	for(i in 1:length(neutralmass)){
		massError<-range(ppmDev(neutralmass[i], DBsearchMS))
		idx<-which(met[,"mass"] > massError[1] & met[,"mass"] < massError[2])
		if(length(idx) >1){
			idxError<-ppm(met[idx,"mass"], neutralmass[i])
			ans[i,"mass"]<-met[idx,][which.min(idxError),]$mass
			ans[i,"name"]<-met[idx,][which.min(idxError),]$name
		} else if(length(idx) == 1){
			ans[i,"mass"]<-met[idx,]$mass
			ans[i,"name"]<-met[idx,]$name
		} else{
			ans[i,]<-c("","")
		}
	}
	return(ans)
}