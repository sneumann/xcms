## On the long run it would be nice to have all generics in here.
if(!isGeneric("getXcmsRaw"))
    setGeneric("getXcmsRaw", function(object, ...) standardGeneric("getXcmsRaw"))

if(!isGeneric("mslevel"))
    setGeneric("mslevel", function(object, ...) standardGeneric("mslevel"))
if(!isGeneric("mslevel<-"))
    setGeneric("mslevel<-", function(object, value) standardGeneric("mslevel<-"))

if(!isGeneric("profinfo"))
    setGeneric("profinfo", function(object) standardGeneric("profinfo"))

if(!isGeneric("profinfo<-"))
    setGeneric("profinfo<-", function(object, value) standardGeneric("profinfo<-"))

if(!isGeneric("scanrange"))
    setGeneric("scanrange", function(object, ...) standardGeneric("scanrange"))
if(!isGeneric("scanrange<-"))
    setGeneric("scanrange<-", function(object, value) standardGeneric("scanrange<-"))

