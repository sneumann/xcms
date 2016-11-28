## Methods for ProcessHistory and XProcessHistory.

setMethod("initialize", "ProcessHistory", function(.Object, ...) {
    classVersion(.Object)["ProcessHistory"] <- "0.0.2"
    callNextMethod(.Object, ...)
})
setMethod("initialize", "XProcessHistory", function(.Object, ...) {
    classVersion(.Object)["XProcessHistory"] <- "0.0.1"
    callNextMethod(.Object, ...)
})

