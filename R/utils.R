# some utilities

capitalize <- function(str) {
    substring(str, 1, 1) <- toupper(substring(str, 1, 1))
    str
}
uncapitalize <- function(str) {
    ## Dont't capitalize ALL CAPS, e.g. abbreviations
    if (str != toupper(str)) {
        substring(str, 1, 1) <- tolower(substring(str, 1, 1))
    }
    str
}
