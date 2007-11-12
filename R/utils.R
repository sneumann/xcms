
# some utilities
capitalize <- function(str) {
  substring(str, 1, 1) <- toupper(substring(str, 1, 1))
  str
}
uncapitalize <- function(str) {
  substring(str, 1, 1) <- tolower(substring(str, 1, 1))
  str
}
