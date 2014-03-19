setClass("xcmsEIC", representation(eic = "list", mzrange = "matrix",
                                   rtrange = "matrix", rt = "character",
                                   groupnames = "character"),
         prototype(eic = list(), mzrange = matrix(nrow = 0, ncol = 0),
                   rtrange = matrix(nrow = 0, ncol = 0),
                   rt = character(0), groupnames = character(0)))
