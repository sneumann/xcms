profBin <- function(x, y, num, xstart = min(x), xend = max(x), 
                    param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfBin",
       x,
       y,
       as.integer(length(x)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = double(num),
       DUP = FALSE, PACKAGE = "xcms")$out
}

profBinM <- function(x, y, zidx, num, xstart = min(x), xend = max(x), 
                     NAOK = FALSE, param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfBinM",
       x,
       y,
       as.integer(length(x)),
       as.integer(zidx),
       as.integer(length(zidx)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = doubleMatrix(num, length(zidx)),
       NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")$out
}

profBinLin <- function(x, y, num, xstart = min(x), xend = max(x), 
                       param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfBinLin",
       x,
       y,
       as.integer(length(x)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = double(num),
       DUP = FALSE, PACKAGE = "xcms")$out
}

profBinLinM <- function(x, y, zidx, num, xstart = min(x), xend = max(x), 
                        NAOK = FALSE, param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfBinLinM",
       x,
       y,
       as.integer(length(x)),
       as.integer(zidx),
       as.integer(length(zidx)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = doubleMatrix(num, length(zidx)),
       NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")$out
}

profBinLinBase <- function(x, y, num, xstart = min(x), xend = max(x), 
                            param = list()) {

    if (is.null(param$baselevel))
       baselevel <- min(y)/2
    else
       baselevel <- param$baselevel
    if (is.null(param$basespace))
       basespace <- 0.075
    else
       basespace <- param$basespace
    
    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfBinLinBaseM",
       x,
       y,
       as.integer(length(x)),
       as.double(baselevel),
       as.double(basespace),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = doubleMatrix(num, length(zidx)),
       NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")$out
}

profBinLinBaseM <- function(x, y, zidx, num, xstart = min(x), xend = max(x), 
                            NAOK = FALSE, param = list()) {

    if (is.null(param$baselevel))
       baselevel <- min(y)/2
    else
       baselevel <- param$baselevel
    if (is.null(param$basespace))
       basespace <- 0.075
    else
       basespace <- param$basespace
    
    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfBinLinBaseM",
       x,
       y,
       as.integer(length(x)),
       as.integer(zidx),
       as.integer(length(zidx)),
       as.double(baselevel),
       as.double(basespace),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = doubleMatrix(num, length(zidx)),
       NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")$out
}

profIntLin <- function(x, y, num, xstart = min(x), xend = max(x), 
                       param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfIntLin",
       x,
       y,
       as.integer(length(x)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = double(num),
       DUP = FALSE, PACKAGE = "xcms")$out
}

profIntLinM <- function(x, y, zidx, num, xstart = min(x), xend = max(x), 
                        NAOK = FALSE, param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfIntLinM",
       x,
       y,
       as.integer(length(x)),
       as.integer(zidx),
       as.integer(length(zidx)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = doubleMatrix(num, length(zidx)),
       NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")$out
}

medianFilter <- function(x, mrad, nrad) {

    if (!is.double(x)) mat <- as.double(x)
    .C("MedianFilter",
       x,
       as.integer(dim(x)[1]),
       as.integer(dim(x)[2]),
       as.integer(mrad),
       as.integer(nrad),
       out = doubleMatrix(dim(x)[1], dim(x)[2]),
       DUP = FALSE, PACKAGE = "xcms")$out
}

descendZero <- function(y, istart = max(y)) {

    if (!is.double(y)) y <- as.double(y)
    unlist(.C("DescendZero",
              y,
              length(y),
              as.integer(istart-1),
              ilower = integer(1),
              iupper = integer(1),
              DUP = FALSE, PACKAGE = "xcms")[4:5]) + 1
}

descendValue <- function(y, value, istart = which.max(y)) {

    if (!is.double(y)) y <- as.double(y)
    unlist(.C("DescendValue",
              y,
              length(y),
              as.integer(istart-1),
              as.double(value),
              ilower = integer(1),
              iupper = integer(1),
              DUP = FALSE, PACKAGE = "xcms")[5:6]) + 1
}

descendMin <- function(y, istart = max(y)) {

    if (!is.double(y)) y <- as.double(y)
    unlist(.C("DescendMin",
              y,
              length(y),
              as.integer(istart-1),
              ilower = integer(1),
              iupper = integer(1),
              DUP = FALSE, PACKAGE = "xcms")[4:5]) + 1
}

findEqualGreaterM <- function(x, values) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(values)) values <- as.double(values)
    .C("FindEqualGreaterM",
       x,
       length(x),
       values,
       length(values),
       index = integer(length(values)),
       DUP = FALSE, PACKAGE = "xcms")$index + 1
}

findEqualGreater <- function(x, value) {

    if (!is.double(x)) x <- as.double(x)
    .C("FindEqualGreater",
       x,
       length(x),
       as.double(value),
       index = integer(1),
       DUP = FALSE, PACKAGE = "xcms")$index + 1
}

findEqualLess <- function(x, value) {

    if (!is.double(x)) x <- as.double(x)
    .C("FindEqualLess",
       x,
       length(x),
       as.double(value),
       index = integer(1),
       DUP = FALSE, PACKAGE = "xcms")$index + 1
}

findRange <- function(x, values, NAOK = FALSE) {

    if (!is.double(x)) x <- as.double(x)
    start <- .C("FindEqualGreater",
                x,
                length(x),
                as.double(values[1]),
                integer(1),
                NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")[[4]]
    end <- .C("FindEqualLess",
              x,
              length(x),
              as.double(values[2]),
              integer(1),
              NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")[[4]]
    c(start, end) + 1
}

colMax <- function (x, na.rm = FALSE, dims = 1) {

    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (!is.array(x) || length(dn <- dim(x)) < 2) 
        stop("`x' must be an array of at least two dimensions")
    if (dims < 1 || dims > length(dn) - 1) 
        stop("invalid `dims'")
    n <- prod(dn[1:dims])
    dn <- dn[-(1:dims)]
    if (!is.double(x)) x <- as.double(x)
    z <- .C("ColMax",
            x, 
            as.integer(n), 
            as.integer(prod(dn)),
            double(prod(dn)),
            DUP = FALSE, PACKAGE = "xcms")[[4]]
    if (length(dn) > 1) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[-(1:dims)]
    }
    else names(z) <- dimnames(x)[[dims + 1]]
    z
}

rowMax <- function (x, na.rm = FALSE, dims = 1) {

    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (!is.array(x) || length(dn <- dim(x)) < 2) 
        stop("`x' must be an array of at least two dimensions")
    if (dims < 1 || dims > length(dn) - 1) 
        stop("invalid `dims'")
    p <- prod(dn[-(1:dims)])
    dn <- dn[1:dims]
    if (!is.double(x)) x <- as.double(x)
    z <- .C("RowMax",
            x, 
            as.integer(prod(dn)), 
            as.integer(p),
            double(prod(dn)),
            DUP = FALSE, PACKAGE = "xcms")[[4]]
    if (length(dn) > 1) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[1:dims]
    }
    else names(z) <- dimnames(x)[[1]]
    z
}

doubleMatrix <- function(nrow = 0, ncol = 0) {

    .Call("DoubleMatrix", 
          as.integer(nrow), 
          as.integer(ncol),
          PACKAGE = "xcms")
}

logicalMatrix <- function(nrow = 0, ncol = 0) {

    .Call("LogicalMatrix", 
          as.integer(nrow), 
          as.integer(ncol),
          PACKAGE = "xcms")
}
