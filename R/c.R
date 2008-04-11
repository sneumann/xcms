profMaxIdx <- function(x, y, num, xstart = min(x), xend = max(x),
                       param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfMaxIdx",
       x,
       y,
       as.integer(length(x)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = integer(num),
       DUP = FALSE, PACKAGE = "xcms")$out
}

profMaxIdxM <- function(x, y, zidx, num, xstart = min(x), xend = max(x),
                        NAOK = FALSE, param = list()) {

    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    .C("ProfMaxIdxM",
       x,
       y,
       as.integer(length(x)),
       as.integer(zidx),
       as.integer(length(zidx)),
       as.double(xstart),
       as.double(xend),
       as.integer(num),
       out = integerMatrix(num, length(zidx)),
       NAOK = NAOK, DUP = FALSE, PACKAGE = "xcms")$out
}

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
       out = double(num),
       DUP = FALSE, PACKAGE = "xcms")$out
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

    dimx <- dim(x)
    if (!is.double(x)) x <- as.double(x)
    .C("MedianFilter",
       x,
       as.integer(dimx[1]),
       as.integer(dimx[2]),
       as.integer(mrad),
       as.integer(nrad),
       out = doubleMatrix(dimx[1], dimx[2]),
       DUP = FALSE, PACKAGE = "xcms")$out
}

descendZero <- function(y, istart = which.max(y)) {

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

descendMin <- function(y, istart = which.max(y)) {

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

which.colMax <- function (x, na.rm = FALSE, dims = 1) {

    if (is.data.frame(x))
        x <- as.matrix(x)
    if (!is.array(x) || length(dn <- dim(x)) < 2)
        stop("`x' must be an array of at least two dimensions")
    if (dims < 1 || dims > length(dn) - 1)
        stop("invalid `dims'")
    n <- prod(dn[1:dims])
    dn <- dn[-(1:dims)]
    if (!is.double(x)) x <- as.double(x)
    z <- .C("WhichColMax",
            x,
            as.integer(n),
            as.integer(prod(dn)),
            integer(prod(dn)),
            DUP = FALSE, PACKAGE = "xcms")[[4]]
    if (length(dn) > 1) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[-(1:dims)]
    }
    else names(z) <- dimnames(x)[[dims + 1]]
    z
}

which.rowMax <- function (x, na.rm = FALSE, dims = 1) {

    if (is.data.frame(x))
        x <- as.matrix(x)
    if (!is.array(x) || length(dn <- dim(x)) < 2)
        stop("`x' must be an array of at least two dimensions")
    if (dims < 1 || dims > length(dn) - 1)
        stop("invalid `dims'")
    p <- prod(dn[-(1:dims)])
    dn <- dn[1:dims]
    if (!is.double(x)) x <- as.double(x)
    z <- .C("WhichRowMax",
            x,
            as.integer(prod(dn)),
            as.integer(p),
            integer(prod(dn)),
            DUP = FALSE, PACKAGE = "xcms")[[4]]
    if (length(dn) > 1) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[1:dims]
    }
    else names(z) <- dimnames(x)[[1]]
    z
}

rectUnique <- function(m, order = seq(length = nrow(m)), xdiff = 0, ydiff = 0) {

    nr <- nrow(m)
    nc <- ncol(m)
    if (!is.double(m))
        m <- as.double(m)
    .C("RectUnique",
       m,
       as.integer(order-1),
       nr,
       nc,
       as.double(xdiff),
       as.double(ydiff),
       logical(nrow(m)),
       DUP = FALSE, PACKAGE = "xcms")[[7]]
}

doubleMatrix <- function(nrow = 0, ncol = 0) {

    .Call("DoubleMatrix",
          as.integer(nrow),
          as.integer(ncol),
          PACKAGE = "xcms")
}

integerMatrix <- function(nrow = 0, ncol = 0) {

    .Call("IntegerMatrix",
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

continuousPtsAboveThreshold <- function(y, threshold, num, istart = 1) {
    if (!is.double(y)) y <- as.double(y)
    if (.C("continuousPtsAboveThreshold",
              y,
              as.integer(istart-1),
              length(y),
              threshold = as.double(threshold),
              num = as.integer(num),
              n = integer(1),
              DUP = FALSE, PACKAGE = "xcms")$n > 0) TRUE else FALSE
}
