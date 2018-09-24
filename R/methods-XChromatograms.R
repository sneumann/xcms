setMethod("show", "XChromatograms", function(object) {
    callNextMethod()
    ## TODO number of peaks per ...
    
    nr <- nrow(object)
    nc <- ncol(object)
    cat(class(object), " with ",
        nr, ifelse(nr == 1, " row and ", " rows and "),
        nc, ifelse(nc == 1, " column\n", " columns\n"),
        sep = "")
    sumFun <- function(z) {
        paste0("length: ", length(z[[1]]))
    }
    if (nr > 0 && nc > 0) {
        if (nr <= 4) {
            out <- apply(object, MARGIN = c(1, 2), sumFun)
            rownames(out) <- paste0("[", 1:nrow(out), ",]")
        }
        else {
            out <- rbind(
                apply(object[c(1, 2), , drop = FALSE], MARGIN = c(1, 2), sumFun),
                rep(" ... ", ncol(object)),
                apply(object[nrow(object) - c(1, 0), , drop = FALSE],
                      MARGIN = c(1, 2), sumFun)
            )
            rownames(out) <- c("[1,]", "[2,]", "...",
                               paste0("[", c(nrow(object) - c(1, 0)), ",]"))
        }
        rn <- rownames(out)
        out <- rbind(rep("<Chromatogram>", ncol(out)), out)
        rownames(out) <- c("", rn)        
        print(out, quote = FALSE, right = TRUE)
    }
    cat("phenoData with", length(varLabels(object@phenoData)), "variables\n")
    cat("featureData with", length(fvarLabels(object)), "variables\n")
})

## FIX: if there is only a single element -> make a Chromatograms object again!
setMethod("[", "XChromatograms",
          function(x, i, j, drop = FALSE) {
              if (missing(i) & missing(j))
                  return(x)
              if (missing(i))
                  i <- seq_len(nrow(x))
              if (missing(j))
                  j <- seq_len(ncol(x))
              if (is.logical(i))
                  i <- which(i)
              if (is.logical(j))
                  j <- which(j)
              ## Return a single element as a Chromatogram
              if (length(i) == 1 & length(j) == 1)
                  return(x@.Data[i, j, drop = TRUE][[1]])
              pd <- x@phenoData
              fd <- x@featureData
              ## Multiple elements, return type depends on drop.
              x <- x@.Data[i = i, j = j, drop = drop]
              if (!drop) {
                  x <- as(x, "Chromatograms")
                  pd <- pd[j, ]
                  ## Drop levels
                  pData(pd) <- droplevels(pData(pd))
                  x@phenoData <- pd
                  fd <- fd[i, ]
                  pData(fd) <- droplevels(pData(fd))
                  x@featureData <- fd
              }
              if (validObject(x))
                  x
          })
