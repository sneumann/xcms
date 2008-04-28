# functions for generating documentation

setGeneric("prompt",
           function(object, filename = NULL, name = NULL, ...)
           standardGeneric("prompt"))

setMethod("prompt", "missing",
          function(object, filename, name, ...)
          prompt(get(name, parent.frame()), filename, name, ...))

promptStageMethods <- function(object, doc)
  {
    # add methods for dispName, inType and outType
    stageValues <- c(dispName = dispName(object),
                     inType = inType(object), outType = outType(object))
    stageMethods <- paste("Returns \"", stageValues, "\".", sep = "")
    fillMethod <- function(i)
      {
        raw <- doc$`section{Methods}`[i+2]
        filled <- sub("...", stageMethods[i], raw, fixed = TRUE)
        link <- paste("\\link{", names(stageValues)[i], "}", sep = "")
        sub(names(stageValues)[i], link, filled, fixed = TRUE)
      }
    filledMethods <- sapply(seq_along(stageMethods), fillMethod)
    doc$`section{Methods}`[seq_along(filledMethods)+2] <- filledMethods
    doc
  }

setMethod("prompt", "xcmsStage",
          function(object, filename, name, protocols = TRUE, ...)
          {
            protos <- sapply(protocolClasses(object), new)
            methods <- paste(role(object), sapply(protos, method), sep = ".")
            
            ############# STAGE CLASS
            
            cl <- class(object)
            cp <- promptClass(cl, NA, ...)
            
            # set title to display name
            cp$title <- rdmacro("title", paste(dispName(object)))

            # set description to something generic
            methodLink <- rdlink(role(object))
            desc <- paste("A role that a protocol might play in a pipeline. ",
                          "See", methodLink, "for details.")
            
            cp$description <- rdmacro("description", desc)

            # objects-from-the-class can be created using xcmsStage()
            ofc <- paste("Create instances with \\code{\\link{xcmsStage}(\"",
                         role(object), "\")}.", sep = "")
            ofc <- rdmacro("section{Objects from the Class}", ofc)
            cp$`section{Objects from the Class}` <- ofc
            
            # add methods for dispName, inType and outType
            cp <- promptStageMethods(object, cp)
            
            # author - the xcms developers
            cp$author <- "The xcms developers"

            # seealso - the main method for performing protocols
            seealso <- paste(methodLink, "to perform a protocol of this stage.")
            cp$seealso <- rdmacro("seealso", seealso)
            
            # add 'internal' keyword
            internal <- rdmacro("keyword", "internal")
            cp$keywords <- paste(cp$keywords, internal, sep = "\n")

            # stuff we don't want to include
            cp$references <- character(0)
            cp$note <- character(0)
            cp$examples <- character(0)
            
            ############ MAIN METHOD
            
            met <- role(object)
            mp <- promptMethods(met, NA, ...)
            mp$aliases <- c(mp$aliases, rdmacro("alias", met))
            mp$title <- rdmacro("title", met)
            mp$keywords <- mp$keywords[1] # drop dummy keyword
            mp$author <- cp$author
            seealso <- paste(paste(rdlink(methods), collapse=", "),
                             "for details on each (built-in) method.")
            mp$seealso <- rdmacro("seealso", seealso)

            ############ PROTOCOL ACCESSORS

            # produce and fill-in getter docs
            getter <- paste(role(object), "Proto", sep = "")
            gp <- promptMethods(getter, NA)
            gp$aliases <- c(gp$aliases, rdmacro("alias", getter))
            title <- paste("Access", dispName(object), "Protocols")
            gp$title <- rdmacro("title", title)
            desc <- paste("Get or replace a", rdlink(role(object)), "protocol.")
            gp$description <- rdmacro("description", desc)
            getCode <- paste(getter, "(object)", sep = "")
            getBlurb <- paste(rdmacro("\\code", getCode), ": Get the first",
                              rdmacro("\\code", role(object)), "protocol.")
            getterMethods <- gp$`section{Methods}`
            filledMethods <- sub("~~.*? }", paste(getBlurb, "}"), getterMethods)
            gp$keywords[2] <- rdmacro("keyword", "internal")
            
            # merge in setter docs
            setter <- paste(getter, "<-", sep = "")
            sp <- promptMethods(setter, NA)
            gp$aliases <- c(gp$aliases, sp$aliases, rdmacro("alias", setter))
            setCode <- paste(getter, "(object) <- value", sep = "")
            setBlurb <- paste(rdmacro("\\code", setCode), ": Replace the first",
                              rdmacro("\\code", role(object)), "protocol.")
            setterMethods <- tail(head(sp$`section{Methods}`, -1), -1)
            findSig <- function(txt) sub(".*object = \"([^\"]*).*", "\\1", txt)
            setterSig <- findSig(setterMethods)
            getterSig <- findSig(getterMethods)
            mm <- match(setterSig, getterSig)
            settersAdded <- sub(" }", paste("\n", setBlurb, "}"),
                                filledMethods[mm])
            filledMethods[mm] <- settersAdded

            gp$`section{Methods}` <- filledMethods
            
            ############ PROTOCOLS

            if (is.null(filename))
              filename <- getwd()
            
            pp <- NULL
            if (protocols)
              pp <- sapply(protos, prompt, filename = filename, ...)

            ############ OUTPUT
            
            if (!is.na(filename)) {
              rdwrite(cp, filename, cl, "class")
              rdwrite(mp, filename, met, "methods")
              rdwrite(gp, filename, getter, "methods")
              invisible(filename)
            } else list(class = cp, method = mp, accessor = gp, protocols = pp) 
          })

setMethod("prompt", "xcmsProtocol",
          function(object, filename, name, ...)
          {
            print(filename)
            cl <- class(object)
            stage <- stage(object)
            role <- role(stage)
            method <- method(object)

            ############## PROTOCOL CLASS
            
            cp <- promptClass(cl, NA, ...)
            
            # set title to dicsplay name
            cp$title <- rdmacro("title", dispName(object))

            # set description to something generic
            met <- paste(role, method, sep = ".")
            methodLink <- rdlink(met)
            desc <- paste("An object that performs an analysis strategy.",
                          "See", methodLink, "for details.")
            cp$description <- rdmacro("description", desc)

            # objects-from-the-class can be created using xcmsStage()
            ofc <- paste("Create instances with \\code{\\link{xcmsProtocol}(\"",
                         role, "\", \"", method, "\")}.", sep = "")
            ofc <- rdmacro("section{Objects from the Class}", ofc)
            cp$`section{Objects from the Class}` <- ofc

            # note that slots are parameters
            slotMessage <- paste("Slots match parameters in ", methodLink, ".",
                                 sep = "")
            cp$`section{Slots}`[1] <- paste(cp$`section{Slots}`[1], slotMessage)
            cp$`section{Slots}` <- sub(" ~~ ", "", cp$`section{Slots}`)
            
            # add the dispName, inType, outType methods
            cp <- promptStageMethods(object, cp)
            
            # author - the xcms developers
            cp$author <- rdmacro("author", "The xcms developers")

            # seealso - the main method for performing protocols
            seealso <- paste(methodLink, "to perform this protocol.")
            cp$seealso <- rdmacro("seealso", seealso)
            
            # add 'internal' keyword
            internal <- rdmacro("keyword", "internal")
            cp$keywords <- c(cp$keywords, internal)

            # stuff we don't want to include
            cp$references <- character(0)
            cp$note <- character(0)
            cp$examples <- character(0)

            ############ MAIN METHOD

            mp <- promptMethods(met, NA, ...)
            mp$aliases <- c(mp$aliases, rdmacro("alias", met))
            mp$title <- rdmacro("title", dispName(object))
            mp$keywords <- mp$keywords[1] # drop dummy keyword
            mp$author <- cp$author
            seealso <- paste(rdlink(role), "which delegates to this function.")
            mp$seealso <- rdmacro("seealso", seealso)

            ############ OUTPUT

            if (is.null(filename))
              filename <- getwd()
            
            if (!is.na(filename)) {
              rdwrite(cp, filename, cl, "class")
              rdwrite(mp, filename, met, "methods")
              invisible(filename)
            } else list(class = cp, method = mp)
          })

# convenience wrappers

promptStage <- function(role, filename = NULL, protocols = TRUE, ...)
{
  prompt(xcmsStage(role), filename, protocols = protocols, ...)
}

promptProtocol <- function(role, method, filename = NULL, ...)
{
  prompt(xcmsProtocol(role, method), filename, ...)
}

# utils

rdmacro <- function(name, args)
  paste(paste("\\", name, "{", args, "}", sep = ""), collapse = "\n")
rdlink <- function(to) rdmacro("code", rdmacro("link", to))
rdwrite <- function(prompt, dir, name, suffix)
{
  filebase <- paste(paste(name, suffix, sep = "-"), "Rd", sep = ".")
  if (!is.null(dir))
    filepath <- file.path(dir, filebase)
  else filepath <- filebase
  writeLines(unlist(prompt), filepath)
}
