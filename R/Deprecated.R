## ## Put all deprecated methods/functions in here, so they will be
## ## easier to defunct/remove later.

xcmsParallelSetup <- function(nSlaves) {
    .Deprecated(msg = "Use of 'xcmsParallelSetup' is deprecated! Use 'BPPARAM' arguments instead.")
    runParallel <- 0
    parMode <- ""
    snowclust <- NULL

    if (!is.null(nSlaves)) {
        if (nSlaves > 1) {
            ## If MPI is available ...
            rmpi = "Rmpi"
            opt.warn <- options("warn")$warn
            options("warn" = -1)

            ## Rmpi does not work on the BioC build machines:
            if ( (Sys.info()["sysname"] != "Windows") && require(rmpi,character.only=TRUE,quietly=TRUE)) {
                if (is.loaded('mpi_initialize')) {
                    mpi.spawn.Rslaves(nslaves=nSlaves, needlog=FALSE)
                    ## If there are multiple slaves AND this process is the master,
                    ## run in parallel.
                    if ((mpi.comm.size() > 2)  && (mpi.comm.rank() == 0)) {
                        runParallel <- 1
                        parMode <- "MPI"
                    }
                }
            } else {
                ## try local sockets using snow package
                snow = "snow"
                if (try(require(snow,character.only=TRUE,quietly=TRUE))) {
                    cat("Starting snow cluster with",nSlaves,"local sockets.\n")
                    snowclust <- makeCluster(nSlaves, type = "SOCK")
                    runParallel <- 1
                    parMode <- "SOCK"
                } else{
                    ## check parallel package... can use the mclapply on local CPUs
                    if(requireNamespace("parallel", quietly=TRUE)){
                        cat("Processing on", nSlaves, "cores.\n")
                        runParallel <- 1
                        parMode <- "parallel"
                        ## setting the number of cores
                        options(mc.cores=nSlaves)
                    }
                }
            }
            options("warn" = opt.warn)
        }
    }
    return (list(runParallel=runParallel,
                 parMode=parMode,
                 snowclust=snowclust))
}


##
## Modified papply function from package papply 0.2 (Duane Currie):
##
## Parts of the slave function were wrapped in try() to make it failsafe
## (if e.g. peak picking in a single file fails - papply would wait forever for the MPI slave to finish ...)
##

"xcmsPapply" <-
    function(arg_sets,papply_action,papply_commondata=list(),
             show_errors=TRUE,do_trace=FALSE,also_trace=c()) {
        .Deprecated(msg = "Use of 'xcmsPapply' is deprecated! Use BiocParallel 'bplapply' instead.")
        ## Check to ensure arguments are of the correct type
        if (!is.list(arg_sets)) {
            print("1st argument to papply must be a list")
            return(NULL)
        }

        if (!is.function(papply_action)) {
            print("2nd argument to papply must be a function")
            return(NULL)
        }

        if (!is.list(papply_commondata)) {
            print("3rd argument to papply must be a list")
            return(NULL)
        }

        papply_also_trace <- also_trace

        ## Default to running serially.  Only run in parallel if Rmpi
        ## is installed AND there's multiple slaves AND papply is called
        ## from the master.
        run_parallel <- 0

        ## Load the MPI Environment if not already there.
        if (!is.loaded('mpi_initialize')) {
            libname <- 'Rmpi'
            try(require(libname, character.only = TRUE, quietly=TRUE))
        }

        ## Now, if Rmpi is loaded, make sure we have a bunch of slaves
        if (is.loaded('mpi_initialize')) {
            ## Spawn as many slaves as possible
            if (mpi.comm.size() < 2) {
                mpi.spawn.Rslaves()
            }

            ## If there's multiple slaves, AND this process is the master,
            ## run in parallel.
            if (mpi.comm.size() > 2) {
                if (mpi.comm.rank() == 0) {
                    run_parallel <- 1
                }
            }
        }

        ## Ideally, by here, we can tell if there's a parallel environment
        ## to run in.  If not, this should work:
        ##    attach(commondata)
        ##    results <- lapply(arg_sets,action)
        ##    detach(commondata)
        ##    return(results)

        if (run_parallel != 1) {
            ## There's either no parallel environment, no parallel environment
            ## worth using, or papply's being called in a slave.  Use the
            ## serial version which just calls lapply
            print("Running serial version of papply\n")
            attach(papply_commondata)
            results <- lapply(arg_sets,papply_action)
            detach(papply_commondata)
            return(results)
        }

        ## Create the driver function for the slave processes.
        ## Essentially, the slaves imports the commondata into the current
        ## namespace, and begins doing tasks.  It's task loop involves
        ## signaling the master that it's ready for a task, receives an argument
        ## set from the master master and applies the given function to it,
        ## and returns the result back to the master.  This continues until
        ## all argument sets have been processed.  Care is taken to preserve
        ## the matching order of argument sets an results
        papply_int_slavefunction <- function() {
            ## Import the environment into the current namespace
            attach(papply_commondata)

            if (get("papply_do_trace")) {
                assign("papply_fn_bodies$papply_action",
                       as.list(body(papply_action)),
                       envir=globalenv())
                trace(papply_action,
                      quote({papply_lineno <- papply_lineno+1 ;
                          cat("papply_action: Line ",papply_lineno, ": ") ;
                          print(get("papply_fn_bodies$papply_action")[[papply_lineno]]) }),
                      quote(cat("\n")),
                      1:length(get("papply_fn_bodies$papply_action")),
                      where=environment(),
                      print=FALSE)

                cat("About to start tracing: ",papply_also_trace,"\n")
                for (fn_name in papply_also_trace) {
                    assign(paste("papply_fn_bodies$",fn_name),
                           as.list(body(fn_name)),
                           envir=globalenv())
                    trace(fn_name,
                          substitute({papply_lineno <- papply_lineno+1 ;
                              cat(fn_name,": Line ",papply_lineno, ": ") ;
                              print(get(paste("papply_fn_bodies$",fn_name))[[papply_lineno]]) },list(fn_name=fn_name)),
                          quote(cat("\n")),
                          1:length(get(paste("papply_fn_bodies$",fn_name))),
                          where=environment(),
                          print=FALSE)
                }
            }

            papply_int_junk <- 0
            papply_int_done <- 0
            while (papply_int_done != 1) {
                ## Signal master that this slave is ready for a task
                mpi.send.Robj(papply_int_junk,0,1)

                ## Receive a task, and get its meta-information
                papply_int_task <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
                papply_int_task_info <- mpi.get.sourcetag()
                papply_int_source <- papply_int_task_info[1]
                papply_int_tag <- papply_int_task_info[2]

                if (papply_int_tag == 1) {
                    ## If a request for processing, runs the user-provided
                    ## function on the current argument set, and compiles
                    ## a results message, indicating the return value, and
                    ## proper order in the results
                    papply_int_seqno <- papply_int_task$seqno
                    ##  Modification to make papply_action() failsafe
                    papply_int_results <- try(papply_action(papply_int_task$data))
                    if (class(papply_int_results) == "try-error") {
                        papply_int_results <- geterrmessage()
                        res_tag <- 4
                    } else
                        res_tag <- 2
                    ##
                    papply_int_result_obj <- list(results=papply_int_results,seqno=papply_int_seqno)

                    ## Send the results to the master
                    mpi.send.Robj(papply_int_result_obj,0,tag=res_tag)
                }
                else if (papply_int_tag == 2) {
                    ## Master says it's all done
                    papply_int_done <- 1
                }
            }

            ## Tell master I'm exiting, and detach from namespace

            mpi.send.Robj(papply_int_junk,0,3)
            untrace(papply_action)
            detach(papply_commondata)
        }

        ## Back in the master.
        ## Send the necessary data to all slaves, and tell them to
        ## run the function just defined above.
        mpi.bcast.Robj2slave(papply_commondata)
        mpi.bcast.Robj2slave(papply_action)
        mpi.bcast.Robj2slave(papply_int_slavefunction)
        mpi.bcast.Robj2slave(papply_also_trace)
        if (show_errors) {
            mpi.bcast.cmd(options(error=quote( {
                                                  cat("Error: ",geterrmessage(),"\n") ;
                                                  assign(".mpi.err", TRUE, env = .GlobalEnv)
                                              })))
        }

        mpi.bcast.cmd(papply_fn_bodies <- list())
        mpi.bcast.cmd(papply_lineno <- 0)
        mpi.bcast.cmd(papply_do_trace <- FALSE)
        if (do_trace) {
            mpi.bcast.cmd(papply_do_trace <- TRUE)
        }

        mpi.bcast.cmd(papply_int_slavefunction())

        ## Prepare for communication with the slaves
        junk <- 0
        n_slaves <- mpi.comm.size() - 1
        exited <- 0
        results <- list()
        current_task <- 1
        tag <- 0

        while ((exited < n_slaves)) {
            ## Get a message from a slave, and get the message's meta-data
            message <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
            ## source=mpi.any.source(); tag=mpi.any.tag() ; comm=1; status=0

            message_info <- mpi.get.sourcetag()
            slave_id <- message_info[1]
            tag <- message_info[2]

            if (tag == 1) {
                ## If the slave is ready for a task, either give it one of the
                ## argument sets as a task, or tell it there's no tasks left.
                if (current_task <= length(arg_sets)) {
                    task <- list(data=arg_sets[[current_task]],seqno=current_task)
                    ## cat("Sending task ##",arg_sets[[current_task]],"\n")
                    cat("Sending task ##",current_task,"\n"); flush.console();
                    mpi.send.Robj(task,slave_id,1)
                    current_task <- current_task + 1
                }  else {
                    mpi.send.Robj(junk,slave_id,2)
                }
            }
            else if (tag == 2) {
                ## The slave gave some results.  Compile it into the results
                ## array
                results[[message$seqno]] <- message$results
            }
            else if (tag == 3) {
                ## A slave exited.
                exited <- exited + 1
            }
            else if (tag == 4) {
                ## An error occured in the slave function
                cat(message$results)
            }
        }

        ## Now all slaves are done doing tasks.
        ## Clean up slaves for future calls.
        ## NOTE:
        ##     If I get real smart about multiple runs, lazy deletions
        ##     could be done, and could use a hash computation to tell
        ##     if data has to be re-sent to the slaves, or it can be
        ##     left as is.  Basically, to avoid sending the same data
        ##     multiple times across papply runs.
        mpi.bcast.cmd(papply_commondata <- NULL)
        mpi.bcast.cmd(papply_action <- NULL)
        mpi.bcast.cmd(papply_int_slavefunction <- NULL)
        mpi.bcast.cmd(papply_lineno <- NULL)
        mpi.bcast.cmd(papply_fn_bodies <- NULL)
        mpi.bcast.cmd(papply_do_trace <- NULL)

        return(results)
    }


## clusterApplyLB / dynamicClusterApply
xcmsClusterApply <- function(cl, x, fun, msgfun=NULL, ...) {
    .Deprecated(msg = "Use of 'xcmsClusterApply' is deprecated! Use 'BPPARAM' arguments instead.")
    argfun <- function(i) c(list(x[[i]]), list(...))
    n <- length(x)

    checkCluster(cl)
    p <- length(cl)
    if (n > 0 && p > 0) {
        submit <- function(node, job) sendCall(cl[[node]], fun,
                                               argfun(job), tag = job)
        for (i in 1:min(n, p)) {
            if (!is.null(msgfun))
                do.call(msgfun,args=list(x=x,i=i));

            submit(i, i)
        }
        val <- vector("list", n)
        for (i in 1:n) {
            d <- recvOneResult(cl)
            j <- i + min(n, p)
            if (j <= n) {
                if (!is.null(msgfun))
                    do.call(msgfun,args=list(x=x,i=j));

                submit(d$node, j)
            }
            val[d$tag] <- list(d$value)
        }
        checkForRemoteErrors(val)
    }

}

## setMethod("extractChromatograms",
##           signature(object = "OnDiskMSnExp"),
##           function(object, rt, mz, aggregationFun = "sum", missing = NA_real_) {
##               .Deprecated(msg = paste0("Use of 'extractChromatograms' is ",
##                                        "deprecated, please use 'chromatogram' ",
##                                        "instead."))
##               chromatogram(object, rt = rt, mz = mz,
##                            aggregationFun = aggregationFun, missing = missing)
##           })

## plotChromatogram <- function(x, rt, col = "#00000060",
##                              lty = 1, type = "l", xlab = "retention time",
##                              ylab = "intensity", main = NULL, ...) {
##     .Deprecated(msg = paste0("Use of 'plotChromatogram' is deprecated, please ",
##                              "use 'plot' instead."))
##     if (!is.list(x) & !is(x, "Chromatogram"))
##         stop("'x' should be a Chromatogram object or a list of Chromatogram",
##              " objects.")
##     if (is(x, "Chromatogram"))
##         x <- list(x)
##     isOK <- lapply(x, function(z) {
##         if (is(z, "Chromatogram")) {
##             return(TRUE)
##         } else {
##             if (is.na(z))
##                 return(TRUE)
##         }
##         FALSE
##     })
##     if (any(!unlist(isOK)))
##         stop("if 'x' is a list it should only contain Chromatogram objects")
##     ## Subset the Chromatogram objects if rt provided.
##     if (!missing(rt)) {
##         rt <- range(rt)
##         x <- lapply(x, function(z) {
##             if (is(z, "Chromatogram"))
##                 filterRt(z, rt = rt)
##         })
##     }
##     if (length(col) != length(x)) {
##         col <- rep(col[1], length(x))
##     }
##     ## If main is NULL use the mz range.
##     if (is.null(main)) {
##         mzr <- range(lapply(x, mz), na.rm = TRUE, finite = TRUE)
##         main <- paste0(format(mzr, digits = 7), collapse = " - ")
##     }
##     ## Number of measurements we've got per chromatogram. This can be different
##     ## between samples, from none (if not a single measurement in the rt/mz)
##     ## to the number of data points that were actually measured.
##     lens <- unique(lengths(x))
##     max_len <- max(lens)
##     max_len_vec <- rep_len(NA, max_len)
##     ## Generate the matrix of rt values, columns are samples, rows retention
##     ## time values. Fill each column with NAs up to the maximum number of values
##     ## we've got in a sample/file.
##     rts <- do.call(cbind, lapply(x, function(z) {
##         cur_len <- length(z)
##         if (cur_len == 0)
##             max_len_vec
##         else {
##             ## max_len_vec[,] <- NA  ## don't need that. get's copied.
##             max_len_vec[seq_len(cur_len)] <- rtime(z)
##             max_len_vec
##         }
##     }))
##     ## Same for the intensities.
##     ints <- do.call(cbind, lapply(x, function(z) {
##         cur_len <- length(z)
##         if (length(z) == 0)
##             max_len_vec
##         else {
##             ## max_len_vec[,] <- NA  ## don't need that. get's copied.
##             max_len_vec[seq_len(cur_len)] <- intensity(z)
##             max_len_vec
##         }
##     }))
##     ## Define the x and y limits
##     x_lim <- c(0, 1)
##     y_lim <- c(0, 1)
##     if (all(is.na(rts)))
##         if (!missing(rt))
##             x_lim <- range(rt)
##     else
##         x_lim <- range(rts, na.rm = TRUE, finite = TRUE)
##     if (!all(is.na(ints)))
##         y_lim <- range(ints, na.rm = TRUE, finite = TRUE)
##     ## Identify columns that have only NAs in either intensity or rt - these
##     ## will not be plotted.
##     keepCol <- which(apply(ints, MARGIN = 2, function(z) any(!is.na(z))) |
##                      apply(rts, MARGIN = 2, function(z) any(!is.na(z))))
##     ## Finally plot the data.
##     if (length(keepCol)) {
##         matplot(x = rts[, keepCol, drop = FALSE],
##                 y = ints[, keepCol, drop = FALSE], type = type, lty = lty,
##                 col = col[keepCol], xlab = xlab, ylab = ylab, main = main,
##                 ...)
##     } else
##         plot(x = 3, y = 3, pch = NA, xlab = xlab, ylab = ylab, main = main,
##              xlim = x_lim, ylim = y_lim)
## }
