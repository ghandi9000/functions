################################################################################
##
##                               Neighborhoods
##
## * Stable version of mnm
##
################################################################################
## Returns a list of neighborhood data for each target in every combination of
## plot/time period.  Targets are determined by parameters (i.e. spec == "ABBA").
## Neighbors of a given target are determined by another set of parameters
## (i.e. neighbor[["height"]] >= target[["height"]]).  Potential neighbors must
## be within a specified neighborhood radius of targets.  The data to report in
## the neighborhood about each neighbor must be supplied to the function.
##
## Parameters
## tPars: target parameters (syntax same as subset())
## nPars: comparisons between targets and potential neighbors, syntax as follows,
##  !is.na(neighbor[["value"]]) & neighbor[["value"]] >= target[["value"]]
## dPars: logical expression applied to data before neighborhood computations,
##  (syntax same as subset())
## nCols: names of columns to keep for neighbors
## dat: data
## nRad: neighborhood radius, the distance around targets to look for neighbors
## pLims: upper and lower bounds for x,y coordinates for targets/neighbors
## cushion: logical value, if TRUE no targets are used within 1 nRad of
##  upper/lower pLim
## parallel: run the computation in parallel (usually slower due to overhead)
mnm <- function(tPars, nPars, dPars, nCols, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE, parallel=FALSE) {
    ## require(plyr)  # causes problems in code where dplyr is loaded
    if (parallel) {
        require(parallel)
        require(doSNOW)
        numCores <- detectCores()
        cl <- makeCluster(numCores, type = "SOCK")
        registerDoSNOW(cl)
    }

    ## Trim data and create targets
    dd <- dat[eval(dPars, dat),]
    dd <- dd[dd[["bqudx"]] <= pLims["xupper"] &
             dd[["bqudx"]] >= pLims["xlower"] &
             dd[["bqudy"]] <= pLims["yupper"] &
             dd[["bqudy"]] >= pLims["ylower"],]
    dd$targ <- eval(tPars, dd)
    if (cushion)
        dd[(dd$bqudx > pLims["xupper"] + 1 - nRad) |
           (dd$bqudx < pLims["xlower"] - 1 + nRad) |
           (dd$bqudy > pLims["yupper"] + 1 - nRad) |
           (dd$bqudy < pLims["ylower"] - 1 + nRad), "targ"] <- FALSE

    ## Pick out columns needed for neighbor compare
    nParsStr <- deparse(nPars)
    matches <- gregexpr('\\[".*?"\\]', nParsStr)
    colStrs <- gsub('"|\\[|\\]', '', unlist(regmatches(nParsStr, matches)))
    compCols <- names(dd)[match(colStrs, names(dd))]
    compCols <- c(compCols, nCols, "bqudx", "bqudy", "id", "time", "pplot", "targ") # always necessary
    compCols <- unique(compCols[!is.na(compCols)])

    ## Drop extra columns
    nn <- dd[, match(compCols, names(dd))]

    ## Main work done on pplot/time subsets
    neighborhood <- dlply(nn, .(pplot, time), .parallel = parallel, function(neb) {
        targs <- which(neb$targ)
        neighbors <- lapply(targs, FUN = function(tt) {
            targ <- as.list(neb[tt, ])
            neb[eval(nPars, list(target = targ, neighbor = neb)) &
                targ[["id"]] != neb[["id"]] &
                targ[["bqudx"]] + nRad > neb[["bqudx"]] &
                targ[["bqudx"]] - nRad < neb[["bqudx"]] &
                targ[["bqudy"]] + nRad > neb[["bqudy"]] &
                targ[["bqudy"]] - nRad < neb[["bqudy"]] &
                targ[["time"]] == neb[["time"]], names(neb) %in% nCols]
        })
        list(neighbors = neighbors, targets = neb[targs,])
    })

    ## Save parameter information
    attr(neighborhood, "radius") <- nRad
    attr(neighborhood, "neighbor_par") <- nPars
    attr(neighborhood, "target_par") <- tPars
    attr(neighborhood, "data_par") <- dPars

    if (parallel)
        stopCluster(cl)

    return( neighborhood )
}

