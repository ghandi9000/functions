################################################################################
##
##                        Convert mnm to matrix form
##
################################################################################
source("~/work/functions/neighborhoods/rewrite/run-mnm.R")

## Converts neighborhoods in list format to neighborhoods in matrix format
## nLst: neighborhoods in list format (output from mnm)
mnm_to_matrix <- function(nLst) {
    ## Number of columns in matrices is the maximum number of neighbors across
    ## all plots/times, number of rows is the sum of all targets
    dims <- sapply(nLst, FUN = function(x) {
        c(neighbors = max(unlist(lapply(x[["neighbors"]], nrow))),
          targets = nrow(x[["targets"]]))
    })
    cols <- max(dims["neighbors",])
    rows <- sum(dims["targets",])
    nCols <- names(nLst[[1]][["neighbors"]][[1]])
    plotTime <- attr(nLst, "split_labels")

    ## Initialize matrices
    nMats <- lapply(nCols, FUN = function(x) {
        matrix(NA, nrow = rows, ncol = cols)
    })
    names(nMats) <- nCols

    ## Fill matrices
    neighbors <- unlist(lapply(nLst, function(x) x[["neighbors"]]), recursive = F)

    for (i in seq_along(neighbors)) {
        for (col in nCols) {
            nMats[[col]][i, ]
        }
    }

}


library(rbenchmark)
benchmark(
    mnm_to_matrix(nLst),
    columns = c("test", "replications", "elapsed", "relative"),
    replications = 1,
    order = "relative"
    )


## ## Question on SO
## ## Example data: list of data.frames
## set.seed(0)
## ns <- c(10,5,3)                                 # number of rows in each data.frame
## ids <- runif(sum(ns), 0, 100)                   # some ids
## cll <- quote(data.frame(var     = rnorm(i),     # call to create data.frames
##                         id      = ids[1:i],
##                         y       = 1:i,
##                         x       = 1:i,
##                         pos     = letters[1:i],
##                         stringsAsFactors = F))
## lstForm <- lapply(ns, function(i)               # create list of data.frames
##                   eval(cll, list(i = i)))

## ## Desired output: list of matrices
## ncols <- max(unlist(lapply(lstForm, nrow)))     # number cols of each matrix
## nrows <- length(lstForm)                        # number of rows for each matrix
## colNames <- names(lstForm[[1]])                 # matrix names
## matForm <- lapply(colNames, FUN = function(x)   # initialize matrices
##                   matrix(NA, nrow = nrows, ncol = ncols))
## names(matForm) <- colNames

## ## Filling the matrices using for loops
## for (i in seq_along(lstForm)) {
##     n <- nrow(lstForm[[i]])                     # number of columns to fill
##     for (column in colNames) {
##         matForm[[column]][i, 1:n] <- lstForm[[i]][[column]]
##     }
## }

