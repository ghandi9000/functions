################################################################################
##
##                        Convert mnm to matrix form
##
################################################################################
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
        n <- nrow(neighbors[[i]])
        if (n > 0) {
            for (col in nCols) {
                nMats[[col]][i, 1:n] <- neighbors[[i]][[col]]
            }
        }
    }

    ## Add additional information about targets and plots
    targets <- do.call(rbind, lapply(nLst, function(x) x[["targets"]]))
    nMats[["targets"]] <- targets
    nMats[["plot"]] <- targets$pplot
    nMats[["time"]] <- targets$time

    ## Save some attributes from neighborhood
    attr(nMats, "radius") <- attr(nLst, "radius")
    attr(nMats, "neighbor_par") <- attr(nLst, "neighbor_par")
    attr(nMats, "target_par") <- attr(nLst, "target_par")
    attr(nMats, "data_par") <- attr(nLst, "data_par")

    return ( nMats )
}


## ## Timing
## library(rbenchmark)
## benchmark(
##     mnm_to_matrix(nLst),
##     columns = c("test", "replications", "elapsed", "relative"),
##     replications = 1,
##     order = "relative"
##     )
