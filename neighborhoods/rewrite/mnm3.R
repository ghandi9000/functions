## Uses is_neb_list()
nPars2 <- quote(
    !is.na(neb[["ba"]]) & neb[["ba"]] >= targ[["ba"]]
    )

nCols <- c("ba","id","bqudx","bqudy","z","x","y","spec")

mnm3 <- function(tPars, nPars, dPars, nCols, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE) {
    require(plyr)

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
    neighborhood <- dlply(nn, .(pplot, time), function(neb) {
        targs <- which(neb$targ)
        neighbors <- lapply(targs, FUN = function(target) {
            targ <- as.list(neb[target, ])
            neb[eval(nPars, list(targ = targ, neb = neb)) &
                targ[["id"]] != neb[["id"]] &
                targ[["bqudx"]] + nRad > neb[["bqudx"]] &
                targ[["bqudx"]] - nRad < neb[["bqudx"]] &
                targ[["bqudy"]] + nRad > neb[["bqudy"]] &
                targ[["bqudy"]] - nRad < neb[["bqudy"]] &
                targ[["time"]] == neb[["time"]], names(neb) %in% nCols]
        })
        list(neighbors = neighbors, targets = neb[targs,],
             plot = unique(neb[["pplot"]]), time = unique(neb[["time"]]))
    })

    return( neighborhood )
}



################################################################################
##
##                               Benchmarking
##
################################################################################
library(rbenchmark)


tst1 <- function(...) {
    sapply(targs, FUN = function(targ) {
        sapply(1:nrow(pp), FUN = function(neb) {
            is_neb_list(as.list(pp[targ, ]), as.list(pp[neb, ]), nRad, nPars2)
        })
    })
}



tst2 <- function(...) {
    lapply(targs, FUN = function(targ) {
        lapply(1:nrow(pp), FUN = function(neb) {
            is_neb_list(as.list(pp[targ, ]), as.list(pp[neb, ]), nRad, nPars2)
        })
    })
}


tst3 <- function(...) {
    sapply(targs, FUN = function(targ) {
        sapply(1:nrow(pp), FUN = function(neb) {
            is_neb(pp[targ, ], pp[neb, ], nRad, nPars)
        })
    })
}

benchmark(
    tst1(), tst2(), tst3(),
    columns = c("test", "replications", "elapsed", "relative"),
    order = "relative", replications = 1
    )


thing <- deparse(nPars2)
## "!is.na(neb[[\"ba\"]]) & neb[[\"ba\"]] >= targ[[\"ba\"]]"
blah <- strsplit(thing, split = "&|\\||>=|<=|<|>|!=")


tstStr <- "!is.na(neb[[\"ba\"]]) & neb[[\"ba\"]] >= targ[[\"ba\"]] & neb[[\"spec\"]] == \"ABBA\""
nPars <- quote(!is.na(neb[["ba"]]) &
               neb[["ba"]] >= targ[["ba"]] &
               neb[["spec"]] == "ABBA")

## Pick out columns that are necessary for neighbor comparisons
m <- gregexpr('".*?"', deparse(nPars))
res <- gsub('"', '', regmatches(tstStr, m)[[1]])
cols <- names(dd)[match(res, names(dd))]
cols <- unique(cols[!is.na(cols)])

res <- mnm3(tPars=tPars, nPars=nPars, dPars=dPars, nCols = nCols, dat=dat, nRad =3)
