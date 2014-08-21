## Uses is_neb_list()
nPars2 <- quote(
    !is.na(neb[["ba"]]) & neb[["ba"]] >= targ[["ba"]]
    )


mnm3 <- function(tPars, nPars, dPars, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE) {
    require(plyr)

    ## Trim data and create targets
    dd <- dat[eval(dPars, dat),]
    dd <- dd[dd$bqudx <= pLims["xupper"] &
             dd$bqudx >= pLims["xlower"] &
             dd$bqudy <= pLims["yupper"] &
             dd$bqudy >= pLims["ylower"],]
    dd$targ <- eval(tPars, dd)
    if (cushion)
        dd[(dd$bqudx > pLims["xupper"] + 1 - nRad) |
           (dd$bqudx < pLims["xlower"] - 1 + nRad) |
           (dd$bqudy > pLims["yupper"] + 1 - nRad) |
           (dd$bqudy < pLims["ylower"] - 1 + nRad), "targ"] <- FALSE

    ## Main work done on pplot/time subsets
    neighbors <- ddply(dd, .(pplot, time), function(pp) {
        targs = which(pp$targ)
        thisMax = numNebs = 0
        idMax = NA
        for (targ in targs) {
            numNebs = 0
            for (neb in 1:nrow(pp)) {
                if (is_neb(targ = pp[targ,], neb = pp[neb,], nRad = nRad, nPars = nPars)) {
                    numNebs = numNebs + 1
                }
            }
            if (numNebs > thisMax) {
                thisMax = numNebs
                idMax = pp[targ, "id"]
            }
        }
        data.frame(maxn = thisMax, id = idMax)
    })

    return( neighbors )
}

tst <- lapply(targs, FUN = function(targ) {
    lapply(1:nrow(pp), FUN = function(neb) {
        is_neb(pp
    })
}
