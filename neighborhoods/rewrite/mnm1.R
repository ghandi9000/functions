## uses doSNOW for parallel backend
## is_neb data.frame version
mnm1 <- function(tPars, nPars, dPars, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE) {
    require(plyr)
    require(parallel)
    require(doSNOW)
    nCores <- detectCores()
    cl <- makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)

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
    neighbors <- dlply(dd, .(pplot, time), .parallel = TRUE, function(pp) {
        source("~/work/functions/neighborhoods/rewrite/is-neb.R")
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
    stopCluster(cl)
    return( neighbors )
}
