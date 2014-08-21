neighbors <- dlply(dd, .(pplot, time), function(pp) {
    ## source("~/work/functions/neighborhoods/rewrite/is-neb.R")
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

