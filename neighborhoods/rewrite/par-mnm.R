## Parallel version of mnm using foreach instead of ddply

## Version to work on Windows: using doSNOW for parallel backend
## compute max neighbor number for each plot
mnm <- function(tPars, nPars, dPars, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE) {
    require(plyr)
    require(snow)
    require(doSNOW)
    require(foreach)
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
    ## Parallel foreach version: surprisingly slow when using 100% CPU
    subs <- expand.grid(pplot = unique(dd$pplot), time = unique(dd$time))
    foreach(i = 1:nrow(subs), .combine = rbind) %dopar% {
        ## source("~/work/functions/neighborhoods/rewrite/is-neb.R")
        pp <- dd[dd$pplot == subs[i, "pplot"] & dd$time == subs[i, "time"], ]
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
        c(pplot = subs[i, "pplot"], time = subs[i, "time"], maxn = thisMax, id = idMax)
    }

    ## ddply with lapply loops
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
    stopCluster(cl)
    return( neighbors )
}

################################################################################
##
##                             Test on real data
##
################################################################################
dat <- read.csv("~/work/data/moose/moose-long.csv")

## Target and neighbor parameters
tPars <- quote(spec == "ABBA" &
               !is.na(ba) &
               ba > 0 &
               !is.na(bagrowth) &
               bagrowth >= 0 &
               time == 98)

nPars <- quote(!is.na(neb$ba) &
               neb$ba >= targ$ba &
               neb$spec == "BECO")

dPars <- quote(!is.na(ba) &
               stat == "ALIVE" &
               !is.na(bqudx) &
               bqudx < 11 &
               bqudx > 0 &
               !is.na(bqudy) &
               bqudy < 11 &
               bqudy > 0 &
               pplot > 3)
nRad <- 3

## Run
tst <- dat[dat$pplot %in% 4:10, ]
mxx <- maxn(tPars=tPars, nPars=nPars, dPars=dPars, dat=dat, nRad=nRad)








################################################################################
##
##                          Number of comparisons
##
################################################################################
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

numComps <- ddply(dd, .(pplot, time), .parallel = T, .fun = function(x) {
    nTargs <- sum(x$targ)
    nComps <- nTargs * nrow(x)
    data.frame(nComps = nComps)
})

nComps <- sum(numComps$nComps)



p4 <- dd[dd$pplot == 4 & dd$time == 98, ]


################################################################################
##
##                                 Profiling
##
################################################################################
tmp <- tempfile()
Rprof(tmp, line.profiling=TRUE)
eval(parse(file = "~/work/functions/neighborhoods/rewrite/test-par.R", keep.source=TRUE))
Rprof(NULL)

summaryRprof(tmp, lines = "show")

## Test is_neb() function
library(rbenchmark)
source("~/work/functions/neighborhoods/rewrite/is-neb.R")
source("~/work/functions/neighborhoods/rewrite/is-neb-list.R")
source("~/work/functions/neighborhoods/rewrite/is-neb-noeval.R")
targ <- p4[which(p4$targ)[[1]], ] # choose the first target in pplot 4, time 98
targ2 <- as.list(targ)
neb <- p4[1,] # potential neighbor is first tree
neb2 <- as.list(neb)

benchmark(is_neb(targ, neb, nRad, nPars),
          is_neb(targ2, neb2, nRad, nPars),
          is_neb_noeval(targ, neb, nRad, nPars),
          is_neb_list(targ2, neb2, nRad, nPars2),
          columns = c("test", "replications", "elapsed", "relative"),
          order = "relative", replications = 50000)

## Note: with eval(nPars) 1.26 times slower
