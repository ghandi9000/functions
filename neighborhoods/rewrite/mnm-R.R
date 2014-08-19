################################################################################
##
##                          Neighbor Matrices in R
##
################################################################################
## Parameters to set:
## Max radius/quadrat number
## Neighbor radius
## target parameters
## neighbor parameters

## Determines if a tree is a legit neighbor of a target tree
is_Neb <- function(targ, neb, nRad, nPars) {
    return(targ$id != neb$id &
           targ$bqudx + nRad > neb$bqudx &
           targ$bqudx - nRad < neb$bqudx &
           targ$bqudy + nRad > neb$bqudy &
           targ$bqudy - nRad < neb$bqudy &
           targ$time == neb$time &
           eval(nPars)
    )
}


is_Neb2 <- function(targ, neb, nRad, nPars) {
    return(t)
}

nPars <- quote(!is.na(neb[["ba"]]) & neb[["ba"]] >= targ[["ba"]])


blah <- sapply(which(pp$targ), FUN = function(targ)
       apply(pp, 1, function(x) findNebs(pp[targ,], x, nRad, nPars))
)

find_nebs <- function(targ, neb, nRad, nPars) {
    return(targ[["id"]] != neb[["id"]] &
           targ[["bqudx"]] + nRad > neb[["bqudx"]] &
           targ[["bqudx"]] - nRad < neb[["bqudx"]] &
           targ[["bqudy"]] + nRad > neb[["bqudy"]] &
           targ[["bqudy"]] - nRad < neb[["bqudy"]] &
           targ[["time"]] == neb[["time"]] &
           eval(nPars))
}

test <- lapply(targs, FUN = function(targ) {
    target <- pp[targ,]
    apply(pp, 1, function(neighbor) {
        target$id != neighbor[["id"]]
    })
})
## Version I: two explicit inner loops
## compute max neighbor number for each plot
maxn <- function(tPars, nPars, dPars, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE) {
    require(plyr)
    require(doMC)
    registerDoMC() # start multicore procs

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
        targs = which(pp$targ)
        thisMax = numNebs = 0
        idMax = NA
        for (targ in targs) {
            numNebs = 0
            for (neb in 1:nrow(pp)) {
                if (isNeb(targ = pp[targ,], neb = pp[neb,], nRad = nRad, nPars = nPars)) {
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
    return( maxes )
}

## Version II: One explicit inner loop
## compute max neighbor number for each plot
maxn2 <- function(tPars, nPars, dPars, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE) {
    require(plyr)
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
    maxes <- ddply(dd, .(pplot), function(pp) {
        targs = which(pp$targ)
        thisMax = numNebs = 0
        idMax = NA
        nebCounts = sapply(targs, FUN = function(targ) {
            numNebs = 0
            for (neb in 1:nrow(pp)) {
                if (isNeb(targ = pp[targ,], neb = pp[neb,], nRad = nRad, nPars = nPars)) {
                    numNebs = numNebs + 1
                }
            }
            c(id = pp[targ, "id"], nnebs = numNebs)
        })

        thisMax = which.max(nebCounts["nnebs",])
        if (nebCounts["nnebs", thisMax] > 0) {
            numNebs = nebCounts["nnebs", thisMax]
            idMax = nebCounts["id", thisMax]
        }
        data.frame(maxn = thisMax, id = idMax)
    })
    return( maxes )
}


## Version III: Store neighbor ids for later use
maxn3 <- function(tPars, nPars, dPars, dat, nRad,
                 pLims=c(xlower=1, xupper=10, ylower=1, yupper=10),
                 cushion=TRUE) {
    require(plyr)
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

    nebs <- ddply(dd, .(pplot, time), function(pp) {
        targs = which(pp$targ)
        thisMax = numNebs = 0
        idMax = NA
        for (targ in targs) {
            numNebs = 0
            for (neb in 1:nrow(pp)) {
                if (isNeb(targ = pp[targ,], neb = pp[neb,], nRad = nRad, nPars = nPars)) {
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
    return( maxes )
}


################################################################################
##
##                                 Test data
##
################################################################################
set.seed(0)
quads <- expand.grid(x = 1:10, y = 1:10)
tst <- data.frame(pplot = rep(1:3, each = 100),
                  bqudx = rep(quads$x, 3),
                  bqudy = rep(quads$y, 3),
                  ba = NA,
                  time = 1,
                  id = rep(1:100, 3),
                  spec = rep("ABBA", 300),
                  stat = rep("ALIVE", 300))

## Plot 1
filled <- sample(100, 30)
tst[tst$pplot == 1 & row.names(tst) %in% filled, "ba"] <- runif(10, 0, 1.5)

p1 <- tst[tst$pplot == 1,]
plot(p1$bqudx, p1$bqudy, type = "n")
abline(v = seq(0.5, 10.5, 1), h = seq(0.5, 10.5, 1))
text(p1[!is.na(p1$ba), "bqudx"],
     p1[!is.na(p1$ba), "bqudy"],
     labels = round(p1[!is.na(p1$ba), "ba"], 2))
text(p1[!is.na(p1$ba), "bqudx"],
     p1[!is.na(p1$ba), "bqudy"] + 0.25,
     labels = p1[!is.na(p1$ba), "id"], col = "red")


################################################################################
##
##                                 Run test
##
################################################################################
## Target and neighbor parameters
tPars <- quote(spec == "ABBA" & pplot %in% c(1) & !is.na(ba) & ba > 0)
nPars <- quote(!is.na(neb$ba) & neb$ba >= targ$ba) # & neb$ht >= targ$ht)
dPars <- quote(!is.na(ba) &
               stat == "ALIVE" &
               !is.na(bqudx) &
               !is.na(bqudy))
nRad <- 3

mxx <- maxn(tPars=tPars, nPars=nPars, dPars=dPars, dat=tst, nRad=nRad)

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

## Sample plot 4
p4 <- dat[dat$pplot == 4, ]
p4$targ <- eval(tPars, p4)

## Visualize
plot(p4$bqudx, p4$bqudy, type = "n")
abline(v = seq(0.5, 10.5, 1), h = seq(0.5, 10.5, 1))
text(jitter(p4[p4$targ, "bqudx"]),
     p4[p4$targ, "bqudy"],
     labels = format(p4[p4$targ, "ba"], sci=T, digits = 2))
text(jitter(p4[p4$targ, "bqudx"]),
     jitter(p4[p4$targ, "bqudy"] + 0.25),
     labels = p4[p4$targ, "id"], col = "red")

ps <- dat[dat$pplot %in% c(4,5,6), ]
mxx <- maxn(tPars=tPars, nPars=nPars, dPars=dPars, dat=dat, nRad=nRad)

mxx2 <- maxn2(tPars=tPars, nPars=nPars, dPars=dPars, dat=ps, nRad=nRad)
