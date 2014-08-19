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


## data
dat <- read.csv("~/work/data/moose/moose-long.csv")

## Target and neighbor parameters
tPars <- quote(spec == "ABBA" & bqudx %in% c(4,5,6))
nPars <- quote(neb$ba >= targ$ba & neb$ht >= targ$ht)

## Determines if a tree is a legit neighbor of a target tree
isNeb <- function(targ, neb, nRad, nPars) {
    return(targ$id != neb$id &
           targ$bqudx + nRad > neb$bqudx &
           targ$bqudx - nRad < neb$bqudx &
           targ$bqudy + nRad > neb$bqudy &
           targ$bqudy - nRad < neb$bqudy &
           targ$time == neb$time &
           eval(nPars)
    )
}


## compute max neighbor number for each plot
maxn <- function(tPars, nPars, dat, nRad) {
    require(plyr)
    dd <- dat[dat$stat == "ALIVE", ]
    dd$targ <- eval(tPars, dd)
    maxes <- ddply(dd, .(pplot), FUN = function(pp) {
        thisMax <- 0
        for (targ in 1:nrow(pp)) {
            if (pp[targ, "targ"]) {
                numNebs <- 0
                for (neb in 1:nrow(pp)) {
                    if (isNeb(targ = pp[targ,], neb = pp[neb,], nRad = nRad, nPars = nPars))
                        numNebs <- numNebs + 1
                }
            }
            print(numNebs)
            if (numNebs > thisMax) thisMax <- numNebs
        }
        data.frame(plot = unique(pp$pplot), maxn = thisMax)
    })
    return( maxes )
}

subset2 <- function(x, condition) {
  condition_call <- substitute(condition)
  r <- eval(condition_call, x)
  x[r, ]
}
subset2(sample_df, a >= 4)



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
filled <- sample(100, 10)
tst[tst$pplot == 1 & row.names(tst) %in% filled, "ba"] <- runif(10, 0, 1.5)

p1 <- tst[tst$pplot == 1,]
plot(p1$bqudx, p1$bqudy, type = "n")
text(p1[!is.na(p1$ba), "bqudx"],
     p1[!is.na(p1$ba), "bqudy"],
     labels = round(p1[!is.na(p1$ba), "ba"], 2))
text(p1[!is.na(p1$ba), "bqudx"],
     p1[!is.na(p1$ba), "bqudy"] + 0.25,
     labels = p1[!is.na(p1$ba), "id"], col = "red")

################################################################################
##
##                                Test isNeb
##
################################################################################
targ <- 37
neb <- 27


################################################################################
##
##                                 Run test
##
################################################################################
## Target and neighbor parameters
tPars <- quote(spec == "ABBA" & pplot %in% c(1) & !is.na(ba) & ba > 0)
nPars <- quote(neb$ba >= targ$ba) # & neb$ht >= targ$ht)

mxx <- maxn(tPars=tPars, nPars=nPars, dat=tst, nRad=2)
