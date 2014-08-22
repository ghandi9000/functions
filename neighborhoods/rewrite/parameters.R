## Target and neighbor parameters
tPars <- quote(!is.na(ba) &
               ba > 0 &
               spec %in% c("ABBA", "BECO", "PIRU"))

nPars <- quote(!is.na(neighbor[["ba"]]) &
               neighbor[["ba"]] >= target[["ba"]])

nCols <- c("ba", "id", "bqudx", "bqudy", "x", "y", "z", "spec", "cpos")

dPars <- quote(!is.na(ba) &
               stat == "ALIVE" &
               !is.na(bqudx) &
               bqudx < 11 &
               bqudx > 0 &
               !is.na(bqudy) &
               bqudy < 11 &
               bqudy > 0 &
               pplot > 3)
nRad <- 2

## Compare parallel and sequential versions
source("~/work/functions/neighborhoods/rewrite/mnm.R")
library(rbenchmark)
dat <- read.csv("~/work/data/moose/moose-long.csv")

benchmark(
    tst1 <- mnm(tPars = tPars, nPars = nPars, dPars = dPars, nCols = nCols,
               nRad = nRad, dat = dat, parallel=T),
    tst2 <- mnm(tPars = tPars, nPars = nPars, dPars = dPars, nCols = nCols,
               nRad = nRad, dat = dat, parallel=F),
    columns = c("test", "replications", "elapsed", "relative"),
    order = "relative", replications = 1
    )

tst <- mnm(tPars = tPars, nPars = nPars, dPars = dPars, nCols = nCols,
           nRad = nRad, dat = dat, parallel=F)

