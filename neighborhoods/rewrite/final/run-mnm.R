source("~/work/functions/neighborhoods/rewrite/parameters.R")
source("~/work/functions/neighborhoods/rewrite/mnm.R")
source("~/work/functions/neighborhoods/rewrite/mnm-to-matrix.R")

## data
dat <- read.csv("~/work/data/moose/moose-long.csv", stringsAsFactors = F)

## Compare parallel and sequential versions
## library(rbenchmark)
## benchmark(
##     tst1 <- mnm(tPars = tPars, nPars = nPars, dPars = dPars, nCols = nCols,
##                nRad = nRad, dat = dat, parallel=T),
##     tst2 <- mnm(tPars = tPars, nPars = nPars, dPars = dPars, nCols = nCols,
##                nRad = nRad, dat = dat, parallel=F),
##     columns = c("test", "replications", "elapsed", "relative"),
##     order = "relative", replications = 1
##     )

nLst <- mnm(tPars = tPars, nPars = nPars, dPars = dPars, nCols = nCols,
           nRad = nRad, dat = dat, parallel=F)
nMat <- mnm_to_matrix(nLst)
