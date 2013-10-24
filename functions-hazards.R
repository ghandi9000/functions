## Functions related to hazard analysis
library(segmented)

bclong1 <- read.csv("~/work/data/data/long-bc-derived.csv")
bc <- subset(bclong1, stat == "ALIVE" & bvgrowth>=0)
bc$sdpclass <- factor(as.numeric(bc$sdpclass))

tst <- subset(bc, install == 1 & plot == 16 & time == 76)
pb <- tst$priorbv
bv <- tst$bvgrowth
bounds <- cont2class(tst, "priorbv", 3)

######################################################################################
##
## rGR functions
##

## Function to select model to calculate rgr.  Goal is to find a well suited
##  model that also predicts bole volume growth that is not correlated to observed
##  bole volume at the end of the prediction period
## Model selection:
## - Choose best model by MSE
## - If best model significantly correlated to bv at end of period, try next best
## - Possible models: linear, power, segmented linear, polynomials
removeCorr <- function(dat, colm) {
    pval <- 0
    segs <- 2
    while (pval < 0.05) { # While the data is correlated
        bounds <- cont2class(dat, colm, segs) # size class boundaries
        ## Fit to each size class
        for (i in 1:nrow(bounds)) {

        }
    }
}

## Function takes a range of polynomials and returns the polynomial best suited
##  (determined by the above conditions)
remCorrPoly <- function(

