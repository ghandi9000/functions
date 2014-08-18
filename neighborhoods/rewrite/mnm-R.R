################################################################################
##

##                          Neighbor Matrices in R
##
################################################################################

## data
dat <- read.csv("~/work/data/moose/moose-long.csv")

## Target parameters
targPs <- list(spec = "ABBA", bqudx %in% c(4,5,6))

## Take row numbers of targets and data.frame of neighbors, sr
maxn <- function(
