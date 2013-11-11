## Generate randomly located data based on bc plots
## Attributes of generated data:
## - same density as real plot
## - same distribution of sizes as real plot
## - CSR for coordinate data

source("~/work/functions/functions.R")

#################################################################################
##
## Functions to generate random spatial data
##
## Function that pulls uniformly samples sizes from given vector, returning
##  a specified number of samples (defaults to the length of the size vector)
randomSize <- function(sizes, numSamp = length(sizes)) {
    samp <- length(sizes)
    sizes[floor(runif(numSamp, 1, samp+1))]
}

## Function to copy and randomize the locations of a given plot
## Assumes that there are columns labeled x and y containing the location data
makeRandomPlot <- function(dat, newSize, sizeVar) {
    ## Get the area of the plot, assume square plots
    dat <- droplevels(dat)
    plotArea = ceiling(max(dat$x, na.rm = TRUE))^2
    newArea = newSize^2
    bound <- c(-newSize, newSize)
    numTrees <- newArea * (nrow(dat) / plotArea)
    newX <- runif(numTrees, bound[1], bound[2])
    newY <- runif(numTrees, bound[1], bound[2])
    newDat <- data.frame(x = newX, y = newY,
                         size = randomSize(dat[,sizeVar], numSamp = numTrees),
                         si = rep(mean(dat$si,na.rm=TRUE),numTrees),
                         sdpclass = rep(mean(dat$sdpclass), numTrees),
                         stat = factor(rep("ALIVE", numTrees), levels = c("ALIVE","DEAD")),
                         time = rep(1, numTrees), tag = 1:numTrees,
                         pplot = rep(1, numTrees),
                         spec = rep("FD", numTrees),
                         install = rep(mean(dat$install), numTrees),
                         plot = rep(mean(dat$plot), numTrees))
    newDat
}

## Function to make a given number of random plots from a dataset containing at least
##  one plot, now it assumes there are colunmns labeled "install" and "plot" -- like
##  bc data
makeManyRandomPlots <- function(dat, newSize, sizeVar, numPlots, plotN1 = "install",
                                plotN2 = "plot") {
    nPlots <- nrow(unique(dat[c(N1, N2)]))

}

## Test plot-generating functions
## bc <- read.csv("~/Allometry/data/long-bc-derived.csv")
## Variables
## newSize <- 50
## sizeVar <- "ba"
## numPlots <- 100

## tst <- subset(bc, install == 1 & plot == 16 & time == 76)
## tst2 <- makeRandomPlot(tst, newSize = newSize, sizeVar = sizeVar)
## ggplot(tst, aes(x,y,size = ba)) + geom_point()
## windows()
## ggplot(tst2, aes(x,y,size = size)) + geom_point()

#################################################################################
##
## Functions to grow the forest
##

## Alter the fit.MLE.functions() to create neighbor matrices for the generated data
## plotSize = 1/2 the length of one side (plots are created as squares around origin)
## nRadius = neighbor radius
## ind.var = attribute to measure neighbors (defaults to basal area)
## realdist = TRUE if exact locations are available
neighborMatrices <- function(dat, nRadius, plotSize, ind.var = "size", realdist = TRUE,
                             bigger = FALSE) {
    ## define the targets (excludes the neighbors in the buffer area)
    if(realdist == TRUE) {
        targets <<- subset(dat, abs(x) < (plotSize-nRadius) &
                           abs(y) < (plotSize-nRadius) & stat=="ALIVE")
        neighbors <<- subset(dat, abs(x) <= plotSize & abs(y) <= plotSize & stat=="ALIVE")
    }
    ## Determine the dimensions of the matrices by tree with maximum neighbors
    ifelse(realdist==FALSE,
           max.neighbors <- maxneighbors(targets, neighbors, nRadius),
           max.neighbors <- mNeighbors(targets, neighbors, nRadius))
                                        # initialize matrices
    distances <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    bas <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    species <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
                                        # populate matrices
    if(realdist==TRUE) {
        for(i in 1:nrow(targets)) {
            ifelse(bigger==TRUE,
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] &
                          neighdist(x,y,targets$x[i], targets$y[i]) <= nRadius &
                          get(ind.var) >= targets[,ind.var][i] &
                          time == targets$time[i]),
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] &
                          neighdist(x,y,targets$x[i], targets$y[i]) <= nRadius &
                          get(ind.var) >= targets[,ind.var][i] &
                          time == targets$time[i]))
            if (nrow(nebs) > 0) {
                distances[i,1:nrow(nebs)] <<-
                    neighdist(targets$x[i],targets$y[i],
                              nebs$x, nebs$y, addifsame = FALSE)
                bas[i,1:nrow(nebs)] <<- nebs[,ind.var]
                species[i,1:nrow(nebs)] <<- nebs$spec
            }
        }
    }
}



# function to find the maximum neighbors for a target given a neighborhood radius in
#  real distance
mNeighbors <- function(targets, neighbors, nRadius, ind.var = "size") {
    max.neighbors <- 0
    for(i in 1:nrow(targets)) {
        nebs <- subset(neighbors, pplot == targets$pplot[i] & tag!=targets$tag[i] &
                       neighdist(x,y,targets$x[i], targets$y[i]) <= nRadius &
                       get(ind.var) >= targets[,ind.var][i] & stat == "ALIVE" &
                       time==targets$time[i])
        max.neighbors <- max(max.neighbors, nrow(nebs), na.rm = TRUE)
    }
    max.neighbors
}
