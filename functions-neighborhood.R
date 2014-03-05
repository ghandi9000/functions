## Functions used to create neighbor matrices and other similar tasks

# function to find the maximum neighbors for a target given a neighborhood radius in
#  real distance
maxneighbors2 <- function(targets, neighbors, sr) {
    max.neighbors <- 0
    for(i in 1:nrow(targets)) {
        nebs <- subset(neighbors, pplot == targets$pplot[i] & tag!=targets$tag[i] &
                       neighdist(x,y,targets$x[i], targets$y[i]) <= sr &
                       get(ind.var) >= targets[,ind.var][i] & stat == "ALIVE" &
                       time==targets$time[i])
        max.neighbors <- max(max.neighbors, nrow(nebs), na.rm = TRUE)
    }
    max.neighbors
}

# function to calculate max neighbors, faster than apply, ~ 1.5 mins for sr2
maxneighbors <- function(targets, neighbors, sr) {
    max.neighbors <- 0
    for(i in 1:nrow(targets)) {
        nebs <- subset(neighbors, pplot==targets$pplot[[i]] & tag!=targets$tag[[i]] &
                       stat=="ALIVE" & bqudx < targets$bqudx[[i]]+sr & bqudx >
                       targets$bqudx[[i]]-sr & bqudy < targets$bqudy[[i]]+sr & bqudy >
                       targets$bqudy[[i]]-sr & time==targets$time[[i]])
        max.neighbors <- max(max.neighbors, nrow(nebs), na.rm = TRUE)
    }
    max.neighbors
}

# distance between neighbors
#  if neighbor is in same quadrat as target give it a small distance between them (0.5)
neighdist<-function(targetx, targety, neighborx, neighbory, addifsame=FALSE) {
  if(addifsame==TRUE) {
      ifelse(sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)==0,
             return(0.5),
             return(sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)))
  }
  if(addifsame==FALSE) sqrt((targetx-neighborx)^2 + (targety-neighbory)^2)
}

## Function to create neighbor matrices.
## Creates (in the calling environment) the following matrices with each row corresponding to a target
##  and the columns to neighbors (number of columns = the maximum number of neighbors for a single target):
## ** Matrices:
## - distances: distance between targets and neighbors
## - species: species of neighbors
## - variable: the variable chosen to quantify neighbors (i.e. size)
## - direction: direction from target to neighbors
##
## ** Parameters:
## - targets: subset of data containing the targets (i.e minus those that would fall in edge shadow)
## - neighbors: subset of data including those in the edge shadow
## - sr: the radius around targets to search for neighbors
## - bigger: TRUE indicates that only neighbors larger than the target (by 'ind.var') should be counted
## - ind.var: variable to measure neighbors
## - realdist: FALSE indicates x,y-coordinates are in quadrats (location assumed to be center of quadrat),
##    TRUE indates exact x,y-coordinates (euclidean distances)
make.neighbor.matrices <- function(targets, neighbors, sr, bigger=FALSE, ind.var="ba",
                                   realdist=FALSE) {
    ifelse(realdist==FALSE,
           max.neighbors <- maxneighbors(targets, neighbors, sr),
           max.neighbors <- maxneighbors2(targets, neighbors, sr))
                                        # initialize matrices
    distances <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    bas <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    species <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    direction <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
                                        # populate matrices
    if(realdist==FALSE) {
        for(i in 1:nrow(targets)) {
            ifelse(bigger==TRUE,
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] & bqudx < targets$bqudx[i]+sr &
                          bqudx > targets$bqudx[i]-sr & bqudy < targets$bqudy[i]+sr &
                          bqudy > targets$bqudy[i]-sr &
                          get(ind.var) >= targets[,ind.var][i] &
                          time==targets$time[i]),
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] & bqudx < targets$bqudx[i]+sr &
                          bqudx > targets$bqudx[i]-sr & bqudy < targets$bqudy[i]+sr &
                          bqudy > targets$bqudy[i]-sr & time==targets$time[i]))
            if (nrow(nebs) > 0) {
                distances[i,1:nrow(nebs)] <<-
                neighdist(targets$bqudx[i]*2,targets$bqudy[i]*2,
                          nebs$bqudx*2, nebs$bqudy*2, addifsame = TRUE)
                bas[i,1:nrow(nebs)] <<- nebs[,ind.var]
                species[i,1:nrow(nebs)] <<- nebs$spec
            }
        }
    }
    if(realdist==TRUE) {
        for(i in 1:nrow(targets)) {
            ifelse(bigger==TRUE,
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] &
                          neighdist(x,y,targets$x[i], targets$y[i]) <= sr &
                          get(ind.var) >= targets[,ind.var][i] &
                          time == targets$time[i]),
                   nebs <-
                   subset(neighbors, pplot == targets$pplot[i] &
                          tag!=targets$tag[i] &
                          neighdist(x,y,targets$x[i], targets$y[i]) <= sr &
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

## Function to calculate the direction to a neighbor
direction <- function()
