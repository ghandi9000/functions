## Functions used to create neighbor matrices and other similar tasks
source("~/work/functions/neighborhoods/maxneighbors.R")

# function to find the maximum neighbors for a target given a neighborhood radius in
#  real distance
maxneighbors_real<- function(targets, neighbors, sr) {
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
maxneighbors_disc <- function(targets, neighbors, sr) {
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
## Returns list of the following matrices with each row corresponding to a target
##  and the columns to neighbors (number of columns = the maximum number of neighbors for a single target):
## ** Matrices:
## - distances: distance between targets and neighbors
## - species: species of neighbors
## - variable: the variable chosen to quantify neighbors (i.e. size)
## - direction_x: direction from target to neighbors in x-direction
## - direction_y: direction from target to neighbors in y-direction
##
## ** Parameters:
## - dat: data (targets will be created by removing those in the edge shadow)
##  * targets: subset of data containing the targets (only includes "ALIVE")
##  * neighbors: subset of data including those in the edge shadow (but not outside range)
## - sr: the radius around targets to search for neighbors
## - bigger: TRUE indicates that only neighbors larger than the target (by 'ind.var') should be counted
## - ind.var: variable to measure neighbors
## - realdist: FALSE indicates x,y-coordinates are in quadrats (location assumed to be center of quadrat),
##    TRUE indates exact x,y-coordinates (euclidean distances)
make.neighbor.matrices <- function(targets, neighbors, sr, ind.var="ba", range=12,
                                   realdist=FALSE, bigger=FALSE) {
    ## define targets and neighbors
    ## targets <- subset(dat, bqudx < (12-sr) & bqudx > (-1 + sr) & bqudy < (12 - sr) &
    ##                   bqudy > (-1 + sr) & stat=="ALIVE")
    ## neighbors <- subset(dat, bqudx < 11 & bqudx > 0 & bqudy < 11 &
    ##                     bqudy > 0 & stat=="ALIVE")
    ifelse(realdist,
           { max.neighbors <- maxneighbors_real(targets, neighbors, sr) },
           { max.neighbors <- maxneighbors_disc(targets, neighbors, sr) })
    ## initialize matrices
    distances <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    bas <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    species <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    direction_x <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    direction_y <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    ## populate matrices
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
                distances[i,1:nrow(nebs)] <-
                neighdist(targets$bqudx[i]*2,targets$bqudy[i]*2,
                          nebs$bqudx*2, nebs$bqudy*2, addifsame = TRUE)
                bas[i,1:nrow(nebs)] <- nebs[,ind.var]
                species[i,1:nrow(nebs)] <- nebs$spec
                direction_x[i,1:nrow(nebs)] <- targets$bqudx[i] - nebs$bqudx
                direction_y[i,1:nrow(nebs)] <- targets$bqudy[i] - nebs$bqudy
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
                distances[i,1:nrow(nebs)] <-
                    neighdist(targets$x[i],targets$y[i],
                              nebs$x, nebs$y, addifsame = FALSE)
                bas[i,1:nrow(nebs)] <- nebs[,ind.var]
                species[i,1:nrow(nebs)] <- nebs$spec
                direction_x[i,1:nrow(nebs)] <- targets$x[i] - nebs$x
                direction_y[i,1:nrow(nebs)] <- targets$y[i] - nebs$y
            }
        }
    }
    return( list(distances = distances, variable = bas, species = species,
                 direction_x = direction_x, direction_y = direction_y) )
}

## Function to find targets for neighborhood analysis
findTargets <- function(sr) {

}

## Function to get number of quadrats in neighborhood not including targets' quadrat
## **Expand to work for realdist later
surround <- function(sr) {
 (sr + 1)^2 - 1
}

## Function to extract the percentage of quadrats containing neighbors from
##  list of neighbor matrices
## - sr: square radius
## - targets: dataframe containing targets (edge effects removed)
## - nmatrices: result from running make.neighbor.matrices on targets
##  - use mats[["direction_x"]]  and mats[["direction_y"]] from nmatrices
percSurround <- function(sr, targets, nmatrices) {
    nebsize <- surround(sr)
    crowd <- sapply(1:nrow(targets), FUN = function(i) {
        rows <- unique(data.frame(
            row_x = mats[["direction_x"]][i,],
            row_y = mats[["direction_y"]][i,]))
        rows <- rows[complete.cases(rows),]
        rows <- rows[!(rows[,1] == 0 & rows[,2] == 0),]
        nrow(rows) / nebsize
    })
}

#####################################################################################
## Testing
mnm <- function(targets, neighbors, sr, ind.var="ba", range=12,
                                   realdist=FALSE, bigger=FALSE) {
    ## define targets and neighbors
    ## targets <- subset(dat, bqudx < (12-sr) & bqudx > (-1 + sr) & bqudy < (12 - sr) &
    ##                   bqudy > (-1 + sr) & stat=="ALIVE")
    ## neighbors <- subset(dat, bqudx < 11 & bqudx > 0 & bqudy < 11 &
    ##                     bqudy > 0 & stat=="ALIVE")
    ifelse(realdist,
           { max.neighbors <- maxneighbors_real(targets, neighbors, sr) },
           { max.neighbors <- maxnebs_disc(targets$tag, neighbors$tag,
                                           neighbors$pplot, neighbors$stat,
                                           neighbors$bqudx, neighbors$bqudy,
                                           neighbors$time, sr) })
    ## initialize matrices
    distances <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    bas <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    species <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    direction_x <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    direction_y <- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    ## populate matrices
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
                distances[i,1:nrow(nebs)] <-
                neighdist(targets$bqudx[i]*2,targets$bqudy[i]*2,
                          nebs$bqudx*2, nebs$bqudy*2, addifsame = TRUE)
                bas[i,1:nrow(nebs)] <- nebs[,ind.var]
                species[i,1:nrow(nebs)] <- nebs$spec
                direction_x[i,1:nrow(nebs)] <- targets$bqudx[i] - nebs$bqudx
                direction_y[i,1:nrow(nebs)] <- targets$bqudy[i] - nebs$bqudy
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
                distances[i,1:nrow(nebs)] <-
                    neighdist(targets$x[i],targets$y[i],
                              nebs$x, nebs$y, addifsame = FALSE)
                bas[i,1:nrow(nebs)] <- nebs[,ind.var]
                species[i,1:nrow(nebs)] <- nebs$spec
                direction_x[i,1:nrow(nebs)] <- targets$x[i] - nebs$x
                direction_y[i,1:nrow(nebs)] <- targets$y[i] - nebs$y
            }
        }
    }
    return( list(distances = distances, variable = bas, species = species,
                 direction_x = direction_x, direction_y = direction_y) )
}
