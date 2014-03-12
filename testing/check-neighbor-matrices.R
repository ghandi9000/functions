##
pp <- read.csv("~/work/data/data/dynamicallometry/moose-long-canopies.csv")
tst <- subset(pp, pplot == 8 & time == 98 & stat == "ALIVE")
## make.neighbor.matrices(tst, tst, sr=2)

library(ggplot2)
ggplot(tst, aes(bqudx, bqudy, col = factor(splot))) + geom_point()


## Make test matrices
# define targets and neighbors, for moose data
sr <- 2
spec = "ABBA"
dep.var <- "bagrowth"
ind.var <- "ba"

targets <- subset(tst, bqudx < (12-sr) & bqudx > (-1 + sr) & bqudy < (12 - sr) &
                  bqudy > (-1 + sr) & stat=="ALIVE")
neighbors <- subset(tst, bqudx < 11 & bqudx > 0 & bqudy < 11 &
                    bqudy > 0 & stat=="ALIVE")

# remove trees that dont satisfy certain conditions
grew <- which(!is.na(targets[,dep.var]) & targets$spec==spec & targets[,dep.var]>0)
targets <- targets[grew,]

## Check neigbors and targets visually
i <- 29
nebs <- subset(neighbors, pplot == targets$pplot[i] &
       tag!=targets$tag[i] & bqudx < targets$bqudx[i]+sr &
       bqudx > targets$bqudx[i]-sr & bqudy < targets$bqudy[i]+sr &
       bqudy > targets$bqudy[i]-sr &
       get(ind.var) >= targets[,ind.var][i] &
       time==targets$time[i])

ggplot(nebs, aes(bqudx-.5, bqudy-.5)) + geom_point(col = "Red") + xlim(0,10) + ylim(0,10) +
    geom_hline(yintercept = seq(0,10)) + geom_vline(xintercept = seq(0,10)) +
    geom_point(data = data.frame(x = targets[i,]$bqudx, y = targets[i,]$bqudy),
               aes(x-.5,y-.5), col = "Blue")
nebs$bqudx

source("~/work/functions/functions-neighborhood.R")
foo <- make.neighbor.matrices(targets, neighbors, sr, ind.var = ind.var, bigger = TRUE)
dim(foo[[1]])
names(foo)

## check distance_x and distance_y calculations
i <- 29
tst <- data.frame(x = c(targets[i,]$bqudx, targets[i,]$bqudx - foo[[4]][i,]),
                  y = c(targets[i,]$bqudy, targets[i,]$bqudy - foo[[5]][i,]),
                  type = c("target", rep("neighbor", length(foo[[1]][1,]))))

dev.new()
ggplot(tst, aes(x - 0.5, y - 0.5, col = type)) + geom_jitter() + xlim(0,10) + ylim(0,10) +
    geom_hline(yintercept = seq(0,10)) + geom_vline(xintercept = seq(0,10))

bqudx < targets$bqudx[i]+sr &
       bqudx > targets$bqudx[i]-sr

## ## Max neighbors
## max.neighbors <- 0
## for(i in 1:nrow(targets)) {
##     nebs <- subset(neighbors, pplot==targets$pplot[[i]] & tag!=targets$tag[[i]] &
##                    stat=="ALIVE" & bqudx < targets$bqudx[[i]]+sr & bqudx >
##                    targets$bqudx[[i]]-sr & bqudy < targets$bqudy[[i]]+sr & bqudy >
##                    targets$bqudy[[i]]-sr & time==targets$time[[i]])
##     max.neighbors <- max(max.neighbors, nrow(nebs), na.rm = TRUE)
##     if (max.neighbors == 36)
##         print (i);
## }
## max.neighbors

## Check make.neighbor.matrices with realdist
tst <- targets;
names(tst)[which(names(tst) %in% c("bqudx", "bqudy"))] <- c("x", "y")
tst.nebs <- neighbors
names(tst.nebs)[which(names(tst.nebs) %in% c("bqudx", "bqudy"))] <- c("x", "y")
foo <- make.neighbor.matrices(tst, tst.nebs, sr, ind.var = ind.var, bigger = TRUE,
                              realdist = TRUE)
