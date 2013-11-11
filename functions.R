### Functions and libraries for Allometry Project
### load packages, functions, settings...
library(ggplot2)
library(plyr)
library(lattice)
##library(spatstat)
##library(odfWeave)
library(grid)
library(xtable)
##library(hydroGOF)
library(bbmle)

## Gives string to match data directory on windows or linux
getDataDir <- function(project) {
    ifelse(Sys.info()["sysname"] == "Linux",
           paste0("~/",project,"/data/"),
           paste0("C:/R/Current/",project,"/data/"))
}

## Just capitalizes the first letter of each string in a character
##  vectors
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

# Function to calculate compound growth, returns either grow rate
#  of actual growth depending on est.yrs parameter, if single values are
#  given for size and prior an exact value is returned, if vectors are given,
#  r is estimated by least squares fit to log-transformed data.
#  y = y0(1 + r)^t
#  y^(1/6) = y0(1+r) # 6 compounded growth periods
#  ln(y) = A + tB + err; where A = ln(y0), B = ln(1 + r)
#  r = exp(B) - 1
get.compr <- function(size, prior, obs.yrs, est.yrs = NULL) {
    ifelse(length(size)>1,
       { comp.fit <- lm(log(size^(1/obs.yrs)) ~ log(prior))
         r <- exp(coef(comp.fit)[[1]]) - 1 },
           r <- (size/prior)^(1/obs.yrs) - 1)
    ifelse(!missing(est.yrs),
           return(prior*(1 + r)^est.yrs),
           return(r))
}

# update the start values to be the mean of other model predictions for corresponding
#  value
update.start <- function(pars) {
    pars[pars$model=="start",unlist(lapply(pars, is.numeric))] <<-
        apply(pars[,unlist(lapply(pars, is.numeric))], 2, function(d) {
            median(d, na.rm=TRUE) })
}

# function to add new parameters to parameters file, takes a named vector of numbers
#  and currentmodel name as inputs.  Finds parameter file, checks for previous versions
#  of model fit and overwrites them or adds new row with new params
add.params <- function(sr, spec, ind.var, dep.var, newpars, currentmodel) {
    parfile <- paste0("parameters.csv")
    newrow <- as.data.frame(t(newpars))
    newrow[,c("msr","mspec","model","mdep.var","mind.var")] <-
        c(sr,spec,currentmodel,dep.var,ind.var)
    if(file.exists(parfile)) {
        pars <- read.csv(parfile)
        parrow <- which(pars$model==currentmodel & pars$mspec==spec & pars$msr==sr &
                        pars$mdep.var==dep.var & pars$mind.var==ind.var)
        if(length(parrow)>0) pars <- pars[-parrow,]
        pars <- rbind.fill(pars, newrow)
        update.start(pars)
        write.csv(pars, parfile, row.names = FALSE)
    }
}

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

# function to create neighbor matrices. Three matrices initialized and filled:
#  neighbor heights, bas, species
make.neighbor.matrices <- function(targets, neighbors, sr, bigger=FALSE, ind.var="ba",
                                   realdist=FALSE) {
    ifelse(realdist==FALSE,
           max.neighbors <- maxneighbors(targets, neighbors, sr),
           max.neighbors <- maxneighbors2(targets, neighbors, sr))
                                        # initialize matrices
    distances <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    bas <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
    species <<- matrix(NA, nrow=nrow(targets), ncol=max.neighbors)
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

# function to extract model AICs
getAIC <- function(list_of_models) {
    ldply(list_of_models,
          function(model) {
              if(!is.null(model)) {
                  c(
                      aic = extractAIC(model)
                      )
              }
          })
}

# function to get mse
getMSE <- function(list_of_models) {
    ldply(list_of_models,
          function(model) {
              if(!is.null(model)) {
                  c(
                      mse = mean(residuals(model)^2)
                      )
              }
          })
}

# function to return dataset of pvals for 2nd degree parameter and signs,
#  split by factor column
split.pvals <- function(indat, factor.col) {
    vals <- ddply(indat, factor.col, function(x) {
        if(nrow(x)>3) {
            c(summary(lm(x[["ht"]] ~ poly(x[["dbh"]],2)))$coefficients[3,][4],
              sign(summary(lm(x[["ht"]] ~ poly(x[["dbh"]],2)))$coefficients[3,][1]))
        }
    })
    names(vals) <- c("class","pval","sign")
    vals
}

# make columns of length l from list of named factors
make.cols <- function(factors,l) {
    as.data.frame(sapply(factors, function(x) rep(x,l)))
}

# function to subset a dataset by a list of factors, list of factors has names
#   matching columns from dataset
make.subset <- function(factors, dat) {
    f <- factors[!is.na(factors)]
    new <- dat
    for(i in 1:length(f)) {
        new <- new[new[,names(f)[i]]==f[[i]],]
    }
    new
}

## return all combinations of vector up to maximum length n
multicombn <- function(dat, n) {
    unlist(lapply(1:n, function(x) combn(dat, x, simplify=F)), recursive=F)
}

# Calculate higher order derivatives
DD <- function(expr, name, order=1) {
    if(order < 1) stop("'order' must be >= 1")
    if(order==1) D(expr, name)
    else DD(D(expr, name), name, order - 1)
}

## function to create classes based on quantiles
quantclass <- function(var,numbreaks,smallest=NULL) {
    if(!missing(smallest)) {
        cut(var, breaks = quantile(var[var>smallest],
                 probs = seq(0,1,1/numbreaks),
                 na.rm=TRUE))
    } else {
        cut(var, breaks = quantile(var, probs = seq(0,1,1/numbreaks),
                 na.rm=TRUE), include.lowest = TRUE)
    }
}

## make size even size classes, but largest is openended
sizeclasses <- function(var, nbreaks, nozero=FALSE) {
    if(nozero) {
        cut(var, breaks = seq(0.00001, max(var,na.rm=TRUE)+nbreaks, nbreaks))
    } else {
        cut(var, breaks = seq(0, max(var,na.rm=TRUE), max(var,na.rm=TRUE)/nbreaks),
            include.lowest = TRUE)
    }
}


### GGPLOT functions
#
## Make scatterplot graphs with 2nd degree polynomial fits: can be
#   colored by 'fit.by'.  If 'fit.by' is categorical, polynomial fits will be fit
#   to each category.  Output can be faceted by 'facet.by'
# input classes: 'indat'=your data, use strings for ('inx', 'iny', 'fit.by'),
#   and formula for 'facet.by'
make.poly <- function(indat, inx=NULL, iny=NULL, fit.by=NULL, facet.by=NULL) {
    d <- ggplot(data=indat, aes_string(x=inx,y=iny)) + geom_point(alpha=0.3,size=2)
    if(!missing(fit.by)) {
        d <- ggplot(data=indat, aes_string(x=inx, y=iny, color=fit.by)) +
            geom_point(alpha=0.3,size=2)
    }
    if(!missing(facet.by)) {
        d <- d + facet_wrap(facet.by)
    }
    d + geom_smooth(method="lm", formula = y ~ poly(x,2), se=FALSE, lwd=1)
}

## Give p-values for and direction for 2nd degree term from polynomial fits
# variables: 'ind'=independent, 'dep'= dependent, 'ind2' and 'ind3'
# are categorical variables by which to subset data for fits.
# 'ind' and 'dep' must be complete for polynomial fitting
poly.ps <- function(indat, ind=NULL, dep=NULL, ind2=NULL, ind3=NULL, ind4=NULL) {
    stopifnot(complete.cases(subset(indat, select=c(ind, dep))))
    # splitting by two categorical variables
    if(!missing(ind2) & !missing(ind3) & missing(ind4)) {
        stuff <- ddply(indat, .(get(ind2),get(ind3)),
                       .fun = function(d) {
                           c(nrow(d), mean(d[[dep]],na.rm=TRUE),
                             mean(d[[ind]],na.rm=TRUE) ) })
        names(stuff) <- c(as.name(ind2),as.name(ind3),"sample size",
                          paste("mean",dep),paste("mean",ind))
        pvals <- ddply(indat, .(get(ind2),get(ind3)),
                       .fun = function(d) {
                           summary(lm(d[[dep]]~poly(d[[ind]],2)))$coefficients[3,][4]
                       })
        names(pvals) <- c(as.name(ind2),as.name(ind3),"P.2nd")
        out.table <- merge(stuff,pvals, by = c(names(pvals)[1:2]))
        return(out.table)
    }
    # splitting by three variables
    if(!missing(ind2) & !missing(ind3) & !missing(ind4)) {
        stuff <- ddply(indat, .(get(ind2),get(ind3),get(ind4)),
                       .fun = function(d) {
                           c(nrow(d), mean(d[[dep]],na.rm=TRUE),
                             mean(d[[ind]],na.rm=TRUE) ) })
        names(stuff) <- c(as.name(ind2),as.name(ind3),as.name(ind4),"sample size",
                          paste("mean",dep),paste("mean",ind))
        pvals <- ddply(indat, .(get(ind2),get(ind3),get(ind4)),
                       .fun = function(d) {
                           summary(lm(d[[dep]]~poly(d[[ind]],2)))$coefficients[3,][4]
                       })
        names(pvals) <- c(as.name(ind2),as.name(ind3),as.name(ind4),"P.2nd")
        out.table <- merge(stuff,pvals, by = c(names(pvals)[1:3]))
        return(out.table)
    }
}


# Function to process data by an input function that will be applied to
#  subsets of dataset by all combinations of categorical variables
#  inputs: indat = data.frame, variable names = list, func = function
#  NOTE: output from function must be atomic or data.frame
allsubs <- function(indat, vars, func=NULL, out.name=NULL) {
    results <- data.frame()
    nvars <- rev(multicombn(vars,length(vars)))
    for(i in 1:length(nvars)) {
        results <-
            rbind.fill(results, ddply(indat, unlist(nvars[i]), func))
    }
    if(!missing(out.name)) names(results)[length(vars)+1] <- out.name
    results
}

# helper function with allsubs()
#  -returns pvals for 2nd degree polynomial fits
allsubs.poly <- function(x) { if(nrow(x) > 3) {
    summary(lm(x$ht ~ poly(x$dbh, 2)))$coefficients[3,][4] }
else NA }

### Like all subs, but outputs a list, so can store models
dlallsubs <- function(indat, vars, func=NULL, out.name=NULL) {
    results <- list()
    nvars <- rev(multicombn(vars,length(vars)))
    for(i in 1:length(nvars)) {
        results <-
            c(results, dlply(indat, unlist(nvars[i]), func))
    }
    results
}

# Function to coerce model output from dlallsubs to data.frame with useful
#  fitting statistics
modelstats <- function(list_of_models) {
    ldply(list_of_models,
          function(model) {
              c(
              #    aic = extractAIC(model),
                  deviance = deviance(model),
                  logLik = logLik(model),
                  confint = confint(model),
                  coef = coef(model)
                  )
          })
}

# Function to copy files from origin folder to backup folder, defaults to
#  current project folder to dropbox backup folder, also can removes files in destination
#  folder that arent in origin folder
backup <- function(origin = "C:/R/Current",
                   dest = "C:/Users/Noah/Dropbox/Backup/R/Allometry",
                   overwrite = TRUE,
                   pattern = "*.R$",
                   removeold = FALSE) {
    if(removeold == TRUE) {
        old = list.files(path = dest, pattern = pattern)
        lapply(old, function(d) {
            file.remove(file.path(dest,d))
        })
    }
    filenames = list.files(path = origin, pattern = pattern)
    lapply(filenames, function(d) {
        ifelse(overwrite == TRUE, file.copy(d, file.path(dest,d), overwrite = TRUE),
               file.copy(d, file.path(dest,d), overwrite = FALSE))
    })
}

### Automated fitting of neighborhood models by MLE
fit.MLE.models <- function(dat, sr, spec, ind.var, dep.var, models=NULL, bigger=TRUE,
                           method="Nelder-Mead", maxit=1000, savefits="currentfits.Rda",
                           realdist = FALSE) {
    srt <- max(sr) # if multiple sr, targets are those in all neighborhoods
    fits <- c()
    if(realdist == FALSE) {
        neighbors <<- subset(dat, bqudx < 11 & bqudx > 0 & bqudy < 11 &
                             bqudy > 0 & stat=="ALIVE")
        targets <<- subset(dat, bqudx < (12-srt) & bqudx > (-1 + srt) &
                           bqudy < (12 - srt) & bqudy > (-1 + srt) & stat=="ALIVE")
    }
    if(realdist == TRUE) {
        targets <<- subset(dat, abs(x) < (11-sr) & abs(y) < (11-sr) & stat=="ALIVE")
        neighbors <<- subset(dat, abs(x) <= 11 & abs(y) <= 11 & stat=="ALIVE")
    }
                                        # remove trees that dont satisfy certain conditions
    grew <- which(!is.na(targets[,dep.var]) & targets$spec==spec & targets[,dep.var]>0 &
                  targets[,ind.var]>0)
    targets <<- targets[grew,]
    for(i in sr) {  # make neighbor matrices
        print(paste("Making neighbor matrices for *", i, "* sized neighborhoods..."))
        make.neighbor.matrices(targets, neighbors, i, ind.var=ind.var, bigger=bigger,
                               realdist = realdist)
                                        # assign matrices in global for later access
        assign(paste0("species",i), species, envir = .GlobalEnv)
        assign(paste0("bas",i), bas, envir = .GlobalEnv)
        assign(paste0("distances",i), distances, envir = .GlobalEnv)
                                        # fit models
        if(!missing(models)) {
            print(paste("Fitting models with sr =", i))
            fits1 <- sapply(models, FUN=function(d) {
                print(paste("Model:", d))
                currentmodel <<- d
                ps <- get.params(sr = i, spec, ind.var, dep.var, d)
                print("Starting Parameters:"); print(unlist(ps,recursive = FALSE))
                parnames(normNLL) <<- c(names(ps))
                fit2 <- mle2(normNLL,
                             start = unlist(ps,recursive = FALSE),
                             data = list(x = targets[,dep.var]),
                             method = method,
                             control = list(maxit = maxit))
                add.params(sr=i, spec, ind.var, dep.var, newpars = coef(fit2), d)
                                        # add fit to current fits saved file
                tmp.env <- new.env() # environment to save fits in
                load(savefits, envir = tmp.env)
                assign(paste(d,sr,spec,sep = "."),fit2,envir=tmp.env)
                save(list=ls(all.names=TRUE, pos=tmp.env),
                     envir=tmp.env, file=savefits)
                rm(tmp.env)
                fit2
            })
            names(fits1) <- paste0(models,i)
            fits <- c(fits, fits1)
        }
    }
    fits
}

# log.likelihood function, normal distribution of residuals, goes with fit.MLE.models
normNLL <- function(params, x, currentmodel=NULL) {
    if(missing(currentmodel)) { currentmodel <- get.model() }
    sd = params[["sd"]]
    ind.var <- get.ind.var()
    mu = do.call(currentmodel, list(params,ind.var))
    -sum(dnorm(x, mean = mu, sd = sd, log = TRUE))
}
