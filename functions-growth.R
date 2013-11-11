## Functions for fitting growth models

### keep track of parameter estimates from different growth models
# function to retrieve starting parameters, returns as list of params
getParams <- function(sr, spec, ind.var, dep.var, currentmodel, newpars=NULL) {
    parfile <- paste0("~/work/data/data/parameters/parameters-growth.csv")
    if(file.exists(parfile)) {
        pars <- read.csv(parfile)
        parrow <- which(pars$model==currentmodel & pars$mspec==spec & pars$msr==sr &
                        pars$mdep.var==dep.var & pars$mind.var==ind.var)
        ifelse(length(parrow)>0,
            ps <- as.list(pars[parrow,-grep("^m",names(pars))]),
            ps <- as.list(pars[pars$model == "start",-grep("^m",names(pars))]))
    ps <- ps[!is.na(ps)]
    }
    ps
}
