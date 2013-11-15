## Functions to help with data transformations

####################################################################################
##
## Data check
##

## check that columns exist, takes names of columns and year suffixes
checkCols <- function(dat, ns, times) {
    cnames <- unlist(lapply(ns, function (x) { paste0(x, times) } ))
    return( cnames[!cnames %in% names(dat)] )
}

####################################################################################
##
## Wide/long transform
##
## Find columns that can be transformed to long (i.e they have data from each time in
##  times)
colTrans <- function(dat, times) {
    toTrans = c()
    trans = grep("[[:digit:]]", names(dat), value=TRUE)
    n = gsub("[[:digit:]]", "", trans)
    for (i in 1:length(n)) {
        if (!n[i] %in% toTrans) {
            if (length(grep(n[i], n)) >= length(times)) {
                toTrans = c(toTrans, n[i])
            }
        }
    }
    return (toTrans)
}

## Given a vector of column names and times, drop columns that vary by time except
##  for those specified, return dataframe in long format
makeLong <- function(dat, ns, times) {
    vnames <- lapply(ns, function (x) { paste0(x, times) } )
    trans = grep("[[:digit:]]", names(dat), value=TRUE)
    noTrans = names(dat)[!names(dat) %in% trans] ## save cols that dont vary by year
    outdf <- dat[,c(noTrans, unlist(vnames))]
    long <- reshape(outdf, varying = vnames,
                    v.names = ns,
                    times = times,
                    direction = "long")
    return (long)
}

## Makes dummy columns of NAs for long transform
addDummies <- function(dat, ns) {
    outdf <- dat
    for (name in ns) {
        outdf[, name] <- rep(NA, nrow(dat))
    }
    return (outdf)
}

## Remove rows that have all NA values for given columns/yrs
removeEmpty <- function(dat, ns, yrs) {
    indices <- which (names(dat) %in% ns)
    outdf <- dat
    toKeep <- apply(dat, 1, function(x) { any(!is.na(x[indices])) })
    return ( outdf[toKeep,] )
}

