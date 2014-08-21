
## Determines if a tree is a legit neighbor of a target tree
is_neb_noeval <- function(targ, neb, nRad, nPars) {
    return(targ$id != neb$id &
           targ$bqudx + nRad > neb$bqudx &
           targ$bqudx - nRad < neb$bqudx &
           targ$bqudy + nRad > neb$bqudy &
           targ$bqudy - nRad < neb$bqudy &
           targ$time == neb$time
           ## eval(nPars)
    )
}
