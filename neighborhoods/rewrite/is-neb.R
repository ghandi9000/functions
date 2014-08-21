## Functions/parameters to source into parallel neighbor matrices
## Target and neighbor parameters
tPars <- quote(spec == "ABBA" &
               !is.na(ba) &
               ba > 0)

nPars <- quote(!is.na(neb$ba) &
               neb$ba >= targ$ba)

dPars <- quote(!is.na(ba) &
               stat == "ALIVE" &
               !is.na(bqudx) &
               bqudx < 11 &
               bqudx > 0 &
               !is.na(bqudy) &
               bqudy < 11 &
               bqudy > 0 &
               pplot > 3)
nRad <- 3

## Determines if a tree is a legit neighbor of a target tree
is_neb <- function(targ, neb, nRad, nPars) {
    return(targ$id != neb$id &
           targ$bqudx + nRad > neb$bqudx &
           targ$bqudx - nRad < neb$bqudx &
           targ$bqudy + nRad > neb$bqudy &
           targ$bqudy - nRad < neb$bqudy &
           targ$time == neb$time &
           eval(nPars)
    )
}
