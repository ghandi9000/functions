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
