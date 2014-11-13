## Target and neighbor parameters
tPars <- quote(!is.na(ba) &
               ba > 0 &
               spec %in% c("ABBA", "BECO"))

nPars <- quote(!is.na(neighbor[["ba"]]) &
               neighbor[["ba"]] >= target[["ba"]])

nCols <- c("ba", "id", "bqudx", "bqudy", "x", "y", "z", "spec", "cpos")

dPars <- quote(!is.na(ba) &
               stat == "ALIVE" &
               !is.na(bqudx) &
               bqudx < 11 &
               bqudx > 0 &
               !is.na(bqudy) &
               bqudy < 11 &
               bqudy > 0 &
               pplot > 3)
nRad <- 2



