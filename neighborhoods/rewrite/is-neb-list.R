## is_neb takes lists instead of data.frames

nPars2 <- quote(
    !is.na(neb[["ba"]]) & neb[["ba"]] >= targ[["ba"]]
    )

is_neb_list <- function(targ, neb, nRad, nPars) {
    return(targ[["id"]] != neb[["id"]] &
           targ[["bqudx"]] + nRad > neb[["bqudx"]] &
           targ[["bqudx"]] - nRad < neb[["bqudx"]] &
           targ[["bqudy"]] + nRad > neb[["bqudy"]] &
           targ[["bqudy"]] - nRad < neb[["bqudy"]] &
           targ[["time"]] == neb[["time"]] &
           eval(nPars)
    )
}
