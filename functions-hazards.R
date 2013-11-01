## Hazard models and other related functions

## Hazard model, as in Scott's manuscript
## lambda = log(N/N0) for each plot, where N = number of trees
hazard <- function(lambda, b0, b1, b2, b3, b4, b5) {
    lambda = pars["lambda"]
    b0 = pars["b0"]
    b0 * lambda * LRS^(b1 + b3 * SDP + b4 * SI + b5 * SDP * SI) * rGR ^ b2
}
