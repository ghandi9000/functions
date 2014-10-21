################################################################################
##
##             Functions having to do with coordinate tranforms
##
################################################################################
## Polar to cartesian
## deg: if input theta is in degrees
pol2cart <- function(r, theta, deg = FALSE, recycle = FALSE) {
    if (deg) theta <- theta * pi/180
    if (length(r) > 1 && length(r) != length(theta) && !recycle)
        stop("'r' vector different length than theta, if recycling 'r' values is desired 'recycle' must be TRUE")
    xx <- r * cos(theta)
    yy <- r * sin(theta)
    ## Account for machine error in trig functions
    xx[abs(xx) < 2e-15] <- 0
    yy[abs(yy) < 2e-15] <- 0
    return( cbind(x = xx, y = yy) )
}

## Cartesian to polar
## deg: returns theta in degrees
cart2pol <- function(x, y, deg = FALSE) {
    r <- sqrt(x^2 + y^2)
    theta <- atan(y / x)
    theta[x < 0] <- theta[x < 0] + pi
    theta[x >= 0 & y < 0] <- theta[x >= 0 & y < 0] + 2*pi
    if (deg) theta <- theta * 180/pi
    return( cbind(r = r, theta = theta) )
}

## Rotate point about origin by angle
## - x,y are initial point
## - angle is angle to rotate in degrees (counterclockwise)
rotate_point <- function(x, y, theta_r, deg = FALSE) {
    if (deg) theta_r <- theta_r * pi/180 # theta in radians
    rot <- matrix(c(cos(theta_r), sin(theta_r), -sin(theta_r), cos(theta_r)),
                  nrow = 2, ncol=2)
    p <- matrix(c(x,y))
    return ( rot %*% p )
}

## Euclidean distance
## euc <- function(pxi, pyi, pxs, pys, recycle=FALSE) {
##     if (!recycle && length(pxi) > 1 && length(pxi) != length(pxs))
##         stop("Should initial points be recycled?")
##     return( sqrt((pxi - pxs)^2 + (pyi - pys)^2) )
## }
euc <- function(a, b) {
    a <- as.matrix(a)
    b <- as.matrix(b)
    return ( sqrt(colSums((a-b)^2)) )
}
