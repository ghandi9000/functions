################################################################################
##
##             Functions having to do with coordinate tranforms
##
################################################################################
## Polar to cartesian
## deg: if theta is in degrees
pol2cart <- function(r, theta, deg = FALSE) {
    if (deg) theta <- theta * pi/180
    xx <- r * cos(theta)
    yy <- r * sin(theta)
    ## Account for machine error in trig functions
    if (abs(xx) < 2e-16) xx <- 0
    if (abs(yy) < 2e-16) yy <- 0
    return( c(x = xx, y = yy) )
}

## Cartesian to polar
## deg: if theta is in degrees
cart2pol <- function(x, y, deg = FALSE) {
    r <- sqrt(x^2 + y^2)
    theta <- atan(y / x)
    if (x < 0) {
        theta <- theta + pi
    } else if (y < 0) {
        theta <- theta + 2*pi
    }
    if (deg) theta <- theta * 180/pi
    return( c(r = r, theta = theta) )
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
