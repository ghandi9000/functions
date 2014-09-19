################################################################################
##
##                                 Geometry
##
################################################################################

################################################################################
##
##                               Solid angles
##
################################################################################

## http://en.wikipedia.org/wiki/Subtended_angle
## Solid angle subtended by circular cone of opening angle, theta
sr_cone <- function(theta, deg=FALSE) {
    if (deg) theta <- theta * pi/180
    2*pi*(1 - cos(theta/2))
}

library(geometry)
require(rgl)

fd = function(p, ...) sqrt((p^2)%*%c(1,1,1)) - 1
                                        # also predefined as `mesh.dsphere'
fh = function(p,...)  rep(1,nrow(p))
                                        # also predefined as `mesh.hunif'
bbox = matrix(c(-1,1),2,3)
p = distmeshnd(fd,fh,0.3,bbox, maxiter=100)
                                        # this may take a while:
                                        # press Esc to get result of current iteration
## End(Not run)

