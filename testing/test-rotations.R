### test-rotations.R --- 
## Filename: test-rotations.R
## Description: Test functions-coordinates.R rotations
## Author: Noah Peart
## Created: Tue Feb 10 15:50:52 2015 (-0500)
## Last-Updated: Tue Feb 10 15:51:30 2015 (-0500)
##           By: Noah Peart
######################################################################

source("~/work/functions/function-coordinates.R")

## Test point (1,0,0,1)
p <- matrix(c(1,0,0,1))

################################################################################
##
##                          Rotations about z-axis
##
################################################################################
## 45 degrees => (sqrt(2)/2, sqrt(2)/2, 0, 1)
rz(pi/4) %*% p

## 135 degrees => (-sqrt(2)/2, sqrt(2)/2, 0, 1)
rz(3*pi/4) %*% p

## 180 degrees => (-1, 0, 0, 1)
rz(pi) %*% p

## 225 degrees => (-sqrt(2)/2, -sqrt(2)/2, 0, 1)
rz(5*pi/4) %*% p

## 270 degrees => (0, -1, 0, 1)
rz(6*pi/4) %*% p

################################################################################
##
##                          Rotations about y-axis
##
################################################################################
## 45 degrees => (sqrt(2)/2, 0, -sqrt(2)/2, 1)
ry(pi/4) %*% p

## 135 degrees => (-sqrt(2)/2, 0, -sqrt(2)/2, 1)
ry(3*pi/4) %*% p

## 180 degrees => (-1, 0, 0, 1)
ry(pi) %*% p

## 225 degrees => (sqrt(2)/2, 0, -sqrt(2)/2, 1)
ry(5*pi/4) %*% p

## 270 degrees => (0, 0, 1, 1)
ry(6*pi/4) %*% p

################################################################################
##
##                          Rotations about x-axis
##
################################################################################
p <- c(0, 1, 0, 1)
## 45 degrees => (0, sqrt(2)/2, sqrt(2)/2, 1)
rx(pi/4) %*% p

## 135 degrees => (0, -sqrt(2)/2, sqrt(2)/2, 0, 1)
rx(3*pi/4) %*% p

## 180 degrees => (0, -1, 0, 1)
rx(pi) %*% p

## 225 degrees => (sqrt(2)/2, 0, -sqrt(2)/2, 1)
ry(5*pi/4) %*% p

## 270 degrees => (0, 0, -1, 1)
rx(6*pi/4) %*% p
