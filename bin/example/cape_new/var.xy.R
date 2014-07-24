## function that computes the variance matrix of the product of 2 matrices

## input
##  x = x matrix
##  y = y matrix
##  sigma.x = covariance matrix of x (square of standard error) same dimensions as x
##  sigma.y = covariance matrix of y, same dimensions as y

## output
##  returns the covariance matrix of x*y based on Goodman's expression

var.xy <- function(x,y,sigma.x,sigma.y) {
  goodman <- function(xx,yy,sx,sy) {
    xx^2*sy + yy^2*sx + sx*sy
  }
  sigma.xy <- matrix(0,nrow=nrow(x),ncol=ncol(y))
  for (i in 1:nrow(x)) {
    xvec <- x[i,]
    sigx.vec <- sigma.x[i,]
    for (j in 1:ncol(y)) {
      yvec <- y[,j]
      sigy.vec <- sigma.y[,j]
      for (k in 1:ncol(x)) {
        sigma.xy[i,j] <- sigma.xy[i,j] + goodman(xvec[k],yvec[k],sigx.vec[k],sigy.vec[k])
      }
    }
  }
  sigma.xy
}