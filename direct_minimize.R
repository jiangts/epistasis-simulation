# Direct search method (directly regress against values of "delta")

################## OPTIMIZATION METHOD SUMMARY ##################
# Method "BFGS" is a quasi-Newton method (also known as a variable metric algorithm),
# specifically that published simultaneously in 1970 by Broyden, Fletcher, Goldfarb 
# and Shanno. This uses function values and gradients to build up a picture of the 
# surface to be optimized.
#
# Method "SANN" is by default a variant of simulated annealing given 
# in Belisle (1992).
#
# convergence codes:
# 1: indicates that the iteration limit maxit had been reached.
# 10: indicates degeneracy of the Nelder–Mead simplex.
#################################################################
library("reshape")
library("MASS")
library("minpack.lm")

# calculate cost of unrolled delta matrix given matrix of beta values
#   This is similar to solving for deltas from the betas. In fact it is 
#   the exact same if we let B = A
cost_fn <- function(unr_A) #, n_v, betas)
{ 
  #shape into matrix and insert diagonal of 0's
  A <- reshape.with.diag(unr_A, n_v)
  B <- A + A%^%2
  
  betas <- betas[-1,] #removes beta_0's
  me <- betas[1:n_v,] #main effect matrix
  betas <- betas[(1+n_v):dim(betas)[1],] #use interaction betas
  combos <- t(combn(1:n_v, 2))
  
  x_i <- matrix(0, dim(betas)[1], 2)
  total_err = 0
  for(i in 1:dim(betas)[1])
  {
    pair <- combos[i,]
    y_i <- matrix(c(betas[i,])) #pair effect. "correct" solns
    wgts <- t(me[pair,]) #needs to be transposed! (to have pheno going horiz)
    
    x_i[,1] = B[pair[2],pair[1]] #21
    x_i[,2] = B[pair[1],pair[2]] #12
    hypo <- as.matrix(rowSums(wgts * x_i))
    sq_err <- sum((hypo - y_i)^2)
    total_err <- total_err + sq_err
  }
  return(total_err)
}
#####
#reshape an unrolled vector into a matrix with 0s on the main diagonal
reshape.with.diag <- function(unr, s)
{
  ind <- s*0:(s-1)
  val <- c(unr, rep(0, length(ind)))
  id <- c(seq_along(unr), ind+0.5)
  return(matrix(val[order(id)], s, s))
}
#####
#TODO: calculate gradient as well.
#################################################################

bfgs_opt <- function()
{
  res <- optim(rep(.5,6), cost_fn, NULL, method = "BFGS")
  return(res)
}

nelder_mead_opt <- function()
{
  res <- optim(c(.1,.2,.3,.4,.5,.6), cost_fn, method = "Nelder-Mead",
               control = list(maxit = 20000))
  return(res)
}

sim_anneal_opt <- function()
{
  res <- optim(rep(.5,6), cost_fn, method = "SANN",
               control = list(maxit = 20000, temp = 20, parscale = rep(.5,6)))
  return(res)
}

#TODO levenberg_marquardt method...
LMA_opt <- function(betas)
{
  
}
